"""Build traceable preparation report, source inventory and companion notebook.
Inputs: completed stage 00-06 records and immutable raw/reference files.
Outputs: outputs/audit, preparation report, README, environment lock inventory.
No scientific analysis is executed; reproducibility seed is study config.
"""
import csv,datetime as dt,importlib.metadata,json,logging,platform,subprocess,sys
from pathlib import Path
import nbformat
from pipeline_utils import ROOT,CONFIG,path,sha256,write_json,run_standard_module,now

def read(name):return json.loads((path("audit")/name).read_text(encoding="utf-8"))

def main():
    required=["00_environment_audit","01_acquire_public_data","02_validate_inputs","03_prepare_r_environment","04_prepare_lincs_matrix","05_check_r_functions","06_acquire_reference_resources"]
    states={n:read(n+"_status.json")["status"] for n in required}
    if any(v!="PASS" for v in states.values()):raise ValueError(f"Unresolved module failures: {states}")
    command=[sys.executable,"-m","unittest","discover","-s",str(ROOT/"tests"),"-v"]
    test=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,encoding="utf-8",errors="replace")
    (path("logs")/"07_input_contract_tests.log").write_text(test.stdout+test.stderr,encoding="utf-8")
    if test.returncode:raise ValueError("Input contract tests failed")
    resources=[]
    for base in [path("raw"),path("reference")]:
        for side in base.rglob("*.provenance.json"):
            meta=json.loads(side.read_text(encoding="utf-8"));target=ROOT/meta["relative_path"]
            if not target.is_file() or target.stat().st_size!=meta["bytes"] or sha256(target)!=meta["sha256"]:
                raise ValueError(f"Source integrity mismatch: {target.name}")
            resources.append(meta)
    fields=["relative_path","source_kind","bytes","sha256","url","retrieved_utc","last_modified"]
    with (path("audit")/"source_manifest.tsv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fields,delimiter="\t",extrasaction="ignore");writer.writeheader();writer.writerows(resources)
    geo=read("geo_input_qc.json");lincs=read("lincs_metadata_qc.json");matrix=read("lincs_matrix_qc.json");r=read("r_functional_checks.json");refs=read("reference_resources_qc.json")
    installed=sorted((d.metadata["Name"],d.version) for d in importlib.metadata.distributions())
    (ROOT/"environment/python_installed_freeze.txt").write_text("\n".join(f"{n}=={v}" for n,v in installed)+"\n",encoding="utf-8")
    inventory=[]
    for folder in ["config","scripts","docs","tests","environment"]:
        for p in (ROOT/folder).rglob("*"):
            if p.is_file() and "R-library" not in p.parts and "__pycache__" not in p.parts:
                inventory.append({"path":p.relative_to(ROOT).as_posix(),"sha256":sha256(p),"bytes":p.stat().st_size})
    write_json(path("audit")/"preparation_code_manifest.json",inventory)
    issues=[
      {"severity":"high","issue":"Human discovery uses fatal sepsis hearts; healthy-donor and cardiomyopathy comparators differ in procurement and clinical background.","resolution":"Retain scope; evaluate metadata; no causal or early-reversible disease claims."},
      {"severity":"medium","issue":"Clinical supplement covers 18 sepsis samples only.","resolution":"No invented full-cohort covariate adjustment."},
      {"severity":"medium","issue":"GSE267388 contains noninteger estimated counts.","resolution":"Preserve values; edgeR input smoke passed; verify quantification method before inference."},
      {"severity":"high","issue":"Single-cell resource is pooled and CD45 selected.","resolution":"Metadata-only optional context; no whole-heart proportion tests or cell-as-animal inference."},
      {"severity":"medium","issue":"LINCS cell systems are mostly noncardiac; quality/activity is not viability.","resolution":"Stratify and retain biological context; perform independent candidate safety review."},
      {"severity":"medium","issue":"Published study overlap exists.","resolution":"Focused novelty audit and independent candidate verification remain required."}]
    write_json(path("audit")/"open_scientific_issues.json",issues)
    size=sum(x["bytes"] for x in resources)
    status={"status":"READY_FOR_G2A_PREPROCESSING","prepared_utc":now(),"protocol_version":CONFIG["protocol_version"],
      "scope":"G0 design and G1 input/environment preparation","formal_analysis_enabled":False,"source_files":len(resources),"source_bytes":size,
      "module_checks":states,"input_contract_tests":"PASS","science_results":"NOT_RUN","wet_lab_validation":"NOT_PERFORMED","submission":"NOT_SUBMITTED"}
    current_path=ROOT/"CURRENT_RELEASE.json"
    current=json.loads(current_path.read_text(encoding="utf-8")) if current_path.exists() else {}
    # Rechecking preparation must not erase a later scientific-analysis release.
    if current.get("status") in (None,"READY_FOR_G2A_PREPROCESSING"):
        write_json(current_path,status)
    write_json(path("audit")/"preparation_readiness.json",status)
    text=f'''# 启动准备审计报告

日期：2026-10-02。状态：**可进入人类心肌预处理与正式质量审查**。本轮没有进行疾病差异分析、化合物排名或论文结果撰写。

## 已完成

- 研究工作方案v{CONFIG['protocol_version']}、数据角色、预设质量门槛与后续任务表。
- Python独立环境、R项目库、Python精确版本锁和R renv锁。
- {len(resources)}个官方来源文件，共{size/1024**3:.2f} GiB，逐文件本地SHA256复核；LINCS与官方SHA512匹配。
- 3个bulk表达矩阵与71个样本完成匹配；另保存10个单细胞文库元数据。
- 原始51份CEL归档核查，1份CEL真实读取测试通过；limma合成方向测试、edgeR估计计数输入测试通过。
- 人类平台注释初步覆盖全部{r['lincs_landmarks_total']}个LINCS landmark genes；此处尚未剔除多对多映射，不代表疾病特征已全部覆盖。
- LINCS矩阵真实尺寸{matrix['full_shape'][0]} × {matrix['full_shape'][1]}；分块提取{matrix['subset_shape'][0]} × {matrix['subset_shape'][1]}子集，所有值有限。
- {refs['one_to_one_pairs']}对操作性一对一人鼠NCBI同源关系、{refs['human_reactome_pathways']}个人类Reactome直接注释通路及{refs['literature_records']}条已核对参考文献。

## 表达输入概况

| 数据 | 特征×样本 | 分组 | 缺失/完全重复样本列 |
|---|---|---|---|
'''
    for x in geo["datasets"]:text+=f"| {x['dataset']} | {x['features']} × {x['samples']} | {x['groups']} | {x['missing_values']} / {len(x['exact_duplicate_sample_pairs'])} |\n"
    text+=f'''
## 化合物资源覆盖（不是筛选结果）

原始化学扰动条目{lincs['chemical_perturbagens']}个。按预设质量、剂量、时间条件，保留{lincs['eligible_signatures']}条signature，涉及{lincs['eligible_compounds_before_multicell_gate']}个化合物；其中{lincs['eligible_compounds_after_multicell_gate']}个达到至少2细胞系的覆盖要求。没有使用疾病表达结果筛选或挑选药物。

## 重要处理与限制

- GSE267388有18,996个小数计数元素，保留原值；不是将FPKM/TPM当原始计数。
- GSE70138没有is_hiq；采用方案中明示的操作性质量过滤。原始um单位、无列名索引及-666哨兵值已按真实格式处理。
- GSE190856包含CD45富集、动物混样、时间和基因型差异，未下载并混合所有单细胞数据用于推断。
- 人类数据来自死亡患者，补充临床信息仅覆盖18例脓毒症样本，不能据此宣称临床疗效或独立人类验证。
- 尚需执行人类RMA、异常样本审查、无歧义注释、冻结疾病signature以及正式的分层化合物筛选。
- Seurat/Scanpy、分子对接和分子动力学均为后续可选工具，当前核心计算无需安装这些大型依赖。

## 审计证据

source_manifest.tsv记录原始来源与SHA256；lincs_publisher_checksums.json为官方SHA512核对；geo_input_qc.json和lincs_matrix_qc.json为真实输入核查；r_functional_checks.json为软件功能测试；07_input_contract_tests.log为输入风险回归测试。失败和修复历史留在各模块追加日志中。

本报告通过的是本地启动准备，不是科学结果验证、湿实验、临床验证、官方格式审查或投稿。
'''
    (path("audit")/"准备审计报告.md").write_text(text,encoding="utf-8")
    nb=nbformat.v4.new_notebook()
    nb.metadata.kernelspec={"display_name":"Sepsis cardiac project","language":"python","name":"sepsis-cardiac-project"}
    nb.cells=[nbformat.v4.new_markdown_cell("# 公开数据启动QC\n本notebook读取已经执行并保存的QC结果。实际计算代码位于scripts/02_validate_inputs.py和04_prepare_lincs_matrix.py；没有新的科学结果。"),
      nbformat.v4.new_code_cell("from pathlib import Path\nimport json, pandas as pd\nroot = Path.cwd()\nif root.name == 'notebooks': root = root.parent\nassert (root / 'config/study_config.json').exists()\ndef load(name): return json.loads((root / 'outputs/audit' / name).read_text(encoding='utf-8'))"),
      nbformat.v4.new_code_cell("pd.DataFrame(load('geo_input_qc.json')['datasets'])"),
      nbformat.v4.new_code_cell("load('lincs_metadata_qc.json')"),
      nbformat.v4.new_code_cell("load('lincs_matrix_qc.json'), load('r_functional_checks.json')"),
      nbformat.v4.new_markdown_cell("需要重新执行输入检查时，在项目目录运行run_preparation.ps1；原始数据不会被覆盖。")]
    nbformat.write(nb,ROOT/"notebooks/01_输入数据审计.ipynb")
    # Safe receipt contains public source labels and reviewed counts, no local paths.
    receipt={"schemaVersion":1,"items":[{"id":"preparation-data","title":"公开数据与本地启动核查","queries":[
      {"id":"bulk-inputs","source":{"label":"NCBI GEO心肌表达数据","links":[{"label":"GSE79962官方数据目录","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79962/"},{"label":"GSE185754官方数据目录","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185754/"},{"label":"GSE267388官方数据目录","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267388/"}],"caveats":["人类心脏组织来自死亡脓毒症患者；小鼠研究提供跨模型支持。"]},"columns":["dataset","features","samples"],"rows":[{k:x[k] for k in ["dataset","features","samples"]} for x in geo["datasets"]]},
      {"id":"lincs-coverage","source":{"label":"GSE70138 LINCS Phase II","links":[{"label":"官方文件目录","url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/"}],"caveats":["7093条合格表达特征和482个跨细胞系覆盖化合物是资源覆盖统计，不是有效药物名单。"]},"summary":"完整扰动矩阵已与官方SHA512核对，并完成基因和signature ID匹配。","columns":["signatures","genes","eligible_signatures","multicell_compounds"],"rows":[{"signatures":lincs["signatures"],"genes":lincs["genes"],"eligible_signatures":lincs["eligible_signatures"],"multicell_compounds":lincs["eligible_compounds_after_multicell_gate"]}]}]}]}
    write_json(path("audit")/"sources_receipt.json",receipt)
    return status
if __name__=="__main__":run_standard_module("07_finalize_preparation",main)
