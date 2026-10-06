"""Packages: numpy/pandas/scipy/matplotlib. Inputs: completed analysis and private QC.
Outputs: report, numerical audit, figures, hashes and CURRENT_RELEASE.json.
Seed: shared configuration. Independently verify contrasts/BH/t P-values and frozen inputs.
"""
import json,logging,subprocess,sys
import numpy as np
import pandas as pd
from scipy.stats import t as student_t, false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT,CONFIG,path,sha256,write_json,now,run_standard_module

def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))

def main():
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    states={n:read(f'outputs/audit/{n}_status.json')['status'] for n in ['08_import_local_rat','09_human_rma','10_preprocessing_qc','11_human_discovery','12_rat_orthology']}
    if any(s!='PASS' for s in states.values()):raise ValueError('Unresolved module state')
    out=path('analysis')/'human'
    matrix=pd.read_csv(ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz',sep='\t',index_col=0);matrix.index=matrix.index.astype(str)
    samples=pd.read_csv(ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv',sep='\t')
    audits=[];tables={}
    for contrast,reference in [('sepsis_vs_nonfailing','nonfailing'),('sepsis_vs_IHD','IHD'),('sepsis_vs_DCM','DCM')]:
        d=pd.read_csv(out/f'primary_{contrast}.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
        means=matrix.loc[:,samples.loc[samples.group=='sepsis','gsm']].mean(axis=1)-matrix.loc[:,samples.loc[samples.group==reference,'gsm']].mean(axis=1)
        contrast_error=float(np.max(np.abs(means.loc[d.index]-d.logFC)))
        if contrast_error>1e-10:raise ValueError('Contrast mismatch')
        if not np.allclose(d['adj.P.Val'],false_discovery_control(d['P.Value'].values),rtol=1e-9,atol=1e-14):raise ValueError('BH mismatch')
        if not np.allclose(d.t,d.logFC/d.moderated_SE,rtol=1e-10,atol=1e-10):raise ValueError('Moderated t/SE mismatch')
        if not np.allclose(d['P.Value'],2*student_t.sf(np.abs(d.t),d.df_total),rtol=1e-9,atol=1e-14):raise ValueError('Moderated P mismatch')
        audits.append({'contrast':contrast,'mean_difference_max_error':contrast_error,'BH_independent_check':'PASS','moderated_t_and_p_check':'PASS'});tables[contrast]=d
    freeze=read('outputs/analysis/human/frozen_signatures/signature_freeze_manifest.json')
    if freeze['study_config_sha256']!=sha256(ROOT/'config/study_config.json'):raise ValueError('Config changed since signature freeze')
    if freeze['human_matrix_sha256']!=sha256(ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz'):raise ValueError('Matrix changed since freeze')
    for entry in freeze['files']:
        if sha256(ROOT/entry['path'])!=entry['sha256']:raise ValueError('Frozen signature mismatch')
    # Meaningful existing input contract checks remain required.
    tests=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (path('logs')/'13_input_contract_tests.log').write_text(tests.stdout+tests.stderr,encoding='utf-8')
    if tests.returncode:raise ValueError('Input contract tests failed')
    write_json(path('audit')/'human_discovery_numeric_verification.json',audits)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for ax,(contrast,d) in zip(axes,tables.items()):
        sig=d['adj.P.Val']<CONFIG['analysis']['fdr'];colors=np.where(sig,np.where(d.logFC>0,'#CC6677','#4477AA'),'#CCCCCC')
        ax.scatter(d.logFC,-np.log10(d['adj.P.Val'].clip(lower=1e-300)),s=5,c=colors,alpha=.6,rasterized=True)
        ax.axhline(-np.log10(.05),color='#555555',lw=.7,ls='--')
        ax.set(title=contrast.replace('_',' '),xlabel='Log2 fold change',ylabel='-log10 BH adjusted P')
    folder=ROOT/'outputs/figures/discovery';folder.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(folder/f'human_contrasts.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    rma=read('outputs/preprocessing/human/rma_annotation_summary.json')
    rat=read('outputs/preprocessing/rat/normalization_qc.json')
    orth=read('outputs/preprocessing/rat/orthology_qc.json')
    stats=read('outputs/analysis/human/discovery_summary.json')
    sensitivity=read('outputs/analysis/human/sensitivity_concordance.json')
    lines=['# 人类发现分析与自有心肌细胞数据接入报告','',
      '日期：2026-10-02。方案版本v0.2.0。已完成真实表达分析；尚未运行药物连接分数、排名或疗效验证。','',
      '## 本轮完成','',
      f'- 人类51份CEL：RMA核心探针标准化，{rma["core_probes"]:,}个探针，去掉无Entrez映射/歧义映射并按全体均值折叠后{rma["unique_genes"]:,}个基因，覆盖{rma["landmarks_covered"]}个LINCS实测基因。',
      '- 自有大鼠心肌细胞：用户确认control 6、LPS 6个独立生物学样本，未上传公共数据库，无功能表型。原始结果只读导入、保存来源和SHA256。',
      f'- 大鼠30,454基因完成输入核查；低表达过滤后{rat["retained_genes"]:,}基因用于TMM/logCPM质量检查。未做组间差异、通路效应或候选物验证。',
      f'- 固定NCBI大鼠注释和人鼠（大鼠）同源资源，低表达过滤后与人类芯片共同覆盖{orth["filtered_genes_shared_with_human_array"]:,}个一对一映射基因。未映射基因未用符号猜配。','',
      '## 人类差异表达','',
      '四组limma标准eBayes；全部51个样本保留，每个比较分别对18,866基因作BH校正。FDR<0.05，无额外倍数阈值。正logFC表示脓毒症组更高。','',
      '| 比较 | 上调 | 下调 | 总计 |','|---|---:|---:|---:|']
    for contrast,v in stats['primary'].items():lines.append(f'| {contrast} | {v["fdr_up"]:,} | {v["fdr_down"]:,} | {v["fdr_up"]+v["fdr_down"]:,} |')
    lines+=['','不能据差异基因数量推断脓毒症特异性、治疗靶点或药物疗效；其他心肌病直接比较也不能消除临床/取材混杂。','',
      '## 质量标记与敏感性','',
      'GSM2109172全局相关性偏低，但原始强度不在样本范围外，无充分技术失败证据。主分析保留；另做去掉此样本的敏感性。此决定记录在差异分析之前。','',
      '| 敏感性 | 全基因logFC秩相关 | 主分析显著基因方向一致率 | 同方向且仍显著比例 |','|---|---:|---:|---:|']
    for v in sensitivity:lines.append(f'| {v["sensitivity"]} | {v["all_gene_logfc_spearman"]:.4f} | {v["direction_agreement_primary_fdr_genes"]:.1%} | {v["same_direction_and_fdr_fraction"]:.1%} |')
    lines+=['','方向一致率分母为主比较3,105个FDR显著基因。作者矩阵为相同样本的处理路线敏感性，不是独立队列验证。',
      '','大鼠PCA显示两组在PC1上分离，但C1/C2的深度和PC2位置不同于其他对照，尚缺批次元数据；保留全部样本，不做未知批次校正。仅PCA分离不能证明模型成立或所有差异都来自LPS。',
      '','## 冻结的筛药输入','',
      '主特征为100个上调、90个下调的FDR合格LINCS实测基因；下调不足100时不填入不显著基因。50版本为50/50，150版本为150/90。每方向覆盖均超过预设15基因门槛。',
      '','冻结文件、选取规则、配置/矩阵哈希位于outputs/analysis/human/frozen_signatures。下一阶段进行自定义双向加权富集、化合物聚合、零分布及重采样稳定性；尚无候选药物名单。',
      '','## 已执行核验','',
      '- Python重新核对三个比较的组均值差、BH校正和moderated t的P值，全部通过。',
      '- 冻结特征和源矩阵哈希复核通过；既有6项输入契约测试通过。',
      '- 人类/大鼠QC图已渲染；预处理图已逐图查看。历史失败及修复保留在追加日志。',
      '','## 仍需补充但不阻塞人类筛选的事项','',
      '大鼠细胞的具体来源/分离方式、LPS浓度和处理时长、实验及建库批次、StringTie计数矩阵生成方法。当前按细胞模型支持定位；没有心超并不阻止计算分析，但不能声称证实整心功能恢复。',
      '','人类发现是死亡患者心肌，部分临床协变量缺失；相应结论仍限于组织表达关联与候选机制。大鼠留出疾病效应和GSE267388留出效应均未用于选择人类特征。']
    report=path('audit')/'人类发现与自有数据接入报告.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    old=ROOT/'CURRENT_RELEASE.json'
    previous=path('audit')/'pre_discovery_release_snapshot.json'
    if old.exists() and not previous.exists():previous.write_bytes(old.read_bytes())
    status={'status':'HUMAN_SIGNATURE_FROZEN_READY_FOR_COMPOUND_SCREENING','updated_utc':now(),'protocol_version':CONFIG['protocol_version'],'formal_analysis_enabled':True,'module_checks':states,'human_samples':51,'private_rat_samples':12,'human_primary_differential_genes':3105,'primary_signature_up':100,'primary_signature_down':90,'compound_ranking':'NOT_RUN','private_rat_disease_testing':'HELD_OUT_NOT_RUN','wet_lab_validation':'NOT_PERFORMED','submission':'NOT_SUBMITTED','report':report.relative_to(ROOT).as_posix()}
    current=json.loads(old.read_text(encoding='utf-8')) if old.exists() else {}
    if not current.get('screening_version'):
        write_json(old,status)
    inventory=[]
    for folder in ['scripts','config','docs','outputs/analysis/human','outputs/preprocessing/rat','outputs/preprocessing/human']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts and 'cel' not in p.parts:
                inventory.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size})
    write_json(path('audit')/'discovery_release_manifest.json',inventory)
    return status
if __name__=='__main__':run_standard_module('13_release_discovery',main)
