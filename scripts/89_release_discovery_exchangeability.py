"""Packages standard library/pandas and shared requests utilities.
Input: central metadata/scan configs and completed source/model audits.
Output: Chinese decision report, manuscript insertions, gate, source/asset manifest.
Reproducibility: verify raw CEL dates with independent byte extraction and frozen
inputs; no random tests. Pipeline audit -> report -> separate release registration.
"""
import gzip,json,struct,math,re
import pandas as pd
from pipeline_utils import ROOT,sha256,now,write_json,run_standard_module

def main():
 cfg=json.loads((ROOT/'config/discovery_exchangeability_v1.json').read_text(encoding='utf8'));out=ROOT/cfg['output'];pub=ROOT/cfg['publication'];pub.mkdir(parents=True,exist_ok=True)
 meta=pd.read_csv(ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv',sep='\t').set_index('gsm');scan=pd.read_csv(out/'cel_scan_metadata.tsv',sep='\t')
 key='affymetrix-scan-date'.encode('utf-16-be')
 for row in scan.itertuples(index=False):
  p=ROOT/meta.loc[row.gsm,'cel_path'];assert sha256(p)==meta.loc[row.gsm,'sha256']
  with gzip.open(p,'rb') as f:b=f.read(100000)
  start=b.index(key)+len(key);length=struct.unpack('>i',b[start:start+4])[0]
  date=b[start+4:start+4+length].decode('utf-16-be').rstrip('\0')
  assert date==getattr(row,'_2'),row.gsm
 for freeze in [out/'scope_freeze.json',out/'scan_sensitivity/amendment_freeze.json']:
  for r in json.loads(freeze.read_text(encoding='utf8'))['files']:assert sha256(ROOT/r['path'])==r['sha256']
 old=ROOT/'outputs/publication/robustness_amendment_v1/amendment_release_manifest.json'
 for r in json.loads(old.read_text(encoding='utf8'))['files']:assert sha256(ROOT/r['path'])==r['sha256']
 refs=json.loads((out/'reference_identity_checks.json').read_text(encoding='utf8'));assert sum(r['status']=='VERIFIED_IDENTITY' for r in refs)==11
 gate={'assessed_utc':now(),'status':'UNRESTRICTED_DISEASE_LABEL_PERMUTATION_NOT_JUSTIFIED','executed':False,'unrestricted_assignments':math.comb(31,20),'scan_day_restricted_assignments':math.comb(23,16)*math.comb(5,4),'scan_day_restriction_sufficient':False,'clinical_covariates':'18/20 sepsis records; 0/11 nonfailing; cannot fit complete-case age/sex adjustment across groups','procurement':'Sepsis rapid-autopsy vs nonfailing brain-dead donor; disease and procurement are structurally linked','within_procurement_label_assignments':1,'scan_covariate_sensitivity_completed':True,'next_route':'Retain results as exploratory; strengthen independent disease-context evidence; do not run unrestricted labels as confirmatory support'}
 write_json(out/'permutation_gate.json',gate)
 report='''# 发现队列临床资料与扫描日期审计

本轮结论：发现队列存在可利用的临床补充表和原始芯片扫描日期，已据此完成两项技术敏感性分析；现有资料仍不支持把未经限制的疾病标签置换解释为可靠的疾病特异性检验。正式支持药物仍为0，本轮没有重新筛药。

## 来源与临床覆盖

核对原始GEO临床工作簿及其下载校验值后，18条临床记录全部唯一匹配至20例脓毒症样本。缺失的2例为GSM2109167、GSM2109180。11例nonfailing、11例IHD和9例DCM均未在该表中获得逐例临床记录。

| 临床字段 | 脓毒症组可用数 | 非衰竭对照可用数 |
|---|---:|---:|
| 年龄与性别 | 18/20 | 0/11 |
| 列出的五种血管活性药物记录 | 各17/20 | 各0/11 |
| LVEF | 11/20 | 0/11 |
| LVEDD和LVESD | 各9/20 | 各0/11 |

18例已知年龄42–93岁，已知男9例、女9例；原文Table 1报告全20例男10例、女10例，不能据此把性别分配给两个具体缺失样本。逐例LVEF范围25%–75%，11例中4例低于50%。该阈值仅用于描述，不用于新建心肌病亚组或重新筛药。

原文报告12例有超声资料、LVEF范围20%–70%，超声至死亡平均间隔5.6±2.5天；与公开表可逐例匹配的11个LVEF值不完全一致。两个来源分别保留，没有插补缺失值或修改原文数字。H9c2本地数据不因此获得心超信息。

上轮“GEO characteristics只含condition和tissue”对那个字段本身正确，但不足以概括全部GEO资料。本轮补充临床工作簿审计，并不改写历史记录。

## 取材方式与交换性

原文Methods说明：脓毒症病例来自ICU死亡后的快速尸检，非衰竭对照来自脑死亡器官供者，IHD/DCM来自移植时切除的心脏。总体取材在死亡或移植后5–180分钟进行；没有可与全部GSM对应的个体取材时长。

疾病组与取材来源在本队列中相互关联。即使有扫描日期，也不能通过技术调整分离脓毒症、供者状态和取材条件的全部影响。缺失全部非衰竭对照的逐例年龄/性别意味着完整病例调整将失去对照组；不能使用均值插补、由表达推断性别或任意替代变量来宣称解决交换性。

## 原始芯片扫描日期

51个CEL文件逐一匹配原SHA256，读取Calvin通用头后，用独立字节定位再核对扫描日期。全部51个文件有日期。

| 扫描日期 | 脓毒症 | 非衰竭对照 | IHD | DCM |
|---|---:|---:|---:|---:|
| 2011-11-10 | 16 | 7 | 9 | 9 |
| 2011-11-11 | 4 | 1 | 2 | 0 |
| 2012-04-17 | 0 | 3 | 0 | 0 |

扫描日期是已测量技术变量，不等同于完整提取、建库或处理批次；与患者取材日期也不同。混合日期允许估计日期调整后的组间差异，但不证明临床交换性。

## 预先记录规则后的敏感性分析

在新模型计算前记录结果后方案修订，保留原51例RMA预处理和18,866基因映射，采用原limma经验贝叶斯设置。只拟合sepsis减nonfailing的比较，不重新选择筛药查询，不重新排序药物。

| 分析 | 样本数 | BH显著上调/下调 | 总DEG数 | 与原分析效应ρ | 原190基因方向一致 | 原190基因同方向且BH显著 |
|---|---:|---:|---:|---:|---:|---:|
| 原四组模型重现 | 51 | 1643/1462 | 3105 | 1 | 190 | 190 |
| 四组模型加扫描日期 | 51 | 1016/1322 | 2338 | 0.9269 | 190 | 144 |
| 限制共同扫描日期并调整日期 | 48 | 1001/1316 | 2317 | 0.9269 | 190 | 144 |

原3,105个显著基因在两项敏感性中均保留方向，仍同方向且显著的分别为2,084和2,067个。190个查询基因中46个不再达到BH阈值，不能将方向稳定写成显著性完全稳健。

两项调整分析的logFC在数值精度范围内相同：原来只含3例对照的单独日期，在全样本模型中已由该日期项吸收；删除这些样本主要改变残差和经验贝叶斯推断。这不是两个独立复现证据。限制样本分析沿用全样本RMA，是固定预处理下的敏感性，并未重新归一化。

## 标签置换的执行决定

31例病例与对照不限制置换共有84,672,315种分配；按扫描日期固定组数时，共有1,225,785种分配。但按扫描日期限制并没有解决取材来源的系统差异。若要求同取材来源内交换病例标签，则当前分组没有非平凡可交换方案（仅原分配1种）。因此本轮不运行疾病标签置换，也不以其可能得到的P值支持候选药物。

若未来得到可匹配的患者、取材时长和完整处理批次，应另立方案：明确条件零假设，比较协变量重叠与模型可识别性；协变量调整后的残差置换也需相应误差交换性假设，不能机械替代标签打乱。任何全流程置换应重拟合差异表达、BH、基因选择和药物聚合；不删除空查询，不临时补足显著基因；仅当检验统计量对所有分配有事前定义时才能执行。次数、精度和检验家族需在执行前固定。本轮没有生成置换P值或运行时间估计。

## 核验与论文定位

原模型的logFC、t、P和BH结果复现误差为0。两项新模型的logFC以Python普通最小二乘独立复算，最大误差分别约1.95×10⁻¹⁴和2.29×10⁻¹⁴。上一版本108项文件校验全部通过，未覆盖原分析和Word文件。

现有11条期刊文献的标题与所列PMID匹配，文中列出的DOI也匹配；第11条是GEO数据集引用，未算作期刊文献核验。更正文章标题的HTML标签已在比对时去除，未修改来源。这是标识符一致性核验，不代表本轮逐项重新审查了所有文献的论证内容。第一条正确PMID为28067713；会话中的数字误读已澄清，实际稿件没有该错误。

文章宜继续定位于脓毒症相关心肌表达改变和药物反转证据的稳健性。不能将全部发现队列称为经统一心功能标准确诊的脓毒症心肌病队列。扫描日期敏感性支持表达方向的稳定，同时暴露显著性强度与取材背景的限制。下一步优先将本轮证据纳入稿件并进行投稿前内容审查；若希望升级药物支持程度，需要独立、疾病相关的扰动证据，继续变换现有评分不会创造独立验证。

## 原始来源

- GEO GSE79962及Clinical_metadata.xls.gz：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79962
- Matkovich等原文：https://pmc.ncbi.nlm.nih.gov/articles/PMC5315660/ ，PMID 28067713，DOI 10.1097/CCM.0000000000002207。
- 原文关键位置：Methods Human hearts、Measurement of gene expression；Results Clinical characteristics；Table 1。对应本轮保存的article_text.txt段20、24、26、30、83。
'''
 (pub/'发现队列临床与扫描日期审计报告.md').write_text(report,encoding='utf8')
 insert='''# Manuscript additions from the discovery metadata audit

These paragraphs are prepared for the next manuscript revision. They are not yet incorporated into the previously frozen Word manuscript. Source reference [1] is already present.

## Methods addition

In a post-result sensitivity analysis specified before the new model fits, scan dates were extracted from the original CEL headers and verified for all 51 arrays. We refitted the four-group limma model with scan day as a categorical covariate, and repeated this analysis after restricting samples to dates containing both sepsis and nonfailing hearts. Both analyses retained the original 18,866-gene mapping and full-cohort RMA preprocessing. We assessed differential-expression counts and the direction and adjusted significance of the original frozen 100-UP/90-DOWN query, without replacing that query or reranking compounds. Scan day represented an observed technical characteristic, rather than a complete laboratory-batch designation.

## Results addition

Scan-date adjustment yielded 2,338 BH-significant genes, compared with 3,105 in the original analysis, and genome-wide log-fold changes remained correlated with the original effects (Spearman rho = 0.927). All 190 frozen query genes retained their original direction, while 144 retained both direction and BH-adjusted significance. Restriction to the two shared scan dates removed three nonfailing controls, leaving 48 samples; this analysis yielded 2,317 significant genes and the same 190 directional and 144 direction-plus-significance retentions. The adjusted effect estimates were numerically equivalent across these two analyses because the control-only scan date was already absorbed by the date term in the full model; their residual and moderated inference differed. These are related technical sensitivities, not independent replications.

## Limitations addition

The discovery cohort contrasted rapid-autopsy sepsis hearts with nonfailing organ-donor hearts, while ischemic and dilated cardiomyopathy samples were obtained at transplantation [1]. Public patient-level clinical annotations covered 18 of the 20 sepsis samples but none of the 11 nonfailing controls, preventing complete-case adjustment for age and sex across the primary comparison. Eleven individual LVEF values could be linked from the public clinical workbook, whereas the source publication described echocardiographic data for 12 patients; these sources were not treated as interchangeable. The cohort therefore represents sepsis-associated myocardial expression rather than a uniformly phenotyped population with established myocardial dysfunction. Scan-date adjustment preserved query direction but reduced adjusted-significance retention and could not resolve procurement-related confounding. A disease-label permutation was not performed because clinical exchangeability was not established; restriction by scan day alone would not address the group-specific tissue sources.
'''
 (pub/'正文补充段落_待并入.md').write_text(insert,encoding='utf8')
 write_json(out/'final_audit.json',{'reviewed_utc':now(),'CEL_independent_date_extractions':51,'previous_release_108_hashes_match':True,'frozen_scope_and_scan_inputs_match':True,'bibliographic_identity_verified':11,'disease_permutation_executed':False,'new_drug_ranking':False,'clinical_data_invented_or_imputed':False})
 pointer=ROOT/'DISCOVERY_EXCHANGEABILITY_RELEASE.json';write_json(pointer,{'version':cfg['version'],'date':now(),'report':(pub/'发现队列临床与扫描日期审计报告.md').relative_to(ROOT).as_posix(),'manuscript_additions':(pub/'正文补充段落_待并入.md').relative_to(ROOT).as_posix(),'status':'METADATA_AND_SCAN_SENSITIVITY_COMPLETE_PERMUTATION_NOT_JUSTIFIED','current_word_unchanged':True,'formal_drug_hits':0})
 files=[pointer,ROOT/'config/discovery_exchangeability_v1.json',ROOT/'config/discovery_scan_sensitivity_v1.json']
 for n in range(86,90):files+=list((ROOT/'scripts').glob(f'{n}_*'))
 manifest=pub/'release_manifest.json'
 files+=[p for d in [out,pub] for p in d.rglob('*') if p.is_file() and p!=manifest]
 rows=[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in sorted(set(files))];write_json(manifest,{'created_utc':now(),'files':rows})
 for r in rows:assert sha256(ROOT/r['path'])==r['sha256']
 return {'manifest_files':len(rows),'status':gate['status'],'report':str(pub/'发现队列临床与扫描日期审计报告.md')}

if __name__=='__main__':run_standard_module('89_release_discovery_exchangeability',main)
