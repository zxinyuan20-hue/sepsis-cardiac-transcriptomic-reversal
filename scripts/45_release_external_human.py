"""Packages numpy/pandas/scipy. Independently audit cohort scores, input integrity,
all-gene mean contrasts and sample roles. Outputs Chinese report, readiness plan,
source manifests/current release. No extra significance analyses or drug re-ranking.
"""
import json,itertools,gzip,re
import numpy as np,pandas as pd
from scipy.stats import rankdata
from pipeline_utils import ROOT,sha256,write_json,now,run_standard_module
OUT=ROOT/'outputs/analysis/external_human';AUDIT=ROOT/'outputs/audit'
def j(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 if j(ROOT/'CURRENT_RELEASE.json').get('raw_heart_version'):
  return {'status':'HISTORICAL_RELEASE_PRESERVED','reason':'Newer raw-heart release is current; do not regenerate historical outputs or downgrade current release'}
 for name in ['43_prepare_external_human','44_external_human_analysis']:assert j(AUDIT/(name+'_status.json'))['status']=='PASS'
 for f in j(OUT/'input_freeze.json')['files']:assert sha256(ROOT/f['path'])==f['sha256']
 prep=j(OUT/'preparation_summary.json');summary=j(OUT/'external_human_summary.json')
 scores=pd.read_csv(OUT/'GSE141864_patient_signature_scores.tsv',sep='\t');genes=pd.read_csv(OUT/'GSE141864_gene_expression.tsv.gz',sep='\t',index_col=0);genes.index=genes.index.astype(str)
 bg=pd.read_csv(OUT/'shared_gene_background.tsv',sep='\t',dtype=str).entrez_id.tolist();x=genes.loc[bg,scores.gsm]
 rank=(rankdata(x.to_numpy(),axis=0,method='average')-.5)/len(x);u=pd.read_csv(ROOT/'outputs/analysis/human/frozen_signatures/landmarks_100_up.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id;d=pd.read_csv(ROOT/'outputs/analysis/human/frozen_signatures/landmarks_100_down.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id
 vv=rank[x.index.isin(u)].mean(axis=0)-rank[x.index.isin(d)].mean(axis=0);assert np.max(abs(vv-scores.primary_disease_rank_score))<1e-12
 iscase=scores.group.eq('sepsis').to_numpy();observed=vv[iscase].mean()-vv[~iscase].mean();null=[]
 for controls in itertools.combinations(range(7),2):
  keep=np.ones(7,bool);keep[list(controls)]=False;null.append(vv[keep].mean()-vv[~keep].mean())
 p=np.mean(np.asarray(null)>=observed-1e-14);assert p==summary['exact_one_sided_p'];assert np.isclose(min(2*p,1),summary['planned_two_cohort_BH_q'])
 genestats=pd.read_csv(OUT/'GSE141864_all_gene_results.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id');actual=genes[scores.loc[iscase,'gsm']].mean(axis=1)-genes[scores.loc[~iscase,'gsm']].mean(axis=1);assert np.max(abs(actual.loc[genestats.index]-genestats.logFC))<1e-10
 df=pd.read_csv(OUT/'GSE237861_individual_descriptive_concordance.tsv',sep='\t');assert len(df)==7;assert int(df.both_query_sets_majority_agree.sum())==summary['GSE237861_descriptive_patients_both_query_majorities']
 coverage=pd.read_csv(OUT/'GSE237861_count_coverage.tsv',sep='\t');assert len(coverage)==14 and coverage.NCBI_count_available.sum()==3
 runs=pd.read_csv(OUT/'GSE237861_heart_run_manifest.tsv',sep='\t');assert len(runs)==14 and runs.title.nunique()==14;total=sum(sum(int(n) for n in s.split(';')) for s in runs.fastq_bytes.astype(str))
 nextplan={'next_stage':'GSE237861_FROM_RAW_QUANTIFICATION','status':'NOT_STARTED','heart_samples':14,'compressed_fastq_bytes':total,'compressed_fastq_GB_decimal':total/1e9,'raw_downloaded':False,'reason':'NCBI counts cover2cases+1control; complete7+7needed','requirements':['Use all14verified heart runs, same reference annotation and gene-level pipeline','Verify read quality, adapters, 3prime stranded library orientation, alignment and gene assignment rates','Recover sex/clinical covariates from paper; preserve individuals and predefined sepsis/control groups','Freeze reconstruction reference versions, aligner/count options and technical QC before new disease results','Do not combine existing partial NCBI counts with a differently quantified remainder','Compare fixed original190gene direction first; do not choose signatures after seeing cohort results'],'stop_or_redirect':'If complete cohort also lacks reproducible direction, stop escalating current three-drug protective mechanism claims; retain transferability/heterogeneity analysis'}
 write_json(OUT/'GSE237861_raw_reconstruction_readiness.json',nextplan)
 # Existing frozen numerical results remain unchanged; editable docs/scripts are new versions.
 for f in j(AUDIT/'independent_evidence_release_manifest.json'):
  if f['path'].startswith('outputs/analysis/independent_evidence/') and not f['path'].endswith('sources/GSE141864_samples.soft'):assert sha256(ROOT/f['path'])==f['sha256']
 for f in j(AUDIT/'specificity_release_manifest.json'):
  if f['path'].startswith('outputs/analysis/specificity/'):assert sha256(ROOT/f['path'])==f['sha256']
 sfreeze=j(ROOT/'outputs/analysis/screening/screening_result_freeze.json');assert sha256(ROOT/'outputs/analysis/screening/complete_screening_results.tsv')==sfreeze['complete_results_sha256']
 sources=[]
 for p in (OUT/'sources').iterdir():
  if p.is_file():sources.append({'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha256(p)})
 write_json(OUT/'source_file_manifest.json',sources)
 lod=pd.read_csv(OUT/'GSE141864_leave_one_patient_out.tsv',sep='\t');pathways=pd.read_csv(OUT/'GSE141864_fixed_pathway_results.tsv',sep='\t');primary=pathways[pathways['mode'].eq('estimated_residual_correlation')]
 lines=['# 独立人类心肌复现：结果与下一步决定','','2026-10-02｜external-human-1.0。',
 '', '**结论：已完成GSE141864的独立患者层面检查，但未获得稳定复现支持。GSE237861取得了部分标准化计数及全部7位病例的作者差异结果，完整7+7组间检验仍需统一重建计数。原药物结论保持不变，不升级单药保护机制。**',
 '', '## 本轮新增了什么', '',
 'GSE141864：从作者提交的全部样本表达表中提取10份心肌记录，按患者与保存方式固定FFPE5例病例和2例对照为主要分析。另3份冷冻心肌来自已计入的病例，只列入技术背景，不重复计算患者。使用官方GPL17586注释，剔除多Entrez歧义映射，保留22,203基因，与发现队列共有17,948基因。原100上调基因全部覆盖，90下调基因覆盖87个。表达值为作者提交的RMA结果，没有声称重新处理CEL。',
 '', 'GSE237861：核实14份心肌样本和14个ENA测序运行。NCBI标准化计数下载接口实际可用，但文件只涵盖20份各组织样本，其中仅2份脓毒症心肌和1份对照心肌。因此纠正此前“未找到计数”的表述为“找到部分计数，完整心肌矩阵尚不齐全”。同时取得七个病例的作者差异结果。',
 '', '## GSE141864：主要疾病特征没有稳定复现', '',
 '| 指标 | 结果 |', '|---|---:|',
 f'| 独立主要样本 | 5病例＋2对照 |',f'| 病例减对照的疾病特征分数 | {100*observed:.3f}秩百分点 |',f'| 单侧精确P（全部21种患者标签分配） | {summary["exact_one_sided_p"]:.6f} |',f'| 预定两队列家族BH q | {summary["planned_two_cohort_BH_q"]:.6f} |',f'| 17,948共有基因logFC描述性Spearman | {summary["all_gene_logFC_spearman"]:.4f} |',f'| 全基因FDR<0.05基因数 | {summary["external_gene_FDR05"]} |',
 '', '疾病特征分数是在样本内计算原上调基因平均百分位秩减去原下调基因平均百分位秩，再比较病例与对照。秩百分点不是表达量变化百分比，更不是心功能改善。正均值只是方向；当前P未达0.05，不能称验证成功。',
 '', '预定家族含两个队列。GSE237861未完成，P=1仅为保守多重校正占位，不是该队列实际检验结果。即便仅看GSE141864未校正P=0.286，结论也不显著。观察性标签可交换假设不能消除病原、年龄、尸检、保存或临床状态的混杂。',
 '', f'原50/50基因敏感性平均差为{100*summary["query50_difference"]:.3f}秩百分点；转录簇中位数合并为{100*summary["median_collapse_difference"]:.3f}秩百分点，平均方向一致。但逐一去掉患者后，差值范围为{100*lod.effect.min():.3f}至{100*lod.effect.max():.3f}秩百分点，出现方向翻转。',
 '', '七位患者的分数如下，直接展示仅有两个对照的影响。分数为样本内基因秩差，绝对大小不跨不同基因背景比较。', '', '| 样本 | 患者 | 组别 | 分数（秩百分点） |','|---|---|---|---:|']
 for r in scores.itertuples():lines.append(f'| {r.gsm} | {r.patient_id} | {r.group} | {100*r.primary_disease_rank_score:.3f} |')
 lines += ['', '其中一位对照的分数高于全部病例，另一位明显较低；小样本对照差异直接影响总体均值。不能因其影响结论而删除任何对照。本轮所有预定患者均保留。',
 '', '## 五条既定机制通路同样没有提供稳健支持', '',
 f'沿用此前固定的五个集合，主要camera使用残差估计的基因间相关性，{int((primary.planned_10_test_BH_q<.05).sum())}条通过预定十检验家族BH。呼吸电子传递的竞争性集合方向为Up、P=0.272，与原发现中Down方向不一致；该集合估计相关性约0.704。',
 '', '将基因相关性固定为0.01的敏感性会得到很小的P值，但不能覆盖预先固定的主要分析，也不能忽视集合内高度相关。保留两种计算结果，避免用假定低相关的显著性替代稳健机制证据。camera的方向是相对背景的集合变化，不等于测得通路活性。',
 '', '## GSE237861：可见异质性，尚不能完成完整队列推断', '',
 '作者原表为每位病例与同性别共享对照的edgeR结果；treated/non_treated列是组名文字，不是样本表达值。各表含约1.6—1.8万基因，并保留大量不显著基因，不是仅列显著基因的小名单。逐病例描述性结果如下：', '', '| 病例 | 共有基因logFC相关 | 原上调方向一致 | 原下调方向一致 | 两部分均过半 |','|---|---:|---:|---:|---|']
 for r in df.itertuples():lines.append(f'| {r.patient} | {r.logFC_spearman_descriptive:.3f} | {r.up_same_direction_fraction:.1%} | {r.down_same_direction_fraction:.1%} | {"是" if r.both_query_sets_majority_agree else "否"} |')
 lines += ['', '只有3/7病例同时呈现两部分多数基因方向一致。这是描述性结果，不能把3/7解释为验证成功率，也不能因病例共享对照而合并P值或把七张差异表当成七个独立队列。不能据此断言GSE237861正式验证失败，因为完整统一分析尚未完成。',
 '', '## 现在的项目判断', '',
 '**目前的证据不足以支持一个稳定、跨人类队列的疾病特征，更不足以解释三药的保护机制。** 新分析进一步表明，简单把一个发现队列的特征用于化合物排名，会受到患者异质性、对照背景和样本处理条件影响。它不证明三个药物无效，也不证明所有脓毒症心肌机制均不存在。',
 '', '下一步若继续投入，最有价值且范围明确的工作是统一重建GSE237861全部14份心肌样本计数，并按固定规则完成7+7分析。其原始压缩FASTQ总量约4.144 GB，样本和下载链接已整理；本轮没有开始下载或重建。需要记录参考版本、链特异性、比对/计数方法及质量审核，不能把现有NCBI三样本计数与另一流程生成的剩余样本直接拼接。',
 '', '在这一步完成前，暂停扩大三药筛选、单药对接及保护机制写作。若完整队列仍无稳定疾病方向，应停止升级当前单药保护论文路线，保留疾病异质性和跨背景可迁移性评估。当前也不保证该比较路线具有足够的发表创新性。',
 '', '## 输出与核验', '',
 '结果目录：outputs/analysis/external_human。包含患者样本表、表达矩阵、注释选择、完整基因差异结果、所有标签排列、逐患者敏感性、五集合结果、逐病例描述性表、计数覆盖以及原始测序运行清单。',
 '', '独立复核使用SciPy重新生成秩分数、按对照组合重新穷举21次分配，并核对R的所有基因均值差、样本数和冻结输入。原始药物筛选、疾病特异性和上一阶段证据矩阵未改动。',
 '', '本轮入口页曾遇到浏览器验证，改用公开GEO文本、FTP及计数下载接口。平台全家族归档为7.63GB，下载前大小检查已拒绝，随后取得仅平台注释文件。失败信息保存在日志；不把错误网页当作数据来源。',
 '', '来源：[GSE141864](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141864)、[GSE237861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237861)、[ENA PRJNA996949](https://www.ebi.ac.uk/ena/browser/view/PRJNA996949)。未开展药物外部验证、湿实验或临床疗效验证，未投稿或上传自有数据。']
 report=AUDIT/'独立人类心肌复现与完整计数准备报告.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
 qa={'scipy_rank_reconstruction':'PASS','complementary_exact_permutations':'PASS','all_R_logFC_mean_differences':'PASS','input_freeze':'PASS','original_numerical_results':'UNCHANGED','GSE237861_count_coverage':'3of14; NOT_ANALYZABLE_AS_PLANNED_FULL_COHORT','GSE237861_raw_runs':14,'report_claims':'ONE_FORMAL_SMALL_COHORT_CHECK_PLUS_DESCRIPTIVE_PATIENT_CONTRASTS','submitted_expression_not_raw_CEL':True}
 write_json(AUDIT/'external_human_numeric_qa.json',qa)
 cp=ROOT/'CURRENT_RELEASE.json';snap=AUDIT/'pre_external_human_release_snapshot.json'
 if not snap.exists():snap.write_bytes(cp.read_bytes())
 current=j(cp);current.update(status='EXTERNAL_HUMAN_PARTIAL_CHECK_NO_STABLE_REPLICATION_FULL_COHORT_RECONSTRUCTION_PENDING',updated_utc=now(),external_human_version=summary['version'],external_human_report=report.relative_to(ROOT).as_posix(),report=report.relative_to(ROOT).as_posix(),continuation_decision='PAUSE_SINGLE_DRUG_MECHANISM_ESCALATION_PENDING_COMPLETE_INDEPENDENT_DISEASE_MATRIX',supported_candidates=0)
 current['module_checks'].update({'43_prepare_external_human':'PASS','44_external_human_analysis':'PASS'});write_json(cp,current)
 paths=[p for p in OUT.rglob('*') if p.is_file()]+[report,cp,ROOT/'README.md',AUDIT/'external_human_numeric_qa.json']
 for folder in ['scripts','config','docs']:paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
 write_json(AUDIT/'external_human_release_manifest.json',[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size} for p in paths])
 return {**summary,'QA':qa,'next':nextplan}
if __name__=='__main__':run_standard_module('45_release_external_human',main)
