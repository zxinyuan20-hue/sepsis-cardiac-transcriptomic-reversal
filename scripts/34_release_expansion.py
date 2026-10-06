"""Packages pandas/numpy/matplotlib. Input completed expansion analyses and frozen
common background. Outputs decision report, full comparison figure/source tables,
independent ES/BH audit, release hashes and current state. Deterministic; no new tests.
"""
import json,logging
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module

def readj(p):return json.loads(p.read_text(encoding='utf-8'))
def reference_es(v,hit):
 order=np.argsort(-v,kind='stable');h=hit[order];w=np.abs(v[order]);n=h.sum()
 inc=np.full(len(v),-1/(len(v)-n));total=w[h].sum();inc[h]=w[h]/total if total>0 else 1/n
 walk=np.r_[0,np.cumsum(inc)];return float(walk.max() if walk.max()>=-walk.min() else walk.min())

def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 out=path('analysis')/'expansion';audit=path('audit');figs=ROOT/'outputs/figures/expansion';figs.mkdir(exist_ok=True)
 states={n:readj(audit/f'{n}_status.json')['status'] for n in ['31_expansion_coverage','32_prepare_expansion','33_expanded_cardiac','35_drug_class_sensitivity']};assert all(v=='PASS' for v in states.values())
 for item in readj(out/'expansion_plan_freeze.json')['files']:assert sha256(ROOT/item['path'])==item['sha256'],item['path']
 sf=readj(ROOT/'outputs/analysis/screening/screening_result_freeze.json')
 for f,k in [('complete_screening_results.tsv','complete_results_sha256'),('supported_shortlist.tsv','supported_shortlist_sha256'),('exploratory_top10.tsv','exploratory_shortlist_sha256')]:assert sha256(ROOT/'outputs/analysis/screening'/f)==sf[k]
 summary=readj(out/'expanded_results_summary.json');feas=readj(out/'feasibility_summary.json');prep=readj(out/'expansion_preparation_summary.json')
 comp=pd.read_csv(out/'cross_context_comparison.tsv',sep='\t');end=pd.read_csv(out/'expanded_endpoint_results.tsv',sep='\t');loo=pd.read_csv(out/'cross_context_leave_one_out.tsv',sep='\t')
 classes=pd.read_csv(out/'drug_class_sensitivity.tsv',sep='\t');kinase_rho=float(classes.loc[classes.scope=='within_class','harmonized_spearman'].iloc[0])
 # Independent BH by sorted ranks and reverse cumulative minimum.
 for size,g in end.groupby('size'):
  p=g.exact_one_sided_p.to_numpy();order=np.argsort(p,kind='stable');q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];result=np.empty_like(q);result[order]=np.minimum(q,1)
  assert np.allclose(result,g.family_BH_q,atol=1e-14)
 # Independent cumulative-walk ES for every real donor vector.
 common=pd.read_csv(out/'common_landmark_background.tsv',sep='\t',dtype=str)
 query={d:set(pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_100_{d}.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id) for d in ['up','down']}
 masks={d:common.ENTREZID.isin(ids).to_numpy() for d,ids in query.items()}
 card=pd.read_csv(out/'harmonized_cardiac_donor_scores.tsv',sep='\t');errors=[]
 for drug,g in card.groupby('state'):
  x=pd.read_csv(out/f'{drug}_donor_logFC.tsv.gz',sep='\t',index_col=0).loc[common.SYMBOL]
  for r in g.itertuples():
   u=reference_es(x[r.donor].to_numpy(),masks['up']);d=reference_es(x[r.donor].to_numpy(),masks['down']);score=(d-u)/2 if d*u<0 else 0
   errors.extend([abs(u-r.es_up),abs(d-r.es_down),abs(score-r.score)])
 assert max(errors)<1e-12
 check={'family_BH_independent_sorted_implementation':'PASS','common_ES_independent_walk':'PASS','real_donor_vectors_checked':len(card),'maximum_ES_error':max(errors),'all_frozen_expansion_inputs':'PASS','original_482_entry_results':'UNCHANGED','VEM_regression':readj(audit/'expanded_numeric_checks.json')['VEM_frozen_regression']}
 write_json(audit/'expansion_release_checks.json',check)
 contract={'conclusion':'Partial cross-context score ordering coexists with no FDR-supported cardiac reversal drug.','archetype':'quantitative comparison','backend':'Python matplotlib','dimensions_mm':[183,165],'panels':{'a':'All 21 donor-level mean scores with unadjusted t95% CIs; color denotes two set directions, not significance','b':'Same weighted ES on 897 landmarks in both contexts; 21 nonindependent drugs, descriptive rho only'},'sources':['cross_context_comparison.tsv'],'exports':['PDF','SVG','PNG','TIFF'],'inference':'4-6 donors per drug, shared controls; family BH=21; no clinical claims'}
 write_json(audit/'expansion_figure_contract.json',contract)
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
 fig,(ax,bx)=plt.subplots(1,2,figsize=(183/25.4,165/25.4),gridspec_kw={'width_ratios':[1.2,1]},layout='constrained')
 display=comp.sort_values('mean_reversal',ascending=False).reset_index(drop=True)
 for i,r in enumerate(display.itertuples()):
  color='#276B8C' if r.both_directions_opposing else '#888888'
  ax.errorbar(100*r.mean_reversal,i,xerr=[[100*(r.mean_reversal-r.t_ci_low)],[100*(r.t_ci_high-r.mean_reversal)]],fmt='o',color=color,capsize=2,ms=3,lw=.8)
 ax.axvline(0,ls='--',color='#555555',lw=.7);ax.set_yticks(range(len(display)),display.drug_name,fontsize=7);ax.invert_yaxis()
 ax.set_xlabel('Cardiac reversal (rank percentage points)');ax.set_title('a  All eligible drugs',loc='left',fontweight='bold',pad=12)
 ax.plot([],[],color='#276B8C',marker='o',ls='',ms=3,label='Both gene sets opposing');ax.plot([],[],color='#888888',marker='o',ls='',ms=3,label='Other set directions');ax.legend(loc='lower right',fontsize=6)
 bx.axhline(0,color='#999999',lw=.6,ls='--');bx.axvline(0,color='#999999',lw=.6,ls='--')
 for r in comp.itertuples():bx.scatter(r.harmonized_lincs_score,r.harmonized_cardiac_score,s=22,color='#276B8C' if r.both_directions_opposing else '#888888',edgecolor='white',lw=.3)
 offsets={'REG':(5,7),'VEM':(-7,9),'SOR':(-5,8),'CAR':(7,-12)}
 for code,(dx,dy) in offsets.items():
  r=comp[comp.state==code].iloc[0];bx.annotate(r.drug_name,(r.harmonized_lincs_score,r.harmonized_cardiac_score),xytext=(dx,dy),textcoords='offset points',ha='right' if dx<0 else 'left',fontsize=6)
 bx.set_xlabel('LINCS common-background ES');bx.set_ylabel('Cardiac common-background ES');bx.set_title('b  Harmonized scoring',loc='left',fontweight='bold',pad=12)
 bx.set_xlim(-.47,.31);bx.set_ylim(-.49,.49);bx.set_box_aspect(1.35)
 bx.text(.03,.97,f'Descriptive Spearman ρ = {summary["comparison"]["harmonized_ES_spearman"]:.3f}\n21 drugs; 897 common landmarks',transform=bx.transAxes,va='top',fontsize=6.5)
 for ext in ['png','pdf','svg','tiff']:fig.savefig(figs/f'multidrug_cross_context.{ext}',dpi=600,facecolor='white')
 plt.close(fig);comp.to_csv(figs/'figure_source_data.tsv',sep='\t',index=False)
 (figs/'figure_caption.md').write_text('''# 多药物跨背景比较图说明

图：21药物存在部分跨背景评分一致性，但没有药物达到本阶段主要多重检验阈值。

a，纳入全部21个结构和实验覆盖合格药物，按平均心肌逆转分数排序。点为供者等权均值，误差线为基于各药物4—6供者的未作多重校正的t分布95%区间。蓝色表示疾病UP集合下降且DOWN集合上升，灰色为其他组合；颜色不代表统计显著。主检验为供者层面穷举单侧符号翻转，21药物统一BH，全部q≥0.164。秩百分点是相对表达位置变化，不是表达量百分比或疗效大小。a中分数不能与b的ES直接比较绝对大小。

b，在所有药物共同通过过滤的897个实测landmark上，用相同93 UP/84 DOWN查询及加权ES、相反方向门槛评分。LINCS使用扰动z值，先细胞内再细胞间中位数；心肌使用供者logFC，汇总为供者中位数。未满足相反方向门槛记为0；0不等于没有表达变化。供者和对照在药物间共享，药物机制亦相关，Spearman ρ仅为描述量，不提供独立药物P或普通bootstrap区间。剂量、处理时间、平台仍有差异，不能将差异归因于细胞类型本身。

总计421份不同记录（360药物、61对照）；同一对照参与多药物配对，不能加总为独立样本。新增比较在已知VEM和原LINCS排名后固定，属于探索扩展。图表由Python生成，源数据为figure_source_data.tsv；PDF/SVG保留矢量文本，PNG/TIFF为600dpi，未确认具体期刊投稿格式。
''',encoding='utf-8')
 selected=comp[comp.exact_one_sided_p<.05].sort_values('mean_reversal',ascending=False)
 positive=comp[comp.both_directions_opposing & (comp.mean_reversal>0)].sort_values('mean_reversal',ascending=False)
 decision={'decision':'CONTINUE_COMPARATIVE_RESEARCH_WITH_BOUNDED_SCOPE_NO_THERAPEUTIC_ADVANCEMENT','basis':['21 identity/dose/repeat-eligible drugs enable comparison beyond a single compound','Partial descriptive rank agreement remains under shared gene background and metric','Six drugs have average bidirectional pattern but no family-corrected support','Top10 coverage limited to VEM; current resource cannot settle other top candidates'],'not_recommended':['Write a single-drug treatment-mechanism paper from these results','Add docking or tune significance thresholds to create positive results','Unlimited additional screens without new independent evidence'],'next_value_gate':'Disease specificity against non-sepsis failing myocardium and/or a genuinely independent human cardiac resource; exploratory signals need independent drug exposure and safety evidence','candidate_evidence_level':'hypothesis only','publishability':'not established by these results alone'}
 decision['post_result_class_caution']={'kinase_inhibitor_drugs':16,'within_class_spearman':kinase_rho,'interpretation':'Overall agreement partly reflects class composition; not reliable within-class prioritization','analysis_timing':'POST_RESULT'}
 write_json(out/'continuation_decision.json',decision)
 lines=['# 是否继续：21药物扩展比较与项目决策','','2026-10-02；cardiac-expansion-1.0。',
 '', '**判断：这一轮继续有价值，已从单药复核推进到21药物比较。项目可保留为跨背景可迁移性与候选证据评估研究；当前证据仍不足以推进单药治疗机制论文或认定有效药物。**',
 '', '## 为什么值得做这一轮', '',
 'VEM单药阴性不能回答整个筛选方法能否迁移。本轮先核对全部54种药物状态，完成47个非生物制剂名称的PubChem结构查询，再依据结构、剂量及配对覆盖确定全集。21药物满足门槛，包含原评分为正、零与负的条目，没有挑选最有利药物。扩展方案在查看其他药物效应前冻结；VEM与原排名已经已知，明确记录为结果后扩展。',
 '', '## 真实覆盖与纳排', '',
 '- 1,171个GEO样本与计数设计键全部对应；22种药物与原库有完整结构键支持。',
 '- 其中vandetanib存在333 nM/333 µM标注差异，暂不进入主要分析，剩余21药物。PAZ、VRP、EDN也有不同数值剂量，未自动按笔误修正。',
 '- Afatinib有两个同名原条目，结构只支持其中一个；使用该条目的原分数，没有按较高排名选择。Azacytidine与azacitidine按结构对应。Daunorubicin、epirubicin、tofacitinib的同名条目立体化学不一致，保留为身份未解决。',
 '- 最初482个筛选单位更准确称“扰动化合物条目”，并不保证全部对应不同化学实体；原482项多重检验和全部排名保持原样。',
 '- 冻结前10名仍只有VEM被本资源覆盖，因此没有验证其余9个，也没有证明其无效。',
 '- 合计421份不同记录：360份药物、61份对照；对照跨药物重复使用。每药物4—6供者，至少每供者2个同培养实验配对。实际独立单位仍是供者，不是421个独立个体。',
 '', '## 新结果', '',
 f'21种药物中，{summary["positive_mean_drugs"]}种平均相对秩分数为正，{summary["positive_and_both_directions"]}种同时呈UP下降、DOWN上升的平均双向模式；{summary["nominal_primary_p_lt_05"]}种未校正单侧P<0.05，但**0种达到21药物BH q<0.05，最小q={comp.family_BH_q.min():.6f}**。',
 '', '| 未校正P<0.05药物（全部列出） | 供者数 | 正分数供者 | 平均逆转（秩百分点） | 单侧P | 21药物BH q | 平均双向模式 |', '|---|---:|---:|---:|---:|---:|---|']
 for r in selected.itertuples():lines.append(f'| {r.drug_name} | {r.donors} | {r.positive_donors}/{r.donors} | {100*r.mean_reversal:.3f} | {r.exact_one_sided_p:.5f} | {r.family_BH_q:.6f} | {"是" if r.both_directions_opposing else "否"} |')
 lines += ['', '这四行不是显著候选名单。其余两个具有平均双向模式的药物为imatinib与rosiglitazone，另外ponatinib虽无未校正P<0.05也有平均双向模式；完整六个双向药物为：'+', '.join(positive.drug_name)+ '。判断方向、未校正P和校正q回答不同问题，不能混为同一种证据。',
 '', 'Regorafenib、sorafenib等可作为下一轮独立验证的问题；它们的原LINCS评分并非显著阳性，因此不能写成“由原筛选成功预测并验证”。用药后表达相反不自动意味着保护心脏，更不证明这些药物适用于脓毒症。当前培养物来自健康供者，尚无LPS/脓毒症细胞被救治的证据。',
 '', '4个药物只有4供者，精确单侧P最小为0.0625；5供者最小0.03125，6供者最小0.015625。小样本与多重比较共同限制检出能力，不能把不显著写成确定无效，也不能通过改检验或阈值制造阳性。',
 '', '## 原筛选信息能否迁移', '',
 f'- 原LINCS分数与心肌主要秩分数的描述性Spearman相关为{summary["comparison"]["original_LINCS_vs_cardiac_rank_spearman"]:.3f}。',
 f'- 为减少评分方法和背景差异，在897个共同实测landmark、93 UP/84 DOWN基因上，重新用同一加权ES计算两端分数，相关为**{summary["comparison"]["harmonized_ES_spearman"]:.3f}**。',
 f'- 去掉已知VEM后为{summary["comparison"]["harmonized_ES_spearman_without_VEM"]:.3f}；逐一去掉供者范围{loo[loo.omission_type=="donor"].spearman.min():.3f}—{loo[loo.omission_type=="donor"].spearman.max():.3f}；逐一去掉药物为{loo[loo.omission_type=="drug"].spearman.min():.3f}—{loo[loo.omission_type=="drug"].spearman.max():.3f}。',
 '- 两端ES均非零的药物只有9个，其中5个同号（55.6%）。其余ES为0常表示两个集合未达到相反方向门槛，不是没有表达变化。',
 '', '因此，存在部分整体排序一致性，但个别药物方向不能由原排名可靠保证。差异还包含剂量、时间、平台和原细胞构成，不能归因于纯粹的“心肌特异效应”。21药物共享供者、对照与作用类别，不按独立药物计算相关P或普通bootstrap置信区间。',
 '', '![多药物比较](../figures/expansion/multidrug_cross_context.png)',
 '', '## 敏感性与核验', '',
 '主要100/90查询之外，50/50、150/90查询同样没有BH阳性。排除已知VEM后重算20药物BH，也没有阳性。供者层面的既有VEM分数与前阶段15条查询×供者结果一致，未用新方法替换旧阴性。',
 '', '各药物内TMM/logCPM与供者配对基因模型保留全部结果；基因BH仅是单药内部校正，不是21药物×全部基因的确证发现。本阶段没有依据这些次要基因结果新选机制。',
 '', '对21药物逐一独立复核样本秩、基因BH、moderated P、供者平均logFC及符号翻转P；新增纯累计和实现复算全部真实供者的ES，并独立排序复算家族BH。冻结输入、原482条目全表及名单校验通过。',
 '', '## 明确的继续/停止边界', '',
 '**继续的内容：** 聚焦“疾病特征驱动筛选的跨背景可迁移性”。当前可用工作题目为《脓毒症相关心肌转录特征的跨细胞背景药物逆转评估》。本轮新增了20个药物比较、化学身份约束、统一评分以及供者层面推断，已经比单药负结果更有信息量。',
 '', '**暂不推进的内容：** 单药治疗机制稿件、依据当前信号安排药物采购或功能结论，以及不限定范围地继续筛选。对接与分子动力学不能补足疾病逆转和疗效证据。',
 '', '**下一步只有能增加独立信息的工作才值得投入：** 核查这些表达逆转是否针对脓毒症特征，还是一般心衰/应激特征；或取得真正独立的人类心肌/药物资源。已有数据内的非脓毒症心肌比较只能作特异性检查，不能冒充新外部验证。若特异性或独立支持仍不足，应收束为方法/资源评估，或暂停“发现有效候选药物”的目标。当前结果本身尚不保证能形成可发表完整论文。',
 '', '来源：[GSE217421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421)、[原始研究全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11390749/)、PubChem逐药物结构查询、冻结GSE79962/GSE70138分析。原始研究包含细胞实验与多尺度分析，其实验疗效或安全性结论不能直接转移给本次二次分析。自有大鼠数据未外传，本轮未新增湿实验或临床验证。',
 '', '主要输出：all_drug_feasibility.tsv（54药物纳排）、compound_overlap_identity.tsv（结构）、expanded_endpoint_results.tsv（所有终点）、expanded_donor_scores.tsv（供者）、cross_context_comparison.tsv（统一比较）、continuation_decision.json（研究去向）。']
 classlines=['','## 结果后药物类别敏感性（必须共同解释）','',f'在21药物总体结果可见后，追加完整类别剔除及类内描述检查，见docs/12。16种激酶抑制剂内部的统一ES相关仅为 **{kinase_rho:.3f}**；去掉两种蛋白酶体抑制剂后，总体相关降至0.348。总体0.494不能被解释为可靠的同类候选排序能力，也不能作为治疗证据。完整五类别剔除结果全部保留，没有筛选有利类别。','', '这些为结果后描述性敏感性，不改变原ρ、主终点或BH结论。按MoA粗分组不等同于确定的毒性机制，也没有功能安全验证。']
 insert=lines.index('## 敏感性与核验');lines[insert:insert]=classlines+['']
 with (figs/'figure_caption.md').open('a',encoding='utf-8') as f:f.write(f'\n结果后药物类别检查：16种激酶抑制剂内部ρ={kinase_rho:.3f}，总体相关部分依赖类别构成；该追加分析不是效应前预设。来源drug_class_sensitivity.tsv。\n')
 # Keep the six-name statement exact without miscounting nominal-significant rows.
 lines=[line.replace('其余两个具有平均双向模式的药物为imatinib与rosiglitazone，另外ponatinib虽无未校正P<0.05也有平均双向模式；','未校正P未达0.05但具有平均双向模式的药物还包括ponatinib、imatinib和rosiglitazone；') for line in lines]
 report=audit/'21药物扩展比较与项目继续价值评估.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
 receipt={'schemaVersion':1,'items':[{'id':'expanded','title':'21药物跨背景比较：来源与结果','queries':[{'id':'full-family','source':{'label':'GSE217421＋冻结GSE70138评分；全部21结构与设计合格药物','links':[{'label':'心肌扰动数据','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421'},{'label':'原始论文全文','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11390749/'}],'caveats':['每药物4—6供者；共享对照；0个通过21药物BH。表中未校正P不能替代校正q。']},'columns':['drug_name','donors','mean_reversal','exact_one_sided_p','family_BH_q'],'rows':selected[['drug_name','donors','mean_reversal','exact_one_sided_p','family_BH_q']].to_dict('records')},{'id':'transfer','source':{'label':'相同897个landmark与加权ES的描述比较','links':[],'caveats':['药物间非独立；无相关性推断P；不能分离细胞、剂量、时间和平台影响。']},'columns':['drugs','spearman','both_nonzero','same_sign'],'rows':[{'drugs':21,'spearman':summary['comparison']['harmonized_ES_spearman'],'both_nonzero':9,'same_sign':5}]}]}]}
 receipt['items'][0]['queries'][1]['source']['caveats'].append(f'结果后检查：16种激酶抑制剂内部ρ={kinase_rho:.3f}，整体相关不代表同类药物排序可靠。')
 write_json(audit/'expansion_sources_receipt.json',receipt)
 old=ROOT/'CURRENT_RELEASE.json';archive=audit/'pre_expansion_release_snapshot.json'
 if not archive.exists():archive.write_bytes(old.read_bytes())
 status=readj(old);preserve_specificity=bool(status.get('specificity_version'));status.update(status='MULTIDRUG_COMPARISON_COMPLETE_PARTIAL_TRANSFER_NO_FDR_SUPPORTED_HIT',updated_utc=now(),expansion_version=summary['version'],expanded_drugs=21,expanded_BH_supported_drugs=0,expanded_harmonized_spearman=summary['comparison']['harmonized_ES_spearman'],supported_candidates=0,continuation_decision=decision['decision'],report=report.relative_to(ROOT).as_posix());status['module_checks'].update(states)
 if not preserve_specificity:write_json(old,status)
 status['post_result_kinase_class_spearman']=kinase_rho
 if not preserve_specificity:write_json(old,status)
 paths=[]
 for folder in ['scripts','config','docs','outputs/analysis/expansion','outputs/figures/expansion','data/reference/cardiac_expansion']:
  paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
 paths += [ROOT/'README.md',old,report,audit/'expansion_release_checks.json',audit/'expansion_sources_receipt.json',audit/'expansion_figure_contract.json']
 write_json(audit/'expansion_release_manifest.json',[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size} for p in paths])
 return status

if __name__=='__main__':run_standard_module('34_release_expansion',main)
