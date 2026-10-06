"""Packages pandas/numpy/scipy/matplotlib. Input frozen specificity results.
Outputs independent numerical QA, decision report, figures/source data, state/manifest.
Reproducibility: no stochastic new inference; verify bootstrap draws and original hashes.
All paths from ROOT/config and all results in outputs; log via run_standard_module.
"""
import itertools,json,logging
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module

def readj(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 out=path('analysis')/'specificity';ex=path('analysis')/'expansion';audit=path('audit');figs=ROOT/'outputs/figures/specificity';figs.mkdir(exist_ok=True)
 cfg=readj(ROOT/'config/specificity_v1.json');summary=readj(out/'specificity_summary.json')
 states={n:readj(audit/f'{n}_status.json')['status'] for n in ['36_prepare_specificity','37_disease_specificity']};assert all(v=='PASS' for v in states.values())
 for f in readj(out/'specificity_plan_freeze.json')['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
 sf=readj(ROOT/'outputs/analysis/screening/screening_result_freeze.json')
 for f,k in [('complete_screening_results.tsv','complete_results_sha256'),('supported_shortlist.tsv','supported_shortlist_sha256'),('exploratory_top10.tsv','exploratory_shortlist_sha256')]:assert sha256(ROOT/'outputs/analysis/screening'/f)==sf[k]
 result=pd.read_csv(out/'specificity_results.tsv',sep='\t');donor=pd.read_csv(out/'donor_profile_scores.tsv',sep='\t');test=pd.read_csv(out/'component_exact_tests.tsv',sep='\t');boot=pd.read_csv(out/'human_bootstrap_endpoint_values.tsv.gz',sep='\t')
 common=pd.read_csv(out/'specificity_common_genes.tsv',sep='\t',dtype=str);profiles=pd.read_csv(out/'human_continuous_profiles.tsv',sep='\t',dtype={'ENTREZID':str})
 # Algebra independently checks both new comparator mean profiles against old contrasts.
 for target in ['IHD','DCM']:
  assert np.allclose(profiles[f'logFC_{target}_vs_nonfailing'],profiles.logFC_sepsis_vs_nonfailing-profiles[f'logFC_sepsis_vs_{target}'],atol=1e-12)
 # Weighted design-matrix reconstruction independently verifies all sample->donor vectors.
 pieces=[]
 for chunk in pd.read_csv(ex/'expanded_counts.tsv.gz',sep='\t',index_col=0,chunksize=1000):pieces.append(chunk.loc[chunk.index.isin(common.SYMBOL)])
 x=pd.concat(pieces).loc[common.SYMBOL];rank=(x.rank(axis=0,method='average')-.5)/len(x)
 saved=pd.read_csv(out/'donor_rank_change_vectors.tsv.gz',sep='\t',index_col=[0,1]);d=pd.read_csv(ex/'expanded_analysis_design.tsv',sep='\t')
 errors=[]
 for (drug,person),g in d.groupby(['drug_state','Cell']):
  c=pd.Series(0.,index=x.columns);n=g.block.nunique()
  for block,b in g.groupby('block'):
   for condition,sign in [('DRUG',1),('CTRL',-1)]:
    cols=b.loc[b.condition==condition,'column_id'];c.loc[cols]+=sign/(n*len(cols))
  v=rank.to_numpy()@c.to_numpy();err=float(np.max(abs(v-saved.loc[(drug,person),common.SYMBOL].to_numpy())));errors.append(err);assert err<1e-12
 # Independently enumerate subsets of negative donor signs for all three tests.
 checks=[]
 for r in result.itertuples():
  mat=donor[donor.drug_state==r.drug_state].pivot(index='donor',columns='profile',values='reversal')
  arrays=[mat.sepsis_vs_nonfailing.to_numpy(),(mat.sepsis_vs_nonfailing-mat.IHD_vs_nonfailing).to_numpy(),(mat.sepsis_vs_nonfailing-mat.DCM_vs_nonfailing).to_numpy()]
  ps=[]
  for name,a in zip(cfg['specificity_endpoints'],arrays):
   null=[]
   for k in range(len(a)+1):
    for idx in itertools.combinations(range(len(a)),k):null.append((a.sum()-2*a[list(idx)].sum())/len(a))
   p=float(np.mean(np.asarray(null)>=a.mean()-1e-14));expected=float(test[(test.drug_state==r.drug_state)&(test.endpoint==name)].exact_one_sided_p.iloc[0]);assert p==expected;ps.append(p)
  assert max(ps)==r.iut_p
  b=boot[boot.drug_state==r.drug_state];freq=(b[cfg['specificity_endpoints']]>0).all(axis=1).mean();assert np.isclose(freq,r.human_bootstrap_joint_positive_frequency)
  checks.append({'drug_state':r.drug_state,'all_component_P':'PASS','IUT_max_P':'PASS','bootstrap_frequency':'PASS'})
 p=result.iut_p.to_numpy();order=np.argsort(p,kind='stable');q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q0=np.empty_like(q);q0[order]=np.minimum(q,1)
 assert np.allclose(q0,result.iut_BH_q);assert np.allclose(np.minimum(q0*np.sum(1/np.arange(1,len(p)+1)),1),result.iut_BY_q)
 # Validate patient strata and several actual bootstrap weight reconstructions.
 h=pd.read_csv(ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz',sep='\t',index_col=0);h.index=h.index.astype(str);h=h.loc[common.ENTREZID]
 hs=pd.read_csv(ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv',sep='\t').set_index('gsm');draws=pd.read_csv(out/'human_bootstrap_draws.tsv.gz',sep='\t');bw=np.load(out/'human_bootstrap_weights.npz')['weights']
 assert len(draws)==4*cfg['bootstrap_iterations'] and draws.groupby('iteration').size().eq(4).all()
 for row in draws.itertuples():
  ids=row.sample_ids.split(';');assert len(ids)==int(hs.group.eq(row.group).sum()) and hs.loc[ids,'group'].eq(row.group).all()
 for b in [0,127,511,999]:
  means={r.group:h[r.sample_ids.split(';')].mean(axis=1).to_numpy() for r in draws[draws.iteration==b].itertuples()}
  for k,group in enumerate(['sepsis','IHD','DCM']):
   fc=means[group]-means['nonfailing'];pos=np.maximum(fc,0);neg=np.maximum(-fc,0);w=neg/neg.sum()-pos/pos.sum();assert np.allclose(w,bw[b,k],atol=1e-14)
 qa={'checks':checks,'sample_to_donor_vectors_max_error':max(errors),'comparator_contrast_algebra':'PASS','IUT_BH_BY_independent':'PASS','all_human_draw_strata':'PASS','four_bootstrap_weight_reconstructions':'PASS','original_482_results_unchanged':'PASS','frozen_specificity_inputs':'PASS'}
 write_json(audit/'specificity_numeric_checks.json',qa)
 passing=result[result.exploratory_followup_gate].sort_values('drug_name');passing.to_csv(out/'independent_followup_hypotheses.tsv',sep='\t',index=False)
 contract={'conclusion':'Three exploratory relative-specificity hypotheses remain, with zero formal family-corrected support.','backend':'Python matplotlib','archetype':'quantitative heatmap plus uncertainty comparison','size_mm':[183,170],'panels':{'a':'All21 mean sepsis reversal and advantages over IHD/DCM in rank points','b':'Human-profile bootstrap stability versus donor IUT BH q; different uncertainty sources'},'sources':['specificity_results.tsv'],'no_efficacy_claim':True,'heatmap_color_capping':'±3 rank points; numeric labels preserve full values'}
 write_json(audit/'specificity_figure_contract.json',contract)
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
 fig,(ax,bx)=plt.subplots(1,2,figsize=(183/25.4,170/25.4),gridspec_kw={'width_ratios':[1.15,1]},layout='constrained')
 show=result.sort_values('sepsis_reversal',ascending=False).reset_index(drop=True);cols=['sepsis_reversal','sepsis_minus_IHD','sepsis_minus_DCM'];vals=show[cols].to_numpy()*100
 im=ax.imshow(vals,aspect='auto',cmap='RdBu_r',vmin=-3,vmax=3)
 ax.set_yticks(range(len(show)),show.drug_name,fontsize=6.5);ax.set_xticks(range(3),['Sepsis\nreversal','Advantage\nvs IHD','Advantage\nvs DCM'],fontsize=6.5);ax.set_title('a  Relative disease specificity',loc='left',fontweight='bold',pad=12)
 for i,row in enumerate(show.itertuples()):
  if row.exploratory_followup_gate:ax.get_yticklabels()[i].set_fontweight('bold')
  for j in range(3):ax.text(j,i,f'{vals[i,j]:.2f}',ha='center',va='center',fontsize=6,color='white' if abs(vals[i,j])>1.9 else 'black')
 fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.85,pad=.05,label='Rank points; colors capped at ±3')
 other=result[~result.exploratory_followup_gate]
 bx.scatter(other.human_bootstrap_joint_positive_frequency,other.iut_BH_q,marker='x',s=25,color='#888888',lw=.8,label='Other drugs')
 bx.scatter(passing.human_bootstrap_joint_positive_frequency,passing.iut_BH_q,s=28,color='#276B8C',label='Exploratory follow-up')
 bx.axhline(.05,color='#AA5555',ls='--',lw=.8);bx.text(.02,.065,'Formal threshold q = 0.05',fontsize=6,color='#994444')
 for name,label,offset in [('Dasatinib','Dasatinib / Cabozantinib',(-4,-17)),('Ponatinib','Ponatinib',(-5,-8)),('Regorafenib','Regorafenib',(-5,-8))]:
  r=result[result.drug_name==name].iloc[0];bx.annotate(label,(r.human_bootstrap_joint_positive_frequency,r.iut_BH_q),xytext=offset,textcoords='offset points',ha='right',fontsize=6)
 bx.set(xlim=(-.05,1.08),ylim=(-.02,1.08),xlabel='Human-profile joint-positive frequency',ylabel='Donor IUT BH q (21 drugs)');bx.xaxis.set_major_formatter(PercentFormatter(1));bx.set_box_aspect(1.3)
 bx.set_title('b  Stability is not statistical support',loc='left',fontweight='bold',pad=12);bx.legend(loc='lower left',bbox_to_anchor=(0,.12),fontsize=6)
 for ext in ['png','pdf','svg','tiff']:fig.savefig(figs/f'disease_specificity_overview.{ext}',dpi=600,facecolor='white')
 plt.close(fig);result.to_csv(figs/'figure_source_data.tsv',sep='\t',index=False)
 caption='''# 疾病特异性检查图注

图：3个探索性相对特异性假设保留，但无药物通过正式家族校正。

a，所有21药物在共同876基因背景的供者等权平均分数。三列依次为脓毒症逆转分数、比IHD模式逆转的优势、比DCM模式逆转的优势。数值乘100为秩百分点，不是表达量百分比；颜色截于±3，但数字保留实际值。粗体药物名代表通过预先固定的探索跟进规则，不代表统计显著或治疗有效。该连续权重分数与上阶段100/90截断查询数值不应直接比较。

b，横轴为1000次人类分组患者重采样中，三个终点同时为正的频率；它只反映人类模式稳定性，不是疗效概率。纵轴为4—6个药物供者的三终点交并检验经21药物BH校正的q。两轴反映不同来源的不确定性。Dasatinib与cabozantinib同处横轴100%、q=0.328；后者因不满足此前平均双向模式前提而未进入探索跟进名单。所有q均高于0.05。

人类分组：sepsis20、IHD11、DCM9、nonfailing11；三种疾病共享同一重采样nonfailing对照。心肌药物沿用21药物配对设计，培养实验与孔不是独立供者。无新外部验证、疾病细胞救治或功能性安全证据。源数据figure_source_data.tsv，Python生成PDF/SVG与600dpi PNG/TIFF；当前为本地研究图，非指定期刊格式认证。
'''
 (figs/'figure_caption.md').write_text(caption,encoding='utf-8')
 lines=['# 疾病特异性检查与下一步决策','','2026-10-02；disease-specificity-1.0。',
 '', '**结论：正式特异性检验仍为0个阳性；预先固定的探索性独立跟进规则保留dasatinib、ponatinib、regorafenib三个假设。可以继续独立证据核查，不能升级为已验证的脓毒症治疗候选。**',
 '', '## 本阶段实际检验了什么', '',
 '人类GSE79962含20例脓毒症、11例IHD、9例DCM和11例nonfailing心肌。固定上一阶段共同landmark与人类RMA的交集876基因，全部21药物继续纳入。新评分使用全背景连续logFC权重，正/负权重各归一化为1，避免不同疾病显著基因数量造成截断偏差。它衡量相对表达模式，不比较疾病严重程度，也不是此前100/90名单分数的更新。旧分数与阴性结论保留。',
 '', '对每个药物要求三项同时为正：脓毒症逆转；逆转强于IHD模式；逆转强于DCM模式。每项均以4—6供者作单侧精确符号翻转，交并检验IUT取三个P的最大值，再对21药物BH；BY作保守敏感性。共享患者与药物对照、细胞背景、剂量及时间差异均保留。',
 '', '在这些876基因内，IHD与DCM的疾病logFC模式相关为0.807；脓毒症与IHD为0.333、与DCM为0.390。这支持模式存在差异，但只是描述性相关；取材、终末状态和其他混杂仍可能贡献差异，不能直接称生物学专一性。',
 '', '## 正式结果和探索跟进要分开', '',
 f'21药物中{summary["all_three_point_positive"]}个三个均值均为正；正式IUT BH q<0.05的药物为 **{summary["formal_specificity_supported"]}个**，最低q=**{summary["minimum_iut_BH_q"]:.6f}**。BY同样无阳性。',
 '', '在查看本阶段新终点前，另固定了限制独立跟进范围的规则：原扩展结果已有平均双向模式，三个主要均值均正，所有逐供者剔除、去人类QC标记样本及top50等权敏感性仍均正，且人类重采样联合正向频率至少80%。它不是另一个显著性阈值或治疗证据标准。',
 '', '| 探索跟进假设 | 脓毒症逆转 | 超过IHD的优势 | 超过DCM的优势 | 人类模式联合正向频率 | IUT BH q |', '|---|---:|---:|---:|---:|---:|']
 for r in passing.itertuples():lines.append(f'| {r.drug_name} | {100*r.sepsis_reversal:.3f} | {100*r.sepsis_minus_IHD:.3f} | {100*r.sepsis_minus_DCM:.3f} | {r.human_bootstrap_joint_positive_frequency:.1%} | {r.iut_BH_q:.6f} |')
 lines += ['', '前三项单位为秩百分点。不能把大小解释为功能改善幅度。频率100%只表示在本批人类患者重采样、固定药物供者平均效应的条件下，三个均值均为正；不表示药物100%有效，也不消除5个药物供者的不确定性。三个药物均对IHD/DCM也有正分数，属于相对偏好，不能称只作用于脓毒症。',
 '', '## 为什么其他线索没有保留', '',
 '- Sorafenib：本阶段相对IHD/DCM的优势更小，人类联合正向频率70.1%，且逐供者剔除不能始终保持三个方向；未达预定探索门槛。',
 '- Cabozantinib：本阶段方向与人类重采样较稳定，IUT q与dasatinib相同；但此前未呈完整平均双向模式，按事先规则不进入三药物名单。不能在看到本阶段结果后删除这个前提。',
 '- Imatinib：对IHD/DCM模式的逆转分数高于对脓毒症，未支持相对脓毒症优先。',
 '- Vemurafenib：对脓毒症的分数也低于IHD/DCM，联合正向频率仅1.1%；本轮进一步不支持把它作为脓毒症特异线索。',
 '', '![疾病特异性概览](../figures/specificity/disease_specificity_overview.png)',
 '', '## 稳定性和不确定性', '',
 '人类1000次重采样按四组独立抽取，并在每次三种疾病比较中复用同一nonfailing重采样均值。另检查逐一去掉药物供者、去掉原先标记的GSM2109172、改用同背景前50基因/方向等权评分。三个保留假设均满足既定方向规则；所有未通过者也完整保留。',
 '', 'sepsis/IHD与sepsis/DCM直接对比模式的描述性评分也已输出；三个保留药物在两项直接对比中均为正，但这不是新增独立证据。连续评分、旧截断查询、人类重采样和供者检验回答不同问题，没有择优替换原结果。',
 '', '最重要的限制：同一个人类队列内比较，仍无独立疾病队列复现；健康供者来源心肌培养物的药物反应，不是脓毒症损伤后的救治；药物供者仅4—6个；人类模式归一化比较形状，不能消除临床混杂或推断绝对通路活性。',
 '', '## 下一步的有边界决定', '',
 '**保留项目，但只把三个药物作为独立跟进假设。** 下一步应查找不重复使用这批患者或这批心肌药物样本的证据，核对暴露条件、作用靶点方向和心肌安全性，再决定是否收敛到一个药物或共同机制。尚不支持据此采购药物、提出临床应用或撰写治疗有效的机制结论。',
 '', '若找不到独立支持，或安全性/作用方向与保护假设冲突，应保留跨背景方法比较，停止“脓毒症特异单药机制”解释的升级。本阶段没有改变原482条目筛选、21药物主要筛选均无FDR支持候选的事实。',
 '', '## 核验与输出', '',
 '人类均值差复现原三个limma对比；新IHD/DCM对比满足对比代数关系。另用样本权重矩阵和pandas秩独立重建全部供者向量，子集穷举复算所有63个分量P，核对IUT、BH/BY及bootstrap频率；全部4000组抽样记录核对组别/样本数，并重建抽样权重。原始冻结结果哈希保持一致。',
 '', '输出：specificity_results.tsv为全部21药物；independent_followup_hypotheses.tsv为3个探索假设；component_exact_tests.tsv为63项分量检验；human_continuous_profiles.tsv为全部疾病模式；human_bootstrap_endpoint_values.tsv.gz与human_bootstrap_draws.tsv.gz保留重采样；specificity_sensitivities.tsv保留全量敏感性。',
 '', '来源：[GSE79962](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79962)、[GSE217421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421)、既有冻结的处理矩阵和实验设计。未新增湿实验、临床验证或投稿；自有数据未外传。']
 report=audit/'疾病特异性检查与三项探索假设报告.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
 receipt={'schemaVersion':1,'items':[{'id':'specificity','title':'疾病特异性：3个探索假设，0个正式阳性','queries':[{'id':'followup','source':{'label':'同一GSE79962内三疾病比较＋GSE217421供者层面药物反应','links':[{'label':'人类心肌数据','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79962'},{'label':'心肌药物数据','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421'}],'caveats':['探索门槛与正式显著性不同；所有IUT q≥0.328。','患者bootstrap频率不是疗效概率；同队列特异性不等于外部验证。']},'columns':['drug_name','sepsis_minus_IHD','sepsis_minus_DCM','human_bootstrap_joint_positive_frequency','iut_BH_q'],'rows':passing[['drug_name','sepsis_minus_IHD','sepsis_minus_DCM','human_bootstrap_joint_positive_frequency','iut_BH_q']].to_dict('records')}]}]}
 write_json(audit/'specificity_sources_receipt.json',receipt)
 old=ROOT/'CURRENT_RELEASE.json';archive=audit/'pre_specificity_release_snapshot.json'
 if not archive.exists():archive.write_bytes(old.read_bytes())
 status=readj(old);preserve_independent=bool(status.get('independent_evidence_version'));status.update(status='SPECIFICITY_REVIEW_COMPLETE_THREE_EXPLORATORY_HYPOTHESES_NO_FDR_SUPPORT',updated_utc=now(),specificity_version=cfg['version'],formal_specificity_supported=summary['formal_specificity_supported'],exploratory_specificity_followup=passing.drug_name.tolist(),supported_candidates=0,continuation_decision=summary['decision'],report=report.relative_to(ROOT).as_posix());status['module_checks'].update(states)
 if not preserve_independent:write_json(old,status)
 paths=[]
 for folder in ['scripts','config','docs','outputs/analysis/specificity','outputs/figures/specificity']:
  paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
 paths += [ROOT/'README.md',old,report,audit/'specificity_numeric_checks.json',audit/'specificity_sources_receipt.json',audit/'specificity_figure_contract.json']
 write_json(audit/'specificity_release_manifest.json',[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size} for p in paths])
 return status

if __name__=='__main__':run_standard_module('38_release_specificity',main)
