"""Packages numpy/pandas/scipy/matplotlib. Input frozen cardiac results.
Outputs report, source-backed publication figures, numerical audit and release manifest.
No new scientific fitting; exhaustive deterministic verification; paths from ROOT/config.
"""
import itertools,json,logging,importlib.metadata
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module
def read_json(p):return json.loads(p.read_text(encoding='utf-8'))

def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 out=path('analysis')/'cardiac';audit=path('audit');figs=ROOT/'outputs/figures/cardiac';figs.mkdir(parents=True,exist_ok=True)
 states={n:read_json(audit/f'{n}_status.json')['status'] for n in ['28_cardiac_design','29_cardiac_effects']};assert all(v=='PASS' for v in states.values())
 cfg=read_json(ROOT/'config/cardiac_v1.json');freeze=read_json(out/'cardiac_plan_freeze.json')
 for r in freeze['files']:assert sha256(ROOT/r['path'])==r['sha256'],r['path']
 sf=read_json(ROOT/'outputs/analysis/screening/screening_result_freeze.json')
 for filename,key in [('complete_screening_results.tsv','complete_results_sha256'),('supported_shortlist.tsv','supported_shortlist_sha256'),('exploratory_top10.tsv','exploratory_shortlist_sha256')]:assert sha256(ROOT/'outputs/analysis/screening'/filename)==sf[key]
 endpoint=pd.read_csv(out/'endpoint_results.tsv',sep='\t');donor=pd.read_csv(out/'donor_reversal_scores.tsv',sep='\t');culture=pd.read_csv(out/'culture_reversal_scores.tsv',sep='\t');coverage=pd.read_csv(out/'query_coverage.tsv',sep='\t')
 scores=donor[(donor.scenario=='primary')&(donor.signature_size==100)].copy();p=endpoint[(endpoint.scenario=='primary')&(endpoint.signature_size==100)].iloc[0]
 # Independently enumerate subsets to negate rather than matrix multiplication.
 a=scores.reversal.to_numpy();null=[]
 for k in range(6):
  for subset in itertools.combinations(range(5),k):null.append((a.sum()-2*a[list(subset)].sum())/5)
 exact=float(np.mean(np.asarray(null)>=a.mean()-1e-14));assert exact==p.one_sided_exact_signflip_p
 assert np.allclose(scores.down_change-scores.up_change,scores.reversal,atol=1e-12)
 sample=pd.read_csv(out/'sample_rank_scores.tsv',sep='\t');sample=sample[(sample.scenario=='primary')&(sample.signature_size==100)]
 check=sample.groupby(['donor','block','condition']).disease_rank.mean().unstack();check['reversal']=check.CTRL-check.VEM
 direct=check.reversal.groupby('donor').mean();assert np.allclose(direct.loc[scores.donor],scores.reversal,atol=1e-12)
 design=read_json(out/'design_summary.json');genes=read_json(out/'gene_analysis_summary.json');con=pd.read_csv(out/'human_drug_concordance.tsv',sep='\t');cam=pd.read_csv(out/'primary_pathways.tsv',sep='\t')
 lod=pd.read_csv(out/'leave_one_donor_out.tsv',sep='\t');loc=pd.read_csv(out/'leave_one_culture_out.tsv',sep='\t')
 checks={'primary_exact_signflip_independent_subset_enumeration':'PASS','sample_to_culture_to_donor_endpoint':'PASS','directional_score_identity':'PASS','original_screening_unchanged':'PASS','frozen_cardiac_inputs_unchanged':'PASS','independent_donors':5,'exact_one_sided_p':exact}
 write_json(audit/'cardiac_release_checks.json',checks)
 # Figure contract is saved before drawing: no new data-driven endpoint choice.
 contract={'conclusion':'Small positive aggregate rank score does not establish significant or bidirectional reversal.','archetype':'asymmetric quantitative figure','backend':'Python matplotlib','width_mm':183,'height_mm':140,'panels':{'a':'Donor scores, nested culture scores, equal-donor mean with t95% CI','b':'Frozen query-size sensitivity: mean and donor t95% CI','c':'Human UP and DOWN set rank changes: individual donors and means'},'inference':'n=5 donors, culture points descriptive only; one primary exact signflip test','units':'rank percentage points (100 times rank fraction), not percent expression change','sources':['donor_reversal_scores.tsv','culture_reversal_scores.tsv','endpoint_results.tsv'],'export':['PDF','SVG','PNG','TIFF'],'status':'local result figure; no journal-specific submission validation'}
 write_json(audit/'cardiac_figure_contract.json',contract)
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False})
 fig=plt.figure(figsize=(183/25.4,140/25.4),layout='constrained');grid=fig.add_gridspec(2,2,width_ratios=[1.12,1])
 ax=fig.add_subplot(grid[:,0]);bax=fig.add_subplot(grid[0,1]);cax=fig.add_subplot(grid[1,1])
 cp=culture[culture.signature_size==100]
 for i,r in enumerate(scores.itertuples()):
  vals=cp[cp.donor==r.donor].reversal.to_numpy()*100
  ax.scatter(vals,i+np.linspace(-.13,.13,len(vals)),s=12,color='#B5B5B5',zorder=2)
  ax.scatter(r.reversal*100,i,s=38,color='#336699',zorder=3)
 ax.errorbar(p.mean_reversal*100,6,xerr=[[100*(p.mean_reversal-p.t_ci_low)],[100*(p.t_ci_high-p.mean_reversal)]],fmt='D',color='#111111',capsize=3,ms=4,lw=1.2)
 ax.axvline(0,color='#777777',ls='--',lw=.7);ax.set_yticks([0,1,2,3,4,6],[r.donor.split('-')[0] for r in scores.itertuples()]+['Mean (95% CI)']);ax.set_ylim(6.7,-.7)
 ax.set_xlabel('Reversal score (rank percentage points)');ax.set_title('a  Donors and matched cultures',loc='left',fontweight='bold',pad=12)
 ax.text(.02,.11,f'5 donors; 16 cultures\nOne-sided exact P = {exact:.3f}',transform=ax.transAxes,fontsize=7)
 ax.scatter([],[],s=38,color='#336699',label='Donor mean');ax.scatter([],[],s=12,color='#B5B5B5',label='Culture pair');ax.legend(loc='upper right',fontsize=6)
 sens=endpoint[endpoint.scenario=='primary'].set_index('signature_size').loc[[50,100,150]]
 for i,(size,r) in enumerate(sens.iterrows()):
  bax.errorbar(r.mean_reversal*100,i,xerr=[[100*(r.mean_reversal-r.t_ci_low)],[100*(r.t_ci_high-r.mean_reversal)]],fmt='o',color='#336699' if size==100 else '#777777',capsize=3,ms=4)
 bax.axvline(0,color='#777777',ls='--',lw=.7);bax.set_yticks(range(3),['50 / 50','100 / 90','150 / 90']);bax.set_ylim(2.6,-.6)
 bax.set_ylabel('Frozen human UP / DOWN genes');bax.set_xlabel('Reversal score (rank percentage points)');bax.set_title('b  Query sensitivity',loc='left',fontweight='bold',pad=12)
 for r in scores.itertuples():cax.plot([0,1],[100*r.up_change,100*r.down_change],color='#B5B5B5',marker='o',ms=3,lw=.7)
 cax.plot([0,1],[100*p.mean_up_change,100*p.mean_down_change],color='#336699',marker='D',ms=5,lw=1.5,label='Donor mean')
 cax.axhline(0,color='#777777',ls='--',lw=.7);cax.set_xticks([0,1],['Human UP set','Human DOWN set']);cax.set_xlim(-.35,1.35);cax.set_ylabel('VEM − CTRL rank change (points)')
 cax.set_title('c  Direction of each gene set',loc='left',fontweight='bold',pad=12)
 for ext in ['png','pdf','svg','tiff']:fig.savefig(figs/f'cardiac_reversal_overview.{ext}',dpi=600,facecolor='white')
 plt.close(fig)
 # Figures are derived directly from the full, saved analysis tables.
 for name,df in [('figure_donor_scores.tsv',scores),('figure_culture_scores.tsv',cp),('figure_query_sensitivity.tsv',sens.reset_index())]:df.to_csv(figs/name,sep='\t',index=False)
 caption='''# 人源心肌细胞扰动结果图说明

图：Vemurafenib在健康供者来源心肌细胞培养物中的表达变化未建立显著或双向疾病特征逆转。

a，主分析使用5个供者细胞系、16个同培养实验配对（16药物记录、24对照记录）。灰点为培养实验配对分数，仅显示供者内差异；蓝点为供者均值。黑菱形为五供者等权均值，误差线为基于5供者、4自由度t分布的95%置信区间。培养实验和孔不作为独立供者。主要单侧P来自全部32种供者符号翻转，零假设下要求独立供者效应对称；不是随机化因果检验。

b，预先冻结的50/50、100/90、150/90人类上调/下调基因查询；100/90为本阶段主要查询。点与误差线定义同a。实际可映射且通过表达过滤的基因数量另见query_coverage.tsv。敏感性不用于挑选显著方案。

c，每条灰线为一个供者，蓝线为供者均值。用药后上调集合与下调集合的平均相对秩均下降，未显示完整双向逆转。

分数为相对表达百分位秩之差乘100，单位为秩百分点，不是表达量百分比。正逆转分数表示DOWN集合相对于UP集合变化更高，并不要求两个集合分别向相反方向变化。这里主要分数置信区间跨0，单侧精确P=0.15625。无疗效或功能性心脏安全结论。

仅一个预先固定的主要终点；基因和机制分析另行BH校正。全部图表由Python生成，来源为同目录figure_*.tsv及完整分析表。PDF/SVG为矢量文字，PNG/TIFF为600dpi。当前为本地分析图，未针对指定期刊作投稿格式确认。
'''
 (figs/'figure_caption.md').write_text(caption,encoding='utf-8')
 lines=['# 人源心肌细胞药物扰动复核报告','','2026-10-02；cardiac-1.0。',
 '', '**本轮结论：vemurafenib在人源心肌细胞背景中未显示有统计支持的疾病特征逆转，也没有形成清楚的上调集合下降、下调集合上升的双向模式。原482个化合物没有主要FDR支持候选的结论保持不变。**',
 '', '## 数据与分析单位', '',
 '来源：[GSE217421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421)及[Hansen等，Nature Communications 2024](https://doi.org/10.1038/s41467-024-52145-4)。已获取并核对全文PMC11390749。原论文为包含细胞实验、转录组与多尺度分析的资源研究，不是单纯二次生信分析。培养物来自健康供者iPSC并具有细胞异质性，不等同于纯成熟心室心肌细胞。VEM暴露2 µM、48小时，CTRL为载体对照。',
 '', 'GEO配置将Cell对应Subject、Experiment对应Culture、Dish对应Replicate、Plate对应Measure、Well对应Sample。Experiment可用于培养实验匹配，但更细的分化/测序批次和随机分配情况并不完整；Plate均为0，无法据此估计板效应。',
 '', f'原子集83份记录（20 VEM、63 CTRL）。主要分析保留{design["primary_records"]}份：{design["primary_VEM"]} VEM、{design["primary_CTRL"]} CTRL，形成16个同供者、同培养实验配对，供者数为5。4份VEM缺少同培养实验对照而未进入主要分析；其余被排除记录为无同实验VEM的对照，包括无VEM供者。完整纳排表为analysis_design.tsv。',
 '', '培养实验内多个对照先等权平均，培养实验再在供者内等权平均，五供者最后等权汇总。分化实验或培养孔不扩充供者样本量；全部78份同五供者记录用于预定的非实验配对敏感性。',
 '', '## 固定终点及主要结果', '',
 '分析方案与输入哈希在查看本阶段药物效应前冻结，但VEM选择发生在原LINCS探索排名已知之后，应称外部细胞背景探索复核，不是独立确证性验证。人类主查询仍为100上调/90下调，本数据中覆盖93/85个。过滤后保留16,142个基因，精确一对一映射14,533个。',
 '', '每样本以全部保留基因作百分位秩背景。疾病分数为UP集合平均秩减DOWN集合平均秩；逆转分数为对照疾病分数减用药疾病分数。正数仅表示相对方向更相反，不必然表示每个集合或每个基因均逆转。分数不是CMap tau，也不是药物疗效。',
 '', '| 供者细胞系 | 配对培养实验数 | 逆转分数（秩百分点） |', '|---|---:|---:|']
 for r in scores.itertuples():lines.append(f'| {r.donor} | {design["cultures_per_donor"][r.donor]} | {100*r.reversal:.3f} |')
 lines += ['', f'五供者平均逆转分数为 **{100*p.mean_reversal:.3f}秩百分点**，t分布95%区间 **{100*p.t_ci_low:.3f}至{100*p.t_ci_high:.3f}**。4/5供者分数为正；单侧精确符号翻转P=**{exact:.5f}**，双侧P={p.two_sided_exact_signflip_p:.4f}。仅5个独立单位，采用完整32种符号组合；检验依赖零假设下供者差值对称，不能视作随机化治疗检验。',
 '', f'更重要的是，UP集合平均变化 **{100*p.mean_up_change:.3f}秩百分点**，DOWN集合也下降 **{100*p.mean_down_change:.3f}秩百分点**。因而平均分数略正来自前者下降更多，而不是完整双向逆转。不能把略正均值或4/5供者为正单独包装为有效。',
 '', '![心肌扰动复核](../figures/cardiac/cardiac_reversal_overview.png)',
 '', '## 敏感性及整体表达一致性', '', '| 方案 | 查询UP/DOWN | 平均逆转分数（秩百分点） | 单侧P（敏感性不作确证） |', '|---|---|---:|---:|']
 for r in endpoint.itertuples():lines.append(f'| {"同培养实验配对" if r.scenario=="primary" else "同供者全部培养实验"} | {r.signature_size}/{50 if r.signature_size==50 else 90} | {100*r.mean_reversal:.3f} | {r.one_sided_exact_signflip_p:.5f} |')
 lines += ['', f'逐一去掉供者后，均值范围为{100*lod.mean_reversal.min():.3f}至{100*lod.mean_reversal.max():.3f}秩百分点；逐一去掉培养实验后为{100*loc.mean_reversal.min():.3f}至{100*loc.mean_reversal.max():.3f}。这些为描述性影响分析，没有借去掉样本获得阳性。查询大小与非配对敏感性均未提供显著逆转证据。',
 '', '与人类疾病比较的共同基因有12,589个。药物logFC与人类sepsis/nonfailing logFC的Spearman相关为−0.024；在人类显著且可映射基因中，药物变化相反者为53.26%。sepsis/IHD为0.060，sepsis/DCM为−0.004，未显示明确、特异的整体疾病逆转结构。这里不把相关基因当独立重复计算显著性。',
 '', '## 基因与机制层面', '',
 '供者×条件平均logCPM形成10个表达谱，用供者固定效应+条件的limma-trend分析，原始残差自由度4，BH覆盖全部16,142基因。主分析仅CISH达到基因FDR<0.05（log2FC=0.709，q=0.0317；5供者方向均上调）。使用同五供者全部培养实验后其q=0.589，且该敏感性没有显著基因。单个主要分析阳性不足以推定药物作用机制，尤其在跨方案不稳定时。',
 '', '五个固定机制集合中TLR4只覆盖7个基因，不满足至少10个的门槛，未检验。其余4项采用实际残差基因相关性的camera检验，BH后均未显著：HSP90循环下调q=0.243；RAF/MAPK下调、呼吸电子传递上调、IL-1下调的q均约0.816。呼吸电子传递的药物方向与疾病下调相反，但没有统计支持，不能称“恢复线粒体功能”。',
 '', '## 对项目的决定', '',
 '维持“无有统计支持的治疗候选”的证据级别。当前不宜围绕vemurafenib单药撰写治疗机制论文；也不以分子对接、换富集方法或调低阈值补救主要结果。',
 '', '这一阶段增加了有价值的研究发现：LINCS非心肌背景中的探索排名，不自动转化为人源心肌背景中的稳健双向逆转。项目仍可围绕跨模型可迁移性与筛选可靠性推进，但完整论文需要足够的比较范围和方法信息量，不能仅凭这一个阴性药物复核保证可发表性。',
 '', '后续若扩展，应先审计剩余冻结候选在公开心肌扰动库中的覆盖，或取得另一独立药物扰动发布，固定扩展药物全集与检验家族，再进行比较。不得按已见结果挑选最有利药物或背景。已有自有大鼠数据继续用于疾病模型一致性，并不承担药物疗效验证。',
 '', '## 复核与可复现性', '',
 '- 主查询覆盖、输入哈希和原药物排名冻结校验通过。', '- R差异基因BH与moderated P经Python独立复算；样本→培养实验→供者的表达聚合及logFC均值经独立实现核对。', '- 主要符号翻转P以独立子集穷举重算为5/32，方向分量恒等式与供者分数一致。', '- 所有样本保留预设规则；没有依据药物效应剔除异常值。', '- 运行中修复了Windows默认文字编码及表索引名称读取问题，未改变科学规则；日志保留失败与重跑记录。', '- 没有进行湿实验、疾病细胞救治试验、功能性心脏安全性验证、临床验证或投稿。',
 '', '主要输出位于outputs/analysis/cardiac：analysis_design.tsv、cardiac_plan_freeze.json、endpoint_results.tsv、donor_reversal_scores.tsv、primary_gene_results.tsv、primary_pathways.tsv、human_drug_concordance.tsv。图表及源表位于outputs/figures/cardiac。']
 report=audit/'人源心肌细胞药物扰动复核报告.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
 receipt={'schemaVersion':1,'items':[{'id':'cardiac','title':'人源心肌药物扰动来源与复核','queries':[{'id':'endpoint','source':{'label':'GSE217421同培养实验配对；冻结GSE79962人类查询','links':[{'label':'GSE217421','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421'},{'label':'原始研究全文','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11390749/'}],'caveats':['5供者，16培养实验，40记录；32种供者符号翻转；单侧P=0.15625；非疾病救治实验。']},'columns':['donor','reversal'],'rows':scores[['donor','reversal']].to_dict('records')}]}]}
 write_json(audit/'cardiac_sources_receipt.json',receipt)
 old=ROOT/'CURRENT_RELEASE.json';archive=audit/'pre_cardiac_release_snapshot.json'
 if not archive.exists():archive.write_bytes(old.read_bytes())
 status=read_json(old);status.update(status='CARDIAC_PERTURBATION_REVIEW_COMPLETE_NO_SUPPORTED_REVERSAL',updated_utc=now(),cardiac_version=cfg['version'],cardiac_perturbation='DONOR_AWARE_ANALYSIS_COMPLETE_NO_SUPPORTED_BIDIRECTIONAL_REVERSAL',external_drug_replication='EXPLORATORY_CARDIAC_CONTEXT_TEST_COMPLETED_NOT_CONFIRMATORY',cardiac_donors=5,cardiac_primary_records=40,cardiac_primary_exact_p=exact,supported_candidates=0,clinical_efficacy='NOT_ESTABLISHED',report=report.relative_to(ROOT).as_posix());status['module_checks'].update(states)
 current=read_json(old)
 if not current.get('expansion_version'):
  write_json(old,status)
 lock=ROOT/'environment/cardiac_requirements.lock';lock.write_text('\n'.join(f'{n}=={v}' for n,v in sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions()))+'\n',encoding='utf-8')
 paths=[]
 for folder in ['scripts','config','docs','outputs/analysis/cardiac','outputs/figures/cardiac']:
  paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
 paths += [report,old,lock,audit/'cardiac_release_checks.json',audit/'cardiac_sources_receipt.json',audit/'cardiac_figure_contract.json',ROOT/'README.md']
 write_json(audit/'cardiac_release_manifest.json',[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size} for p in paths])
 return status

if __name__=='__main__':run_standard_module('30_release_cardiac',main)
