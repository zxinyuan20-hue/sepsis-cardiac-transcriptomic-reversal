"""Packages: pandas/numpy/matplotlib. Input fixed screening results and annotations.
Outputs: Chinese report, figures, version manifest, sources receipt, CURRENT_RELEASE.json.
Reports negative results and all missingness; no retuning or strong efficacy claims.
"""
import json,logging,importlib.metadata
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module

def main():
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    states={n:json.loads((path('audit')/f'{n}_status.json').read_text())['status'] for n in ['14_prepare_screening','15_human_bootstrap','16_compound_screen','17_screening_stability','18_candidate_annotation']}
    if any(x!='PASS' for x in states.values()):raise ValueError('Unresolved screening module')
    out=path('analysis')/'screening';freeze=json.loads((out/'screening_result_freeze.json').read_text());summary=json.loads((out/'primary_screen_summary.json').read_text())
    for filename,key in [('complete_screening_results.tsv','complete_results_sha256'),('supported_shortlist.tsv','supported_shortlist_sha256'),('exploratory_top10.tsv','exploratory_shortlist_sha256')]:
        if sha256(out/filename)!=freeze[key]:raise ValueError('Frozen screening results changed')
    rank=pd.read_csv(out/'annotated_complete_results.tsv',sep='\t');top=rank.head(10)
    a=np.load(out/'screening_inputs.npz');null=np.load(out/'matched_null_compound_scores.npy')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(13,8.4),layout='constrained')
    y=np.arange(len(top))
    axes[0,0].barh(y,top.reversal_score,color='#4477AA');axes[0,0].set_yticks(y,top.pert_iname,fontsize=8);axes[0,0].invert_yaxis();axes[0,0].set(title='A  Exploratory top 10: observed score',xlabel='Reversal score (positive = opposing expression)')
    axes[0,1].barh(y,top.bootstrap_top20_frequency,color='#228833');axes[0,1].set_yticks(y,top.pert_iname,fontsize=8);axes[0,1].invert_yaxis();axes[0,1].set(xlim=(0,1),title='B  Human bootstrap ranking stability',xlabel='Fraction in top 20 across 1,000 resamples');axes[0,1].axvline(.5,ls='--',lw=.7,color='gray')
    candidate=rank.sort_values('empirical_p').iloc[0];ix=np.flatnonzero(a['compounds']==candidate.pert_id)[0]
    axes[1,0].hist(null[:,ix],bins=40,color='#AAAAAA',edgecolor='white');axes[1,0].axvline(candidate.reversal_score,color='#CC6677',lw=2,label='Observed')
    axes[1,0].set(title=f'C  Matched null: {candidate.pert_iname}',xlabel='Reversal score',ylabel='Random gene-set draws');axes[1,0].legend(frameon=False);axes[1,0].text(.02,.95,f'P = {candidate.empirical_p:.4f}; BH q = {candidate.empirical_bh_q:.3f}',transform=axes[1,0].transAxes,va='top')
    cols=['score_size_50','score_size_100','score_size_150','score_without_qc_flag','score_submitted']
    limit=max(.35,float(np.abs(top[cols]).max().max()))
    im=axes[1,1].imshow(top[cols],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
    axes[1,1].set_xticks(range(5),['50','100 primary','150','Without QC flag','Submitted'],rotation=35,ha='right',fontsize=8)
    axes[1,1].set_yticks(y,top.pert_iname,fontsize=8);axes[1,1].set_title('D  Query and preprocessing sensitivities');fig.colorbar(im,ax=axes[1,1],label='Reversal score',shrink=.85)
    fig.suptitle(f'First compound screen: {len(rank)} compounds; {freeze["q_lt_05_positive"]} pass empirical FDR < 0.05',fontsize=13)
    figs=ROOT/'outputs/figures/screening';figs.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(figs/f'screening_overview.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    lines=['# 首轮化合物筛选报告','',
      '日期：2026-10-02。方案v0.2＋排名前固定的screening-1.0实施细则。','',
      f'**核心结果：{len(rank)}个化合物中，{summary["positive_compounds"]}个汇总分数支持表达逆转；经过10,000次匹配随机基因集及482项BH校正后，{freeze["q_lt_05_positive"]}个达到q<0.05。通过全部预设统计/稳定性条件的候选为{freeze["supported_gate_pass"]}个。**','',
      '这一结果不证明所有化合物都无效，但当前数据与检验不能支持显著候选药物的结论。探索名单不能直接作为疗效或治疗推荐。没有调低阈值或换方法追求阳性。','',
      '## 已完成的计算','',
      '- 筛选482个具备至少两个细胞系的化合物，使用6,769条合格LINCS表达扰动特征。原资源另有250个仅单细胞系覆盖的化合物，不纳入主要推断。',
      '- 固定人类100上调/90下调查询；加权双向富集，先细胞内中位数、后跨细胞中位数。非相反方向置零，并保留原始分量。',
      '- 10,000次均值×标准差分层随机基因集比较，竞争性单侧P，全部482个化合物一起BH校正。P分辨率约0.0001；不是官方CMap tau或临床疗效P值。',
      f'- 人类四组内样本bootstrap 1,000次，每次重跑完整基因limma和BH并重新选特征；{freeze["bootstrap_valid"]}次满足覆盖门槛，{freeze["bootstrap_failed_coverage"]}次失败。',
      '- 50/100/150查询大小、去掉QC标记样本、作者表达矩阵、时间/剂量分层、逐一去掉细胞系的诊断结果全部保存。',
      '- 自有大鼠细胞及GSE267388疾病效应未用于上述排序。','',
      '## 探索性前10（均未通过主要多重检验）','',
      '| 排名 | 名称 | 逆转分数 | 经验P | BH q | bootstrap前20频率 | 覆盖细胞系 |','|---:|---|---:|---:|---:|---:|---:|']
    for row in top.itertuples():lines.append(f'| {row.rank} | {row.pert_iname} | {row.reversal_score:.3f} | {row.empirical_p:.4f} | {row.empirical_bh_q:.3f} | {row.bootstrap_top20_frequency:.1%} | {row.n_cells} |')
    lines+=['','bootstrap排名稳定性与随机基因集显著性回答不同问题：一个化合物可以经常排名靠前，但分数仍不足以明显超出相同表达背景的随机集合。不能用稳定性替代FDR。',
      '','## 初步机制核查','',
      '采用Broad Drug Repurposing Hub 2020-03-24固定快照作精确名称匹配；这是历史数据库注释，不是当前审批状态，未完成结构身份和心脏毒性的独立核查。前10中6个获得精确名称注释：','',
      '| 化合物 | 数据库作用类别 | 记录靶点 |','|---|---|---|']
    for row in top[top.annotation_match].itertuples():lines.append(f'| {row.pert_iname} | {row.moa} | {str(row.target).replace("|", ", ")} |')
    lines+=['','bisindolylmaleimide、LGX-818、AZD-8330、MLN-2480在该快照中未精确匹配，保留缺失，不按相似名称补填靶点。','',
      '探索排名提示可以核查HSP90及RAF/MEK相关机制，但不构成这些通路应被抑制的证据。人类主比较中HSP90AA1与RAF1的mRNA反而下调；转录丰度不等于蛋白活性，不能将逆转分数拼接成简单的抑制机制。BRAF虽有符号映射，但不在最终主分析基因矩阵中，相应表达效应留空。','',
      '## 对项目的实际影响','',
      '当前应进入机制与外部证据核查，而不是开始围绕某一个“显著药物”写结论。可预先指定HSP90/RAF-MEK等探索问题后，在自有细胞数据及独立模型中检验疾病过程方向；这种检验也不能替代药物处理实验。',
      '','还需独立扰动发布/心肌相关细胞背景、化合物身份核对、功能方向与毒性证据。LINCS Phase I尚未完成，不能声称跨发布验证。对接与分子动力学不能补救目前的主要统计证据不足。',
      '','随机基因零分布保持分层数量但不保持基因模块相关性，适用于本次竞争性比较。它不是患者标签置换检验；这一局限和全部零/负结果应保留在论文中。',
      '','## 核验与可复现性','',
      '- 加权ES经40组独立参考实现比较，正/反方向和细胞平等聚合测试通过。',
      '- 重新用pandas聚合、复算BH和零分布超越计数、逐层核查随机查询上/下基因数，均通过。',
      '- null分批保存并附校验和；实现中的重复解压性能问题修复后复用了同一批10,000次结果，未改变统计规则。早于主模块结束启动的稳定性模块被前置检查阻止，已待主模块完成后重跑通过。',
      '- Numba仅安装在项目Python环境；当前精确版本另存screening_requirements.lock。数据未外传。',
      '','## 主要输出','',
      '- complete_screening_results.tsv：482个化合物的完整分数、经验检验、稳定性与敏感性。',
      '- supported_shortlist.tsv：通过全部条件的名单；本轮为空表，保留列结构。',
      '- exploratory_top10.tsv / annotated_top10.tsv：探索名单及带缺失标记的机制注释。',
      '- candidate_target_expression_context.tsv：记录靶点在人类三个比较中的表达背景。',
      '- screening_result_freeze.json：名单、配置及结果哈希。',
      '- outputs/figures/screening/screening_overview.png与PDF：可视化总览。']
    report=path('audit')/'首轮化合物筛选报告.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    freeze_file=ROOT/'environment/screening_requirements.lock'
    packages=sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions())
    freeze_file.write_text('\n'.join(f'{n}=={v}' for n,v in packages)+'\n',encoding='utf-8')
    archive=path('audit')/'pre_screening_release_snapshot.json'
    if not archive.exists():archive.write_bytes((ROOT/'CURRENT_RELEASE.json').read_bytes())
    status={'status':'PRIMARY_COMPOUND_SCREEN_COMPLETE_NO_FDR_SUPPORTED_HITS','updated_utc':now(),'protocol_version':'0.2.0','screening_version':'screening-1.0','module_checks':states,'compounds':len(rank),'null_iterations':freeze['null_iterations'],'bootstrap_iterations':freeze['bootstrap_valid'],'supported_candidates':freeze['supported_gate_pass'],'exploratory_candidates':len(top),'rat_effects':'HELD_OUT_NOT_RUN','phase1':'NOT_RUN','mechanism_annotation':'DATED_DATABASE_NAME_MATCH_ONLY','clinical_efficacy':'NOT_ESTABLISHED','report':report.relative_to(ROOT).as_posix()}
    current=json.loads((ROOT/'CURRENT_RELEASE.json').read_text(encoding='utf-8'))
    if not current.get('mechanism_version'):
        write_json(ROOT/'CURRENT_RELEASE.json',status)
    inventory=[]
    for folder in ['scripts','config','docs','outputs/analysis/screening']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts:inventory.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size})
    write_json(path('audit')/'screening_release_manifest.json',inventory)
    receipt={'schemaVersion':1,'items':[{'id':'screen','title':'本轮化合物筛选来源与核验','queries':[{'id':'screen-evidence','source':{'label':'GSE70138与冻结的人类GSE79962表达特征；本地独立统计复核','links':[{'label':'LINCS原始数据','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138'}],'caveats':['482个化合物；10,000随机查询与1,000人类bootstrap；0个通过主要FDR。探索排名不是疗效结论。']},'columns':['compound','score','P','BH_q','bootstrap_top20'],'rows':[{'compound':r.pert_iname,'score':round(r.reversal_score,4),'P':round(r.empirical_p,5),'BH_q':r.empirical_bh_q,'bootstrap_top20':r.bootstrap_top20_frequency} for r in top.head(5).itertuples()]},{'id':'mechanism','source':{'label':'Broad Drug Repurposing Hub 2020-03-24','links':[{'label':'固定注释快照','url':'https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt'}],'caveats':['仅精确名称匹配；未核验结构身份；临床阶段字段非当前审批状态。']},'columns':['exact_name_matched_top10','total_top10'],'rows':[{'exact_name_matched_top10':int(top.annotation_match.sum()),'total_top10':len(top)}]}]}]}
    write_json(path('audit')/'screening_sources_receipt.json',receipt)
    return status
if __name__=='__main__':run_standard_module('19_release_screening',main)
