"""Packages: pandas, numpy, matplotlib. Inputs: completed modules 20-26.
Paths: centralized project ROOT and configured analysis/audit directories.
Reproducibility: no new inference or random sampling; verify frozen results.
Outputs: Chinese report, figure PNG/PDF, release manifest, source receipt, log.
"""
import json
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT, path, sha256, write_json, now, run_standard_module

def main():
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    audit=path('audit'); source=path('analysis')/'mechanism'
    names=['20_prepare_mechanism','21_cross_model','22_mechanism_sources','23_verify_mechanism_evidence','24_cardiac_resource_audit','25_prepare_cardiac_perturbation','26_chemical_standardization']
    states={n:json.loads((audit/f'{n}_status.json').read_text(encoding='utf-8'))['status'] for n in names}
    assert all(v=='PASS' for v in states.values())
    fixed=pd.read_csv(source/'pathway_evidence.tsv',sep='\t')
    before=pd.read_csv(source/'pathway_evidence_before_correlation_sensitivity.tsv',sep='\t')
    common=[c for c in before.columns if c in fixed.columns]
    pd.testing.assert_frame_equal(fixed[common],before[common])
    estimated=pd.read_csv(source/'pathway_correlation_sensitivity.tsv',sep='\t')
    assert len(estimated[estimated.primary & estimated.global_primary_FDR.notna()])==19
    assert (estimated.loc[estimated.primary,'global_primary_FDR']<.05).sum()==0
    screen=path('analysis')/'screening'
    freeze=json.loads((screen/'screening_result_freeze.json').read_text())
    for filename,key in [('complete_screening_results.tsv','complete_results_sha256'),('supported_shortlist.tsv','supported_shortlist_sha256'),('exploratory_top10.tsv','exploratory_shortlist_sha256')]:
        assert sha256(screen/filename)==freeze[key], filename
    checks=json.loads((audit/'cross_model_numeric_checks.json').read_text())
    assert all(x['gene_BH']=='PASS' and x['moderated_P']=='PASS' for x in checks)
    concordance=pd.read_csv(source/'cross_model_concordance.tsv',sep='\t')
    summary=json.loads((source/'cross_model_summary.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    cohorts=['human','GSE185754','GSE267388','LOCAL_RAT_LPS_2025']
    labels=['Human myocardium','Mouse GSE185754','Mouse GSE267388','Rat cardiomyocytes']
    pathways=['R-HSA-3371497','R-HSA-5673001','R-HSA-166016','R-HSA-611105','R-HSA-9020702']
    plabels=['HSP90 receptor cycle','RAF / MAPK','TLR4 cascade','Respiratory electron transport','IL-1 signaling']
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
    for ax,df,title in zip(axes,[fixed,estimated],['A  Prespecified correlation = 0.01','B  Estimated residual correlation (post-result)']):
        matrix=np.full((5,4),np.nan)
        for i,p in enumerate(pathways):
            for j,c in enumerate(cohorts):
                row=df[(df.pathway_id==p)&(df.cohort==c)].iloc[0]
                if pd.notna(row.global_primary_FDR):matrix[i,j]=(-1 if row.Direction=='Down' else 1)*min(4,-np.log10(row.global_primary_FDR))
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee')
        im=ax.imshow(matrix,cmap=cmap,vmin=-4,vmax=4,aspect='auto')
        ax.set_xticks(range(4),labels,rotation=35,ha='right',fontsize=9)
        ax.set_yticks(range(5),plabels);ax.set_title(title,fontsize=11)
        for i,p in enumerate(pathways):
            for j,c in enumerate(cohorts):
                row=df[(df.pathway_id==p)&(df.cohort==c)].iloc[0];q=row.global_primary_FDR
                text='Not eligible' if pd.isna(q) else f'{row.Direction}\nq={q:.2g}'
                ax.text(j,i,text,ha='center',va='center',fontsize=8,color='white' if np.isfinite(matrix[i,j]) and abs(matrix[i,j])>2.3 else 'black')
    fig.colorbar(im,ax=axes,shrink=.65,label='Direction × −log10(global BH q), capped at ±4')
    fig.suptitle('Pathway evidence depends on within-set gene correlation',fontsize=14)
    figs=ROOT/'outputs/figures/mechanism';figs.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(figs/f'pathway_correlation_comparison.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    lines=['# 跨模型机制与候选核查报告','','日期：2026-10-02；mechanism-1.0。',
    '', '**结论：本轮完成跨模型疾病表达比较、五组机制集合检验、探索候选身份与文献核查，以及人源心肌细胞扰动数据准备。没有候选升级为有统计支持的治疗药物，也没有获得稳健的跨模型机制验证。**',
    '', '## 1. 自有数据提供了什么', '',
    '大鼠心肌细胞为独立样本，LPS与对照各6份。主要分析保留全部12份，发现3,077个FDR<0.05差异基因（上调1,694，下调1,383）。这说明两组转录状态存在差异，但不能由此推断其完整复现人类脓毒症心肌病。',
    '', '| 与人类比较的模型 | 共同基因数 | logFC Spearman相关 | 人类显著基因同向比例 |', '|---|---:|---:|---:|']
    for r in concordance.itertuples():lines.append(f'| {r.cohort} | {r.shared_genes:,} | {r.logfc_spearman:.3f} | {r.same_direction_fraction_in_human_fdr:.2%} |')
    lines += ['', '相关性为描述性结果，基因不是相互独立的生物学重复。大鼠主分析相关性仅0.062；去掉预先标记的C1/C2后为0.051，未显著改善整体一致性。细胞来源、LPS浓度/时长、批次与StringTie计数生成细节尚待补齐；自有计数分析暂定位为探索性支持，没有心超或功能终点。',
    '', '## 2. 机制方向与统计稳健性', '',
    '在查看留出疾病效应前固定HSP90受体伴侣循环、RAF/MAPK、TLR4、呼吸电子传递及IL-1五组Reactome直接注释。各研究单独处理，计数采用TMM/voom/limma，保留小数计数；基因BH在各完整保留基因集合内完成，再作一对一同源映射。通路至少10个基因且覆盖原集合25%；大鼠TLR4仅8个，未进入检验。四个主要研究共19项通路检验统一BH。',
    '', '首轮camera使用预设相关系数0.01：人类HSP90循环、人类及大鼠呼吸电子传递下调，以及两个小鼠研究的IL-1上调达到校正阈值。',
    '', '**在看到首轮结果后，统一追加了估计残差相关性的敏感性分析，并明确记录为结果后分析。19项主要通路中没有任何一项仍达到FDR<0.05。** 人类呼吸电子传递q从1.62×10⁻³²变为0.145，大鼠从2.79×10⁻⁵变为0.897；其残差相关性分别约0.430、0.418。大鼠去掉C1/C2后的该通路q为0.829（单独敏感性检验家族）。',
    '', '因此可以保留“人类与大鼠呼吸电子传递相关表达同向下调”的线索，不能写成稳健的线粒体机制验证。mRNA变化也不等于蛋白活性、呼吸功能或治疗方向。',
    '', '![相关性假设比较](../figures/mechanism/pathway_correlation_comparison.png)',
    '', '## 3. 候选身份和心脏证据', '',
    '冻结的482化合物主要筛选结果未改变：0个通过主要FDR；前10仅为探索名单。10条源SMILES均可复算为对应LINCS InChIKey。六个名称与PubChem完整结构键一致，另外四个需要以下解释：',
    '', '- bisindolylmaleimide是含糊名称；源结构对应bisindolylmaleimide I / GF109203X，而不是仅凭通用名称检索得到的arcyriarubin A。',
    '- NVP-AUY922与luminespib的原始结构键不同，但经RDKit规范互变异构体处理后相同；应记录为互变异构等价，不能误判为错药。',
    '- KU-60019与LGX-818仅连接性键一致，立体化学尚未解决；不能自动替换为数据库中的指定异构体或获批制剂。',
    '', '结构计算核对不验证实际实验批次、纯度或制剂。',
    '', '**bisindolylmaleimide I存在与心肌场景直接相关的安全性和选择性问题。** [PMID12679140](https://pubmed.ncbi.nlm.nih.gov/12679140/)报告动作电位延长；[PMID21120453](https://pubmed.ncbi.nlm.nih.gov/21120453/)报告直接钾电流/hERG阻断（文献摘要中hERG EC50为0.76±0.04 µM）。这些为体外电生理证据，不能直接换算成人体风险，但足以阻止将探索排名第一解释为最值得治疗开发的药物。',
    '', '**HSP90作用有细胞背景和干预方式差异。** [2013年研究，全文已读](https://pubmed.ncbi.nlm.nih.gov/23425388/)提示Hsp90/Akt稳定性与心肌细胞保护相关；[2025年研究，摘要核查](https://pubmed.ncbi.nlm.nih.gov/40967330/)报告心脏微血管内皮细胞Hsp90aa1下调的保护效应。两者的细胞、模型与干预不同，不能合并成“全局抑制HSP90治疗脓毒症心肌病”。[17-DMAG研究，全文已读](https://pubmed.ncbi.nlm.nih.gov/27224288/)也不能替代luminespib本身的疗效证据。RAF/MEK类文献提供心脏安全性背景，不能外推每个候选具有相同毒性。',
    '', '本轮为PubMed/CrossRef靶向检索，非系统综述；9篇选定文献完成DOI/题名核验，取得4篇全文，其余使用摘要。完整检索式、检索失败及来源记录均保留。',
    '', '## 4. 已获取的人源心肌细胞数据', '',
    '[GSE217421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421)原始发布计数矩阵已下载并核对：38,478基因×1,171列，设计表与矩阵列完全匹配，无重复基因/列、非有限值、负值或小数计数。',
    '', '已提取vemurafenib 20份及CTRL 63份记录；六个供者来源细胞系中，只有五个具有药物处理，另外一个只有5份对照。药物暴露为2 µM、48小时，不同于LINCS的6/24小时。83份记录不等于83个独立供者。尚未实施药物差异分析；必须先明确供者、分化/实验批次、培养皿及孔的重复结构。',
    '', 'GSE174773主要为相关基线对照，不能当作独立药物重复验证；GSE253490涉及心肌/内皮共培养及其他药物。PXD014791为蛋白组资源，不与转录组混作同一种重复证据。',
    '', '## 5. 下一阶段与论文定位', '',
    '先审查GSE217421的重复结构，固定供者/实验层面的分析方案，再研究vemurafenib在心肌细胞中的表达方向与潜在损伤背景。选择它是因为冻结探索名单与可用心肌资源存在交集，不是因为它已通过药物筛选。匹配供者主分析、批次处理和多重检验应在查看药物效应前固定；健康细胞药物扰动也不能直接证明逆转患病细胞。',
    '', '目前更合适的工作题目是：**脓毒症相关心肌转录特征的跨模型可迁移性与候选化合物逆转证据评估**。能否形成完整论文，取决于心肌扰动复核的结果与新增信息量；当前不宜写“某药通过某通路治疗脓毒症心肌病”。',
    '', '## 6. 复核及输出', '',
    '模块20—26运行通过；独立复核各模型基因BH和moderated P通过。新增敏感性前后原始通路结果一致；冻结药物全表、支持名单及探索名单哈希一致。原始/自有数据保留且未外传。本次未进行湿实验、临床疗效验证或投稿。',
    '', '- pathway_evidence.tsv：原始固定相关性通路结果。', '- pathway_correlation_sensitivity.tsv：结果后实际相关性敏感性，必须与前表共同解释。', '- cross_model_concordance.tsv：跨模型整体方向一致性。', '- reviewed_compound_identity.tsv：最终结构身份解释。', '- verified_mechanism_literature.json：选定文献与来源核验。', '- outputs/preparation/GSE217421/：已准备的药物/对照子集及设计表。', '- mechanism_release_manifest.json：代码、配置、文档、分析和图表文件校验和。']
    report=audit/'跨模型机制与候选核查报告.md';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    old=ROOT/'CURRENT_RELEASE.json';archive=audit/'pre_mechanism_release_snapshot.json'
    if not archive.exists():archive.write_bytes(old.read_bytes())
    status={'status':'CROSS_MODEL_MECHANISM_REVIEW_COMPLETE_NO_ROBUST_THERAPEUTIC_HIT','updated_utc':now(),'protocol_version':'0.2.0','screening_version':'screening-1.0','mechanism_version':'mechanism-1.0','module_checks':states,'supported_candidates':0,'correlation_robust_primary_pathways':0,'rat_effects':'EXPLORATORY_ANALYSIS_COMPLETE_METADATA_INCOMPLETE','cardiac_perturbation':'INPUT_QC_PASSED_REPEAT_STRUCTURE_PENDING','external_drug_replication':'NOT_PERFORMED','clinical_efficacy':'NOT_ESTABLISHED','report':report.relative_to(ROOT).as_posix()}
    current=json.loads(old.read_text(encoding='utf-8'))
    if not current.get('cardiac_version'):
        write_json(old,status)
    write_json(audit/'mechanism_release_checks.json',{'primary_pathway_table_unchanged':True,'frozen_screening_hashes_unchanged':True,'gene_statistics_independent_checks':'PASS','module_checks':states})
    receipt={'schemaVersion':1,'items':[{'id':'mechanism','title':'跨模型分析及心肌扰动资源','queries':[{'id':'cross-model','source':{'label':'本地模型分析；原始固定相关性与结果后敏感性完整保留','links':[{'label':'人类发现数据GSE79962','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79962'}],'caveats':['自有大鼠数据为探索支持；估计基因相关性后0项主要通路通过FDR。']},'columns':['cohort','logfc_spearman','same_direction_fraction_in_human_fdr'],'rows':concordance[['cohort','logfc_spearman','same_direction_fraction_in_human_fdr']].to_dict('records')},{'id':'cardiac','source':{'label':'GSE217421官方计数与实验设计','links':[{'label':'GEO数据','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421'}],'caveats':['仅完成输入准备；尚未完成药物效应分析。样本列不等于独立供者。']},'columns':['VEM_records','CTRL_records','matched_donor_lines'],'rows':[{'VEM_records':20,'CTRL_records':63,'matched_donor_lines':5}]}]}]}
    write_json(audit/'mechanism_sources_receipt.json',receipt)
    inventory=[]
    for folder in ['scripts','config','docs','outputs/analysis/mechanism','outputs/figures/mechanism']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts:inventory.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size})
    for p in [report,old,audit/'mechanism_release_checks.json',audit/'mechanism_sources_receipt.json']:
        inventory.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size})
    write_json(audit/'mechanism_release_manifest.json',inventory)
    return status

if __name__=='__main__':run_standard_module('27_release_mechanism',main)
