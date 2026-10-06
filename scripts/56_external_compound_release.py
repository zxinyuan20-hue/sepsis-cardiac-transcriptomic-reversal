"""Packages: pandas, Pillow; standard-library XML/JSON/hash utilities.
Inputs: centralized external_compound_v1 config, frozen results and audit records.
Outputs: exposure table, export/visual QA, Chinese report, release manifest.
Reproducibility: shared seeded wrapper; no additional biological tests or selection.
"""
import json
import xml.etree.ElementTree as ET
import pandas as pd
from PIL import Image
from pipeline_utils import ROOT, now, sha256, write_json, run_standard_module

def readj(p):
    return json.loads(p.read_text(encoding='utf8'))

def main():
    cfg=readj(ROOT/'config/external_compound_v1.json')
    out=ROOT/cfg['out']; aud=ROOT/'outputs/audit'; fig=out/'figures'
    modules=['52_prepare_external_compound','53_external_compound_analysis','54_audit_external_compound','55_external_compound_figures']
    for name in modules:
        assert readj(aud/(name+'_status.json'))['status']=='PASS'
    for item in readj(out/'plan_freeze.json')['files']:
        assert sha256(ROOT/item['path'])==item['sha256'], item['path']
    summary=readj(out/'summary.json'); numeric=readj(out/'numeric_qa.json')
    assert numeric['status']=='PASS' and summary['primary_formal_hits']==0 and summary['escalation_gate_hits']==0
    results=pd.read_csv(out/'all_drug_results.tsv',sep='\t')
    cross=pd.read_csv(out/'cross_cohort_decision.tsv',sep='\t').set_index('drug_state')
    design=pd.read_csv(ROOT/'outputs/analysis/expansion/expanded_analysis_design.tsv',sep='\t')
    exposure=design.loc[design.condition.eq('DRUG'),['drug_state','drug_name','concentration','dose_normalized','hours']].drop_duplicates().sort_values('drug_state')
    assert len(exposure)==21 and exposure.drug_state.is_unique and exposure.hours.eq(48).all()
    exposure.to_csv(out/'drug_exposure_audit.tsv',sep='\t',index=False)
    exports=[]
    for stem in ['figure1_disease_replication','figure2_all_drug_contexts','figureS1_gene_and_pathway_evidence']:
        for ext in ['png','tiff']:
            p=fig/(stem+'.'+ext)
            with Image.open(p) as im:
                dpi=im.info['dpi']; assert all(abs(float(x)-600)<.1 for x in dpi)
                assert abs(im.width/float(dpi[0])*25.4-182.88)<.2
                exports.append({'file':p.relative_to(ROOT).as_posix(),'pixels':list(im.size),'dpi':[float(x) for x in dpi]})
        svg=ET.parse(fig/(stem+'.svg'))
        text_count=len(svg.findall('.//{http://www.w3.org/2000/svg}text'))
        assert text_count>10
        assert (fig/(stem+'.pdf')).read_bytes().startswith(b'%PDF')
        exports.append({'stem':stem,'SVG_text_elements':text_count,'PDF_header':'PASS'})
    qa={'status':'PASS','reviewed_utc':now(),'exports':exports,
        'visual_review':'All three final PNGs visually inspected during module55 review; text, legends, pathway table, cohort labels and panel alignment checked.',
        'closed_issue':'Figure2 colorbar originally reduced heatmap height; dedicated colorbar axes restored alignment with component panel; corrected export visually reviewed.',
        'numerical_audit':'numeric_qa.json PASS; source values and tests audited separately',
        'limitations':'PDF header and SVG text checked; no independent PDF raster rendering or target-journal format certification. No clinical validation.',
        'display_only_addition':'Supplementary fixed gene/pathway display introduces no new biological test or threshold.',
        'prior_run_fixes':'Module53 Windows UTF-8 decoding and R character Entrez indexing repaired; failed attempts retained in logs; successful result independently audited.'}
    write_json(out/'figure_release_qa.json',qa)
    primary=results[results.cohort.eq('GSE237861')].set_index('drug_state')
    secondary=results[results.cohort.eq('GSE141864')].set_index('drug_state')
    dloo=pd.read_csv(out/'donor_loo.tsv',sep='\t'); ploo=pd.read_csv(out/'patient_loo.tsv',sep='\t')
    def fraction(frame,drug,cohort):
        rows=frame[frame.drug_state.eq(drug)&frame.cohort.eq(cohort)]
        return f'{int(rows.both_positive.sum())}/{len(rows)}'
    lines=['# 独立疾病背景下21种心肌扰动药物复核与审计报告','',
        '2026-10-05｜external-compound-1.0｜主要分析、独立数值复核、图表检查均已完成。','',
        '**核心结论：21种药物中没有药物通过主要多重检验，也没有药物通过机制升级门槛。Ponatinib与Regorafenib的转录方向比Dasatinib稳定，但仍属于探索性信号；当前证据不能支持单药保护作用论文。**','',
        '## 研究问题与本轮范围','',
        '主题保持为脓毒症心肌转录组的化合物筛选与机制解释。本轮用独立疾病队列重新检验此前固定的全部21种心肌扰动药物，未只检验原先挑出的三药。疾病背景改变，健康人iPSC心肌细胞药物处理实验及其供者被复用，因此不能称为独立药效实验验证。此前筛选结果已知；本轮方案在计算本轮分数前冻结，属于筛选后的扩展分析。','',
        '## 数据与统计方法','',
        '- 主要疾病背景：GSE237861，7例脓毒症与7例非感染性危重症对照，完整原始心肌测序重建，edgeR模型调整性别。次要背景：GSE141864固定FFPE子集，5例脓毒症与2例对照；GSE79962作为发现背景描述性参照。',
        '- 三个疾病背景与原心肌扰动数据共有783个基因，未依据差异P值或效应大小筛选。药物表达的样本内秩仍以原876基因为背景；正负疾病权重分别在783基因交集中归一化。此加权药物分析与原100上调/90下调基因的疾病特征复现分数是不同分析。',
        '- 培养孔先在匹配培养区组内平均，再在供者内平均；供者是药物统计单位，共重建103条药物–供者对比向量，并非103位不同受试者。每药4–6名供者。全部药物为48小时处理，每药一个剂量；详见drug_exposure_audit.tsv。',
        '- 分开检验疾病上调基因抑制与下调基因恢复。主要IUT的P值为两个单侧穷尽供者符号翻转P值的最大值，随后对全部21药进行BH校正；BY为敏感性分析。检验依赖独立供者与零假设下对比的对称性。净逆转为两分量之和，不能替代双向终点。',
        '- GSE141864的21药家族及两个外部队列四分量IUT的21药家族分别列为次要分析。逐供者剔除、GSE237861的14次患者剔除重拟合、GSE141864的7次患者剔除均保留。等权每方向前50基因作为敏感性分析，不替代主要结果。','',
        '## 主要结果','',
        f'21药中，{summary["primary_bidirectional_means"]}药在GSE237861呈平均双向逆转，{summary["both_external_bidirectional"]}药在两个外部疾病背景均呈平均双向逆转；主要正式命中0，最小BH q={summary["primary_min_q"]:.6f}，机制升级命中0。','',
        '| 原探索药物 | 主要q（GSE237861） | 次要q（GSE141864） | 跨队列四分量q | 供者剔除双向：237/141 | 患者剔除双向：237/141 |','|---|---:|---:|---:|---|---|']
    for drug in ['DAS','PON','REG']:
        a=primary.loc[drug]; b=secondary.loc[drug]
        lines.append(f'| {a.drug_name} | {a.iut_BH_q:.6f} | {b.iut_BH_q:.6f} | {cross.loc[drug,"secondary_BH_q"]:.6f} | {fraction(dloo,drug,"GSE237861")} / {fraction(dloo,drug,"GSE141864")} | {fraction(ploo,drug,"GSE237861")} / {fraction(ploo,drug,"GSE141864")} |')
    lines+=['','原三药均为5名供者。Ponatinib与Regorafenib主要未校正P均为0.03125，恰为5名供者穷尽单侧检验可达到的最小P；经21药校正后均为q=0.328125。样本量小与P值离散限制检验能力，因此不能解释成“证明药物无效”，也不能据此降低门槛。','',
        'Cabozantinib也在两个外部背景的全部供者及患者剔除分析中保持双向方向；原三药不是唯一方向稳定的药物。Cabozantinib同样未通过主要门槛，不新增为已支持的候选治疗药物。','',
        '### 全部21药主要结果','',
        '| 药物 | 供者数 | UP抑制 | DOWN恢复 | 净逆转 | IUT P | BH q |','|---|---:|---:|---:|---:|---:|---:|']
    for r in primary.sort_values('drug_name').itertuples():
        lines.append(f'| {r.drug_name} | {r.donors} | {r.UP_suppression:.6f} | {r.DOWN_restoration:.6f} | {r.net_reversal:.6f} | {r.iut_p:.6f} | {r.iut_BH_q:.6f} |')
    lines+=['','数值单位为加权样本内秩变化，不能换算为心功能改善。完整三背景的63条结果及区间见[all_drug_results.tsv](../analysis/external_compound_v1/all_drug_results.tsv)。','',
        '## 疾病复现与机制解释的边界','',
        '完整GSE237861固定疾病特征组间差为+9.070853秩百分点，患者标签精确P=0.0250583，两队列家族BH q=0.0501166，严格未达到q<0.05；性别分层置换P=0.0321429。全部患者剔除方向为正。GSE141864对应差为+3.82537秩百分点，P和本轮完整两队列家族q均为0.2857143，患者剔除方向不稳定。历史仅完成一队列时的q保留在历史报告，当前完整家族结果以raw-heart-1.0为准。','',
        '与发现背景的全基因logFC相关分别为0.362和0.099；固定通路稳健检验没有命中。GSE237861的TLR4通路仅覆盖8基因、42.1%，不符合检验覆盖要求，列为未检验；用于校正的P=1占位不得解释成实际检验值。没有新的通路或靶点机制被确认。','',
        'GSE237861对照包括心衰、心肌梗死、心源性休克等危重症，不能称健康对照，也没有统一心超标准定义的脓毒症心肌病表型。全部14份归档FASTQ的质量字符为恒定编码Q30，不能代表原仪器碱基质量；实际读长76bp与论文75bp差异已保留。技术处理通过不能消除此来源限制。','',
        '健康iPSC心肌细胞的单剂量48小时扰动只能反映表达方向，缺乏脓毒症细胞救援、剂量反应、心肌功能与临床结局。已有心脏不良作用及疾病特异性限制仍然成立。自有大鼠LPS/对照数据继续作为探索背景，不赋予缺失的心超或药物干预终点。','',
        '## 审计与图表交付','',
        f'独立重建{numeric["reconstructed_cardio_donor_vectors"]}条药物–供者向量，替代实现复算秩、权重、63个IUT、BH家族及供者/患者剔除结果。最大数值差{numeric["vector_max_abs_error"]:.3g}，冻结来源哈希保持一致；详见[numeric_qa.json](../analysis/external_compound_v1/numeric_qa.json)。',
        '', '三幅图均有PDF、可编辑文本SVG以及600dpi PNG/TIFF，附来源表与英文图注。图2色标造成行错位的问题已修复并复查。数值审计、图像视觉检查、文件导出检查分开记录；未完成指定期刊版式认证或临床验证。',
        '', '- [图1：独立疾病复现](../analysis/external_compound_v1/figures/figure1_disease_replication.pdf)',
        '- [图2：全部21药跨疾病背景及双分量结果](../analysis/external_compound_v1/figures/figure2_all_drug_contexts.pdf)',
        '- [补充图：全基因与固定通路证据](../analysis/external_compound_v1/figures/figureS1_gene_and_pathway_evidence.pdf)',
        '- [导出及视觉检查记录](../analysis/external_compound_v1/figure_release_qa.json)','',
        '## 下一阶段的研究决定','',
        '**结束本轮计算扩展，停止将当前任一单药升级为保护机制主线。下一阶段应先评估论文定位与文献增量，再决定是否写作。** 推荐评估题目：“脓毒症心肌转录组药物重定位信号的跨队列稳健性与证据分层”。保留化合物筛选主题，将跨患者、跨疾病背景、供者稳定性、阴性结果和机制证据边界作为中心。该定位的发表价值尚需与既有研究比较，不能保证创新性或录用。','',
        '具体先形成主张–证据表与主图对应关系，核实最接近的已发表研究，明确现有数据能支持的新增知识；若缺乏足够增量，应保留为方法/探索报告。若未来坚持单药保护机制主线，需要新增独立疾病细胞药物扰动或功能救援证据，当前表达方向不能替代。不为获取阳性结果继续添加对接、挑选通路或修改阈值。','']
    report=aud/'独立疾病背景21药物复核与审计报告.md'
    report.write_text('\n'.join(lines),encoding='utf8')
    release_path=ROOT/'CURRENT_RELEASE.json'; release=readj(release_path)
    if release.get('external_compound_version')!=cfg['version']:
        snapshot=aud/('CURRENT_RELEASE_before_external_compound_'+now().replace(':','-')+'.json')
        snapshot.write_bytes(release_path.read_bytes())
    release.update(status='COMPLETE_EXTERNAL_COMPOUND_AUDITED_NO_FORMAL_HITS',updated_utc=now(),
        external_compound_version=cfg['version'],external_compound_summary=(out/'summary.json').relative_to(ROOT).as_posix(),
        external_compound_report=report.relative_to(ROOT).as_posix(),report=report.relative_to(ROOT).as_posix(),
        supported_candidates=0,external_compound_primary_min_BH_q=summary['primary_min_q'],
        continuation_decision='STOP_SINGLE_DRUG_MECHANISM_ESCALATION_ASSESS_CROSS_COHORT_ROBUSTNESS_MANUSCRIPT',
        external_drug_replication='NEW_DISEASE_CONTEXT_REUSED_CARDIAC_DONORS_NO_FORMAL_HITS')
    release['module_checks'].update({m:'PASS' for m in modules})
    write_json(release_path,release)
    return {'report':report.relative_to(ROOT).as_posix(),'status':'PASS','supported_candidates':0,'figures':3}

if __name__=='__main__':
    run_standard_module('56_external_compound_release',main)
    # Publish module56 status only after the shared wrapper has saved completion.
    release=readj(ROOT/'CURRENT_RELEASE.json')
    release['module_checks']['56_external_compound_release']='PASS'
    write_json(ROOT/'CURRENT_RELEASE.json',release)
    cfg=readj(ROOT/'config/external_compound_v1.json')
    files=list((ROOT/cfg['out']).rglob('*'))
    files += [ROOT/'config/external_compound_v1.json',ROOT/'docs/17_独立疾病背景下化合物复核方案.md',ROOT/'CURRENT_RELEASE.json',ROOT/'README.md',ROOT/release['report']]
    files += [p for p in (ROOT/'scripts').glob('*') if p.name[:2] in ['50','52','53','54','55','56']]
    files += [p for p in (ROOT/'outputs/audit').glob('*_status.json') if p.name[:2] in ['52','53','54','55','56']]
    manifest=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in sorted(set(files)) if p.is_file()]
    write_json(ROOT/'outputs/audit/external_compound_release_manifest.json',manifest)
