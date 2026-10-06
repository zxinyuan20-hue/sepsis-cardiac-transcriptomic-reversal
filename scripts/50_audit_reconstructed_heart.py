"""Independent numerical audit; pandas,numpy,scipy. Seed fixed in shared wrapper.
Inputs reconstructed counts+scores+prior release. Outputs audit/report/manifests.
No automatic efficacy upgrade. Holds are reported separately from analysis success.
"""
import itertools,json
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
from pipeline_utils import ROOT,sha256,write_json,run_standard_module,now
OUT=ROOT/'outputs/analysis/raw_heart_v1';AUD=ROOT/'outputs/audit'
def main():
    j=lambda p:json.loads(p.read_text(encoding='utf8'))
    if j(ROOT/'CURRENT_RELEASE.json').get('external_compound_version'):
        return {'status':'HISTORICAL_RELEASE_PRESERVED','reason':'Newer external-compound release is current; historical results and current release preserved'}
    quant=j(OUT/'quantification_audit.json');raw=j(OUT/'raw_acquisition_audit.json')
    assert raw['verified']==14 and not raw['failures'] and quant['samples']==14
    counts=pd.read_csv(OUT/'gene_counts_14.tsv.gz',sep='\t',index_col=0)
    qc=pd.read_csv(OUT/'technical_qc.tsv',sep='\t');assert len(qc)==14 and qc.title.is_unique
    for r in qc.itertuples():
        assert counts[r.title].sum()==r.assigned_s1
        c=j(OUT/'checksums'/f'{r.run}.json');assert sha256(OUT/'sources/fastq'/f'{r.run}.fastq.gz')==c['sha256']
        a=j(OUT/'raw_qc'/f'{r.run}.json');assert a['reads']==r.raw_reads and a['gzip_crc_and_records_pass']
    oldmanifest=j(AUD/'external_human_release_manifest.json');preserved=0
    for f in oldmanifest:
        if f['path'].startswith('outputs/analysis/external_human/'):
            assert sha256(ROOT/f['path'])==f['sha256'];preserved+=1
    summary=None;numeric={}
    if quant['disease_inference_ready']:
        summary=j(OUT/'disease_replication_summary.json')
        s=pd.read_csv(OUT/'patient_signature_scores.tsv',sep='\t');cases=s.group.eq('sepsis').to_numpy();v=s.primary_disease_rank_score.to_numpy()
        x=pd.read_csv(OUT/'entrez_counts_14.tsv.gz',sep='\t',index_col=0);x.index=x.index.astype(str)
        retained=pd.read_csv(OUT/'expression_filter.tsv',sep='\t',dtype={'entrez_id':str});bg=set(retained.loc[retained.retained.astype(str).str.lower().eq('true'),'entrez_id'])
        h=pd.read_csv(ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str});bg=sorted(bg&set(h.entrez_id))
        a=x.loc[bg,s.title].to_numpy();rank=np.apply_along_axis(rankdata,0,a);rank=(rank-.5)/len(bg)
        sets={d:set(pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_100_{d}.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id) for d in ['up','down']}
        rebuilt=rank[[g in sets['up'] for g in bg],:].mean(axis=0)-rank[[g in sets['down'] for g in bg],:].mean(axis=0)
        assert np.allclose(rebuilt,v,atol=1e-12,rtol=0)
        observed=v[cases].sum()/7-v[~cases].sum()/7;null=[]
        # Complement control assignments and total-sum formula, unlike module49.
        for control in itertools.combinations(range(14),7):
            z=v[list(control)].sum();null.append((v.sum()-2*z)/7)
        p=np.mean(np.asarray(null)>=observed-1e-14)
        assert abs(p-summary['primary_one_sided_exact_p'])<1e-12
        assert abs(observed-summary['primary_effect'])<1e-12
        female=s.sex.eq('F').to_numpy();vx=v.copy()
        for g in [female,~female]:vx[g]-=vx[g].mean()
        case_res=cases.astype(float)
        for g in [female,~female]:case_res[g]-=case_res[g].mean()
        adjusted=case_res@vx/(case_res@case_res);sn=[]
        for cm in itertools.combinations(np.where(~female)[0],4):
            for cf in itertools.combinations(np.where(female)[0],3):
                use=np.ones(14);use[list(cm)+list(cf)]=0
                for g in [female,~female]:use[g]-=use[g].mean()
                sn.append(use@vx/(use@use))
        sp=np.mean(np.asarray(sn)>=adjusted-1e-14)
        assert len(sn)==840 and abs(sp-summary['sex_stratified_sensitivity_p'])<1e-12
        family=pd.read_csv(OUT/'two_cohort_primary_family.tsv',sep='\t')
        assert np.allclose(false_discovery_control(family.p),family.BH_q,atol=1e-12)
        pathways=pd.read_csv(OUT/'two_cohort_pathway_family.tsv',sep='\t')
        for _,t in pathways.groupby('mode'):assert len(t)==10 and np.allclose(false_discovery_control(t.PValue.fillna(1)),t.planned_10_test_BH_q,atol=1e-12)
        numeric={'rank_from_raw_counts':'PASS','complement_3432_permutations':'PASS','residualized_840_sex_permutations':'PASS','BH_scipy':'PASS'}
    status='COMPLETE_COUNTS_QC_HOLD' if summary is None else 'COMPLETE_DISEASE_RECONSTRUCTION_AUDITED_NO_DRUG_EFFICACY_CLAIM'
    audit={'status':status,'samples':14,'count_sums_vs_featureCounts':'PASS','fastq_sha256':'PASS','archive_record_counts':'PASS','prior_external_analysis_files_unchanged':preserved,'numeric':numeric,'qc_hold_samples':quant['qc_hold_samples'],'supported_therapeutic_candidates':0,'verified_utc':now()}
    write_json(AUD/'raw_heart_final_qa.json',audit)
    lines=['# 完整人类心肌重建与疾病复现审计','','2026-10-03｜raw-heart-1.0。', '',f'本轮状态：**{status}**。完整计数包含脓毒症7例和非感染性危重死亡对照7例。研究主题仍是为化合物筛选检验疾病特征的可靠性，没有新增药物有效性证据。', '', '## 数据与方法', '', '全部14份ENA FASTQ已完成官方MD5、文件大小、本地SHA256、gzip CRC与记录数检查。GENCODE v47主组装、Subjunc 2.1.1和featureCounts 2.1.1统一处理全部样本。方法与原作者HISAT2/Trimmomatic有差异；本轮是统一原始读段重分析，不是逐参数重现。', '', '病例6男1女、对照4男3女；信息来自论文Table S1。对照含心衰、心源性休克、心梗等，不能称健康对照，也不能据此定义超声证实的脓毒症心肌病。部分病例已有心肌病，终末期与既往疾病混杂仍存在。', '', f'共有{int(qc.constant_archive_quality.sum())}/14份FASTQ的编码质量值在全文件恒定。不能用其Q30比例证明原始测序质量；质量修剪不提供原始仪器低质量筛除证据。实际读长76 bp，论文记75 bp，ENA碱基总数与76 bp相符。上述差异已在计数前记录。', '', '## 技术质量', '', '| 样本 | 读段保留比例 | 正链基因归属比例 | 正链方向比例 | 审查项 |','|---|---:|---:|---:|---|']
    for r in qc.itertuples():lines.append(f'| {r.title} | {r.retained_read_fraction:.1%} | {r.assigned_forward_fraction:.1%} | {r.forward_share_of_forward_plus_reverse:.1%} | {r.qc_hold if pd.notna(r.qc_hold) else "无预定阈值触发"} |')
    if summary is None:
        lines+=['','## 疾病推断尚未放行','','存在预定技术审查项，未运行正式疾病差异检验。所有样本与计数均保留，不自动删掉不利样本。需要明确技术原因并记录修订，再决定是否可以推断。原药物证据维持原结论。']
    else:
        lines+=['','## 固定疾病特征的结果','','| 指标 | 结果 |','|---|---:|',f'| 原100UP/90DOWN疾病评分病例减对照 | {summary["rank_percentage_points"]:.4f}秩百分点 |',f'| 全3432次排列单侧P | {summary["primary_one_sided_exact_p"]:.6f} |',f'| 两队列家族BH q | {summary["two_cohort_BH_q"]:.6f} |',f'| 性别分层840次敏感性P | {summary["sex_stratified_sensitivity_p"]:.6f} |',f'| 全部逐患者去除后方向为正 | {summary["all_loo_effects_positive"]} |',f'| 共有基因logFC描述性Spearman | {summary["gene_logFC_spearman_descriptive"]:.4f} |',f'| 性别调整edgeR全基因FDR<0.05数 | {summary["genes_FDR05"]} |','','主P和敏感性P分开保留，未择小报告；原GSE141864结果未变，两队列校正更新仅在新版本保存。秩百分点不是表达量百分比或心功能改善。基因相关仅作描述，不将基因当独立患者。', '', '完整结果、固定通路和两队列家族见本目录结果文件。疾病表达方向即使一致，也不能自动证明药物保护作用；若主检验不支持或患者方向不稳定，应停止升级当前三药保护机制解释，保留疾病背景差异的评估。']
    lines+=['','## 审计范围','','重新从原始计数构造秩分数、互补组合计算主排列、残差化计算性别分层排列，并以独立BH实现核对。技术审查未通过时不进行这些疾病数值审计。历史外部分析文件逐一核对哈希保持不变。没有新增药物筛选、分子对接、湿实验、临床疗效验证或正式投稿。','','方案见docs/16_完整心肌原始数据重建方案.md；完整数据及日志见outputs/analysis/raw_heart_v1。']
    report=AUD/'完整人类心肌重建与疾病复现审计.md';report.write_text('\n'.join(lines)+'\n',encoding='utf8')
    cp=ROOT/'CURRENT_RELEASE.json';snap=AUD/'pre_raw_heart_release_snapshot.json'
    if not snap.exists():snap.write_bytes(cp.read_bytes())
    current=j(cp);current.update(status=status,updated_utc=now(),raw_heart_version='raw-heart-1.0',raw_heart_report=report.relative_to(ROOT).as_posix(),report=report.relative_to(ROOT).as_posix(),supported_candidates=0,continuation_decision='REVIEW_FULL_DISEASE_REPLICATION_BEFORE_ANY_DRUG_MECHANISM_ESCALATION')
    write_json(cp,current)
    paths=[p for p in OUT.rglob('*') if p.is_file() and not p.name.endswith('.part')]+[report,AUD/'raw_heart_final_qa.json',cp]
    write_json(AUD/'raw_heart_release_manifest.json',[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size} for p in paths])
    return audit
if __name__=='__main__':run_standard_module('50_audit_reconstructed_heart',main)
