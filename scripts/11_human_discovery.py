"""Packages: pandas/numpy/scipy, R limma. Config input; outputs/analysis/human.
Run prespecified human contrasts and two sensitivities, freeze landmark signatures and hashes.
Seed from config. No LINCS scoring or held-out model disease tests.
"""
import json,subprocess
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from pipeline_utils import ROOT,CONFIG,RUNTIME,path,sha256,write_json,now,r_environment,run_standard_module

def main():
    if not CONFIG['analysis']['formal_analysis_enabled']:raise RuntimeError('Analysis gate disabled')
    qc=json.loads((path('audit')/'preprocessing_qc.json').read_text())
    if set(qc['human']['review_flags'])!=set(CONFIG['analysis']['qc_flagged_samples']):raise ValueError('Unreviewed QC flags')
    with (path('logs')/'11_R_discovery_details.log').open('a',encoding='utf-8') as log:
        p=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/11_human_discovery.R')],cwd=ROOT,env=r_environment(),stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError('Discovery R failed')
    out=path('analysis')/'human'
    primary=pd.read_csv(out/'primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
    if primary.index.duplicated().any() or not np.isfinite(primary[['logFC','t','P.Value','adj.P.Val']].values).all():raise ValueError('Invalid discovery output')
    lm=pd.read_csv(path('preparation')/'lincs_gene_dictionary.tsv',sep='\t',dtype={'pr_gene_id':str})
    universe=set(lm.loc[lm.pr_is_lm==1,'pr_gene_id'])
    ranked=primary.loc[primary.index.intersection(universe)].copy()
    frozen=out/'frozen_signatures';frozen.mkdir(exist_ok=True)
    rows=[]
    for size in CONFIG['analysis']['signature_sizes_per_direction']:
        for direction in ['up','down']:
            subset=ranked[(ranked['adj.P.Val']<CONFIG['analysis']['fdr']) & ((ranked.logFC>0) if direction=='up' else (ranked.logFC<0))].copy()
            subset['abs_t']=subset.t.abs();subset=subset.reset_index().sort_values(['abs_t','entrez_id'],ascending=[False,True]).head(size)
            target=frozen/f'landmarks_{size}_{direction}.tsv';subset.to_csv(target,sep='\t',index=False)
            rows.append({'requested_per_direction':size,'direction':direction,'actual_genes':len(subset),'coverage_pass':len(subset)>=CONFIG['analysis']['minimum_mapped_genes_per_direction'],'path':target.relative_to(ROOT).as_posix(),'sha256':sha256(target)})
    pd.DataFrame(rows).to_csv(frozen/'signature_coverage.tsv',sep='\t',index=False)
    robustness=[]
    for label in ['without_qc_flag','submitted']:
        other=pd.read_csv(out/f'{label}_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id').loc[primary.index]
        significant=primary['adj.P.Val']<CONFIG['analysis']['fdr']
        robustness.append({'sensitivity':label,'all_gene_logfc_spearman':float(spearmanr(primary.logFC,other.logFC).statistic),'primary_fdr_genes':int(significant.sum()),'direction_agreement_primary_fdr_genes':float((np.sign(primary.loc[significant,'logFC'])==np.sign(other.loc[significant,'logFC'])).mean()),'same_direction_and_fdr_fraction':float(((np.sign(primary.loc[significant,'logFC'])==np.sign(other.loc[significant,'logFC']))&(other.loc[significant,'adj.P.Val']<CONFIG['analysis']['fdr'])).mean())})
    write_json(out/'sensitivity_concordance.json',robustness)
    manifest={'frozen_utc':now(),'protocol_version':CONFIG['protocol_version'],'query_background':'948 mapped measured landmarks; significance controlled across all 18866 annotated genes per contrast','primary_size':CONFIG['analysis']['primary_signature_size_per_direction'],'selection':'BH FDR < 0.05, separate sign, decreasing abs moderated t, Entrez lexical tie break','human_matrix_sha256':sha256(ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz'),'study_config_sha256':sha256(ROOT/'config/study_config.json'),'primary_results_sha256':sha256(out/'primary_sepsis_vs_nonfailing.tsv'),'files':rows,'compound_ranking_performed':False,'held_out_disease_tests_performed':False}
    write_json(frozen/'signature_freeze_manifest.json',manifest)
    return {'signatures':rows,'sensitivity':robustness,'compound_ranking_performed':False}
if __name__=='__main__':run_standard_module('11_human_discovery',main)
