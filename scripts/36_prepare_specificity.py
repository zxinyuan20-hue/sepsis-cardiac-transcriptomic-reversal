"""Packages pandas. Pre-specificity contrast/universe freeze; no new drug scores.
Input frozen RMA, expanded counts/design, mapping and configs. Output common876 map,
plan manifest and input audit. Paths from ROOT/config, deterministic pipeline seed/log.
"""
import json
import pandas as pd
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module

def main():
 out=path('analysis')/'specificity';out.mkdir(exist_ok=True)
 if (out/'specificity_results.tsv').exists():raise RuntimeError('Specificity effects already exist; do not rewrite pre-score freeze')
 ex=path('analysis')/'expansion';freeze=json.loads((ex/'expansion_plan_freeze.json').read_text(encoding='utf-8'))
 for f in freeze['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
 h=ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz';ids=pd.read_csv(h,sep='\t',usecols=[0],dtype=str).iloc[:,0]
 m=pd.read_csv(ex/'common_landmark_background.tsv',sep='\t',dtype=str);m=m[m.ENTREZID.isin(ids)].sort_values('ENTREZID')
 assert len(m)==876 and m.SYMBOL.is_unique and m.ENTREZID.is_unique
 m.to_csv(out/'specificity_common_genes.tsv',sep='\t',index=False)
 samples=ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv';s=pd.read_csv(samples,sep='\t')
 assert s.group.value_counts().to_dict()=={'sepsis':20,'IHD':11,'nonfailing':11,'DCM':9}
 humanfreeze=json.loads((ROOT/'outputs/analysis/human/frozen_signatures/signature_freeze_manifest.json').read_text(encoding='utf-8'));assert sha256(h)==humanfreeze['human_matrix_sha256']
 inputs=[ROOT/'config/specificity_v1.json',ROOT/'docs/13_疾病特异性分析方案_比较前.md',h,samples,out/'specificity_common_genes.tsv',ex/'expanded_counts.tsv.gz',ex/'expanded_analysis_design.tsv',ex/'cross_context_comparison.tsv',ex/'expanded_donor_scores.tsv',ROOT/'config/study_config.json']
 for name in ['nonfailing','IHD','DCM']:inputs.append(ROOT/f'outputs/analysis/human/primary_sepsis_vs_{name}.tsv')
 write_json(out/'specificity_plan_freeze.json',{'frozen_utc':now(),'prior_all21_drug_results_known':True,'new_specificity_endpoints_not_computed':True,'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in inputs]})
 result={'human_samples':s.group.value_counts().to_dict(),'common_genes':len(m),'drug_family':21,'known_human_QC_flag':['GSM2109172'],'scope':'same-cohort disease specificity, not external validation'}
 write_json(out/'specificity_preparation_summary.json',result);return result

if __name__=='__main__':run_standard_module('36_prepare_specificity',main)
