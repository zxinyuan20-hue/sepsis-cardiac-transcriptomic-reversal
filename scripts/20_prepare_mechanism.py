"""Packages pandas; R AnnotationDbi. Inputs config/frozen ranked list and mapped raw gene IDs.
Outputs outputs/analysis/mechanism: fixed five pathways, crosswalks and analysis manifests.
Fixed seed inherited. Does not inspect held-out disease effects.
"""
import json,subprocess
import pandas as pd
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,now,r_environment,run_standard_module
def main():
    cfg=json.loads((ROOT/'config/mechanism_v1.json').read_text());out=path('analysis')/'mechanism';out.mkdir(exist_ok=True)
    freeze=json.loads((path('analysis')/'screening/screening_result_freeze.json').read_text())
    if sha256(path('analysis')/'screening/exploratory_top10.tsv')!=freeze['exploratory_shortlist_sha256']:raise ValueError('Exploratory list changed')
    with (path('logs')/'20_R_mapping_details.log').open('a',encoding='utf-8') as log:
        p=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/20_prepare_mechanism.R')],cwd=ROOT,env=r_environment(),stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError('Mouse mapping failed')
    rx=pd.read_csv(path('preparation')/'reactome_human_direct_membership.tsv.gz',sep='\t',dtype=str);rx=rx[rx.pathway_id.isin(cfg['pathways'])]
    if set(rx.pathway_id)!=set(cfg['pathways']):raise ValueError('Missing predefined pathways')
    rx.to_csv(out/'fixed_pathways.tsv',sep='\t',index=False)
    mo=pd.read_csv(path('preparation')/'human_mouse_ortholog_one_to_one.tsv',sep='\t',dtype=str)
    en=pd.read_csv(out/'mouse_ensembl_entrez_one_to_one.tsv',sep='\t',dtype=str)
    one=en.merge(mo,left_on='ENTREZID',right_on='mouse_entrez',validate='one_to_one').rename(columns={'ENSEMBL':'feature_id'})
    one[['feature_id','human_entrez']].to_csv(out/'GSE185754_mapping.tsv',sep='\t',index=False)
    mo.rename(columns={'mouse_entrez':'feature_id'})[['feature_id','human_entrez']].to_csv(out/'GSE267388_mapping.tsv',sep='\t',index=False)
    rat=pd.read_csv(ROOT/'outputs/preprocessing/rat/rat_to_human_mapping.tsv',sep='\t',dtype=str)
    rat.rename(columns={'rat_ensembl':'feature_id'})[['feature_id','human_entrez']].to_csv(out/'LOCAL_RAT_LPS_2025_mapping.tsv',sep='\t',index=False)
    manifest=[]
    for dataset in ['GSE185754','GSE267388','LOCAL_RAT_LPS_2025']:
        if dataset.startswith('GSE'):
            counts=path('preparation')/f'{dataset}_counts.tsv.gz';s=pd.read_csv(path('preparation')/f'{dataset}_sample_map.tsv',sep='\t').rename(columns={'gsm':'sample_id'})
        else:counts=ROOT/'outputs/preprocessing/rat/rat_counts.tsv.gz';s=pd.read_csv(ROOT/'outputs/preprocessing/rat/sample_manifest.tsv',sep='\t')
        s=s[['sample_id','group']];s['group']=s.group.str.lower();s.to_csv(out/f'{dataset}_samples.tsv',sep='\t',index=False)
        manifest.append({'dataset':dataset,'counts':counts.relative_to(ROOT).as_posix(),'counts_sha256':sha256(counts),'samples':f'outputs/analysis/mechanism/{dataset}_samples.tsv','mapping':f'outputs/analysis/mechanism/{dataset}_mapping.tsv'})
    write_json(out/'cohort_inputs.json',manifest)
    write_json(out/'mechanism_plan_freeze.json',{'utc':now(),'config_sha256':sha256(ROOT/'config/mechanism_v1.json'),'pathways_sha256':sha256(out/'fixed_pathways.tsv'),'cohorts':manifest,'screening_freeze':freeze['complete_results_sha256'],'stage':'before held-out effects'})
    return {'mouse_ensembl_human_pairs':len(one),'pathways':rx.groupby('pathway_id').size().to_dict()}
if __name__=='__main__':run_standard_module('20_prepare_mechanism',main)
