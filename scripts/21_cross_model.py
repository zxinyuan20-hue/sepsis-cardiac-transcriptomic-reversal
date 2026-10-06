"""Packages pandas/numpy/scipy; R edgeR/limma. Inputs frozen cohort/pathway contracts.
Outputs cohort effects, global pathway BH, concordance and target expression contexts.
Fixed config and seed; keep primary versus sensitivity families separate.
"""
import json,subprocess
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control,spearmanr,t
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,r_environment,run_standard_module
def main():
 out=path('analysis')/'mechanism';cfg=json.loads((ROOT/'config/mechanism_v1.json').read_text());freeze=json.loads((out/'mechanism_plan_freeze.json').read_text())
 if sha256(ROOT/'config/mechanism_v1.json')!=freeze['config_sha256']:raise ValueError('Mechanism config changed')
 for cohort in freeze['cohorts']:
  if sha256(ROOT/cohort['counts'])!=cohort['counts_sha256']:raise ValueError('Cohort counts changed')
 with (path('logs')/'21_R_cross_model_details.log').open('a',encoding='utf-8') as log:
  p=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/21_cross_model.R')],cwd=ROOT,env=r_environment(),stdout=log,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError('Cross-model R failed')
 cam=pd.read_csv(out/'camera_all_results.tsv',sep='\t');cam['global_primary_FDR']=np.nan;cam['sensitivity_FDR']=np.nan
 primary=cam.primary&cam.PValue.notna();cam.loc[primary,'global_primary_FDR']=false_discovery_control(cam.loc[primary,'PValue'].to_numpy())
 sens=~cam.primary&cam.PValue.notna();cam.loc[sens,'sensitivity_FDR']=false_discovery_control(cam.loc[sens,'PValue'].to_numpy())
 names=pd.read_csv(out/'fixed_pathways.tsv',sep='\t',dtype=str).drop_duplicates('pathway_id')[['pathway_id','pathway_name']]
 cam=cam.merge(names,on='pathway_id',validate='many_to_one');cam.to_csv(out/'pathway_evidence.tsv',sep='\t',index=False)
 alt=pd.read_csv(out/'camera_estimated_correlation_results.tsv',sep='\t');alt['global_primary_FDR']=np.nan;alt['sensitivity_FDR']=np.nan
 prim=alt.primary&alt.PValue.notna();sec=~alt.primary&alt.PValue.notna()
 alt.loc[prim,'global_primary_FDR']=false_discovery_control(alt.loc[prim,'PValue'].to_numpy())
 alt.loc[sec,'sensitivity_FDR']=false_discovery_control(alt.loc[sec,'PValue'].to_numpy())
 alt.to_csv(out/'pathway_correlation_sensitivity.tsv',sep='\t',index=False)
 human=pd.read_csv(path('analysis')/'human/primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
 annotation=pd.read_csv(ROOT/'outputs/preprocessing/human/probe_mapping_all.tsv',sep='\t',dtype=str).dropna().drop_duplicates(['SYMBOL','ENTREZID'])
 annotation=annotation[annotation.SYMBOL.isin(cfg['target_symbols'])][['SYMBOL','ENTREZID']].drop_duplicates();annotation=annotation[~annotation.SYMBOL.duplicated(keep=False)]
 concord=[];targets=[];checks=[]
 for name in ['human','GSE185754','GSE267388','LOCAL_RAT_LPS_2025','LOCAL_RAT_LPS_2025_without_C1_C2']:
  if name=='human':d=human.copy()
  else:
   full=pd.read_csv(out/f'{name}_all_gene_results.tsv',sep='\t',dtype={'human_entrez':str})
   if not np.allclose(false_discovery_control(full['P.Value'].to_numpy()),full['adj.P.Val'].to_numpy(),rtol=1e-9,atol=1e-14):raise ValueError('Gene BH mismatch')
   if not np.allclose(2*t.sf(abs(full.t),full.df_total),full['P.Value'],rtol=1e-9,atol=1e-14):raise ValueError('Gene moderated p mismatch')
   d=full.dropna(subset=['human_entrez']).set_index('human_entrez')
   if d.index.duplicated().any():raise ValueError('Nonunique mapped effects')
   shared=d.index.intersection(human.index);h=human.loc[shared];other=d.loc[shared];sig=h['adj.P.Val']<.05
   concord.append({'cohort':name,'shared_genes':len(shared),'logfc_spearman':float(spearmanr(h.logFC,other.logFC).statistic),'human_fdr_genes_mapped':int(sig.sum()),'same_direction_fraction_in_human_fdr':float((np.sign(h.loc[sig,'logFC'])==np.sign(other.loc[sig,'logFC'])).mean()),'role':'descriptive; genes not independent replicates'})
   checks.append({'cohort':name,'gene_BH':'PASS','moderated_P':'PASS'})
  for row in annotation.itertuples():
   entry={'cohort':name,'symbol':row.SYMBOL,'human_entrez':row.ENTREZID,'available':row.ENTREZID in d.index}
   if entry['available']:entry.update(logFC=float(d.loc[row.ENTREZID,'logFC']),FDR=float(d.loc[row.ENTREZID,'adj.P.Val']))
   targets.append(entry)
 pd.DataFrame(concord).to_csv(out/'cross_model_concordance.tsv',sep='\t',index=False);pd.DataFrame(targets).to_csv(out/'fixed_target_context.tsv',sep='\t',index=False)
 write_json(path('audit')/'cross_model_numeric_checks.json',checks)
 return {'primary_pathway_tests':int(primary.sum()),'significant_primary_pathway_tests':int((cam.global_primary_FDR<.05).sum()),'scope':'disease associations only; estimated-count cohorts provisional','concordance':concord}
if __name__=='__main__':run_standard_module('21_cross_model',main)
