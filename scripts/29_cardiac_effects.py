"""Packages numpy/pandas/scipy; R edgeR/limma. Inputs frozen cardiac design/config.
Output donor-level exact sign-flip test, all sensitivities, gene/pathway numeric QA.
Paths use ROOT/config; deterministic exhaustive 32 patterns (fixed seed for pipeline).
"""
import itertools,json,subprocess
import numpy as np,pandas as pd
from scipy.stats import rankdata,t,false_discovery_control,spearmanr
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,r_environment,run_standard_module

def test_scores(scores):
 a=np.asarray(scores,dtype=float);n=len(a);assert n==5
 null=np.asarray(list(itertools.product([-1,1],repeat=n)))@a/n
 obs=float(a.mean());se=float(a.std(ddof=1)/np.sqrt(n));crit=t.ppf(.975,n-1)
 return {'donors':n,'mean_reversal':obs,'sd_across_donors':float(a.std(ddof=1)),'t_ci_low':float(obs-crit*se),'t_ci_high':float(obs+crit*se),'positive_donors':int((a>0).sum()),'one_sided_exact_signflip_p':float(np.mean(null>=obs-1e-14)),'two_sided_exact_signflip_p':float(np.mean(np.abs(null)>=abs(obs)-1e-14)),'sign_patterns':len(null)},null

def main():
 out=path('analysis')/'cardiac';cfg=json.loads((ROOT/'config/cardiac_v1.json').read_text(encoding='utf-8'))
 freeze=json.loads((out/'cardiac_plan_freeze.json').read_text(encoding='utf-8'))
 for row in freeze['files']:assert sha256(ROOT/row['path'])==row['sha256'],row['path']
 subprocess.run([RUNTIME['rscript'],str(ROOT/'scripts/29_cardiac_effects.R')],cwd=ROOT,env=r_environment(),check=True)
 d=pd.read_csv(out/'analysis_design.tsv',sep='\t',dtype={'Experiment':str})
 mapping=pd.read_csv(out/'symbol_entrez_one_to_one.tsv',sep='\t',dtype=str)
 mp=dict(zip(mapping.SYMBOL,mapping.ENTREZID))
 scores=[];blocks=[];endpoints=[];coverage=[];sample_rows=[];checks=[];lodo=[];loco=[];nulls=[]
 for scenario in ['primary','all_same_donors']:
  e=pd.read_csv(out/f'{scenario}_sample_logCPM.tsv.gz',sep='\t',index_col=0)
  assert np.isfinite(e.to_numpy()).all() and e.columns.is_unique and e.index.is_unique
  ranks=(rankdata(e.to_numpy(),axis=0,method='average')-.5)/len(e)
  ids=np.array([mp.get(x,'') for x in e.index]);s=d.set_index('column_id').loc[e.columns].rename_axis('column_id').reset_index()
  for size in [cfg['primary_signature_size']]+cfg['sensitivity_signature_sizes']:
   sets={};eligible=True
   for direction in ['up','down']:
    query=pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_{size}_{direction}.tsv',sep='\t',dtype={'entrez_id':str})
    mask=np.isin(ids,query.entrez_id);sets[direction]=mask
    covered=int(mask.sum());ratio=covered/len(query)
    coverage.append({'scenario':scenario,'signature_size':size,'direction':direction,'source_genes':len(query),'covered_genes':covered,'coverage':ratio})
    eligible &= covered>=cfg['minimum_signature_genes_per_direction'] and ratio>=cfg['minimum_signature_coverage_per_direction']
   if not eligible:raise ValueError(f'Insufficient frozen query coverage {scenario} {size}')
   s['up_rank']=ranks[sets['up']].mean(axis=0);s['down_rank']=ranks[sets['down']].mean(axis=0)
   s['disease_rank']=s.up_rank-s.down_rank
   for r in s.itertuples():sample_rows.append({'scenario':scenario,'signature_size':size,'column_id':r.column_id,'donor':r.Cell,'condition':r.State,'block':r.block,'up_rank':r.up_rank,'down_rank':r.down_rank,'disease_rank':r.disease_rank})
   culture=s.groupby(['Cell','block','State'])[['up_rank','down_rank','disease_rank']].mean().reset_index()
   donor_rows=[]
   if scenario=='primary':
    for (donor,block),g in culture.groupby(['Cell','block']):
     g=g.set_index('State')[['up_rank','down_rank','disease_rank']];v=g.loc['VEM']-g.loc['CTRL']
     blocks.append({'signature_size':size,'donor':donor,'block':block,'up_change':float(v.up_rank),'down_change':float(v.down_rank),'reversal':float(-v.disease_rank)})
    b=pd.DataFrame([x for x in blocks if x['signature_size']==size]);agg=b.groupby('donor')[['up_change','down_change','reversal']].mean()
   else:
    means=culture.groupby(['Cell','State'])[['up_rank','down_rank','disease_rank']].mean()
    for donor in sorted(s.Cell.unique()):
     v=means.loc[(donor,'VEM')]-means.loc[(donor,'CTRL')]
     donor_rows.append({'donor':donor,'up_change':float(v.up_rank),'down_change':float(v.down_rank),'reversal':float(-v.disease_rank)})
    agg=pd.DataFrame(donor_rows).set_index('donor')
   for donor,r in agg.iterrows():scores.append({'scenario':scenario,'signature_size':size,'donor':donor,**r.to_dict()})
   result,null=test_scores(agg.reversal);result.update(scenario=scenario,signature_size=size,mean_up_change=float(agg.up_change.mean()),mean_down_change=float(agg.down_change.mean()),both_directions_opposing=bool(agg.up_change.mean()<0 and agg.down_change.mean()>0),role='primary planned exploratory endpoint' if scenario=='primary' and size==100 else 'sensitivity; do not select best P')
   endpoints.append(result)
   for i,value in enumerate(null):nulls.append({'scenario':scenario,'signature_size':size,'pattern_index':i,'null_mean':float(value)})
   if scenario=='primary' and size==100:
    for donor in agg.index:lodo.append({'excluded_donor':donor,'remaining_donors':4,'mean_reversal':float(agg.drop(index=donor).reversal.mean())})
    for row in b.itertuples():
     reduced=b[b.block!=row.block].groupby('donor').reversal.mean()
     assert len(reduced)==5
     loco.append({'excluded_block':row.block,'donor':row.donor,'mean_reversal':float(reduced.mean())})
  genes=pd.read_csv(out/f'{scenario}_gene_results.tsv',sep='\t')
  bh=false_discovery_control(genes['P.Value']);assert np.allclose(bh,genes['adj.P.Val'],rtol=1e-9,atol=1e-12)
  pv=2*t.sf(abs(genes['t']),genes.df_total);assert np.allclose(pv,genes['P.Value'],rtol=1e-8,atol=1e-12)
  delta=pd.read_csv(out/f'{scenario}_donor_logFC.tsv.gz',sep='\t',index_col=0)
  assert np.allclose(delta.mean(axis=1),genes.logFC,atol=1e-10)
  assert (genes.df_residual==4).all()
  # Independent aggregation directly from sample expression, rather than model fit.
  for donor in sorted(s.Cell.unique()):
   means={}
   for cond in ['CTRL','VEM']:
    sub=s[(s.Cell==donor)&(s.State==cond)]
    means[cond]=np.mean([e[g.column_id].mean(axis=1).to_numpy() for _,g in sub.groupby('block')],axis=0)
   assert np.allclose(means['VEM']-means['CTRL'],delta[donor],atol=1e-10)
  checks.append({'scenario':scenario,'gene_BH':'PASS','moderated_P':'PASS','gene_logFC_equal_donor_mean':'PASS','sample_to_donor_aggregation':'PASS','original_residual_df':4})
 for filename,rows in [('donor_reversal_scores.tsv',scores),('culture_reversal_scores.tsv',blocks),('endpoint_results.tsv',endpoints),('query_coverage.tsv',coverage),('sample_rank_scores.tsv',sample_rows),('leave_one_donor_out.tsv',lodo),('leave_one_culture_out.tsv',loco),('exact_signflip_distribution.tsv',nulls)]:
  pd.DataFrame(rows).to_csv(out/filename,sep='\t',index=False)
 cam=pd.read_csv(out/'primary_pathways.tsv',sep='\t');valid=cam.PValue.notna()
 assert np.allclose(false_discovery_control(cam.loc[valid,'PValue']),cam.loc[valid,'FDR'],atol=1e-12)
 concordance=[];genes=pd.read_csv(out/'primary_gene_results.tsv',sep='\t',dtype={'entrez_id':str})
 for contrast in ['nonfailing','IHD','DCM']:
  human=pd.read_csv(ROOT/f'outputs/analysis/human/primary_sepsis_vs_{contrast}.tsv',sep='\t',dtype={'entrez_id':str})
  merged=genes.dropna(subset=['entrez_id']).merge(human,on='entrez_id',suffixes=('_drug','_human'),validate='one_to_one')
  assert len(merged)>1000
  selected=merged['adj.P.Val_human']<.05
  concordance.append({'human_contrast':f'sepsis_vs_{contrast}','shared_genes':len(merged),'spearman_logFC':float(spearmanr(merged.logFC_drug,merged.logFC_human).statistic),'human_significant_mapped_genes':int(selected.sum()),'opposite_fraction_in_human_significant':float((merged.loc[selected,'logFC_drug']*merged.loc[selected,'logFC_human']<0).mean()),'role':'descriptive; correlated genes not independent replicates'})
  if contrast=='nonfailing':merged.to_csv(out/'human_drug_gene_comparison.tsv',sep='\t',index=False)
 pd.DataFrame(concordance).to_csv(out/'human_drug_concordance.tsv',sep='\t',index=False)
 write_json(path('audit')/'cardiac_numeric_checks.json',{'gene_checks':checks,'pathway_BH':'PASS','frozen_inputs':'PASS'})
 primary=next(x for x in endpoints if x['scenario']=='primary' and x['signature_size']==100)
 result={'version':cfg['version'],'primary':primary,'all_endpoint_results':endpoints,'original_screen_supported_candidates':0,'scope':'healthy cardiac culture expression response; no disease rescue or safety proof'}
 write_json(out/'cardiac_effect_summary.json',result)
 return result

if __name__=='__main__':run_standard_module('29_cardiac_effects',main)
