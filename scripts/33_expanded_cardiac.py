"""Packages numpy/pandas/scipy/numba; R edgeR/limma. Frozen 21-drug family.
Outputs donor-aware endpoints/BH and harmonized common-landmark ES comparisons.
Paths/config from ROOT; exhaustive deterministic tests; logs through standard module.
"""
import itertools,json,subprocess,logging
import numpy as np,pandas as pd
from scipy.stats import rankdata,t,false_discovery_control,spearmanr
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,r_environment,run_standard_module
from screening_core import signature_es

def readj(p):return json.loads(p.read_text(encoding='utf-8'))
def exact(a):
 a=np.asarray(a);n=len(a);patterns=np.asarray(list(itertools.product([-1,1],repeat=n)));means=patterns@a/n;obs=a.mean();se=a.std(ddof=1)/np.sqrt(n);margin=t.ppf(.975,n-1)*se
 return {'donors':n,'mean_reversal':obs,'t_ci_low':obs-margin,'t_ci_high':obs+margin,'positive_donors':int((a>0).sum()),'exact_one_sided_p':float(np.mean(means>=obs-1e-14)),'minimum_possible_p':1/(2**n)}
def es(z,q):
 order=np.argsort(-z,kind='stable');u,d=signature_es(q,order,np.abs(z[order]));return u,d,(d-u)/2 if u*d<0 else 0.
def correlation(df,x,y):return float(spearmanr(df[x],df[y]).statistic) if df[x].nunique()>1 and df[y].nunique()>1 else None

def main():
 out=path('analysis')/'expansion';cfg=readj(ROOT/'config/expansion_v1.json')
 for row in readj(out/'expansion_plan_freeze.json')['files']:assert sha256(ROOT/row['path'])==row['sha256'],row['path']
 subprocess.run([RUNTIME['rscript'],str(ROOT/'scripts/33_expanded_cardiac.R')],cwd=ROOT,env=r_environment(),check=True)
 sample=pd.read_csv(out/'expanded_sample_rank_scores.tsv',sep='\t');cohort=pd.read_csv(out/'all_drug_feasibility.tsv',sep='\t');cohort=cohort[cohort.eligible]
 donors=[];blocks=[];endpoints=[];leave=[];checks=[]
 for (drug,size),sub in sample.groupby(['drug_state','size']):
  g=sub.groupby(['donor','block','condition'])[['up_rank','down_rank']].mean().unstack('condition')
  g['up_change']=g[('up_rank','DRUG')]-g[('up_rank','CTRL')];g['down_change']=g[('down_rank','DRUG')]-g[('down_rank','CTRL')]
  changes=pd.DataFrame({'up_change':g.up_change,'down_change':g.down_change});changes['reversal']=changes.down_change-changes.up_change
  for (donor,block),r in changes.iterrows():blocks.append({'drug_state':drug,'size':size,'donor':donor,'block':block,**r.to_dict()})
  agg=changes.groupby('donor').mean()
  for donor,r in agg.iterrows():donors.append({'drug_state':drug,'size':size,'donor':donor,**r.to_dict()})
  result=exact(agg.reversal);result.update(drug_state=drug,size=size,up_change=agg.up_change.mean(),down_change=agg.down_change.mean(),both_directions_opposing=bool(agg.up_change.mean()<0 and agg.down_change.mean()>0));endpoints.append(result)
  if size==100:
   for donor in agg.index:leave.append({'drug_state':drug,'excluded_donor':donor,'mean_reversal':agg.drop(index=donor).reversal.mean()})
 allend=pd.DataFrame(endpoints);allend['family_BH_q']=np.nan;allend['nonVEM_sensitivity_BH_q']=np.nan
 for size in [100,50,150]:
  mask=allend['size'].eq(size);assert mask.sum()==len(cohort)
  allend.loc[mask,'family_BH_q']=false_discovery_control(allend.loc[mask,'exact_one_sided_p'])
  m=mask&allend.drug_state.ne('VEM');allend.loc[m,'nonVEM_sensitivity_BH_q']=false_discovery_control(allend.loc[m,'exact_one_sided_p'])
 allend=allend.merge(cohort[['state','drug_name','doses','frozen_lincs_score','lincs_ranks']],left_on='drug_state',right_on='state',validate='many_to_one')
 for name,df in [('expanded_endpoint_results.tsv',allend),('expanded_donor_scores.tsv',pd.DataFrame(donors)),('expanded_culture_scores.tsv',pd.DataFrame(blocks)),('expanded_leave_one_donor_out.tsv',pd.DataFrame(leave))]:df.to_csv(out/name,sep='\t',index=False)
 # Verify sample percentile ranks independently from raw counts (normalization
 # and positive logCPM transform preserve within-sample ranks).
 x=pd.read_csv(out/'expanded_counts.tsv.gz',sep='\t',index_col=0);filters=pd.read_csv(out/'expanded_gene_filters.tsv',sep='\t',index_col=0)
 mapping=pd.read_csv(ROOT/'outputs/analysis/cardiac/symbol_entrez_one_to_one.tsv',sep='\t',dtype=str).set_index('SYMBOL').ENTREZID
 queries={di:set(pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_100_{di}.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id) for di in ['up','down']}
 for drug in sorted(cohort.state):
  sub=sample[(sample.drug_state==drug)&sample['size'].eq(100)];v=x.loc[filters[drug],sub.column_id]
  ranks=(rankdata(v.to_numpy(),axis=0)-.5)/len(v);ids=mapping.reindex(v.index).fillna('')
  for di in ['up','down']:assert np.allclose(ranks[ids.isin(queries[di])].mean(axis=0),sub[f'{di}_rank'],atol=1e-12)
  gt=pd.read_csv(out/f'{drug}_gene_results.tsv',sep='\t');delta=pd.read_csv(out/f'{drug}_donor_logFC.tsv.gz',sep='\t',index_col=0)
  assert np.allclose(false_discovery_control(gt['P.Value']),gt['adj.P.Val'],atol=1e-12)
  assert np.allclose(2*t.sf(abs(gt.t),gt.df_total),gt['P.Value'],atol=1e-12)
  assert np.allclose(gt.logFC,delta.mean(axis=1),atol=1e-10)
  a=pd.DataFrame(donors).query('drug_state==@drug and size==100').reversal.to_numpy();obs=a.mean();values=[]
  for k in range(len(a)+1):
   for subset in itertools.combinations(range(len(a)),k):values.append((a.sum()-2*a[list(subset)].sum())/len(a))
  assert np.isclose(np.mean(np.asarray(values)>=obs-1e-14),float(allend.loc[allend.drug_state.eq(drug)&allend['size'].eq(100),'exact_one_sided_p'].iloc[0]))
  checks.append({'drug_state':drug,'sample_rank':'PASS','gene_BH':'PASS','moderated_P':'PASS','donor_mean_logFC':'PASS','exact_P_subset_enumeration':'PASS'})
 # Mandatory regression: VEM same design/filter/scoring as previous frozen stage.
 previous=pd.read_csv(ROOT/'outputs/analysis/cardiac/donor_reversal_scores.tsv',sep='\t');previous=previous[previous.scenario=='primary']
 current=pd.DataFrame(donors);current=current[current.drug_state=='VEM'].merge(previous,left_on=['size','donor'],right_on=['signature_size','donor'],suffixes=('_new','_old'))
 assert len(current)==15 and np.allclose(current.reversal_new,current.reversal_old,atol=1e-12)
 # Common-background, same-ES descriptive comparison; never overwrites prior scores.
 common=pd.read_csv(out/'common_landmark_background.tsv',sep='\t',dtype=str);geneids=common.ENTREZID.to_numpy();q=np.zeros(len(common),dtype=np.int8)
 q[np.isin(geneids,list(queries['up']))]=1;q[np.isin(geneids,list(queries['down']))]=-1
 a=np.load(ROOT/'outputs/analysis/screening/screening_inputs.npz');meta=pd.read_csv(ROOT/'outputs/analysis/screening/signature_metadata.tsv.gz',sep='\t')
 match=pd.read_csv(out/'compound_overlap_identity.tsv',sep='\t');match=match[match.structure_supported&match.state.isin(cohort.state)]
 index=[list(a['genes']).index(g) for g in geneids];z=a['z'];lincs=[]
 for row in match.itertuples():
  for ix,m in meta[meta.pert_id==row.pert_id].iterrows():
   u,d,s=es(z[ix,index],q);lincs.append({'state':row.state,'pert_id':row.pert_id,'sig_id':m.sig_id,'cell_id':m.cell_id,'es_up':u,'es_down':d,'score':s})
 lincs=pd.DataFrame(lincs);lincs.to_csv(out/'harmonized_lincs_signature_scores.tsv',sep='\t',index=False)
 cell=lincs.groupby(['state','pert_id','cell_id']).score.median().reset_index();cell.to_csv(out/'harmonized_lincs_cell_scores.tsv',sep='\t',index=False)
 drugl=cell.groupby(['state','pert_id']).score.median().groupby('state').median()
 cardiac=[]
 for drug in sorted(cohort.state):
  delta=pd.read_csv(out/f'{drug}_donor_logFC.tsv.gz',sep='\t',index_col=0).loc[common.SYMBOL]
  for donor in delta.columns:
   u,d,s=es(delta[donor].to_numpy(),q);cardiac.append({'state':drug,'donor':donor,'es_up':u,'es_down':d,'score':s,'direction':'reversal' if u<0 and d>0 else ('mimic' if u>0 and d<0 else 'non_opposed')})
 cardiac=pd.DataFrame(cardiac);cardiac.to_csv(out/'harmonized_cardiac_donor_scores.tsv',sep='\t',index=False)
 comp=allend[allend['size']==100].copy().set_index('state');comp['harmonized_lincs_score']=drugl;comp['harmonized_cardiac_score']=cardiac.groupby('state').score.median();comp['harmonized_reversing_donor_fraction']=cardiac.groupby('state').direction.apply(lambda v:(v=='reversal').mean())
 comp=comp.reset_index();comp.to_csv(out/'cross_context_comparison.tsv',sep='\t',index=False)
 loo=[]
 for donor in sorted(cardiac.donor.unique()):
  alt=comp.copy();sc=cardiac[cardiac.donor!=donor].groupby('state').score.median();alt['harmonized_cardiac_score']=alt.state.map(sc)
  loo.append({'omission_type':'donor','omitted':donor,'spearman':correlation(alt,'harmonized_lincs_score','harmonized_cardiac_score')})
 for drug in sorted(comp.state):loo.append({'omission_type':'drug','omitted':drug,'spearman':correlation(comp[comp.state!=drug],'harmonized_lincs_score','harmonized_cardiac_score')})
 pd.DataFrame(loo).to_csv(out/'cross_context_leave_one_out.tsv',sep='\t',index=False)
 nz=(comp.harmonized_lincs_score!=0)&(comp.harmonized_cardiac_score!=0)
 report={'version':cfg['version'],'drugs':len(comp),'new_drugs_beyond_VEM':len(comp)-1,'BH_supported_primary_drugs':int((comp.family_BH_q<.05).sum()),'nominal_primary_p_lt_05':int((comp.exact_one_sided_p<.05).sum()),'positive_mean_drugs':int((comp.mean_reversal>0).sum()),'positive_and_both_directions':int(((comp.mean_reversal>0)&comp.both_directions_opposing).sum()),'comparison':{'original_LINCS_vs_cardiac_rank_spearman':correlation(comp,'frozen_lincs_score','mean_reversal'),'harmonized_ES_spearman':correlation(comp,'harmonized_lincs_score','harmonized_cardiac_score'),'harmonized_ES_spearman_without_VEM':correlation(comp[comp.state!='VEM'],'harmonized_lincs_score','harmonized_cardiac_score'),'both_nonzero_drugs':int(nz.sum()),'same_sign_fraction_among_both_nonzero':float((np.sign(comp.loc[nz,'harmonized_lincs_score'])==np.sign(comp.loc[nz,'harmonized_cardiac_score'])).mean()) if nz.any() else None,'no_independence_based_P':True},'four_donor_drugs':int((comp.donors==4).sum()),'common_landmarks':len(common),'shared_control_warning':True,'clinical_efficacy':'NOT_ESTABLISHED'}
 write_json(out/'expanded_results_summary.json',report)
 write_json(path('audit')/'expanded_numeric_checks.json',{'checks':checks,'VEM_frozen_regression':'PASS','frozen_inputs':'PASS'})
 logging.info('Expanded result: %s',report)
 return report

if __name__=='__main__':run_standard_module('33_expanded_cardiac',main)
