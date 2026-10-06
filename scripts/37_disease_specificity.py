"""Packages numpy/pandas/scipy. Inputs frozen human RMA and 21-drug cardiac design.
Outputs continuous-profile specificity, all21 IUT/BH/BY tests, donor/human stability.
Central config/ROOT paths, seed20261002, exhaustive donor flips,1000 human bootstraps.
No existing ranks or therapeutic claims modified; all outputs under analysis/specificity.
"""
import itertools,json,logging
import numpy as np,pandas as pd
from scipy.stats import rankdata,t,false_discovery_control,spearmanr
from pipeline_utils import ROOT,path,sha256,write_json,run_standard_module

def readj(p):return json.loads(p.read_text(encoding='utf-8'))
def weights(fc,top=None):
 pos=np.clip(fc,0,None);neg=np.clip(-fc,0,None)
 if top is not None:
  ip=np.flatnonzero(pos>0);im=np.flatnonzero(neg>0);assert len(ip)>=top and len(im)>=top
  p=np.zeros(len(pos));m=np.zeros(len(neg));p[ip[np.argsort(-pos[ip],kind='stable')[:top]]]=1;m[im[np.argsort(-neg[im],kind='stable')[:top]]]=1;pos,neg=p,m
 assert pos.sum()>0 and neg.sum()>0
 pos=pos/pos.sum();neg=neg/neg.sum();return neg-pos,pos,neg
def infer(a):
 a=np.asarray(a,dtype=float);n=len(a);obs=float(a.mean());null=np.asarray(list(itertools.product([-1,1],repeat=n)))@a/n
 se=a.std(ddof=1)/np.sqrt(n);margin=t.ppf(.975,n-1)*se
 return {'mean':obs,'t_ci_low':float(obs-margin),'t_ci_high':float(obs+margin),'exact_one_sided_p':float(np.mean(null>=obs-1e-14)),'positive_donors':int((a>0).sum()),'donors':n}
def endpoints(scores):return np.column_stack([scores[:,0],scores[:,0]-scores[:,1],scores[:,0]-scores[:,2]])

def main():
 out=path('analysis')/'specificity';ex=path('analysis')/'expansion';cfg=readj(ROOT/'config/specificity_v1.json')
 for item in readj(out/'specificity_plan_freeze.json')['files']:assert sha256(ROOT/item['path'])==item['sha256'],item['path']
 common=pd.read_csv(out/'specificity_common_genes.tsv',sep='\t',dtype=str)
 h=pd.read_csv(ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz',sep='\t',index_col=0);h.index=h.index.astype(str)
 hs=pd.read_csv(ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv',sep='\t').set_index('gsm').loc[h.columns]
 hv=h.loc[common.ENTREZID].to_numpy();groups=['nonfailing','sepsis','IHD','DCM'];indices={g:np.flatnonzero(hs.group.to_numpy()==g) for g in groups}
 profile_names=cfg['primary_human_profiles']+cfg['secondary_human_profiles'];pairs=[p.split('_vs_') for p in profile_names]
 def make_profiles(use=None):
  mean={g:hv[:,ix if use is None else ix[use[ix]]].mean(axis=1) for g,ix in indices.items()}
  return np.column_stack([mean[a]-mean[b] for a,b in pairs])
 fc=make_profiles();ws=[weights(fc[:,i]) for i in range(5)];w=np.column_stack([a[0] for a in ws]);up=np.column_stack([a[1] for a in ws]);down=np.column_stack([a[2] for a in ws])
 assert np.allclose(w.sum(axis=0),0,atol=1e-14) and np.allclose(up.sum(axis=0),1) and np.allclose(down.sum(axis=0),1)
 table=common.copy()
 for k,name in enumerate(profile_names):table[f'logFC_{name}']=fc[:,k];table[f'weight_{name}']=w[:,k];table[f'UP_weight_{name}']=up[:,k];table[f'DOWN_weight_{name}']=down[:,k]
 table.to_csv(out/'human_continuous_profiles.tsv',sep='\t',index=False)
 # Direct mean differences must reproduce frozen limma contrast coefficients.
 errors=[]
 for label in ['nonfailing','IHD','DCM']:
  prev=pd.read_csv(ROOT/f'outputs/analysis/human/primary_sepsis_vs_{label}.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
  e=float(np.max(abs(fc[:,profile_names.index('sepsis_vs_'+label)]-prev.loc[common.ENTREZID,'logFC'].to_numpy())));errors.append(e);assert e<1e-10
 pd.DataFrame(fc,columns=profile_names).corr(method='spearman').to_csv(out/'human_profile_correlations.tsv',sep='\t')
 flag=readj(ROOT/'config/study_config.json')['analysis']['qc_flagged_samples'];without=~hs.index.isin(flag)
 fc_qc=make_profiles(without);w_qc=np.column_stack([weights(fc_qc[:,i])[0] for i in range(5)]);w_top=np.column_stack([weights(fc[:,i],top=50)[0] for i in range(5)])
 # Read only needed rows from immutable-derived count copy; no new exclusions.
 pieces=[]
 for chunk in pd.read_csv(ex/'expanded_counts.tsv.gz',sep='\t',index_col=0,chunksize=1000):pieces.append(chunk.loc[chunk.index.isin(common.SYMBOL)])
 x=pd.concat(pieces).loc[common.SYMBOL];assert len(x)==876 and np.isfinite(x.to_numpy()).all()
 ranks=(rankdata(x.to_numpy(),axis=0,method='average')-.5)/len(x);ranks=pd.DataFrame(ranks,index=common.SYMBOL,columns=x.columns)
 design=pd.read_csv(ex/'expanded_analysis_design.tsv',sep='\t');prior=pd.read_csv(ex/'cross_context_comparison.tsv',sep='\t').set_index('state')
 donor_vectors={};profile_scores=[];testrows=[];rowresults=[];loo=[];sensitivity=[]
 for drug,ds in design.groupby('drug_state'):
  dd={}
  for donor,ss in ds.groupby('Cell'):
   contrasts=[]
   for block,b in ss.groupby('block'):
    v=ranks[b.loc[b.condition=='DRUG','column_id']].mean(axis=1)-ranks[b.loc[b.condition=='CTRL','column_id']].mean(axis=1);contrasts.append(v.to_numpy())
   dd[donor]=np.mean(contrasts,axis=0)
  dm=pd.DataFrame(dd,index=common.SYMBOL).T;donor_vectors[drug]=dm
  a=dm.to_numpy()@w;ep=endpoints(a);tests=[infer(ep[:,i]) for i in range(3)]
  record={'drug_state':drug,'drug_name':prior.loc[drug,'drug_name'],'donors':len(dm),'prior_bidirectional':bool(prior.loc[drug,'both_directions_opposing']),'sepsis_reversal':a[:,0].mean(),'IHD_reversal':a[:,1].mean(),'DCM_reversal':a[:,2].mean(),'sepsis_minus_IHD':ep[:,1].mean(),'sepsis_minus_DCM':ep[:,2].mean(),'all_three_positive':bool((ep.mean(axis=0)>0).all()),'iut_p':max(ti['exact_one_sided_p'] for ti in tests)}
  for i,name in enumerate(cfg['specificity_endpoints']):
   record[f'{name}_p']=tests[i]['exact_one_sided_p'];testrows.append({'drug_state':drug,'endpoint':name,**tests[i]})
  for i,donor in enumerate(dm.index):
   for j,name in enumerate(profile_names):profile_scores.append({'drug_state':drug,'donor':donor,'profile':name,'reversal':a[i,j],'UP_change':float(dm.iloc[i].to_numpy()@up[:,j]),'DOWN_change':float(dm.iloc[i].to_numpy()@down[:,j])})
  leave_values=[]
  for i,donor in enumerate(dm.index):
   val=np.delete(ep,i,axis=0).mean(axis=0);leave_values.append(val);loo.append({'drug_state':drug,'excluded_donor':donor,**dict(zip(cfg['specificity_endpoints'],val)),'all_three_positive':bool((val>0).all())})
  record['all_LODO_positive']=bool((np.asarray(leave_values)>0).all())
  for name,sw in [('without_human_QC_flag',w_qc),('equal_weight_top50',w_top)]:
   val=endpoints(dm.to_numpy()@sw).mean(axis=0);sensitivity.append({'drug_state':drug,'sensitivity':name,**dict(zip(cfg['specificity_endpoints'],val)),'all_three_positive':bool((val>0).all())});record[name+'_positive']=bool((val>0).all())
  rowresults.append(record)
 pd.concat(donor_vectors,names=['drug_state','donor']).to_csv(out/'donor_rank_change_vectors.tsv.gz',sep='\t')
 result=pd.DataFrame(rowresults).sort_values('drug_state');result['iut_BH_q']=false_discovery_control(result.iut_p);result['iut_BY_q']=false_discovery_control(result.iut_p,method='by')
 # Patient bootstrap uses one resampled shared-control mean for all comparisons.
 rng=np.random.default_rng(cfg['seed']);B=cfg['bootstrap_iterations'];bootw=np.empty((B,3,len(common)));draws=[]
 for b in range(B):
  means={}
  for g,ix in indices.items():
   use=rng.choice(ix,size=len(ix),replace=True);means[g]=hv[:,use].mean(axis=1);draws.append({'iteration':b,'group':g,'sample_ids':';'.join(h.columns[use])})
  for k,g in enumerate(['sepsis','IHD','DCM']):bootw[b,k]=weights(means[g]-means['nonfailing'])[0]
 np.savez_compressed(out/'human_bootstrap_weights.npz',weights=bootw,entrez=common.ENTREZID.to_numpy(dtype=str))
 pd.DataFrame(draws).to_csv(out/'human_bootstrap_draws.tsv.gz',sep='\t',index=False)
 bootrows=[]
 for idx,r in result.iterrows():
  vector=donor_vectors[r.drug_state].mean(axis=0).to_numpy();scores=bootw@vector;ep=endpoints(scores);freq=float((ep>0).all(axis=1).mean());result.loc[idx,'human_bootstrap_joint_positive_frequency']=freq
  for j,name in enumerate(cfg['specificity_endpoints']):
   result.loc[idx,'bootstrap_'+name+'_p025']=float(np.quantile(ep[:,j],.025));result.loc[idx,'bootstrap_'+name+'_p975']=float(np.quantile(ep[:,j],.975))
  for b,v in enumerate(ep):bootrows.append({'drug_state':r.drug_state,'iteration':b,**dict(zip(cfg['specificity_endpoints'],v))})
 result['formal_specificity_gate']=(result.iut_BH_q<.05)&result.all_three_positive
 result['exploratory_followup_gate']=result.prior_bidirectional & result.all_three_positive & result.all_LODO_positive & result.without_human_QC_flag_positive & result.equal_weight_top50_positive & (result.human_bootstrap_joint_positive_frequency>=cfg['exploratory_gate_threshold'])
 result.to_csv(out/'specificity_results.tsv',sep='\t',index=False)
 for name,rows in [('donor_profile_scores.tsv',profile_scores),('component_exact_tests.tsv',testrows),('leave_one_donor_specificity.tsv',loo),('specificity_sensitivities.tsv',sensitivity),('human_bootstrap_endpoint_values.tsv.gz',bootrows)]:pd.DataFrame(rows).to_csv(out/name,sep='\t',index=False)
 means=pd.DataFrame(profile_scores).groupby(['drug_state','profile']).reversal.mean().unstack();means.to_csv(out/'all_five_profile_mean_scores.tsv',sep='\t')
 human_corr=pd.DataFrame(fc[:,:3],columns=profile_names[:3]).corr(method='spearman')
 summary={'version':cfg['version'],'drugs':len(result),'common_genes':len(common),'formal_specificity_supported':int(result.formal_specificity_gate.sum()),'exploratory_followup_pass':int(result.exploratory_followup_gate.sum()),'exploratory_followup_names':result.loc[result.exploratory_followup_gate,'drug_name'].tolist(),'all_three_point_positive':int(result.all_three_positive.sum()),'minimum_iut_BH_q':float(result.iut_BH_q.min()),'maximum_human_bootstrap_joint_positive_frequency':float(result.human_bootstrap_joint_positive_frequency.max()),'human_bootstrap_iterations':B,'human_profile_spearman':human_corr.to_dict(),'prior_therapeutic_supported_candidates':0,'decision':'ONLY_INDEPENDENT_SPECIFICITY_HYPOTHESIS_FOLLOWUP' if result.exploratory_followup_gate.any() or result.formal_specificity_gate.any() else 'STOP_ESCALATING_CURRENT_SEPSIS_SPECIFIC_SINGLE_DRUG_CLAIM','scope':'conditional same-cohort relative specificity; no external confirmation or clinical efficacy'}
 write_json(out/'specificity_summary.json',summary);write_json(out/'implementation_checks.json',{'old_human_contrast_max_errors':errors,'normalized_weight_sums':'PASS','same_shared_control_bootstrap':'single sample draw per group per iteration','frozen_inputs':'PASS'})
 logging.info('Specificity result: %s',summary);return summary

if __name__=='__main__':run_standard_module('37_disease_specificity',main)
