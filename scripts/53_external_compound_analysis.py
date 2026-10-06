"""Packages numpy,pandas,scipy; R edgeR. Inputs frozen52/config; output external_compound_v1.
Seed20261005, exhaustive donor sign flips. Conditions on disease estimates; no efficacy.
"""
import json,itertools,subprocess
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control,t
from pipeline_utils import ROOT,RUNTIME,sha256,write_json,run_standard_module,r_environment
OUT=ROOT/'outputs/analysis/external_compound_v1'
def weights(v,top=False):
 v=np.asarray(v,float);p=np.clip(v,0,None);n=np.clip(-v,0,None)
 if top:
  assert (p>0).sum()>=50 and (n>0).sum()>=50
  p=np.isin(np.arange(len(v)),np.argsort(-p,kind='stable')[:50]).astype(float)
  n=np.isin(np.arange(len(v)),np.argsort(-n,kind='stable')[:50]).astype(float)
 assert p.sum()>0 and n.sum()>0
 return np.c_[-p/p.sum(),n/n.sum()]
def inference(v):
 v=np.asarray(v);n=len(v);m=v.mean();null=np.array(list(itertools.product([-1,1],repeat=n)))@v/n
 width=t.ppf(.975,n-1)*v.std(ddof=1)/np.sqrt(n)
 return float(np.mean(null>=m-1e-14)),float(m-width),float(m+width)
def main():
 for f in json.loads((OUT/'plan_freeze.json').read_text(encoding='utf8'))['files']:assert sha256(ROOT/f['path'])==f['sha256']
 common=pd.read_csv(OUT/'common_genes.tsv',sep='\t',dtype=str);drugs=pd.read_csv(OUT/'frozen_drugs.tsv',sep='\t')
 raw=pd.read_csv(ROOT/'outputs/analysis/specificity/donor_rank_change_vectors.tsv.gz',sep='\t',index_col=[0,1])
 assert raw.shape[1]==876
 profiles={};paths={'GSE79962':'human/primary_sepsis_vs_nonfailing.tsv','GSE237861':'raw_heart_v1/GSE237861_all_gene_results.tsv','GSE141864':'external_human/GSE141864_all_gene_results.tsv'}
 for k,p in paths.items():profiles[k]=pd.read_csv(ROOT/'outputs/analysis'/p,sep='\t',dtype={'entrez_id':str}).set_index('entrez_id').loc[common.ENTREZID,'logFC'].to_numpy()
 r=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/53_external_compound_loo.R')],cwd=ROOT,env=r_environment(),capture_output=True,text=True,encoding='utf8',errors='replace')
 (OUT/'R_console.txt').write_text(r.stdout+r.stderr,encoding='utf8');assert r.returncode==0,r.stderr[-1800:]
 loo237=pd.read_csv(OUT/'GSE237861_patient_loo_logFC.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id').loc[common.ENTREZID]
 assert np.max(abs(loo237['full'].to_numpy()-profiles['GSE237861']))<1e-8
 x=pd.read_csv(ROOT/'outputs/analysis/external_human/GSE141864_gene_expression.tsv.gz',sep='\t',index_col=0);x.index=x.index.astype(str)
 s=pd.read_csv(ROOT/'outputs/analysis/external_human/GSE141864_heart_design.tsv',sep='\t');s=s[s.primary.astype(str).str.lower().eq('true')]
 x=x.loc[common.ENTREZID,s.gsm];g=s.group.eq('sepsis').to_numpy();loo141={}
 for i,patient in enumerate(s.gsm):
  keep=np.arange(len(s))!=i;loo141[patient]=x.iloc[:,keep].iloc[:,g[keep]].mean(axis=1)-x.iloc[:,keep].iloc[:,~g[keep]].mean(axis=1)
 pd.DataFrame(loo141).to_csv(OUT/'GSE141864_patient_loo_logFC.tsv',sep='\t')
 weights_table=common.copy();rows=[];donorrows=[];lod=[];lop=[];sensitivity=[]
 for cohort,fc in profiles.items():
  w=weights(fc);weights_table[cohort+'_logFC']=fc;weights_table[cohort+'_UP_weight']=-w[:,0];weights_table[cohort+'_DOWN_weight']=w[:,1]
  for drug in drugs.itertuples():
   dm=raw.loc[drug.drug_state,common.SYMBOL];assert len(dm)==drug.donors
   scores=dm.to_numpy()@w;pv=[];rec={'cohort':cohort,'drug_state':drug.drug_state,'drug_name':drug.drug_name,'donors':len(dm),'minimum_attainable_p':2.**(-len(dm))}
   for j,name in enumerate(['UP_suppression','DOWN_restoration','net_reversal']):
    v=scores[:,j] if j<2 else scores.sum(axis=1);p,lo,hi=inference(v);rec.update({name:float(v.mean()),name+'_p':p,name+'_ci_low':lo,name+'_ci_high':hi});pv.append(p)
   rec['iut_p']=max(pv[:2]);rec['bidirectional']=bool((scores.mean(axis=0)>0).all())
   for i,donor in enumerate(dm.index):
    donorrows.append({'cohort':cohort,'drug_state':drug.drug_state,'donor':donor,'UP_suppression':scores[i,0],'DOWN_restoration':scores[i,1],'net_reversal':scores[i].sum()})
    value=np.delete(scores,i,axis=0).mean(axis=0);lod.append({'cohort':cohort,'drug_state':drug.drug_state,'omitted':donor,'UP_suppression':value[0],'DOWN_restoration':value[1],'both_positive':bool((value>0).all())})
   rec['all_donor_loo_bidirectional']=bool(all(z['both_positive'] for z in lod if z['cohort']==cohort and z['drug_state']==drug.drug_state))
   if cohort!='GSE79962':
    human=loo237.drop(columns='full') if cohort=='GSE237861' else pd.DataFrame(loo141)
    flags=[]
    for patient in human:
     value=dm.mean(axis=0).to_numpy()@weights(human[patient].to_numpy());flag=bool((value>0).all());flags.append(flag)
     lop.append({'cohort':cohort,'drug_state':drug.drug_state,'omitted':patient,'UP_suppression':value[0],'DOWN_restoration':value[1],'both_positive':flag})
    rec['patient_loo_bidirectional_fraction']=np.mean(flags);rec['all_patient_loo_bidirectional']=all(flags)
   value=dm.mean(axis=0).to_numpy()@weights(fc,top=True);sensitivity.append({'cohort':cohort,'drug_state':drug.drug_state,'UP_suppression':value[0],'DOWN_restoration':value[1],'both_positive':bool((value>0).all())})
   rows.append(rec)
 result=pd.DataFrame(rows)
 for cohort,ix in result.groupby('cohort').groups.items():
  result.loc[ix,'iut_BH_q']=false_discovery_control(result.loc[ix,'iut_p']);result.loc[ix,'iut_BY_q']=false_discovery_control(result.loc[ix,'iut_p'],method='by')
 result['formal_bidirectional_gate']=(result.iut_BH_q<.05)&result.bidirectional
 a=result[result.cohort.eq('GSE237861')].set_index('drug_state');b=result[result.cohort.eq('GSE141864')].set_index('drug_state').loc[a.index]
 cross=drugs.set_index('drug_state').loc[a.index].copy();cross['four_component_iut_p']=np.maximum(a.iut_p,b.iut_p);cross['secondary_BH_q']=false_discovery_control(cross.four_component_iut_p)
 cross['both_external_bidirectional']=a.bidirectional & b.bidirectional
 cross['all_external_donor_loo']=a.all_donor_loo_bidirectional & b.all_donor_loo_bidirectional
 cross['all_external_patient_loo']=a.all_patient_loo_bidirectional.astype(bool) & b.all_patient_loo_bidirectional.astype(bool)
 cross['escalation_gate']=a.formal_bidirectional_gate & cross.both_external_bidirectional & cross.all_external_donor_loo & cross.all_external_patient_loo
 weights_table.to_csv(OUT/'disease_weights.tsv',sep='\t',index=False);result.to_csv(OUT/'all_drug_results.tsv',sep='\t',index=False);cross.to_csv(OUT/'cross_cohort_decision.tsv',sep='\t')
 for name,items in [('donor_scores',donorrows),('donor_loo',lod),('patient_loo',lop),('equal_top50_sensitivity',sensitivity)]:pd.DataFrame(items).to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
 summary={'version':'external-compound-1.0','common_genes':len(common),'drugs':21,'primary_formal_hits':int(a.formal_bidirectional_gate.sum()),'primary_min_q':float(a.iut_BH_q.min()),'primary_bidirectional_means':int(a.bidirectional.sum()),'both_external_bidirectional':int(cross.both_external_bidirectional.sum()),'escalation_gate_hits':int(cross.escalation_gate.sum()),'escalation_names':cross.loc[cross.escalation_gate,'drug_name'].tolist(),'therapeutic_candidates_supported':0,'patient_loo_fits':{'GSE237861':14,'GSE141864':7},'interpretation':'New disease context, reused cardiac perturbation donors; not independent drug efficacy validation'}
 write_json(OUT/'summary.json',summary);return summary
if __name__=='__main__':run_standard_module('53_external_compound_analysis',main)
