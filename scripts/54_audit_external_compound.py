"""Independent audit from original measured cardiac counts and block design.
Packages pandas,numpy; input frozen52; output numeric_qa.json. No raw alterations.
Fixed seeds not needed for exhaustive checks. Logs via shared wrapper.
"""
import itertools,json
import numpy as np,pandas as pd
from pipeline_utils import ROOT,sha256,write_json,run_standard_module
OUT=ROOT/'outputs/analysis/external_compound_v1'
def bh(p):
 p=np.asarray(p);order=np.argsort(p);ans=np.empty(len(p));ans[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return ans
def prob(v):
 # Enumerate subsets of negative signs via binary masks; total -2*subset sum.
 n=len(v);obs=sum(v)/n;null=[]
 for mask in range(2**n):null.append((sum(v)-2*sum(v[i] for i in range(n) if mask&(1<<i)))/n)
 return sum(z>=obs-1e-14 for z in null)/len(null)
def main():
 for f in json.loads((OUT/'plan_freeze.json').read_text(encoding='utf8'))['files']:assert sha256(ROOT/f['path'])==f['sha256']
 common876=pd.read_csv(ROOT/'outputs/analysis/specificity/specificity_common_genes.tsv',sep='\t',dtype=str)
 pieces=[]
 for chunk in pd.read_csv(ROOT/'outputs/analysis/expansion/expanded_counts.tsv.gz',sep='\t',index_col=0,chunksize=1000):pieces.append(chunk.loc[chunk.index.isin(common876.SYMBOL)])
 x=pd.concat(pieces).loc[common876.SYMBOL];rank=(x.rank(axis=0,method='average')-.5)/len(x)
 design=pd.read_csv(ROOT/'outputs/analysis/expansion/expanded_analysis_design.tsv',sep='\t');previous=pd.read_csv(ROOT/'outputs/analysis/specificity/donor_rank_change_vectors.tsv.gz',sep='\t',index_col=[0,1])
 rebuild={};maxerr=0
 for (drug,donor),s in design.groupby(['drug_state','Cell']):
  block_values=[]
  for block in s.block.unique():
   b=s[s.block.eq(block)];drugcols=b.loc[b.condition.eq('DRUG'),'column_id'].tolist();ctrlcols=b.loc[b.condition.eq('CTRL'),'column_id'].tolist();assert drugcols and ctrlcols
   block_values.append(rank[drugcols].sum(axis=1)/len(drugcols)-rank[ctrlcols].sum(axis=1)/len(ctrlcols))
  vector=pd.concat(block_values,axis=1).mean(axis=1);rebuild[drug,donor]=vector
  maxerr=max(maxerr,float(np.max(abs(vector.loc[previous.columns].values-previous.loc[(drug,donor)].values))))
 assert maxerr<1e-12
 w=pd.read_csv(OUT/'disease_weights.tsv',sep='\t',dtype={'ENTREZID':str});scores=pd.read_csv(OUT/'donor_scores.tsv',sep='\t');results=pd.read_csv(OUT/'all_drug_results.tsv',sep='\t')
 assert len(results)==63 and len(results.drug_state.unique())==21
 error=0
 for row in scores.itertuples():
  z=rebuild[row.drug_state,row.donor].loc[w.SYMBOL].to_numpy();fc=w[row.cohort+'_logFC'].to_numpy();pos=fc>0;neg=fc<0
  u=-np.average(z[pos],weights=fc[pos]);d=np.average(z[neg],weights=-fc[neg]);error=max(error,abs(u-row.UP_suppression),abs(d-row.DOWN_restoration));assert abs(u+d-row.net_reversal)<1e-12
 assert error<1e-12
 for cohort,t in results.groupby('cohort'):
  ps=[]
  for row in t.itertuples():
   d=scores[scores.cohort.eq(cohort)&scores.drug_state.eq(row.drug_state)];p=max(prob(d.UP_suppression.to_numpy()),prob(d.DOWN_restoration.to_numpy()));assert abs(p-row.iut_p)<1e-12;ps.append(p)
   assert len(d)==row.donors and row.minimum_attainable_p==2**(-len(d))
  assert np.allclose(bh(ps),t.iut_BH_q,atol=1e-12)
 for typ in ['donor','patient']:
  loo=pd.read_csv(OUT/(typ+'_loo.tsv'),sep='\t')
  assert (loo.both_positive==((loo.UP_suppression>0)&(loo.DOWN_restoration>0))).all()
  for (cohort,drug),ss in loo.groupby(['cohort','drug_state']):
   row=results[results.cohort.eq(cohort)&results.drug_state.eq(drug)].iloc[0]
   if typ=='donor':
    d=scores[scores.cohort.eq(cohort)&scores.drug_state.eq(drug)]
    for r in ss.itertuples():assert np.allclose(d.loc[d.donor.ne(r.omitted),['UP_suppression','DOWN_restoration']].mean(),[r.UP_suppression,r.DOWN_restoration],atol=1e-12)
   else:
    fc=pd.read_csv(OUT/(cohort+'_patient_loo_logFC.tsv'),sep='\t',index_col=0)
    dm=pd.DataFrame([rebuild[drug,donor].loc[w.SYMBOL] for donor in scores.loc[scores.drug_state.eq(drug),'donor'].unique()]).mean().to_numpy()
    for r in ss.itertuples():
     v=fc[r.omitted].to_numpy();p=v>0;n=v<0;assert np.allclose([-np.average(dm[p],weights=v[p]),np.average(dm[n],weights=-v[n])],[r.UP_suppression,r.DOWN_restoration],atol=1e-12)
   assert bool(ss.both_positive.all())==bool(row['all_'+typ+'_loo_bidirectional'])
 cross=pd.read_csv(OUT/'cross_cohort_decision.tsv',sep='\t');a=results[results.cohort.eq('GSE237861')].set_index('drug_state').loc[cross.drug_state];b=results[results.cohort.eq('GSE141864')].set_index('drug_state').loc[cross.drug_state]
 assert np.allclose(cross.four_component_iut_p,np.maximum(a.iut_p,b.iut_p));assert np.allclose(cross.secondary_BH_q,bh(cross.four_component_iut_p))
 gate=a.formal_bidirectional_gate.to_numpy() & cross.both_external_bidirectional & cross.all_external_donor_loo & cross.all_external_patient_loo
 assert (gate==cross.escalation_gate).all()
 qa={'status':'PASS','reconstructed_cardio_donor_vectors':len(rebuild),'vector_max_abs_error':maxerr,'independent_weighted_score_max_abs_error':error,'IUT_tests_audited':63,'BH_families':3,'four_component_secondary_BH':'PASS','donor_patient_omission_checks':'PASS','upstream_hashes':'UNCHANGED','scope':'Numerical/source audit, not independent wet-lab or clinical validation'}
 write_json(OUT/'numeric_qa.json',qa);return qa
if __name__=='__main__':run_standard_module('54_audit_external_compound',main)
