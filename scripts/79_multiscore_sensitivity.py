"""Packages numpy/pandas/scipy/numba. Inputs fixed signed query and frozen perturbations.
Paths: robustness_amendment_v1 config. Deterministic, no random hypothesis testing.
Pipeline: negative cosine and rank correlation with same query, fixed aggregation.
Outputs complete entry/entity/cardiac method scores and descriptive concordance.
"""
import itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
from pipeline_utils import ROOT,run_standard_module,write_json
from screening_core import aggregate
from robustness_utils import CFG,OUT,verify_inputs,layout

def alternatives(query,z):
 q=query.astype(float);z=np.asarray(z,float)
 den=np.linalg.norm(z,axis=1)*np.linalg.norm(q);assert (den>0).all()
 cosine=-(z@q)/den
 rq=rankdata(q);rq-=rq.mean();rz=rankdata(z,axis=1);rz-=rz.mean(axis=1)[:,None]
 sr=-(rz@rq)/(np.linalg.norm(rz,axis=1)*np.linalg.norm(rq))
 return cosine,sr

def comparisons(table,label):
 methods=['weighted_ES','negative_cosine','negative_Spearman'];rows=[]
 for m,n in itertools.combinations(methods,2):
  a=table[m].to_numpy();b=table[n].to_numpy();k=min(20,len(table));ta=a>=np.sort(a)[-k];tb=b>=np.sort(b)[-k]
  rows.append({'scope':label,'method_a':m,'method_b':n,'n':len(a),'rank_spearman':spearmanr(a,b).statistic,'sign_agreement_including_zero':np.mean(np.sign(a)==np.sign(b)),'positive_both':int(((a>0)&(b>0)).sum()),'topk':k,'tie_inclusive_a_n':int(ta.sum()),'tie_inclusive_b_n':int(tb.sum()),'topk_tie_inclusive_jaccard':(ta&tb).sum()/(ta|tb).sum(),'inference':'descriptive only; no independence-based P'})
 return rows

def main():
 verify_inputs();out=OUT/'multiscore';out.mkdir(exist_ok=True);src=ROOT/'outputs/analysis/screening'
 a=np.load(src/'screening_inputs.npz');q=a['queries'][list(a['query_names']).index('size_100')]
 c,s=alternatives(q,a['z']);original=pd.read_csv(src/'complete_screening_results.tsv',sep='\t').set_index('pert_id').loc[a['compounds']]
 args=[a[k] for k in ['cell_order','cell_starts','compound_starts']]
 entry=pd.DataFrame({'pert_id':a['compounds'],'name':original.pert_iname.values,'weighted_ES':original.reversal_score.values,'negative_cosine':aggregate(c,*args),'negative_Spearman':aggregate(s,*args)})
 entry.to_csv(out/'entry_method_scores.tsv',sep='\t',index=False)
 meta=pd.read_csv(src/'signature_metadata.tsv.gz',sep='\t');mapping=pd.read_csv(OUT/'entities/entry_to_entity.tsv',sep='\t').set_index('pert_id').entity_id
 meta['entity_id']=meta.pert_id.map(mapping);entities,order,starts,comp=layout(meta,'entity_id')
 er=pd.read_csv(OUT/'entities/entity_results.tsv',sep='\t').set_index('entity_id').loc[entities]
 entity=pd.DataFrame({'entity_id':entities,'names':er.names.values,'weighted_ES':er.score.values,'negative_cosine':aggregate(c,order,starts,comp),'negative_Spearman':aggregate(s,order,starts,comp)})
 entity.to_csv(out/'entity_method_scores.tsv',sep='\t',index=False)
 exp=ROOT/'outputs/analysis/expansion';common=pd.read_csv(exp/'common_landmark_background.tsv',sep='\t',dtype=str)
 query=pd.Series(q,index=a['genes']).loc[common.ENTREZID].to_numpy();card=pd.read_csv(exp/'cross_context_comparison.tsv',sep='\t');rows=[];donor=[]
 for r in card.itertuples():
  delta=pd.read_csv(exp/f'{r.state}_donor_logFC.tsv.gz',sep='\t',index_col=0).loc[common.SYMBOL]
  cc,ss=alternatives(query,delta.to_numpy().T)
  rows.append({'drug_state':r.state,'drug_name':r.drug_name,'weighted_ES':r.harmonized_cardiac_score,'negative_cosine':np.median(cc),'negative_Spearman':np.median(ss)})
  donor.extend({'drug_state':r.state,'donor':d,'negative_cosine':v,'negative_Spearman':w} for d,v,w in zip(delta.columns,cc,ss))
 cardiac=pd.DataFrame(rows);cardiac.to_csv(out/'cardiac_method_scores.tsv',sep='\t',index=False);pd.DataFrame(donor).to_csv(out/'cardiac_donor_method_scores.tsv',sep='\t',index=False)
 correlations=comparisons(entry,'482 LINCS entries')+comparisons(entity,'operational entities')+comparisons(cardiac,'21 cardiac drugs')
 pd.DataFrame(correlations).to_csv(out/'method_concordance.tsv',sep='\t',index=False)
 # Independent scipy pairwise calculation confirms vectorized Spearman orientation.
 err=max(abs(s[i]+spearmanr(q,a['z'][i]).statistic) for i in [0,1,10,100,len(s)-1]);assert err<1e-12
 summary={'query_up':int((q==1).sum()),'query_down':int((q==-1).sum()),'LINCS_background':len(q),'cardiac_background':len(query),'cardiac_query_up':int((query==1).sum()),'cardiac_query_down':int((query==-1).sum()),'Spearman_independent_max_error':err,'entity_all_three_positive':int((entity[['weighted_ES','negative_cosine','negative_Spearman']]>0).all(axis=1).sum()),'cardiac_all_three_positive':int((cardiac[['weighted_ES','negative_cosine','negative_Spearman']]>0).all(axis=1).sum()),'scope':'Fixed-query descriptive algorithm sensitivity; no alternative-score null tests or new formal support'}
 write_json(out/'summary.json',summary);return summary

if __name__=='__main__':run_standard_module('79_multiscore_sensitivity',main)
