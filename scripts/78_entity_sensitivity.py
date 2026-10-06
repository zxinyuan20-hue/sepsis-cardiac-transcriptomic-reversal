"""Packages numpy/pandas/scipy/numba/RDKit already in project runtime.
Inputs/outputs: robustness_amendment_v1 config; frozen LINCS arrays, keys, nulls, bootstrap.
Fixed seed and query draws reused; no new chemical eligibility or biological data.
Pipeline: exact verified structure mapping -> signature aggregation -> null/BH/stability.
Saves operational entity mapping, full scores/null/bootstrap, duplicate audit and summary.
"""
import json
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
from numba import set_num_threads
from pipeline_utils import ROOT,run_standard_module,write_json
from screening_core import score_queries
from robustness_utils import CFG,OUT,verify_inputs,layout

def main():
 verify_inputs();set_num_threads(CFG['threads']);out=OUT/'entities';out.mkdir(exist_ok=True)
 src=ROOT/'outputs/analysis/screening';a=np.load(src/'screening_inputs.npz')
 ranking=pd.read_csv(src/'complete_screening_results.tsv',sep='\t')
 dictionary=pd.read_csv(ROOT/'outputs/preparation/lincs_compound_dictionary.tsv',sep='\t',dtype=str)
 keys=pd.read_csv(ROOT/'outputs/analysis/expansion/screened_structure_keys.tsv',sep='\t',dtype=str)
 ident=ranking[['pert_id','pert_iname']].merge(dictionary[['pert_id','inchi_key']],on='pert_id',validate='one_to_one').merge(keys[['pert_id','computed_key','tautomer_key','source_key_agrees']],on='pert_id',how='left',validate='one_to_one')
 valid=ident.inchi_key.str.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]').fillna(False)&ident.inchi_key.eq(ident.computed_key)
 ident['identity_verified']=valid
 ident['entity_id']=np.where(valid,'KEY:'+ident.inchi_key,'UNRESOLVED:'+ident.pert_id)
 ident['mapping_rule']=np.where(valid,'exact_full_key_and_computed_SMILES_agreement','retain_original_unresolved_ID')
 ident.to_csv(out/'entry_to_entity.tsv',sep='\t',index=False)
 meta=pd.read_csv(src/'signature_metadata.tsv.gz',sep='\t');assert meta.sig_id.is_unique
 assert len(meta)==len(a['orders'])
 # Preserve positional alignment with frozen signature arrays.
 mapping=ident.set_index('pert_id').entity_id;meta['entity_id']=meta.pert_id.map(mapping);assert meta.entity_id.notna().all()
 entity_order=sorted(ident.entity_id.unique());oldidx={k:i for i,k in enumerate(a['compounds'])}
 group=ident.groupby('entity_id',sort=True)
 observed_old=ranking.set_index('pert_id').loc[a['compounds'],'reversal_score'].to_numpy()
 oldnull=np.load(src/'matched_null_compound_scores.npy');oldboot=np.load(src/'bootstrap_compound_scores.npy')
 allscores=np.empty(len(entity_order));null=np.empty((len(oldnull),len(entity_order)));boot=np.empty((len(oldboot),len(entity_order)))
 bootq=pd.read_csv(src/'bootstrap_queries.tsv.gz',sep='\t');assert list(bootq.columns)==list(a['genes'])
 quality=pd.read_csv(src/'bootstrap_query_quality.tsv',sep='\t').valid.to_numpy(bool)
 q=np.load(src/'matched_null_queries.npz')['queries'];primary=a['queries'][list(a['query_names']).index('size_100')]
 duplicate_rows=[]
 for j,e in enumerate(entity_order):
  ids=group.get_group(e).pert_id.tolist()
  if len(ids)==1:
   k=oldidx[ids[0]];allscores[j]=observed_old[k];null[:,j]=oldnull[:,k];boot[:,j]=oldboot[:,k]
  else:
   ix=np.flatnonzero(meta.entity_id.eq(e));subset=meta.iloc[ix].reset_index(drop=True)
   entities,order,starts,comp=layout(subset,'entity_id');assert entities==[e]
   args=[a['orders'][ix],a['weights'][ix],order,starts,comp]
   allscores[j]=score_queries(primary[None,:],*args)[0,0]
   null[:,j]=score_queries(q,*args)[:,0]
   boot[:,j]=np.nan;boot[quality,j]=score_queries(bootq.to_numpy(np.int8)[quality],*args)[:,0]
   duplicate_rows.append({'entity_id':e,'pert_ids':'|'.join(ids),'signatures':len(ix),'cells':subset.cell_id.nunique(),'score':allscores[j]})
 exceed=(null>=allscores).sum(axis=0);p=(1+exceed)/(len(null)+1);p[allscores<=0]=1;bh=false_discovery_control(p)
 brank=np.full_like(boot,np.nan);brank[quality]=rankdata(-boot[quality],axis=1,method='average')
 result=pd.DataFrame({'entity_id':entity_order,'score':allscores,'empirical_p':p,'BH_q':bh,'null_exceedances':exceed,'bootstrap_positive_frequency':((boot>0)&np.isfinite(boot)).mean(axis=0),'bootstrap_rank_median':np.nanmedian(brank,axis=0),'bootstrap_rank_p05':np.nanquantile(brank,.05,axis=0),'bootstrap_rank_p95':np.nanquantile(brank,.95,axis=0),'rank_average_ties':rankdata(-allscores),'source_entry_count':[len(group.get_group(e)) for e in entity_order],'names':['|'.join(sorted(set(group.get_group(e).pert_iname))) for e in entity_order]})
 result.sort_values(['rank_average_ties','entity_id']).to_csv(out/'entity_results.tsv',sep='\t',index=False)
 pd.DataFrame(duplicate_rows).to_csv(out/'merged_entities.tsv',sep='\t',index=False)
 np.savez_compressed(out/'entity_distributions.npz',entity_ids=np.array(entity_order),null=null,bootstrap=boot)
 # Independent pandas recomputation of the observed entity score.
 sig=pd.read_csv(src/'primary_signature_scores.tsv.gz',sep='\t');sig['entity_id']=sig.pert_id.map(mapping)
 check=sig.groupby(['entity_id','cell_id']).reversal_score.median().groupby('entity_id').median().loc[entity_order].to_numpy()
 err=float(np.max(abs(check-allscores)));assert err<1e-12
 # No new supported candidate gate is claimed from statistical sensitivity alone.
 summary={'entries':len(ident),'verified_source_structure_entries':int(valid.sum()),'verified_unique_full_keys':int(ident.loc[valid,'entity_id'].nunique()),'unresolved_entries':ident.loc[~valid,'pert_id'].tolist(),'operational_entities_including_unresolved':len(result),'merged_groups':len(duplicate_rows),'positive_entities':int((allscores>0).sum()),'BH_q_lt_05_positive_entities':int(((bh<.05)&(allscores>0)).sum()),'minimum_positive_BH_q':float(bh[allscores>0].min()),'null_queries':len(q),'bootstrap_draws':len(boot),'valid_bootstrap_draws':int(quality.sum()),'independent_aggregation_max_error':err,'scope':'exact-source-structure operational sensitivity; not a fully standardized chemical-parent count or new therapeutic validation'}
 write_json(out/'summary.json',summary);return summary

if __name__=='__main__':run_standard_module('78_entity_sensitivity',main)
