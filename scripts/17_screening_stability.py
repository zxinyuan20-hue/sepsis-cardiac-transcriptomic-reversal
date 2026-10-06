"""Packages: numpy/pandas/numba/scipy. Inputs: completed primary/null and human bootstrap.
Outputs: complete stability table, fixed shortlist, per-cell/condition diagnostics and audit.
No candidate chosen until all prespecified gates applied; fixed config/seed; all outputs local.
"""
import json,logging
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from numba import set_num_threads
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module
from screening_core import score_queries

def main():
    cfg=json.loads((ROOT/'config/screening_v1.json').read_text());set_num_threads(cfg['threads'])
    for n in ['15_human_bootstrap','16_compound_screen']:
        if json.loads((path('audit')/f'{n}_status.json').read_text())['status']!='PASS':raise ValueError(f'{n} not complete')
    out=path('analysis')/'screening';a=np.load(out/'screening_inputs.npz')
    bootstrap_hashes=json.loads((out/'bootstrap_file_hashes.json').read_text())
    for name,digest in bootstrap_hashes.items():
        if sha256(out/name)!=digest:raise ValueError('Bootstrap input hash mismatch')
    q=pd.read_csv(out/'bootstrap_queries.tsv.gz',sep='\t');quality=pd.read_csv(out/'bootstrap_query_quality.tsv',sep='\t')
    if list(q.columns)!=list(a['genes']) or len(q)!=cfg['bootstrap_iterations']:raise ValueError('Bootstrap axis mismatch')
    valid=quality.valid.values
    if not valid.any():raise ValueError('All bootstrap queries failed coverage')
    bscore=np.full((len(q),len(a['compounds'])),np.nan)
    args=[a[k] for k in ['orders','weights','cell_order','cell_starts','compound_starts']]
    for start in range(0,len(q),100):
        ix=np.flatnonzero(valid[start:start+100])+start
        if len(ix):bscore[ix]=score_queries(q.iloc[ix].to_numpy(np.int8),*args)
        logging.info('BOOTSTRAP SCORING %d/%d',min(start+100,len(q)),len(q))
    np.save(out/'bootstrap_compound_scores.npy',bscore)
    rank=pd.read_csv(out/'primary_compound_ranking.tsv',sep='\t').set_index('pert_id').loc[a['compounds']].reset_index()
    top=np.zeros_like(bscore,dtype=bool);ranks=np.full_like(bscore,np.nan)
    for i in np.flatnonzero(valid):
        order=np.lexsort((a['compounds'],-bscore[i]));ranks[i,order]=np.arange(1,len(order)+1)
        take=order[:cfg['top_k_stability']];take=take[bscore[i,take]>0];top[i,take]=True
    rank['bootstrap_top20_frequency']=top.mean(axis=0)
    rank['bootstrap_positive_frequency']=((bscore>0)&np.isfinite(bscore)).mean(axis=0)
    rank['bootstrap_score_p025']=np.nanquantile(bscore,.025,axis=0)
    rank['bootstrap_score_median']=np.nanmedian(bscore,axis=0)
    rank['bootstrap_score_p975']=np.nanquantile(bscore,.975,axis=0)
    rank['bootstrap_rank_median']=np.nanmedian(ranks,axis=0)
    sensitivity=['score_size_50','score_size_150','score_without_qc_flag','score_submitted']
    rank['all_sensitivity_scores_positive']=(rank[sensitivity]>0).all(axis=1)
    rank['supported_candidate']=(rank.empirical_bh_q<.05)&(rank.reversal_score>0)&(rank.positive_cell_fraction>=.5)&(rank.bootstrap_top20_frequency>=.5)&rank.all_sensitivity_scores_positive
    rank=rank.sort_values('rank');rank.to_csv(out/'complete_screening_results.tsv',sep='\t',index=False)
    supported=rank[rank.supported_candidate].head(cfg['shortlist_maximum'])
    exploratory=rank[rank.reversal_score>0].head(cfg['shortlist_maximum'])
    supported.to_csv(out/'supported_shortlist.tsv',sep='\t',index=False)
    exploratory.to_csv(out/'exploratory_top10.tsv',sep='\t',index=False)
    # Independent pandas aggregation, BH, matched-label and count checks.
    sig=pd.read_csv(out/'primary_signature_scores.tsv.gz',sep='\t')
    agg=sig.groupby(['pert_id','cell_id']).reversal_score.median().groupby('pert_id').median()
    diff=float(np.max(np.abs(agg.loc[rank.pert_id].values-rank.reversal_score.values)))
    if diff>1e-12:raise ValueError('Hierarchical score mismatch')
    nullq=np.load(out/'matched_null_queries.npz')['queries'];primary=a['queries'][1]
    for s in np.unique(a['gene_strata']):
        ix=a['gene_strata']==s
        for label in [1,-1]:
            if not np.all((nullq[:,ix]==label).sum(axis=1)==(primary[ix]==label).sum()):raise ValueError('Null matching mismatch')
    if not np.allclose(false_discovery_control(rank.empirical_p.to_numpy()),rank.empirical_bh_q.to_numpy(),atol=1e-12):raise ValueError('BH mismatch')
    null=np.load(out/'matched_null_compound_scores.npy')
    original=rank.set_index('pert_id').loc[a['compounds']]
    exceed=(null>=original.reversal_score.to_numpy()[None,:]).sum(axis=0)
    if not np.array_equal(exceed,original.null_exceedances.to_numpy()):raise ValueError('Null count mismatch')
    # Leave-one-cell diagnostics preserve single-cell remainder labels.
    cells=pd.read_csv(out/'primary_compound_cell_scores.tsv',sep='\t');loo=[]
    for pert,sub in cells.groupby('pert_id'):
        for row in sub.itertuples():
            rest=sub[sub.cell_id!=row.cell_id]
            loo.append({'pert_id':pert,'removed_cell':row.cell_id,'remaining_cells':len(rest),'reversal_score':float(rest.reversal_score.median())})
    pd.DataFrame(loo).to_csv(out/'leave_one_cell_out_scores.tsv',sep='\t',index=False)
    summary={'frozen_utc':now(),'family_size':len(rank),'null_iterations':len(null),'bootstrap_planned':len(q),'bootstrap_valid':int(valid.sum()),'bootstrap_failed_coverage':int((~valid).sum()),'q_lt_05_positive':int(((rank.empirical_bh_q<.05)&(rank.reversal_score>0)).sum()),'supported_gate_pass':int(rank.supported_candidate.sum()),'supported_shortlist_size':len(supported),'exploratory_shortlist_size':len(exploratory),'aggregation_max_error':diff,'independent_BH_check':'PASS','matched_strata_count_check':'PASS','empirical_exceedance_check':'PASS','candidate_gate':cfg['supported_candidate_gate'],'scope':'computational prioritization only; phase I and mechanism/safety review pending','private_rat_disease_effects':'NOT_ACCESSED','config_sha256':sha256(ROOT/'config/screening_v1.json'),'complete_results_sha256':sha256(out/'complete_screening_results.tsv'),'supported_shortlist_sha256':sha256(out/'supported_shortlist.tsv'),'exploratory_shortlist_sha256':sha256(out/'exploratory_top10.tsv')}
    write_json(out/'screening_result_freeze.json',summary)
    return summary
if __name__=='__main__':run_standard_module('17_screening_stability',main)
