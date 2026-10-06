"""Packages: numpy/pandas/scipy/numba. Fixed matched-null screening with checkpoint batches.
Inputs: frozen screening arrays/config. Outputs: complete ranked family, null draws/scores.
No external upload; no held-out rat/group effects. Seed and sampling fixed in config.
"""
import json,logging,time
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from numba import set_num_threads
from pipeline_utils import ROOT,path,sha256,write_json,run_standard_module
from screening_core import score_queries,signature_es

def main():
    cfg=json.loads((ROOT/'config/screening_v1.json').read_text());set_num_threads(cfg['threads'])
    out=path('analysis')/'screening';freeze=json.loads((out/'input_freeze.json').read_text())
    if sha256(ROOT/'config/screening_v1.json')!=freeze['screening_config_sha256'] or sha256(out/'screening_inputs.npz')!=freeze['input_npz_sha256']:raise ValueError('Screening inputs changed')
    a=np.load(out/'screening_inputs.npz');args=[a[k] for k in ['orders','weights','cell_order','cell_starts','compound_starts']]
    scores=score_queries(a['queries'],*args);pd.DataFrame(scores,index=a['query_names'],columns=a['compounds']).to_csv(out/'query_compound_scores.tsv',sep='\t')
    primary=a['queries'][list(a['query_names']).index('size_100')]
    rng=np.random.default_rng(cfg['seed']);nullq=np.repeat(primary[None,:],cfg['null_iterations'],axis=0)
    for b in range(len(nullq)):
        for s in np.unique(a['gene_strata']):
            if s<0:continue
            ix=np.flatnonzero(a['gene_strata']==s);nullq[b,ix]=rng.permutation(primary[ix])
    np.savez_compressed(out/'matched_null_queries.npz',queries=nullq)
    checkpoint=out/'null_batches';checkpoint.mkdir(exist_ok=True)
    allnull=[];start=time.perf_counter()
    for begin in range(0,len(nullq),cfg['null_batch_size']):
        end=min(begin+cfg['null_batch_size'],len(nullq));p=checkpoint/f'null_{begin:05d}_{end:05d}.npy';manifest=p.with_suffix('.json')
        key={'input_sha256':freeze['input_npz_sha256'],'config_sha256':freeze['screening_config_sha256'],'begin':begin,'end':end}
        if p.exists():
            record=json.loads(manifest.read_text())
            if any(record[k]!=v for k,v in key.items()) or record['sha256']!=sha256(p):raise ValueError('Null checkpoint mismatch')
            result=np.load(p)
        else:
            result=score_queries(nullq[begin:end],*args);np.save(p,result);write_json(manifest,{**key,'sha256':sha256(p)})
        allnull.append(result)
        if end%500==0 or begin==0:logging.info('MATCHED NULL %d/%d elapsed %.1fs',end,len(nullq),time.perf_counter()-start)
    null=np.concatenate(allnull);np.save(out/'matched_null_compound_scores.npy',null)
    observed=scores[list(a['query_names']).index('size_100')]
    exceed=(null>=observed[None,:]).sum(axis=0)
    p=(1+exceed)/(len(null)+1);p[observed<=0]=1
    q=false_discovery_control(p)
    meta=pd.read_csv(out/'signature_metadata.tsv.gz',sep='\t')
    names=meta.groupby('pert_id').pert_iname.agg(lambda x:' | '.join(sorted(set(x.astype(str)))))
    result=pd.DataFrame({'pert_id':a['compounds'],'pert_iname':names.loc[a['compounds']].values,'reversal_score':observed,'null_exceedances':exceed,'empirical_p':p,'empirical_bh_q':q,'null_mean':null.mean(axis=0),'null_sd':null.std(axis=0,ddof=1)})
    for i,name in enumerate(a['query_names']):result[f'score_{name}']=scores[i]
    # Preserve both ES components and raw/non-opposed connectivity for each signature.
    # npz keys decompress on access; reuse the already loaded arrays in this loop.
    components=np.array([signature_es(primary,args[0][i],args[1][i]) for i in range(len(meta))])
    meta['es_disease_up']=components[:,0];meta['es_disease_down']=components[:,1]
    meta['raw_reversal']=(components[:,1]-components[:,0])/2
    meta['direction']=np.where((components[:,0]<0)&(components[:,1]>0),'reversal',np.where((components[:,0]>0)&(components[:,1]<0),'mimic','non_opposed'))
    meta['reversal_score']=np.where(components[:,0]*components[:,1]<0,meta.raw_reversal,0)
    meta.to_csv(out/'primary_signature_scores.tsv.gz',sep='\t',index=False)
    cell=meta.groupby(['pert_id','cell_id'],sort=True).agg(reversal_score=('reversal_score','median'),n_signatures=('sig_id','size'),fraction_reversing=('direction',lambda x:float((x=='reversal').mean()))).reset_index()
    cell.to_csv(out/'primary_compound_cell_scores.tsv',sep='\t',index=False)
    per=cell.groupby('pert_id').agg(n_cells=('cell_id','size'),positive_cell_fraction=('reversal_score',lambda x:float((x>0).mean())),cell_score_min=('reversal_score','min'),cell_score_max=('reversal_score','max'))
    result=result.merge(per,on='pert_id',validate='one_to_one')
    result['n_signatures']=result.pert_id.map(meta.groupby('pert_id').size())
    result=result.sort_values(['reversal_score','pert_id'],ascending=[False,True]);result.insert(0,'rank',np.arange(1,len(result)+1))
    result.to_csv(out/'primary_compound_ranking.tsv',sep='\t',index=False)
    for field in ['time_h','dose_um']:
        stage=meta.groupby(['pert_id',field,'cell_id']).reversal_score.median().reset_index()
        stage.groupby(['pert_id',field]).agg(reversal_score=('reversal_score','median'),n_cells=('cell_id','nunique')).to_csv(out/f'{field}_stratified_scores.tsv',sep='\t')
    summary={'family_size':len(result),'null_iterations':len(null),'minimum_possible_p':1/(len(null)+1),'positive_compounds':int((result.reversal_score>0).sum()),'positive_q_lt_05':int(((result.reversal_score>0)&(result.empirical_bh_q<.05)).sum()),'ranking_sha256':sha256(out/'primary_compound_ranking.tsv'),'scope':'competitive matched gene-set null, not clinical efficacy; bootstrap assessment pending'}
    write_json(out/'primary_screen_summary.json',summary)
    return summary
if __name__=='__main__':run_standard_module('16_compound_screen',main)
