"""Packages: numpy/pandas/h5py/numba/scipy. Inputs: frozen human signatures and LINCS subset.
Outputs: outputs/analysis/screening prepared arrays, query definitions, algorithm tests.
Fixed seed/config; checksum/axis validation before any compound ranking.
"""
import json,time
import h5py,numpy as np,pandas as pd
from numba import set_num_threads
from pipeline_utils import ROOT,CONFIG,path,sha256,write_json,run_standard_module
from screening_core import signature_es,score_queries,aggregate

def reference_es(z,query,label):
    order=np.argsort(-z,kind='stable');hit=query[order]==label;w=np.abs(z[order]);total=w[hit].sum()
    increments=np.where(hit,w/total if total else 1/hit.sum(),-1/(len(z)-hit.sum()))
    walk=np.cumsum(increments);hi=max(0,walk.max());lo=min(0,walk.min())
    return hi if hi>=-lo else lo

def main():
    cfg=json.loads((ROOT/'config/screening_v1.json').read_text());set_num_threads(cfg['threads'])
    out=path('analysis')/'screening';out.mkdir(exist_ok=True)
    freeze=json.loads((path('analysis')/'human/frozen_signatures/signature_freeze_manifest.json').read_text())
    if freeze['study_config_sha256']!=sha256(ROOT/'config/study_config.json'):raise ValueError('Changed human config')
    for item in freeze['files']:
        if sha256(ROOT/item['path'])!=item['sha256']:raise ValueError('Changed human query')
    inp=path('preparation')/'lincs/eligible_landmarks.h5'
    audit=json.loads((path('audit')/'lincs_matrix_qc.json').read_text())
    if sha256(inp)!=audit['subset_sha256']:raise ValueError('LINCS subset changed')
    with h5py.File(inp) as h:
        genes=np.array([x.decode() for x in h['gene_id'][:]]);ids=[x.decode() for x in h['sig_id'][:]];z=h['zscore'][:]
    # Lexical gene order makes ties deterministic independent of archive ordering.
    gi=np.argsort(genes,kind='stable');genes=genes[gi];z=z[:,gi].astype(np.float64)
    meta=pd.read_csv(path('preparation')/'lincs_eligible_signatures.tsv.gz',sep='\t').set_index('sig_id').loc[ids].reset_index()
    coverage=meta.groupby('pert_id').agg(n_signatures=('sig_id','size'),n_cells=('cell_id','nunique'))
    compounds=coverage[(coverage.n_signatures>=CONFIG['lincs']['minimum_high_quality_signatures'])&(coverage.n_cells>=CONFIG['lincs']['minimum_distinct_cell_lines'])].index.tolist()
    select=meta.pert_id.isin(compounds).values;meta=meta.loc[select].reset_index(drop=True);z=z[select]
    cell_order=[];cell_starts=[0];compound_starts=[0];cell_records=[]
    for pert in compounds:
        for cell in sorted(meta.loc[meta.pert_id==pert,'cell_id'].unique()):
            index=np.flatnonzero((meta.pert_id==pert)&(meta.cell_id==cell))
            cell_order.extend(index);cell_starts.append(len(cell_order));cell_records.append({'pert_id':pert,'cell_id':cell,'n_signatures':len(index)})
        compound_starts.append(len(cell_records))
    orders=np.argsort(-z,axis=1,kind='stable').astype(np.int32)
    weights=np.take_along_axis(np.abs(z),orders,axis=1)
    names=[];queries=[]
    for size in [50,100,150]:
        q=np.zeros(len(genes),dtype=np.int8)
        for direction,label in [('up',1),('down',-1)]:
            table=pd.read_csv(path('analysis')/f'human/frozen_signatures/landmarks_{size}_{direction}.tsv',sep='\t',dtype={'entrez_id':str})
            q[np.isin(genes,table.entrez_id)]=label
        names.append(f'size_{size}');queries.append(q)
    for label in ['without_qc_flag','submitted']:
        table=pd.read_csv(path('analysis')/f'human/{label}_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str})
        table=table[table.entrez_id.isin(genes)&(table['adj.P.Val']<CONFIG['analysis']['fdr'])].copy();table['abs_t']=table.t.abs()
        q=np.zeros(len(genes),dtype=np.int8)
        for sign in [1,-1]:
            selected=table[table.logFC*sign>0].sort_values(['abs_t','entrez_id'],ascending=[False,True]).head(100)
            q[np.isin(genes,selected.entrez_id)]=sign
        names.append(label);queries.append(q)
    queries=np.array(queries)
    if any(min((q==1).sum(),(q==-1).sum())<15 for q in queries):raise ValueError('Query coverage insufficient')
    # Independent reference implementation, sign, and hierarchical aggregation tests.
    rng=np.random.default_rng(cfg['seed']);errors=[]
    for _ in range(40):
        zz=rng.normal(size=37);q=np.zeros(37,np.int8);q[rng.permutation(37)[:7]]=1;q[rng.choice(np.flatnonzero(q==0),5,replace=False)]=-1
        o=np.argsort(-zz,kind='stable');u,d=signature_es(q,o,np.abs(zz[o]))
        errors.extend([abs(u-reference_es(zz,q,1)),abs(d-reference_es(zz,q,-1))])
    zz=np.arange(20,dtype=float)-10;q=np.zeros(20,np.int8);q[:5]=1;q[-5:]=-1;o=np.argsort(-zz)
    eu,ed=signature_es(q,o,np.abs(zz[o]));assert eu<0 and ed>0
    eu2,ed2=signature_es(-q,o,np.abs(zz[o]));assert eu2>0 and ed2<0
    assert np.allclose(aggregate(np.array([1.,1.,1.,-1.,.4]),np.arange(5),np.array([0,3,4,5]),np.array([0,2,3])),[0.,.4])
    if max(errors)>1e-12:raise ValueError('ES reference mismatch')
    write_json(path('audit')/'screening_algorithm_tests.json',{'reference_tests':40,'max_abs_error':max(errors),'reversal_and_mimic_direction':'PASS','equal_cell_weighting':'PASS'})
    human=pd.read_csv(path('analysis')/'human/primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
    expr=pd.read_csv(ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz',sep='\t',index_col=0);expr.index=expr.index.astype(str)
    shared=[g for g in genes if g in human.index]
    means=expr.loc[shared].mean(axis=1);sd=expr.loc[shared].std(axis=1)
    # Rank-based equal-count bins are deterministic in lexical ID order.
    strata=pd.DataFrame({'gene_id':shared,'mean':means.values,'sd':sd.values})
    strata['mean_bin']=pd.qcut(strata['mean'].rank(method='first'),cfg['null_strata_mean_bins'],labels=False)
    strata['sd_bin']=pd.qcut(strata['sd'].rank(method='first'),cfg['null_strata_sd_bins'],labels=False)
    strata['stratum']=strata.mean_bin*cfg['null_strata_sd_bins']+strata.sd_bin
    strata.to_csv(out/'null_gene_strata.tsv',sep='\t',index=False)
    gene_strata=np.full(len(genes),-1,np.int32)
    for row in strata.itertuples():gene_strata[np.flatnonzero(genes==row.gene_id)[0]]=row.stratum
    np.savez_compressed(out/'screening_inputs.npz',genes=genes,orders=orders,weights=weights,z=z,queries=queries,query_names=np.array(names),cell_order=np.array(cell_order,np.int32),cell_starts=np.array(cell_starts,np.int32),compound_starts=np.array(compound_starts,np.int32),compounds=np.array(compounds),gene_strata=gene_strata)
    pd.DataFrame({'gene_id':genes}).to_csv(out/'gene_order.tsv',sep='\t',index=False)
    meta.to_csv(out/'signature_metadata.tsv.gz',sep='\t',index=False)
    pd.DataFrame(cell_records).to_csv(out/'compound_cell_order.tsv',sep='\t',index=False)
    coverage.to_csv(out/'all_compound_coverage.tsv',sep='\t')
    details=[{'query':n,'up':int((q==1).sum()),'down':int((q==-1).sum())} for n,q in zip(names,queries)]
    write_json(out/'input_freeze.json',{'screening_config_sha256':sha256(ROOT/'config/screening_v1.json'),'human_freeze_sha256':sha256(path('analysis')/'human/frozen_signatures/signature_freeze_manifest.json'),'input_npz_sha256':sha256(out/'screening_inputs.npz'),'signature_count':len(meta),'compounds':len(compounds),'queries':details})
    # Time a small algorithm pilot; no candidate names or scores used for decisions.
    t=time.perf_counter();_ = score_queries(np.repeat(queries[1:2],4,axis=0),orders,weights,np.array(cell_order),np.array(cell_starts),np.array(compound_starts));elapsed=time.perf_counter()-t
    return {'signatures':len(meta),'compounds':len(compounds),'queries':details,'compiled_pilot_seconds':elapsed}
if __name__=='__main__':run_standard_module('14_prepare_screening',main)
