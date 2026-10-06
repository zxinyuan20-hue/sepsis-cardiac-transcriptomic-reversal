"""Numba-accelerated weighted ES and fixed hierarchical aggregation.
Packages: numpy,numba. Pure numerical utilities, no filesystem or stochastic side effects.
Queries use 1=disease up, -1=disease down, 0=unselected in original gene order.
"""
import numpy as np
from numba import njit,prange

@njit(cache=True)
def signature_es(query,order,weights):
    n=len(order);nu=0;nd=0;wu=0.;wd=0.
    for j in range(n):
        lab=query[order[j]]
        if lab==1:nu+=1;wu+=weights[j]
        elif lab==-1:nd+=1;wd+=weights[j]
    if nu==0 or nd==0 or nu==n or nd==n:return np.nan,np.nan
    ru=0.;rd=0.;maxu=0.;minu=0.;maxd=0.;mind=0.
    for j in range(n):
        lab=query[order[j]]
        if lab==1:ru+=(weights[j]/wu if wu>0 else 1./nu)
        else:ru-=1./(n-nu)
        if lab==-1:rd+=(weights[j]/wd if wd>0 else 1./nd)
        else:rd-=1./(n-nd)
        if ru>maxu:maxu=ru
        if ru<minu:minu=ru
        if rd>maxd:maxd=rd
        if rd<mind:mind=rd
    eu=maxu if maxu>=-minu else minu
    ed=maxd if maxd>=-mind else mind
    return eu,ed

@njit(cache=True)
def score_one(query,orders,weights):
    scores=np.empty(len(orders),dtype=np.float64)
    for s in range(len(orders)):
        up,down=signature_es(query,orders[s],weights[s])
        scores[s]=(down-up)/2 if up*down<0 else 0.
    return scores

@njit(cache=True)
def aggregate(scores,cell_order,cell_starts,compound_starts):
    cells=np.empty(len(cell_starts)-1,dtype=np.float64)
    for c in range(len(cells)):
        cells[c]=np.median(scores[cell_order[cell_starts[c]:cell_starts[c+1]]])
    compounds=np.empty(len(compound_starts)-1,dtype=np.float64)
    for c in range(len(compounds)):
        compounds[c]=np.median(cells[compound_starts[c]:compound_starts[c+1]])
    return compounds

@njit(cache=True,parallel=True)
def score_queries(queries,orders,weights,cell_order,cell_starts,compound_starts):
    out=np.empty((len(queries),len(compound_starts)-1),dtype=np.float64)
    for b in prange(len(queries)):
        out[b]=aggregate(score_one(queries[b],orders,weights),cell_order,cell_starts,compound_starts)
    return out
