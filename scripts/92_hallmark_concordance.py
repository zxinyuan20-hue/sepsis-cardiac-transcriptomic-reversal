"""Packages: numpy pandas scipy gseapy 1.3.1; inputs/outputs in frozen JSON.
Reproducibility: seed, input hashes, common universe, complete Hallmark family.
Pipeline: freeze amendment BEFORE scores -> signed z ranks -> preranked enrichment
-> all-pathway concordance -> independent ES/BH checks -> saved TSV/JSON/logs.
"""
import json, logging, itertools
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr, false_discovery_control
import gseapy
from pipeline_utils import ROOT, now, sha256, write_json, run_standard_module

def main():
 cfg=json.loads((ROOT/'config/biology_extension_v1.json').read_text(encoding='utf8'))
 out=ROOT/cfg['analysis'];out.mkdir(parents=True,exist_ok=True)
 inputs=[r['path'] for r in cfg['datasets']]+[cfg['gmt']]
 hashes={p:sha256(ROOT/p) for p in inputs}
 freeze=out/'hallmark_protocol_freeze.json'
 if freeze.exists():
  old=json.loads(freeze.read_text());assert old['config']==cfg and old['input_hashes']==hashes
 else:write_json(freeze,{'recorded_utc':now(),'config':cfg,'input_hashes':hashes,'status':'post-result extension; before new enrichment scores'})
 sets={}
 for line in (ROOT/cfg['gmt']).read_text().splitlines():
  f=line.split('\t');sets[f[0]]=list(dict.fromkeys(f[2:]))
 assert len(sets)==50
 data={};coverage=[]
 for r in cfg['datasets']:
  d=pd.read_csv(ROOT/r['path'],sep='\t',dtype={r['id']:'string'})
  d=d.dropna(subset=[r['id'],r['p'],'logFC']).copy()
  d['gene']=d[r['id']].str.replace(r'\.0$','',regex=True)
  assert not d.gene.duplicated().any(),r['name']
  assert d[r['p']].between(0,1).all()
  d['signed_z']=np.sign(d.logFC)*norm.isf(np.maximum(d[r['p']],cfg['p_floor'])/2)
  data[r['name']]=d.set_index('gene')
 common=set.intersection(*(set(d.index) for d in data.values()))
 pd.Series(sorted(common),name='human_entrez').to_csv(out/'common_genes.tsv',sep='\t',index=False)
 results=[];max_error=0
 for scope in ['common','full']:
  for name,d in data.items():
   ranking=d.loc[d.index.isin(common)].copy() if scope=='common' else d.copy()
   ranking['gene_sort']=ranking.index.astype(int)
   ranking=ranking.sort_values(['signed_z','gene_sort'],ascending=[False,True])
   ranking[['logFC','signed_z']].to_csv(out/f'{scope}_{name}_ranks.tsv',sep='\t')
   for term,genes in sets.items():
    k=len(set(genes)&set(ranking.index));coverage.append({'scope':scope,'dataset':name,'pathway':term,'set_size':len(genes),'covered':k,'background':len(ranking),'eligible':cfg['min_size']<=k<=cfg['max_size']})
   logging.info('ENRICHMENT %s %s genes=%s',scope,name,len(ranking))
   fit=gseapy.prerank(rnk=ranking.signed_z,gene_sets=sets,min_size=cfg['min_size'],max_size=cfg['max_size'],seed=cfg['seed'],threads=1,weight=1,ascending=None,method='multilevel',sample_size=cfg['sample_size'],eps=cfg['eps'],no_plot=True,outdir=None)
   for term,v in fit.results.items():
    # Independent direct running-sum implementation, no library score reuse.
    hit=ranking.index.isin(sets[term]);w=abs(ranking.signed_z.to_numpy());inc=np.where(hit,w/(w[hit].sum()),-1/(len(hit)-hit.sum()));rs=np.cumsum(inc)
    es=rs.max() if abs(rs.max())>abs(rs.min()) else rs.min()
    max_error=max(max_error,abs(es-v['es']))
    results.append({'scope':scope,'dataset':name,'pathway':term,'ES':v['es'],'NES':v['nes'],'P':v['pval'],'log2err':v.get('log2err'),'leading_genes':v['lead_genes']})
 res=pd.DataFrame(coverage).merge(pd.DataFrame(results),on=['scope','dataset','pathway'],how='left',validate='one_to_one')
 for scope,g in res.groupby('scope'):
  res.loc[g.index,'BH_300']=false_discovery_control(g.P.fillna(1))
  for name,h in g.groupby('dataset'):res.loc[h.index,'BH_50']=false_discovery_control(h.P.fillna(1))
 res.loc[~res.eligible,['BH_300','BH_50']]=np.nan
 res.to_csv(out/'hallmark_enrichment.tsv',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(out/'hallmark_coverage.tsv',sep='\t',index=False)
 pairs=[]
 for scope in ['common','full']:
  m=res[res.scope==scope].pivot(index='pathway',columns='dataset',values='NES').reindex(sorted(sets))
  m.to_csv(out/f'{scope}_NES_matrix.tsv',sep='\t')
  for a,b in itertools.combinations(data,2):
   x=data[a];y=data[b];genes=sorted(common if scope=='common' else set(x.index)&set(y.index));p=m[[a,b]].dropna()
   pairs.append({'scope':scope,'dataset_a':a,'dataset_b':b,'genes':len(genes),'pathways':len(p),'gene_logFC_rho':spearmanr(x.loc[genes,'logFC'],y.loc[genes,'logFC']).statistic,'gene_signed_z_rho':spearmanr(x.loc[genes,'signed_z'],y.loc[genes,'signed_z']).statistic,'pathway_NES_rho':spearmanr(p[a],p[b]).statistic,'pathway_direction_agreement':float(np.mean(np.sign(p[a])==np.sign(p[b])))})
 pd.DataFrame(pairs).to_csv(out/'cross_model_concordance.tsv',sep='\t',index=False)
 assert max_error<1e-10,max_error
 assert len(res)==600 and res.eligible.sum()==len(results)
 write_json(out/'hallmark_audit.json',{'common_genes':len(common),'all_50_sets_each_dataset_both_scopes':len(res)==600,'eligible_tests':int(res.eligible.sum()),'exclusions':res.loc[~res.eligible,['scope','dataset','pathway','covered']].to_dict('records'),'independent_ES_max_error':max_error,'gseapy_version':gseapy.__version__,'seed':cfg['seed'],'inference':'conditional competitive gene-set null; not independent biological replication; BH across 300 planned tests separately by scope, unavailable tests assigned P=1 for denominator only and reported NA; correlations descriptive, overlapping pathways dependent','input_hashes':hashes})
 return {'common_genes':len(common),'tests_primary':300,'ES_check':max_error}

if __name__=='__main__':run_standard_module('92_hallmark_concordance',main)
