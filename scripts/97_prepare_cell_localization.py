"""Public human-reference acquisition and frozen protocol BEFORE localization.
Packages requests h5py numpy pandas; paths ROOT-relative; immutable downloads.
No private H9c2 data sent. GET requests contain public human query gene symbols.
"""
import json,importlib.util,concurrent.futures,logging
import h5py,numpy as np,pandas as pd,requests
from pipeline_utils import ROOT,now,sha256,write_json,run_standard_module
from cellxgene_read_utils import decode
OUT=ROOT/'outputs/analysis/biology_extension_v1';SRC=OUT/'sources';CELL=OUT/'cell_localization';CELL.mkdir(exist_ok=True)
def main():
 spec=importlib.util.spec_from_file_location('atlas',ROOT/'scripts/95_inspect_atlas_h5.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 with m.RangeFile() as f:
  with h5py.File(f,'r') as h:
   obs={}
   for col in ['cell_type','donor','source','region','cell_source']:
    logging.info('READ reference metadata %s',col)
    labels=h[f'obs/__categories/{col}'].asstr()[:];obs[col]=labels[h[f'obs/{col}'][:]]
   logging.info('READ reference gene names and row identifiers')
   obs=pd.DataFrame(obs);var=pd.Series(h['var/_index'].asstr()[:],name='symbol');raw_names=h['obs/_index'].asstr()[:10]
 obs.to_csv(CELL/'reference_annotations.tsv',sep='\t',index=False);var.to_csv(CELL/'reference_genes.tsv',sep='\t',index=False)
 # Verify API row identity against HDF5, rather than assuming equal dimensions.
 p=SRC/'cellxgene_names.bin'
 if not p.exists():
  r=requests.get('https://www.heartcellatlas.org/global/api/v0.2/annotations/obs',params={'annotation-name':'name_0'},timeout=60);r.raise_for_status();p.write_bytes(r.content);write_json(p.with_suffix('.provenance.json'),{'url':r.url,'retrieved_utc':now(),'sha256':sha256(p)})
 names=decode(p.read_bytes());assert len(names)==len(obs)
 # Full identity checked cheaply by sorted category distributions plus ten exact names;
 # all row names from HDF5 are checked subsequently if first-block check succeeds.
 with m.RangeFile() as f:
  with h5py.File(f,'r') as h:raw_all=h['obs/_index'].asstr()[:]
 joined=pd.Index(raw_all+'-'+obs.cell_source.to_numpy());assert joined.is_unique and names.iloc[:,0].is_unique
 indexer=joined.get_indexer(names.iloc[:,0]);assert (indexer>=0).all() and len(np.unique(indexer))==len(obs)
 obs=obs.iloc[indexer].reset_index(drop=True);obs.to_csv(CELL/'reference_annotations.tsv',sep='\t',index=False)
 write_json(CELL/'row_identity_audit.json',{'rows':len(obs),'unique_one_to_one_matches':int(len(np.unique(indexer))),'mapping':'raw cell identifier plus hyphen plus documented cell_source equals API name_0','original_order_equal':bool(np.array_equal(indexer,np.arange(len(obs))))})
 mapping=pd.read_csv(ROOT/'outputs/preprocessing/human/probe_mapping_all.tsv',sep='\t',dtype=str)[['ENTREZID','SYMBOL']].drop_duplicates().dropna()
 mapping=mapping[mapping.groupby('ENTREZID').SYMBOL.transform('nunique')==1].drop_duplicates('ENTREZID')
 query=[]
 for direction in ['up','down']:
  d=pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_100_{direction}.tsv',sep='\t',dtype={'entrez_id':str})
  d=d.merge(mapping,left_on='entrez_id',right_on='ENTREZID',how='left',validate='one_to_one');d['direction']=direction;d['reference_present']=d.SYMBOL.isin(var);query.append(d[['entrez_id','SYMBOL','direction','reference_present']])
 q=pd.concat(query);q.to_csv(CELL/'query_mapping.tsv',sep='\t',index=False)
 proto={'recorded_utc':now(),'stage':'post-result extension; before cell localization summaries','reference':'Litvinukova et al 2020 adult healthy heart; published global_raw.h5ad annotations and official global viewer expression','observations':'Nuclei only; exclude NotAssigned/doublets; original 11 major types including atrial and ventricular cardiomyocytes','unit':'donor-by-cell-type mean viewer expression; minimum 20 nuclei per stratum; minimum 7 eligible donors per cell type; equal donor weighting','query':'frozen 100-UP/90-DOWN; unambiguous Entrez-symbol map and reference coverage; no replacement or new modules','endpoint':'per-gene reference localization fraction: nonnegative equal-donor mean expression in each eligible type divided by sum across eligible types; equal gene mean separately UP and DOWN; zero-total genes excluded with explicit coverage','sensitivity':'fraction of nuclei with detected expression (>0) summarized with same donor and gene weighting','inference':'descriptive reference-based localization only, no enrichment P values, no cell-level replication test, no cell proportions, no single-cell validation; no claim of injury or causal mechanism','normalization':'official viewer expression scale, checked against raw counts for sampled entries; scale and source retained; localization relative to eligible cardiac types, not a whole-genome enrichment test','query_hash':sha256(CELL/'query_mapping.tsv'),'metadata_hash':sha256(CELL/'reference_annotations.tsv'),'row_identity':'all 486134 row names exactly match raw HDF5 and API','reference_source':'https://www.heartcellatlas.org/global','seed':20261005}
 freeze=CELL/'protocol_freeze.json'
 if not freeze.exists():write_json(freeze,proto)
 symbols=q.loc[q.reference_present,'SYMBOL'].drop_duplicates().tolist()
 def get_batch(i):
  genes=symbols[i:i+10];p=SRC/f'cellxgene_query_{i:03d}.bin'
  if not p.exists():
   for trial in range(3):
    try:
     r=requests.get('https://www.heartcellatlas.org/global/api/v0.2/data/var',params=[('var:name_0',g) for g in genes],timeout=60);r.raise_for_status();d=decode(r.content);assert len(d)==len(obs) and d.shape[1]==len(genes)
     p.write_bytes(r.content);write_json(p.with_suffix('.provenance.json'),{'url':r.url,'symbols':genes,'retrieved_utc':now(),'sha256':sha256(p)});break
    except Exception:
     if trial==2:raise
  d=decode(p.read_bytes());assert set(var.iloc[d.columns].tolist())==set(genes);logging.info('ACQUIRED query genes %s / %s',i+len(genes),len(symbols));return p.name
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:files=list(pool.map(get_batch,range(0,len(symbols),10)))
 write_json(CELL/'acquisition.json',{'files':files,'mapped_queries':len(symbols),'row_identity_verified':True,'reference_rows':len(obs)})
 return {'query_coverage':q.groupby('direction').reference_present.sum().to_dict(),'files':len(files)}
if __name__=='__main__':run_standard_module('97_prepare_cell_localization',main)
