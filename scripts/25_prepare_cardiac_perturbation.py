"""Packages pandas/numpy. Public GSE217421 processed data acquisition/input-only audit.
Inputs GEO counts/design/config and sample metadata; outputs cardiac resource QC/subset.
No differential drug effect or independence claim; preserve donor/experiment/dish/well.
"""
import gzip,json
import numpy as np,pandas as pd
from pipeline_utils import ROOT,path,acquire,write_json,run_standard_module
def main():
 out=path('analysis')/'mechanism';folder=path('preparation')/'GSE217421';folder.mkdir(exist_ok=True)
 raw=path('raw')/'geo/GSE217421';base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217421/suppl/'
 for suffix in ['Experiment-Design','Configs','Read-Counts']:
  name=f'GSE217421_Conv-RNAseq-{suffix}.tsv.gz';acquire(base+name,raw/name,'GEO human cardiomyocyte drug perturbation processed resource')
 design=pd.read_csv(raw/'GSE217421_Conv-RNAseq-Experiment-Design.tsv.gz',sep='\t',dtype=str)
 parts=['State','Cell','Experiment','Dish','Plate','Well'];design['column_id']=design[parts].agg('.'.join,axis=1)
 if design.column_id.duplicated().any():raise ValueError('Duplicate design compound key')
 cp=raw/'GSE217421_Conv-RNAseq-Read-Counts.tsv.gz'
 with gzip.open(cp,'rt',encoding='utf-8') as f:header=f.readline().rstrip('\n\r').split('\t')
 if len(header)!=len(set(header)) or set(header)!=set(design.column_id):raise ValueError('Count/design column disagreement')
 design=design.set_index('column_id').loc[header].reset_index();design.to_csv(folder/'full_design.tsv',sep='\t',index=False)
 selected=design[design.State.isin(['VEM','CTRL'])].copy();keep=selected.column_id.tolist()
 rows=0;nonfinite=0;negative=0;fractional=0;genes=[];subsets=[]
 for chunk in pd.read_csv(cp,sep='\t',names=['gene_symbol']+header,skiprows=1,chunksize=1000):
  data=chunk[header].to_numpy(dtype=float);rows+=len(chunk);nonfinite+=int((~np.isfinite(data)).sum());negative+=int((data<0).sum());fractional+=int((data!=np.floor(data)).sum());genes.extend(chunk.gene_symbol.astype(str));subsets.append(chunk[['gene_symbol']+keep])
 if nonfinite or negative or len(genes)!=len(set(genes)):raise ValueError('Count value/ID QC failed')
 pd.concat(subsets,ignore_index=True).to_csv(folder/'VEM_CTRL_counts.tsv.gz',sep='\t',index=False)
 selected.to_csv(folder/'VEM_CTRL_design.tsv',sep='\t',index=False)
 cov=selected.groupby(['Cell','State']).size().unstack(fill_value=0);cov.to_csv(folder/'donor_condition_coverage.tsv',sep='\t')
 metadata=pd.read_csv(out/'cardiac_resource_sample_metadata.tsv',sep='\t').fillna('');metadata=metadata[metadata.accession=='GSE217421']
 result={'accession':'GSE217421','genes':rows,'count_columns':len(header),'geo_samples':len(metadata),'design_columns_matched':True,'nonfinite':nonfinite,'negative':negative,'fractional':fractional,'VEM_samples':int((selected.State=='VEM').sum()),'CTRL_samples':int((selected.State=='CTRL').sum()),'cell_lines':selected.Cell.nunique(),'donor_condition_counts':cov.reset_index().to_dict('records'),'readiness':'INPUT_QC_PASSED; donor/experiment biological-repeat contract still needed before drug-effect inference','drug_context':'Vemurafenib 2uM, 48h in GEO sample metadata; not same exposure as LINCS 6/24h','scope':'no drug differential analysis performed'}
 write_json(out/'cardiac_perturbation_input_qc.json',result)
 return result
if __name__=='__main__':run_standard_module('25_prepare_cardiac_perturbation',main)
