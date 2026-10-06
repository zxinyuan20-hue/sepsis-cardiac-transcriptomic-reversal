"""Packages numpy, pandas, scipy, requests; R limma/jsonlite.
Input/output sections in config/discovery_scan_sensitivity_v1.json.
Fixed seed, pre-computation amendment and immutable source checks.
Pipeline: freeze -> R model fits -> independent numerical checks -> summaries.
Saved outputs: per-gene effects, designs, concordance and numerical audit; logs.
"""
import json,subprocess
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from pipeline_utils import ROOT,RUNTIME,sha256,now,write_json,run_standard_module,r_environment,path

def main():
 cp=ROOT/'config/discovery_scan_sensitivity_v1.json';cfg=json.loads(cp.read_text(encoding='utf8'));out=ROOT/cfg['output'];out.mkdir(parents=True,exist_ok=True)
 matrix=ROOT/'outputs/preprocessing/human/human_rma_genes.tsv.gz';original=ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv'
 inputs=[cp,matrix,original,ROOT/'outputs/analysis/discovery_exchangeability_v1/cel_scan_metadata.tsv',ROOT/'scripts/88_discovery_scan_sensitivity.R']
 query=[]
 for direction in ['up','down']:
  p=ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_100_{direction}.tsv';inputs.append(p)
  q=pd.read_csv(p,sep='\t',dtype={'entrez_id':str});query.extend(q.entrez_id)
 freeze=out/'amendment_freeze.json'
 if not freeze.exists():write_json(freeze,{'recorded_utc':now(),'config':cfg,'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in inputs]})
 for row in json.loads(freeze.read_text(encoding='utf8'))['files']:assert sha256(ROOT/row['path'])==row['sha256']
 with (path('logs')/'88_R_scan_sensitivity.log').open('a',encoding='utf8') as log:
  subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/88_discovery_scan_sensitivity.R')],env=r_environment(),cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
 read=lambda p:pd.read_csv(p,sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
 base=read(original);rerun=read(out/'primary_reproduction.tsv').loc[base.index]
 reproduction={c:float(np.max(np.abs(base[c]-rerun[c]))) for c in ['logFC','t','P.Value','adj.P.Val']};assert max(reproduction.values())<1e-8
 x=pd.read_csv(matrix,sep='\t',index_col=0);x.index=x.index.astype(str);x=x.loc[base.index]
 rows=[];audit={}
 for label in ['scan_adjusted','shared_days_adjusted']:
  other=read(out/f'{label}.tsv').loc[base.index];design=pd.read_csv(out/f'{label}_design.tsv',sep='\t',index_col=0)
  b=np.linalg.lstsq(design.values,x[design.index].values.T,rcond=None)[0];effect=b[1]-b[0]
  audit[label]=float(np.max(np.abs(effect-other.logFC.values)));assert audit[label]<1e-8
  sig=base['adj.P.Val']<.05;agree=np.sign(base.logFC)==np.sign(other.logFC)
  qbase=base.loc[query];qother=other.loc[query];qagree=np.sign(qbase.logFC)==np.sign(qother.logFC)
  rows.append({'model':label,'all_gene_logFC_spearman':float(spearmanr(base.logFC,other.logFC).statistic),'original_DEGs':int(sig.sum()),'original_DEGs_same_direction':int((sig&agree).sum()),'original_DEGs_same_direction_BH':int((sig&agree&(other['adj.P.Val']<.05)).sum()),'frozen_query_genes':len(query),'frozen_query_same_direction':int(qagree.sum()),'frozen_query_same_direction_BH':int((qagree&(qother['adj.P.Val']<.05)).sum())})
  qout=qother.copy();qout['original_logFC']=qbase.logFC;qout['same_direction']=qagree;qout.to_csv(out/f'{label}_frozen_query.tsv',sep='\t')
 pd.DataFrame(rows).to_csv(out/'sensitivity_concordance.tsv',sep='\t',index=False)
 result={'original_reproduction_max_errors':reproduction,'independent_numpy_OLS_max_errors':audit,'concordance':rows,'new_drug_scores_computed':False,'disease_label_permutation_executed':False}
 write_json(out/'numeric_audit.json',result);print(json.dumps(result,indent=2));return result

if __name__=='__main__':run_standard_module('88_discovery_scan_sensitivity',main)
