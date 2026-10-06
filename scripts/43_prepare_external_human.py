"""Packages requests/pandas/numpy. Input frozen plan and official GEO resources.
Outputs immutable acquired files, matrix/schema QC and input freeze under outputs/.
No external differential result selection; group-blind mapping and fixed seed/log.
"""
import json,gzip,io,re,logging,itertools
import pandas as pd,numpy as np
from pipeline_utils import ROOT,acquire,sha256,write_json,run_standard_module,now
OUT=ROOT/'outputs/analysis/external_human';SRC=OUT/'sources'
def main():
 SRC.mkdir(parents=True,exist_ok=True);cfg=json.loads((ROOT/'config/external_human_v1.json').read_text(encoding='utf-8'))
 plan=OUT/'analysis_plan_freeze.json'
 if not plan.exists():write_json(plan,{'time_utc':now(),'config_sha256':sha256(ROOT/'config/external_human_v1.json'),'config':cfg})
 else:assert json.loads(plan.read_text(encoding='utf-8'))['config_sha256']==sha256(ROOT/'config/external_human_v1.json')
 if (OUT/'input_freeze.json').exists():
  for f in json.loads((OUT/'input_freeze.json').read_text(encoding='utf-8'))['files']:assert sha256(ROOT/f['path'])==f['sha256']
  return {'reused_frozen_preparation':True,**json.loads((OUT/'preparation_summary.json').read_text(encoding='utf-8'))}
 errors=[]
 # The platform family archive includes all linked samples (7.63GB), not just
 # annotation. Use the official platform-only text response instead.
 pp=SRC/'platform_self_full.txt'
 if not pp.exists():
  acquire('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17586&targ=self&form=text&view=full',pp,'Official GEO platform-only annotation',max_gb=.3)
 else:
  meta=pp.with_suffix(pp.suffix+'.source.json');assert sha256(pp)==json.loads(meta.read_text(encoding='utf-8'))['sha256']
 for i in range(1,8):
  fn=f'GSE237861_edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis{i}.heart.tsv.gz'
  try:acquire('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237861/suppl/'+fn,SRC/fn,'Author individual heart contrast (not sample matrix)',max_gb=.05)
  except Exception as e:errors.append({'source':fn,'error':str(e)})
 write_json(OUT/'acquisition_errors.json',errors)
 text=pp.read_text(encoding='utf-8');part=text.split('!platform_table_begin')[1].split('!platform_table_end')[0].strip()
 annot=pd.read_csv(io.StringIO(part),sep='\t',dtype=str,usecols=['ID','gene_assignment']);annot.to_csv(OUT/'platform_annotation.tsv.gz',sep='\t',index=False)
 del text,part
 # Parse and audit identifiers without using group differences for selection.
 mapping=[]
 for r in annot.itertuples():
  ids=set();symbols=set()
  for block in str(r.gene_assignment).split(' /// '):
   fields=[z.strip() for z in block.split(' // ')]
   if len(fields)>=5 and fields[4].isdigit():ids.add(fields[4]);symbols.add(fields[1])
  mapping.append({'cluster':r.ID,'entrez_id':next(iter(ids)) if len(ids)==1 else None,'n_entrez':len(ids),'symbols':';'.join(sorted(symbols))})
 mp=pd.DataFrame(mapping);mp.to_csv(OUT/'cluster_mapping_audit.tsv',sep='\t',index=False)
 src=ROOT/'outputs/analysis/independent_evidence/sources/GSE141864_samples.soft';meta=ROOT/'outputs/analysis/independent_evidence/independent_sample_metadata.tsv'
 sm=pd.read_csv(meta,sep='\t',dtype=str);sm=sm[sm.accession.eq('GSE141864')&sm.tissue.eq('Heart')].copy()
 sm['primary']=sm['type of storage tissue'].eq('FFPE');sm['group']=np.where(sm.diagnosis.eq('Meningococcal septic shock (MSS)'),'sepsis','control');sm['patient_id']=sm['patients id nr']
 assert sm[sm.primary].patient_id.nunique()==7 and sm[sm.primary].group.eq('sepsis').sum()==5
 values={};gsm=None;active=False;lines=[]
 for line in src.open(encoding='utf-8'):
  if line.startswith('^SAMPLE = '):gsm=line.strip().split(' = ')[1]
  elif line.startswith('!sample_table_begin'):active=True;lines=[]
  elif line.startswith('!sample_table_end'):
   if gsm in set(sm.gsm):
    t=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype={'ID_REF':str});assert t.ID_REF.is_unique;values[gsm]=t.set_index('ID_REF').VALUE
   active=False
  elif active:lines.append(line)
 x=pd.DataFrame(values);assert len(values)==10 and np.isfinite(x.to_numpy()).all()
 sm=sm.set_index('gsm').loc[x.columns].rename_axis('gsm').reset_index();sm.to_csv(OUT/'GSE141864_heart_design.tsv',sep='\t',index=False)
 primary=sm.loc[sm.primary,'gsm'].tolist();v=x[primary];assert v.to_numpy().max()<30 and v.to_numpy().min()>-20,'Unexpected scale: stop before inference'
 usable=mp[(mp.n_entrez==1)&mp.cluster.isin(x.index)].copy();usable['mean_primary']=usable.cluster.map(v.mean(axis=1))
 selected=usable.sort_values(['entrez_id','mean_primary','cluster'],ascending=[True,False,True]).drop_duplicates('entrez_id');selected.to_csv(OUT/'selected_cluster_per_gene.tsv',sep='\t',index=False)
 genes=x.loc[selected.cluster].copy();genes.index=selected.entrez_id;genes.index.name='entrez_id';genes.to_csv(OUT/'GSE141864_gene_expression.tsv.gz',sep='\t')
 alt=x.loc[usable.cluster].copy();alt.index=usable.entrez_id;alt=alt.groupby(level=0).median();alt.to_csv(OUT/'GSE141864_median_gene_expression.tsv.gz',sep='\t')
 qc=pd.DataFrame({'gsm':x.columns,'q01':x.quantile(.01).values,'q50':x.median().values,'q99':x.quantile(.99).values,'missing':x.isna().sum().values});qc.to_csv(OUT/'GSE141864_sample_qc.tsv',sep='\t',index=False)
 x.corr().to_csv(OUT/'GSE141864_sample_correlations.tsv',sep='\t')
 common=set(pd.read_csv(ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id)&set(genes.index)
 coverage=[]
 for size in [50,100]:
  for direction in ['up','down']:
   q=pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_{size}_{direction}.tsv',sep='\t',dtype={'entrez_id':str});n=q.entrez_id.isin(common).sum();coverage.append({'query':size,'direction':direction,'original':len(q),'available':int(n),'fraction':float(n/len(q))})
 assert min(r['fraction'] for r in coverage if r['query']==100)>=cfg['query_coverage_min_fraction'];pd.DataFrame(coverage).to_csv(OUT/'frozen_query_coverage.tsv',sep='\t',index=False)
 ncbi=pd.read_csv(SRC/'counts_api.tsv.gz',sep='\t',nrows=0);heart237=pd.read_csv(meta,sep='\t');heart237=heart237[heart237.accession.eq('GSE237861')&heart237.tissue.eq('heart')][['gsm','title']].copy();heart237['NCBI_count_available']=heart237.gsm.isin(ncbi.columns);heart237.to_csv(OUT/'GSE237861_count_coverage.tsv',sep='\t',index=False)
 ready={'GSE141864_primary_samples':7,'GSE141864_primary_cases':5,'GSE141864_primary_controls':2,'heart_records':10,'unique_patients':7,'mapped_genes':len(genes),'common_discovery_genes':len(common),'query_coverage':coverage,'GSE237861_expected_heart':14,'GSE237861_NCBI_heart_available':int(heart237.NCBI_count_available.sum()),'GSE237861_full_matrix_ready':False,'GSE237861_formal_test':'NOT_RUN_INCOMPLETE_MATRIX','scale':'submitted_RMA_log2','raw_CEL_reprocessed':False}
 write_json(OUT/'preparation_summary.json',ready)
 frozen=[ROOT/'config/external_human_v1.json',src,pp,meta,OUT/'GSE141864_gene_expression.tsv.gz',OUT/'GSE141864_median_gene_expression.tsv.gz',OUT/'GSE141864_heart_design.tsv']+list((ROOT/'outputs/analysis/human/frozen_signatures').glob('*.tsv'))
 write_json(OUT/'input_freeze.json',{'utc':now(),'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in frozen]})
 write_json(OUT/'schema_inventory.json',{'platform_columns':annot.columns.tolist(),'platform_rows':len(annot),'author_contrast_columns':{p.name:pd.read_csv(p,sep='\t',nrows=0).columns.tolist() for p in SRC.glob('*individual*.gz')}})
 return {**ready,'annotation_rows':len(annot),'author_heart_contrasts':len(list(SRC.glob('*individual*.gz'))),'errors':errors,'no_replication_scores_computed':True}
if __name__=='__main__':run_standard_module('43_prepare_external_human',main)
