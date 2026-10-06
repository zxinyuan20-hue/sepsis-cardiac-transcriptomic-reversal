"""Packages stdlib/requests/pandas. Inspect public human cardiac drug perturbation metadata.
Input accession leads explicitly found in verified article full text; outputs raw GEO and audit.
No private data upload; no independent replication claim before matrices/design are evaluated.
"""
import gzip,json,re
import pandas as pd
from pipeline_utils import ROOT,path,acquire,write_json,run_standard_module
def main():
 out=path('analysis')/'mechanism';lead=json.loads((out/'cardiac_perturbation_resource_leads.json').read_text(encoding='utf-8'))
 accessions=sorted(set(g for r in lead for g in r['accessions'] if g.startswith('GSE')))
 summaries=[];samples=[];errors=[]
 for accession in accessions:
  try:
   series=accession[:-3]+'nnn';dest=path('raw')/'geo'/accession/(accession+'_family.soft.gz')
   acquire(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{series}/{accession}/soft/{accession}_family.soft.gz',dest,'NCBI GEO cardiac drug resource metadata')
   current=None;meta={};rows=[]
   with gzip.open(dest,'rt',encoding='utf-8',errors='replace') as f:
    for line in f:
     if line.startswith('^SAMPLE = '):
      if current:rows.append(current)
      current={'accession':accession,'gsm':line.strip().split(' = ',1)[1]}
     elif line.startswith('^'):
      if current:rows.append(current);current=None
     elif line.startswith('!Series_') and ' = ' in line:
      k,v=line.strip().split(' = ',1);meta.setdefault(k,[]).append(v)
     elif current is not None and line.startswith('!Sample_') and ' = ' in line:
      k,v=line.strip().split(' = ',1)
      if any(t in k for t in ['title','source_name','organism','characteristics','treatment_protocol','data_processing','supplementary_file']):current[k]=current.get(k,'')+' | '+v
   if current:rows.append(current)
   samples.extend(rows)
   summaries.append({'accession':accession,'samples':len(rows),'title':meta.get('!Series_title'), 'summary':meta.get('!Series_summary'),'overall_design':meta.get('!Series_overall_design'),'supplementary_files':meta.get('!Series_supplementary_file'),'pubmed_ids':meta.get('!Series_pubmed_id')})
  except Exception as e:errors.append({'accession':accession,'error':str(e)})
 pd.DataFrame(samples).to_csv(out/'cardiac_resource_sample_metadata.tsv',sep='\t',index=False)
 write_json(out/'cardiac_resource_audit.json',{'datasets':summaries,'errors':errors,'scope':'metadata only; not a completed independent compound replication'})
 return {'datasets':len(summaries),'sample_records':len(samples),'errors':errors}
if __name__=='__main__':run_standard_module('24_cardiac_resource_audit',main)
