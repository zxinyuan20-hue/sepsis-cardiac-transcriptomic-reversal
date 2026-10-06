"""Packages pandas and standard library, shared requests logging.
Input/output: central discovery_exchangeability_v1 config; CEL manifest.
Reproducibility: verify each original CEL hash, read header only; no normalization.
Pipeline: parse Calvin generic header recursively, extract documented scan fields,
cross-tabulate scan date and disease group. Save audit JSON and TSV under outputs.
Scan dates are technical metadata, not automatically complete processing batches.
"""
import gzip,struct,json
import pandas as pd
from pipeline_utils import ROOT,sha256,write_json,run_standard_module

def parse_header(f,depth=0):
 assert depth<10
 def integer():return struct.unpack('>i',f.read(4))[0]
 def string(wide=False):
  n=integer();assert 0<=n<100000
  return f.read(n*(2 if wide else 1)).decode('utf-16-be' if wide else 'ascii',errors='replace')
 h={'type':string(),'id':string(),'created':string(True),'locale':string(True),'parameters':{}}
 n=integer();assert 0<=n<10000
 for _ in range(n):
  key=string(True);nbytes=integer();assert 0<=nbytes<100000
  blob=f.read(nbytes);mime=string(True)
  value=blob.decode('utf-16-be',errors='replace').rstrip('\0') if mime=='text/plain' else blob.decode('ascii',errors='replace').rstrip('\0') if mime=='text/ascii' else None
  if value is not None:h['parameters'][key]=value
 n=integer();assert 0<=n<100
 h['parents']=[parse_header(f,depth+1) for _ in range(n)]
 return h

def main():
 cfg=json.loads((ROOT/'config/discovery_exchangeability_v1.json').read_text(encoding='utf8'));out=ROOT/cfg['output']
 meta=pd.read_csv(ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv',sep='\t');rows=[];hierarchies=[]
 for r in meta.itertuples():
  p=ROOT/r.cel_path;assert sha256(p)==r.sha256
  with gzip.open(p,'rb') as f:
   assert f.read(2)==b'\x3b\x01';f.read(8);h=parse_header(f)
  hierarchies.append({'gsm':r.gsm,'header':h})
  def flatten(h):
   d=list(h['parameters'].items())
   for parent in h['parents']:d+=flatten(parent)
   return d
  pairs=flatten(h);row={'gsm':r.gsm,'group':r.group}
  for key in ['affymetrix-scan-date','affymetrix-scanner-type','affymetrix-scanner-serialnumber','affymetrix-array-type']:
   values=sorted(set(v for k,v in pairs if k==key));assert len(values)<=1,(r.gsm,key,values)
   row[key]=values[0] if values else ''
  rows.append(row)
 frame=pd.DataFrame(rows);frame['scan_day']=frame['affymetrix-scan-date'].str[:10]
 frame.to_csv(out/'cel_scan_metadata.tsv',sep='\t',index=False)
 table=pd.crosstab(frame.scan_day,frame.group);table.to_csv(out/'scan_day_by_group.tsv',sep='\t')
 write_json(out/'cel_generic_headers.json',hierarchies)
 relevant=frame[frame.group.isin(['sepsis','nonfailing'])];overlap=set(relevant.loc[relevant.group.eq('sepsis'),'scan_day'])&set(relevant.loc[relevant.group.eq('nonfailing'),'scan_day'])
 summary={'cel_files_verified':len(rows),'scan_dates_available':int(frame.scan_day.ne('').sum()),'unique_scan_days':frame.scan_day.nunique(),'sepsis_nonfailing_shared_scan_days':sorted(overlap),'caution':'Scan date is a documented technical field, not patient procurement time or full laboratory processing batch. Date overlap does not establish clinical exchangeability.'}
 write_json(out/'cel_header_summary.json',summary);print(table.to_string());return summary

if __name__=='__main__':run_standard_module('87_discovery_cel_header_audit',main)
