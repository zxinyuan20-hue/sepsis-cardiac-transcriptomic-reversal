"""Packages: pandas, xlrd, requests; standard library.
Inputs/outputs: ROOT-relative config/discovery_exchangeability_v1.json.
Reproducibility: immutable acquisition and source hashes; no model or permutation.
Pipeline: freeze scope, verify workbook, join GSM metadata, archive paper and references.
Saves sample coverage, missingness, source passages, reference identity and audit files.
"""
import gzip,io,json,re,concurrent.futures,html
import pandas as pd
import requests
from pipeline_utils import ROOT,run_standard_module,write_json,sha256,now,acquire

def main():
 cfgpath=ROOT/'config/discovery_exchangeability_v1.json';cfg=json.loads(cfgpath.read_text(encoding='utf8'))
 out=ROOT/cfg['output'];out.mkdir(parents=True,exist_ok=True)
 inputs=[cfgpath,ROOT/cfg['clinical_workbook'],ROOT/cfg['sample_map'],ROOT/cfg['source_manuscript']]
 freeze=out/'scope_freeze.json'
 if not freeze.exists():write_json(freeze,{'recorded_utc':now(),'scope':cfg['scope'],'known_before_audit':'Prior zero hits known; original workbook contains 18 sepsis records; no score changes authorized by this audit','files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in inputs]})
 for r in json.loads(freeze.read_text(encoding='utf8'))['files']:assert sha256(ROOT/r['path'])==r['sha256']
 raw=ROOT/cfg['clinical_workbook'];side=json.loads(raw.with_name(raw.name+'.provenance.json').read_text(encoding='utf8'));assert sha256(raw)==side['sha256']
 excel=pd.ExcelFile(io.BytesIO(gzip.decompress(raw.read_bytes())),engine='xlrd')
 assert len(excel.sheet_names)==1
 df=pd.read_excel(excel,header=None);df.to_csv(out/'clinical_workbook_exact.tsv',sep='\t',index=False,header=False)
 cols=['array_id','gsm','age_years','male_1_female_0','clinical_details','dobutamine','epinephrine','norepinephrine','vasopressin','phenylephrine','LVEDD_cm','LVESD_cm','LVEF_percent']
 clinical=df.iloc[2:].copy();clinical.columns=cols;clinical=clinical[clinical.gsm.notna()].copy()
 for c in cols[2:]:
  if c!='clinical_details':clinical[c]=pd.to_numeric(clinical[c],errors='coerce')
 assert clinical.gsm.is_unique
 sm=pd.read_csv(ROOT/cfg['sample_map'],sep='\t');assert sm.gsm.is_unique
 assert set(clinical.gsm)<=set(sm.gsm)
 merged=sm[['gsm','title','source','group']].merge(clinical,on='gsm',how='left',validate='one_to_one',indicator=True)
 merged['clinical_record_available']=merged.pop('_merge').eq('both')
 merged.to_csv(out/'sample_clinical_coverage.tsv',sep='\t',index=False)
 rows=[]
 for group,g in merged.groupby('group'):
  d={'group':group,'n':len(g),'clinical_records':int(g.clinical_record_available.sum())}
  d.update({c+'_available':int(g[c].notna().sum()) for c in cols[2:] if c!='clinical_details'})
  rows.append(d)
 pd.DataFrame(rows).to_csv(out/'group_coverage.tsv',sep='\t',index=False)
 sepsis=merged[merged.group.eq('sepsis')];ef=sepsis.LVEF_percent.dropna()
 summary={'source_workbook_hash_verified':True,'sheets':excel.sheet_names,'group_coverage':rows,'sepsis_missing_clinical_gsm':sepsis.loc[~sepsis.clinical_record_available,'gsm'].tolist(),'sepsis_age_range':[float(sepsis.age_years.min()),float(sepsis.age_years.max())],'sepsis_known_male':int(sepsis.male_1_female_0.eq(1).sum()),'sepsis_known_female':int(sepsis.male_1_female_0.eq(0).sum()),'LVEF_available':len(ef),'LVEF_range':[float(ef.min()),float(ef.max())],'LVEF_below_50':int((ef<50).sum()),'clinical_complete_case_sepsis_nonfailing_design_possible':False}
 write_json(out/'metadata_summary.json',summary)
 url='https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC5315660/unicode'
 paper=out/'sources/PMC5315660_BioC.json';acquire(url,paper,source_kind='NCBI PMC BioC full text')
 data=json.loads(paper.read_text(encoding='utf8'));passages=[]
 for collection in data:
  for doc in collection['documents']:
   for idx,p in enumerate(doc['passages']):passages.append({'index':idx,'type':p.get('infons',{}).get('section_type',''),'text':p.get('text','')})
 write_json(out/'article_passages.json',passages)
 (out/'article_text.txt').write_text('\n\n'.join(f"[{p['index']}] {p['type']}\n{p['text']}" for p in passages),encoding='utf8')
 manuscript=(ROOT/cfg['source_manuscript']).read_text(encoding='utf8');refs=re.findall(r'^\[(\d+)\] (.+)$',manuscript,re.M)
 def check(item):
  n,t=item;m=re.search(r'PMID: (\d+)',t)
  if not m:return {'reference':n,'status':'DATASET_REFERENCE','text':t}
  pmid=m[1];dest=out/f'sources/PMID_{pmid}.json'
  if not dest.exists():
   for attempt in range(3):
    try:
     r=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params={'query':f'EXT_ID:{pmid} AND SRC:MED','format':'json'},timeout=30);r.raise_for_status();r.json();dest.write_bytes(r.content);break
    except requests.RequestException:
     if attempt==2:return {'reference':n,'PMID':pmid,'status':'NETWORK_UNVERIFIED','original':t}
  results=json.loads(dest.read_text(encoding='utf8'))['resultList']['result'];assert len(results)==1
  r=results[0];norm=lambda s:re.sub(r'[^a-z0-9]','',re.sub(r'<[^>]+>','',html.unescape(s)).lower())
  title_ok=norm(r['title']) in norm(t);doi=re.search(r'DOI: (\S+)',t);doi_ok=not doi or doi[1].rstrip('.').lower()==r.get('doi','').lower()
  return {'reference':n,'PMID':pmid,'title':r['title'],'doi':r.get('doi'),'year':r.get('pubYear'),'title_matches':title_ok,'doi_matches':doi_ok,'status':'VERIFIED_IDENTITY' if title_ok and doi_ok else 'REVIEW_REQUIRED','source':'Europe PMC MED bibliographic record','original':t}
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:checks=list(pool.map(check,refs))
 write_json(out/'reference_identity_checks.json',checks)
 return summary

if __name__=='__main__':run_standard_module('86_discovery_metadata_audit',main)
