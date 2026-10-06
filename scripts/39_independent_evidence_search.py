"""Packages: requests, pandas (project environment). Public literature/label/GEO queries.
Inputs: config/independent_evidence_v1.json. Outputs: outputs/analysis/independent_evidence.
No private data transmitted, no raw overwrites; deterministic queries, fixed date cutoff.
Execution/provenance logged with pipeline_utils; preserves every response and hash.
"""
import json,time,logging,hashlib,xml.etree.ElementTree as ET
from urllib.parse import urlencode
import requests,pandas as pd
from pipeline_utils import ROOT,write_json,now,run_standard_module

OUT=ROOT/'outputs/analysis/independent_evidence'
RAW=OUT/'sources'
def fetch(key,url):
 RAW.mkdir(parents=True,exist_ok=True);p=RAW/key;meta=p.with_suffix(p.suffix+'.source.json')
 if p.exists():
  m=json.loads(meta.read_text(encoding='utf-8'));assert hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256'];assert m['url']==url
  return p
 for attempt in range(3):
  try:
   r=requests.get(url,timeout=40);r.raise_for_status();p.write_bytes(r.content)
   write_json(meta,{'url':url,'retrieved_utc':now(),'sha256':hashlib.sha256(r.content).hexdigest(),'bytes':len(r.content)});return p
  except Exception:
   if attempt==2:raise
   time.sleep(1+attempt)
def text(e):return '' if e is None else ''.join(e.itertext())
def parse_articles(p,records,qid):
 for a in ET.parse(p).findall('.//PubmedArticle'):
  pmid=a.findtext('.//MedlineCitation/PMID');ids={e.get('IdType'):e.text for e in a.findall('./PubmedData/ArticleIdList/ArticleId')}
  r=records.setdefault(pmid,{'pmid':pmid,'doi':ids.get('doi',''),'pmcid':ids.get('pmc',''),'title':text(a.find('.//ArticleTitle')),'journal':a.findtext('.//Journal/Title'),'year':a.findtext('.//PubDate/Year') or a.findtext('.//PubDate/MedlineDate'),'abstract':' '.join(text(e) for e in a.findall('.//Abstract/AbstractText')),'publication_types':[text(e) for e in a.findall('.//PublicationType')],'query_ids':[]})
  r['query_ids'].append(qid)
def main():
 OUT.mkdir(parents=True,exist_ok=True);cfg=json.loads((ROOT/'config/independent_evidence_v1.json').read_text(encoding='utf-8'));queries=[];records={};errors=[]
 tasks=[]
 for drug in cfg['drugs']:
  for domain,terms in [('disease',cfg['disease_terms']),('cardiac',cfg['cardiac_terms'])]:tasks.append((drug+'_'+domain,f'({cfg["synonyms"][drug]}) AND ({terms})'))
 tasks.append(('human_disease_data','(sepsis OR septic) AND (myocardium OR myocardial OR cardiac OR heart) AND (transcriptome OR transcriptomic OR microarray OR RNA-seq) AND (human OR patient)'))
 for qid,query in tasks:
  query=f'({query}) AND ("1900/01/01"[Date - Publication] : "2026/10/02"[Date - Publication])'
  try:
   u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urlencode({'db':'pubmed','term':query,'retmode':'json','retmax':cfg['literature_max_per_query'],'sort':'relevance'})
   x=json.loads(fetch(qid+'_search.json',u).read_text(encoding='utf-8'))['esearchresult'];ids=x['idlist'];queries.append({'id':qid,'query':query,'total_hits':int(x['count']),'retrieved':len(ids),'truncated':int(x['count'])>len(ids),'pmids':ids});time.sleep(.4)
   if ids:
    u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml'});parse_articles(fetch(qid+'_abstracts.xml',u),records,qid);time.sleep(.4)
   logging.info('%s hits=%s fetched=%s',qid,x['count'],len(ids))
  except Exception as e:errors.append({'id':qid,'error':str(e)})
 write_json(OUT/'pubmed_queries.json',queries);write_json(OUT/'literature_records.json',list(records.values()));pd.DataFrame(records.values()).to_csv(OUT/'literature_records.tsv',sep='\t',index=False)
 # Europe PMC title/abstract search adds independent discovery coverage (records can overlap).
 for drug in cfg['drugs']:
  try:
   q=f'TITLE_ABS:{drug} AND (TITLE_ABS:sepsis OR TITLE_ABS:septic OR TITLE_ABS:lipopolysaccharide OR TITLE_ABS:endotoxemia) AND FIRST_PDATE:[1900-01-01 TO 2026-10-02]'
   fetch(drug+'_europepmc.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode({'query':q,'format':'json','pageSize':60,'resultType':'core'}))
  except Exception as e:errors.append({'id':drug+'_europepmc','error':str(e)})
  try:
   fetch(drug+'_fda_label.json','https://api.fda.gov/drug/label.json?'+urlencode({'search':f'openfda.brand_name.exact:"{cfg["brands"][drug]}"','limit':1}))
  except Exception as e:errors.append({'id':drug+'_fda','error':str(e)})
 # GEO search is a lead inventory, not proof of independent validation.
 geo_tasks=[('human_heart_sepsis','(sepsis OR septic OR endotoxemia OR lipopolysaccharide) AND (heart OR cardiac OR myocardium OR cardiomyocyte) AND "Homo sapiens"[Organism] AND gse[Entry Type]'),('human_heart_sepsis_broad','(sepsis OR septic) AND "Homo sapiens"[Organism] AND gse[Entry Type]'),('three_drug_cardiac','(dasatinib OR ponatinib OR regorafenib) AND (heart OR cardiac OR cardiomyocyte) AND gse[Entry Type]')]
 geo=[]
 for qid,q in geo_tasks:
  try:
   u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urlencode({'db':'gds','term':q,'retmode':'json','retmax':cfg['geo_max_per_query']});x=json.loads(fetch(qid+'_geo_search.json',u).read_text(encoding='utf-8'))['esearchresult'];ids=x['idlist'];geo.append({'id':qid,'query':q,'count':int(x['count']),'retrieved':len(ids),'truncated':int(x['count'])>len(ids)});time.sleep(.4)
   if ids:fetch(qid+'_geo_summary.json','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urlencode({'db':'gds','id':','.join(ids),'retmode':'json'}));time.sleep(.4)
  except Exception as e:errors.append({'id':qid,'error':str(e)})
 write_json(OUT/'geo_searches.json',geo);write_json(OUT/'search_errors.json',errors)
 return {'records':len(records),'pubmed_queries':len(queries),'geo_queries':len(geo),'errors':errors,'scope':'targeted review; capped searches explicitly recorded'}
if __name__=='__main__':run_standard_module('39_independent_evidence_search',main)
