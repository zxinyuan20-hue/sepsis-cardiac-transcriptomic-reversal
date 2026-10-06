"""Packages requests/pandas; public PubMed, CrossRef and PubChem APIs.
Inputs predeclared search config and public LINCS chemical keys; no private expression uploaded.
Outputs local raw responses, identity and literature tables; preserve errors explicitly.
"""
import json,time,xml.etree.ElementTree as ET
from urllib.parse import urlencode,quote
import pandas as pd,requests
from pipeline_utils import ROOT,path,acquire,write_json,run_standard_module
def main():
 cfg=json.loads((ROOT/'config/mechanism_search.json').read_text());out=path('analysis')/'mechanism';raw=path('reference')/'mechanism';raw.mkdir(exist_ok=True)
 errors=[];records={};queries=[]
 for item in cfg['pubmed_queries']:
  try:
   url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urlencode({'db':'pubmed','term':item['query'],'retmode':'json','retmax':cfg['max_per_query'],'sort':'relevance'})
   dest=raw/(item['id']+'_search.json');acquire(url,dest,'PubMed targeted search');data=json.loads(dest.read_text());ids=data['esearchresult']['idlist']
   queries.append({**item,'total_hits':int(data['esearchresult']['count']),'retrieved_pmids':ids});time.sleep(.4)
   if not ids:continue
   url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml'})
   dest=raw/(item['id']+'_abstracts.xml');acquire(url,dest,'PubMed indexed metadata and abstracts');tree=ET.parse(dest)
   for article in tree.findall('.//PubmedArticle'):
    pmid=article.findtext('.//MedlineCitation/PMID');title=''.join(article.find('.//ArticleTitle').itertext());abstract=' '.join(''.join(e.itertext()) for e in article.findall('.//Abstract/AbstractText'))
    main_ids=article.findall('./PubmedData/ArticleIdList/ArticleId')
    doi=next((e.text for e in main_ids if e.attrib.get('IdType')=='doi'),'')
    pmc=next((e.text for e in main_ids if e.attrib.get('IdType')=='pmc'),'')
    rec=records.setdefault(pmid,{'pmid':pmid,'doi':doi,'pmcid':pmc,'title':title,'journal':article.findtext('.//Journal/Title'),'year':article.findtext('.//PubDate/Year') or article.findtext('.//PubDate/MedlineDate'),'abstract':abstract,'query_ids':[],'evidence_scope':'PubMed metadata/abstract; full text not yet reviewed'})
    rec['query_ids'].append(item['id'])
   time.sleep(.4)
  except Exception as e:errors.append({'stage':'pubmed','id':item['id'],'error':str(e)})
 write_json(out/'literature_queries.json',queries);write_json(out/'literature_records.json',list(records.values()))
 pd.DataFrame(records.values()).to_csv(out/'literature_records.tsv',sep='\t',index=False)
 top=pd.read_csv(path('analysis')/'screening/exploratory_top10.tsv',sep='\t')
 dictionary=pd.read_csv(path('preparation')/'lincs_compound_dictionary.tsv',sep='\t',dtype=str)
 results=[]
 for row in top.itertuples():
  source=dictionary[dictionary.pert_id==row.pert_id]
  rec={'rank':row.rank,'pert_id':row.pert_id,'lincs_name':row.pert_iname,'lincs_inchikey':None,'status':'unresolved'}
  try:
   keys=source.inchi_key.dropna().unique();keys=[k for k in keys if len(k)==27]
   if len(keys)!=1:raise ValueError('No unique valid full InChIKey')
   key=keys[0];rec['lincs_inchikey']=key
   url=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{quote(key)}/property/Title,InChIKey/JSON'
   dest=raw/f'{row.pert_id}_pubchem_properties.json';acquire(url,dest,'PubChem full InChIKey identity lookup')
   props=json.loads(dest.read_text())['PropertyTable']['Properties'];exact=[p for p in props if p.get('InChIKey')==key]
   if len(exact)!=1:raise ValueError(f'Expected one exact chemical structure; got {len(exact)}')
   p=exact[0];rec.update(status='full_inchikey_match',pubchem_cid=p['CID'],pubchem_title=p.get('Title'))
   url=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{p["CID"]}/synonyms/JSON';dest=raw/f'{row.pert_id}_pubchem_synonyms.json';acquire(url,dest,'PubChem synonym evidence')
   rec['synonyms']=json.loads(dest.read_text())['InformationList']['Information'][0].get('Synonym',[])
  except Exception as e:rec['error']=str(e);errors.append({'stage':'identity','id':row.pert_id,'error':str(e)})
  results.append(rec)
 write_json(out/'compound_identity_records.json',results)
 pd.DataFrame([{k:v for k,v in r.items() if k!='synonyms'} for r in results]).to_csv(out/'compound_identity.tsv',sep='\t',index=False)
 write_json(out/'source_search_errors.json',errors)
 return {'literature_records':len(records),'queries_completed':len(queries),'identities_matched':sum(r['status']=='full_inchikey_match' for r in results),'errors':len(errors)}
if __name__=='__main__':run_standard_module('22_mechanism_sources',main)
