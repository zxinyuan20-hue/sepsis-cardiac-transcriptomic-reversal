"""Packages requests/pandas; input stage39 metadata and explicit reviewed IDs below.
Output immutable full texts, GEO metadata, label fallback and Crossref verification.
Deterministic evidence retrieval, no new drug ranking. All outputs under outputs/.
"""
import importlib,json,concurrent.futures,xml.etree.ElementTree as ET
from urllib.parse import urlencode,quote
from pipeline_utils import ROOT,write_json,run_standard_module
search=importlib.import_module('39_independent_evidence_search');fetch=search.fetch;OUT=search.OUT
PMIDS=['29760707','32565338','32492702','26320862','38323474','38406775','33166663','39348200','38189000','42607350']
GEO=['GSE237861','GSE141864','GSE314561','GSE146096','GSE244740','GSE221595','GSE289264','GSE220121','GSE217423']
def main():
 records=json.loads((OUT/'literature_records.json').read_text(encoding='utf-8'));errors=[];tasks=[]
 for r in records:
  if r['pmid'] not in PMIDS:continue
  if r['pmcid']:tasks.append((r['pmcid']+'_fulltext.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/'+r['pmcid']+'/fullTextXML'))
  if r['doi']:tasks.append((r['pmid']+'_crossref.json','https://api.crossref.org/works/'+quote(r['doi'],safe='')))
 for drug in ['dasatinib','ponatinib','regorafenib']:
  tasks.append((drug+'_fda_generic.json','https://api.fda.gov/drug/label.json?'+urlencode({'search':f'openfda.generic_name:"{drug}"','limit':1})))
 for acc in GEO:tasks.append((acc+'_series.soft','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'+urlencode({'acc':acc,'targ':'self','form':'text','view':'full'})))
 # Adaptive follow-up discovered from the initial searches; these sources were
 # retrieved in the same review and are explicitly reproducible on rerun.
 for acc in ['GSE237861','GSE314561','GSE289264']:
  tasks.append((acc+'_samples.soft','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'+urlencode({'acc':acc,'targ':'gsm','form':'text','view':'full'})))
 tasks.append(('GSE141864_samples_brief.soft','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'+urlencode({'acc':'GSE141864','targ':'gsm','form':'text','view':'brief'})))
 for pmc in ['PMC7896277','PMC10940206','PMC7483871']:
  tasks.append((pmc+'_ncbi.xml','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode({'db':'pmc','id':pmc[3:],'retmode':'xml'})))
 for pmc in ['PMC10667505','PMC10568675','PMC7045056','PMC12408628']:
  tasks.append((pmc+'_fulltext.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML'))
 tasks.append(('additional_pubmed.xml','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode({'db':'pubmed','id':'37731199,32154187,37996554,37069550,40918936,32968055','retmode':'xml'})))
 def get(t):
  k,u=t
  try:fetch(k,u);return {'file':k,'status':'PASS'}
  except Exception as e:return {'file':k,'status':'FAIL','error':str(e)}
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:status=list(pool.map(get,tasks))
 for p in (OUT/'sources').glob('*_fulltext.xml'):
  try:
   tree=ET.parse(p);parts=[]
   for e in tree.findall('.//body//sec'):
    for c in e:
     if c.tag in ['title','p','table-wrap']:parts.append(' '.join(''.join(c.itertext()).split()))
   (OUT/(p.stem+'.txt')).write_text('\n\n'.join(parts),encoding='utf-8')
  except Exception as e:errors.append({'file':p.name,'error':str(e)})
 write_json(OUT/'detail_retrieval_status.json',status);write_json(OUT/'fulltext_parse_errors.json',errors)
 return {'requests':len(status),'failures':[s for s in status if s['status']=='FAIL'],'parse_errors':errors}
if __name__=='__main__':run_standard_module('40_independent_source_details',main)
