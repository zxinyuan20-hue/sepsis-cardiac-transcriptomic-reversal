"""Packages requests/pandas; configuration project_review_v1.json; shared seeded log wrapper.
Input: public search terms and explicit literature identifiers, no private data.
Output: immutable public API responses, provenance, query receipt, deduplicated records.
No biological analyses. T1 PubMed and Crossref direct APIs replace unavailable MCP.
"""
import json, time, logging, re, xml.etree.ElementTree as ET
from urllib.parse import urlencode, quote
import requests, pandas as pd
from pipeline_utils import ROOT, now, sha256, write_json, run_standard_module
CFG=json.loads((ROOT/'config/project_review_v1.json').read_text(encoding='utf8'))
OUT=ROOT/CFG['search_output']; RAW=OUT/'sources'

def fetch(key,url):
    RAW.mkdir(parents=True,exist_ok=True)
    p=RAW/key; meta=p.with_suffix(p.suffix+'.source.json')
    if p.exists():
        m=json.loads(meta.read_text(encoding='utf8')); assert m['url']==url and sha256(p)==m['sha256']
        return p
    for attempt in range(2):
        try:
            r=requests.get(url,timeout=35,headers={'User-Agent':'Academic evidence review/1.0'});r.raise_for_status()
            p.write_bytes(r.content);write_json(meta,{'url':url,'resolved_url':r.url,'retrieved_utc':now(),'sha256':sha256(p),'bytes':len(r.content)})
            return p
        except Exception:
            if attempt: raise
            time.sleep(1)

def xt(e): return '' if e is None else ''.join(e.itertext())

def pubmed(p):
    rows=[]
    for a in ET.parse(p).findall('.//PubmedArticle'):
        ids={x.get('IdType'):xt(x) for x in a.findall('./PubmedData/ArticleIdList/ArticleId')}
        au=a.findall('.//AuthorList/Author'); first=au[0].findtext('LastName','') if au else ''
        rows.append({'pmid':a.findtext('.//MedlineCitation/PMID'),'doi':ids.get('doi','').lower(),'pmcid':ids.get('pmc',''),
            'title':xt(a.find('.//ArticleTitle')),'journal':a.findtext('.//Journal/Title'),
            'year':a.findtext('.//PubDate/Year') or a.findtext('.//PubDate/MedlineDate'),
            'first_author':first,'authors':'; '.join(' '.join([x.findtext('LastName',''),x.findtext('Initials','')]).strip() for x in au),
            'abstract':' '.join(xt(x) for x in a.findall('.//Abstract/AbstractText')),
            'publication_types':'; '.join(xt(x) for x in a.findall('.//PublicationType')),
            'corrections':json.dumps([{'type':x.get('RefType'),'pmid':x.findtext('PMID')} for x in a.findall('.//CommentsCorrections')]),
            'sources':['PubMed'],'source_files':[p.relative_to(ROOT).as_posix()]})
    return rows

def crossref(x,p):
    dates=x.get('published',{}).get('date-parts',[[]]); au=x.get('author',[])
    return {'pmid':'','doi':x.get('DOI','').lower(),'pmcid':'','title':next(iter(x.get('title',[])),''),
        'journal':next(iter(x.get('container-title',[])),''),'year':str(dates[0][0]) if dates and dates[0] else '',
        'first_author':au[0].get('family','') if au else '',
        'authors':'; '.join((a.get('family','')+' '+a.get('given','')).strip() for a in au),
        'abstract':re.sub('<[^>]+>',' ',x.get('abstract','')),'publication_types':x.get('type',''),
        'corrections':json.dumps(x.get('update-to',[])),'sources':['Crossref'],'source_files':[p.relative_to(ROOT).as_posix()]}

def main():
    OUT.mkdir(parents=True,exist_ok=True); records=[]; receipts=[]; errors=[]
    base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    # First real queries also serve as API reachability preflight; failures are recorded.
    for qid,q in CFG['pubmed_queries'].items():
        q=f'({q}) AND ("1900/01/01"[Date - Publication] : "{CFG["cutoff"]}"[Date - Publication])'
        try:
            p=fetch(qid+'_search.json',base+'esearch.fcgi?'+urlencode({'db':'pubmed','term':q,'retmode':'json','retmax':CFG['pubmed_limit'],'sort':'relevance'}))
            x=json.loads(p.read_text(encoding='utf8'))['esearchresult'];ids=x['idlist'];time.sleep(.4)
            receipts.append({'source':'PubMed','id':qid,'query':q,'total':int(x['count']),'retrieved':len(ids),'truncated':int(x['count'])>len(ids)})
            if ids:
                p=fetch(qid+'_articles.xml',base+'efetch.fcgi?'+urlencode({'db':'pubmed','id':','.join(ids),'retmode':'xml'}))
                for r in pubmed(p): r['query_ids']=[qid];records.append(r)
                time.sleep(.4)
            logging.info('%s hits=%s retrieved=%s',qid,x['count'],len(ids))
        except Exception as e: errors.append({'query':qid,'error':str(e)})
    try:
        p=fetch('seed_articles.xml',base+'efetch.fcgi?'+urlencode({'db':'pubmed','id':','.join(CFG['seed_pmids']),'retmode':'xml'}))
        for r in pubmed(p):r['query_ids']=['prior_source_seed'];records.append(r)
    except Exception as e:errors.append({'query':'seeds','error':str(e)})
    for i,q in enumerate(CFG['crossref_queries']):
        try:
            p=fetch(f'crossref_{i}.json','https://api.crossref.org/works?'+urlencode({'query.bibliographic':q,'rows':CFG['crossref_limit'],'filter':'until-pub-date:'+CFG['cutoff'].replace('/','-')}))
            x=json.loads(p.read_text(encoding='utf8'))['message']
            receipts.append({'source':'Crossref','id':f'crossref_{i}','query':q,'total':x['total-results'],'retrieved':len(x['items']),'truncated':x['total-results']>len(x['items'])})
            for it in x['items']:
                r=crossref(it,p);r['query_ids']=[f'crossref_{i}'];records.append(r)
        except Exception as e:errors.append({'query':f'crossref_{i}','error':str(e)})
    for i,doi in enumerate(CFG['seed_dois']):
        try:
            p=fetch(f'seed_doi_{i}.json','https://api.crossref.org/works/'+quote(doi,safe=''))
            r=crossref(json.loads(p.read_text(encoding='utf8'))['message'],p);r['query_ids']=['prior_source_seed'];records.append(r)
        except Exception as e:errors.append({'query':doi,'error':str(e)})
    merged=[];duplicates=0
    def tokens(s):return set(re.findall(r'\w+',s.lower()))-set('a an the in of for on to and with by at'.split())
    for r in records:
        match=None
        for prior in merged:
            same_doi=bool(r['doi'] and r['doi']==prior['doi'])
            same_pmid=bool(r['pmid'] and r['pmid']==prior['pmid'])
            a=tokens(r['title']);b=tokens(prior['title'])
            fallback=(not r['doi'] or not prior['doi']) and r['first_author'] and r['first_author'].lower()==prior['first_author'].lower() and len(a&b)/max(1,len(a|b))>=.9
            if same_doi or same_pmid or fallback:match=prior;break
        if match is None:merged.append(r)
        else:
            duplicates+=1
            for k in ['sources','query_ids','source_files']:match[k]=sorted(set(match[k]+r[k]))
            for k in ['doi','pmid','pmcid','abstract','authors']:
                if not match[k] and r[k]:match[k]=r[k]
    write_json(OUT/'search_queries.json',receipts);write_json(OUT/'search_errors.json',errors);write_json(OUT/'literature_records.json',merged)
    pd.DataFrame(merged).to_csv(OUT/'literature_records.tsv',sep='\t',index=False)
    result={'status':'COMPLETE_WITH_SOURCE_ERRORS' if errors else 'PASS','retrieved_records':len(records),'deduplicated_records':len(merged),'duplicates_merged':duplicates,'queries':len(receipts),'errors':errors,'scope':CFG['scope'],'api_route':'Direct PubMed and Crossref; MCP unavailable'}
    write_json(OUT/'search_summary.json',result);return result

if __name__=='__main__':run_standard_module('57_review_literature',main)
