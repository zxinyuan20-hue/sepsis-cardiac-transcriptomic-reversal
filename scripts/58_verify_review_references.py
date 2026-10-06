"""Packages requests/pandas; independent public metadata/full-text retrieval.
Input module57 records and explicitly selected closest literature; output source receipts.
Shared wrapper fixes seed/log path; no statistics or biological selection changes.
"""
import importlib,json,re,concurrent.futures,xml.etree.ElementTree as ET
from urllib.parse import quote
from difflib import SequenceMatcher
import pandas as pd
from pipeline_utils import ROOT,write_json,run_standard_module,sha256
search=importlib.import_module('57_review_literature');OUT=search.OUT;fetch=search.fetch
SELECTED=['28067713','37996554','39861106','40732369','34475469','39261481','42674821','42123636','41863100','37072707','38323474','32492702']

def main():
    allrows=json.loads((OUT/'literature_records.json').read_text(encoding='utf8'))
    rows=[r for r in allrows if r['pmid'] in SELECTED];assert len(rows)==len(SELECTED)
    def process(r):
        row={k:r[k] for k in ['pmid','doi','pmcid','title','journal','year','authors','publication_types','corrections']}; errors=[]
        try:
            p=fetch(r['pmid']+'_crossref.json','https://api.crossref.org/works/'+quote(r['doi'],safe=''))
            x=json.loads(p.read_text(encoding='utf8'))['message']
            norm=lambda t:re.sub(r'[^a-z0-9]','',re.sub('<[^>]+>','',t).lower())
            similarity=SequenceMatcher(None,norm(r['title']),norm(x['title'][0])).ratio()
            row.update(title_similarity=similarity,crossref_status='VERIFIED' if similarity>=.9 and x['DOI'].lower()==r['doi'] else 'MANUAL_NEEDED',crossref_file=p.relative_to(ROOT).as_posix())
        except Exception as e: row['crossref_status']='FAILED';errors.append(str(e))
        row['fulltext_status']='NOT_AVAILABLE';row['fulltext_file']=''
        if r['pmcid']:
            try:
                pmc=r['pmcid']
                p=fetch(pmc+'_fulltext.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML')
                tree=ET.parse(p);body=tree.find('.//body')
                if body is None:raise ValueError('No article body')
                paragraphs=[' '.join(''.join(e.itertext()).split()) for e in body.iter() if e.tag in ['p','title','table']]
                text='\n\n'.join(paragraphs);assert len(text)>100
                path=OUT/(pmc+'_body.txt');path.write_text(text,encoding='utf8')
                row.update(fulltext_status='ACQUIRED_NOT_AUTOMATICALLY_REVIEWED',fulltext_file=path.relative_to(ROOT).as_posix(),fulltext_sha256=sha256(path))
            except Exception as e: errors.append(str(e))
        row['retrieval_errors']='; '.join(errors);return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(process,rows))
    pd.DataFrame(results).to_csv(OUT/'selected_reference_verification.tsv',sep='\t',index=False)
    write_json(OUT/'selected_reference_verification.json',results)
    return {'selected':len(results),'crossref_verified':sum(r['crossref_status']=='VERIFIED' for r in results),'fulltexts_acquired':sum(bool(r['fulltext_file']) for r in results),'errors':[{'pmid':r['pmid'],'error':r['retrieval_errors']} for r in results if r['retrieval_errors']]}

if __name__=='__main__':run_standard_module('58_verify_review_references',main)
