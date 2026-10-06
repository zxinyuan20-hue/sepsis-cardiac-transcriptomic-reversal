"""Packages requests/pandas. Public chemical name-structure and DOI verification.
Inputs fixed public candidate identifiers and selected PubMed records; output audited evidence.
Preserve full-key vs connectivity vs mismatch; retrieve available PMC article text and data links.
"""
import json,re,xml.etree.ElementTree as ET
from difflib import SequenceMatcher
from urllib.parse import quote
import pandas as pd
from pipeline_utils import ROOT,path,acquire,write_json,run_standard_module
def main():
 out=path('analysis')/'mechanism';raw=path('reference')/'mechanism';errors=[]
 identities=json.loads((out/'compound_identity_records.json').read_text(encoding='utf-8'));results=[]
 for row in identities:
  rec={k:v for k,v in row.items() if k!='synonyms'}
  try:
   name=row['lincs_name'];url=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(name,safe="")}/property/Title,InChIKey/JSON'
   dest=raw/f'{row["pert_id"]}_name_properties.json';acquire(url,dest,'PubChem name-to-structure cross-check')
   props=json.loads(dest.read_text(encoding='utf-8'))['PropertyTable']['Properties'];keys=[p.get('InChIKey','') for p in props];key=row.get('lincs_inchikey','')
   rec['name_lookup_keys']=';'.join(keys);rec['name_lookup_titles']=';'.join(str(p.get('Title','')) for p in props)
   rec['name_structure_status']='full_key_agreement' if key in keys else 'connectivity_only_agreement' if any(k[:14]==key[:14] for k in keys) else 'name_structure_disagreement'
  except Exception as e:rec['name_structure_status']='lookup_unresolved';rec['name_lookup_error']=str(e);errors.append({'stage':'name_identity','id':row['pert_id'],'error':str(e)})
  results.append(rec)
 pd.DataFrame(results).to_csv(out/'compound_identity_crosscheck.tsv',sep='\t',index=False)
 records=json.loads((out/'literature_records.json').read_text(encoding='utf-8'))
 selected=['40967330','23425388','27224288','12679140','21120453','15821757','35192138','35058449','41627811']
 evidence=[];accessions=[]
 for r in records:
  if r['pmid'] not in selected:continue
  item=r.copy()
  try:
   dest=raw/f'PMID{r["pmid"]}_crossref.json';acquire('https://api.crossref.org/works/'+quote(r['doi'],safe=''),dest,'CrossRef DOI metadata verification')
   d=json.loads(dest.read_text(encoding='utf-8'))['message'];item['crossref_doi_match']=d.get('DOI','').casefold()==r['doi'].casefold();item['crossref_title']=d.get('title',[''])[0]
   normalize=lambda s:re.sub(r'[^a-z0-9]','',s.casefold())
   item['crossref_title_similarity']=SequenceMatcher(None,normalize(r['title']),normalize(item['crossref_title'])).ratio()
   if item['crossref_title_similarity']<.85:raise ValueError('CrossRef title mismatch')
  except Exception as e:item['crossref_error']=str(e);errors.append({'stage':'crossref','id':r['pmid'],'error':str(e)})
  if r.get('pmcid'):
   try:
    dest=raw/f'{r["pmcid"]}_fulltext.xml';acquire(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{r["pmcid"]}/fullTextXML',dest,'Europe PMC available full article XML')
    tree=ET.parse(dest);text=' '.join(tree.getroot().itertext());item['fulltext_local']=dest.relative_to(ROOT).as_posix()
    if r['pmid'] in ['35058449','41627811']:
     ids=sorted(set(re.findall(r'\b(?:GSE\d+|PXD\d+|PRJNA\d+)\b',text)))
     links=[e.attrib.get('{http://www.w3.org/1999/xlink}href','') for e in tree.findall('.//ext-link')]
     accessions.append({'pmid':r['pmid'],'accessions':ids,'links':sorted(set(l for l in links if any(k in l.lower() for k in ['geo','synapse','lincs','predictox','pride','figshare','zenodo'])))})
   except Exception as e:item['fulltext_error']=str(e);errors.append({'stage':'fulltext','id':r['pmid'],'error':str(e)})
  evidence.append(item)
 write_json(out/'verified_mechanism_literature.json',evidence);write_json(out/'cardiac_perturbation_resource_leads.json',accessions);write_json(out/'mechanism_verification_errors.json',errors)
 return {'name_structure_statuses':pd.DataFrame(results).name_structure_status.value_counts().to_dict(),'selected_papers':len(evidence),'DOI_verified':sum(r.get('crossref_doi_match',False) for r in evidence),'fulltexts_acquired':sum('fulltext_local' in r for r in evidence),'errors':len(errors)}
if __name__=='__main__':run_standard_module('23_verify_mechanism_evidence',main)
