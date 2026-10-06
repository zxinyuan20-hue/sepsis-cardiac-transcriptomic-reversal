"""Public source acquisition only; no enrichment calculations. requests + stdlib.
ROOT-relative output, immutable success files, source URLs and SHA256.
"""
import requests,re,zipfile,io,concurrent.futures
from pipeline_utils import ROOT,now,write_json,sha256,run_standard_module
OUT=ROOT/'outputs/analysis/biology_extension_v1/sources'
def fetch(item):
 name,url=item
 try:
  p=OUT/name
  if p.exists():return {'file':name,'reused':True,'bytes':p.stat().st_size}
  r=requests.get(url,timeout=40)
  result={'url':url,'status':r.status_code,'bytes':len(r.content),'resolved_url':r.url}
  if r.status_code==200:
   p.write_bytes(r.content);result.update(sha256=sha256(p),retrieved_utc=now());write_json(p.with_name(p.name+'.provenance.json'),result)
   if name.endswith('.html'):result['links']=[x for x in re.findall(r'href="([^"]+)"',r.text) if 'MOESM4' in x]
   if zipfile.is_zipfile(p):result['zip_entries']=zipfile.ZipFile(p).namelist()[:60]
  return result
 except Exception as e:return {'url':url,'error':str(e)[:250]}
def main():
 urls=[('heart_pmc.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC7681775/'),('heart_static.zip','https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2797-4/MediaObjects/41586_2020_2797_MOESM4_ESM.zip'),('heart_media2.zip','https://media.springernature.com/supplementary/10.1038/s41586-020-2797-4/41586_2020_2797_MOESM4_ESM.zip')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:res=list(pool.map(fetch,urls))
 write_json(OUT/'reference_acquisition_attempts.json',res);print(res);return res
if __name__=='__main__':run_standard_module('93_acquire_heart_reference',main)
