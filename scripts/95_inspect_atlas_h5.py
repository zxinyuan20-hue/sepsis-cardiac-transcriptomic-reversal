"""Read publicly hosted HDF5 metadata with verified HTTP byte ranges.
No full matrix download and no pickle execution. requests h5py; cached blocks
and provenance under outputs. This inspection precedes localization protocol.
"""
import io,requests,h5py,json
from pipeline_utils import ROOT,write_json,sha256,now
OUT=ROOT/'outputs/analysis/biology_extension_v1/sources/atlas_ranges';OUT.mkdir(parents=True,exist_ok=True)
URL='https://cellgeni.cog.sanger.ac.uk/heartcellatlas/data/global_raw.h5ad'
class RangeFile(io.RawIOBase):
 def __init__(self):self.pos=0;self.size=5068859764;self.s=requests.Session();self.requests=[]
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,off,whence=0):self.pos=off if whence==0 else self.pos+off if whence==1 else self.size+off;return self.pos
 def readinto(self,b):
  data=self.read(len(b));b[:len(data)]=data;return len(data)
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  if n>100*1024**2:raise ValueError('Metadata inspection budget exceeded')
  result=bytearray();block=1048576
  while n:
   start=self.pos//block*block;end=min(start+block,self.size)-1;p=OUT/f'{start}_{end}.bin'
   if not p.exists():
    for attempt in range(3):
     try:
      request_url=URL if attempt==0 else URL+f'?reference_range={start}_{end}_{attempt}'
      r=requests.get(request_url,headers={'Range':f'bytes={start}-{end}'},timeout=30);break
     except requests.RequestException:
      if attempt==2:raise
    assert r.status_code==206 and r.headers['Content-Range']==f'bytes {start}-{end}/{self.size}' and len(r.content)==end-start+1
    tmp=p.with_suffix('.reader.part');tmp.write_bytes(r.content)
    if not p.exists():
     try:tmp.rename(p)
     except FileExistsError:assert p.read_bytes()==r.content
    else:assert p.read_bytes()==r.content
    write_json(p.with_suffix('.json'),{'url':URL,'range':r.headers['Content-Range'],'etag':r.headers.get('ETag'),'sha256':sha256(p),'retrieved_utc':now()})
   content=p.read_bytes();part=content[self.pos-start:self.pos-start+n];result.extend(part);self.pos+=len(part);n-=len(part)
   if not part:break
  return bytes(result)
if __name__=='__main__':
 with RangeFile() as f:
  with h5py.File(f,'r') as h:
   print('ROOT',list(h));print('UNS',list(h['uns']) if 'uns' in h else [])
   for name in ['obs','var','uns']:
    if name in h:
     def show(path,obj):
      if isinstance(obj,h5py.Dataset):print(name+'/'+path,obj.shape,obj.dtype)
     h[name].visititems(show)
