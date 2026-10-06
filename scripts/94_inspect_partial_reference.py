"""Read-only ZIP member CRC verification; recover only fully transferred members.
Incomplete public EPMC HTTP archives are retained. No results are computed.
"""
import struct,zlib,zipfile
from pipeline_utils import ROOT,write_json,sha256,now
p=ROOT/'outputs/analysis/biology_extension_v1/sources/heart_epmc_supplements.zip.attempt2.part'
b=p.read_bytes();i=0;rows=[]
while b[i:i+4]==b'PK\x03\x04':
 vals=struct.unpack_from('<IHHHHHIIIHH',b,i);_,v,flags,method,tm,dt,crc,cs,us,n,e=vals
 name=b[i+30:i+30+n].decode('utf8');start=i+30+n+e
 if method==8:
  obj=zlib.decompressobj(-15)
  try:raw=obj.decompress(b[start:])
  except zlib.error as err:print(name,err);break
  if not obj.eof:print('INCOMPLETE',name,len(raw));break
  end=len(b)-len(obj.unused_data)
 elif method==0 and cs:raw=b[start:start+cs];end=start+cs
 else:print('UNSUPPORTED',name,flags,method,cs);break
 if flags&8:
  j=end+4 if b[end:end+4]==b'PK\x07\x08' else end
  if len(b)<j+12:break
  crc,cs,us=struct.unpack_from('<III',b,j);end=j+12
 ok=zlib.crc32(raw)==crc and len(raw)==us
 row={'name':name,'uncompressed':len(raw),'crc_verified':ok};rows.append(row);print(row)
 if ok and 'MOESM4' in name:
  dest=p.parent/'heart_tables_recovered_crc.zip'
  assert not dest.exists();dest.write_bytes(raw)
  with zipfile.ZipFile(dest) as z:assert z.testzip() is None;print(z.namelist())
  write_json(dest.with_suffix('.provenance.json'),{'source':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7681775/supplementaryFiles','retrieved_utc':now(),'archive_incomplete':True,'complete_member':name,'member_CRC32_verified':True,'nested_archive_all_member_CRCs_verified':True,'sha256':sha256(dest),'partial_sha256':sha256(p)})
 i=end
write_json(p.parent/'partial_member_audit.json',rows)
