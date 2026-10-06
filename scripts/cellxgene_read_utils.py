"""Read-only flatbuffer decoder matching cellxgene 1.3.0 NetEncoding schema.
Only numeric vectors and JSON arrays decoded; no executable data formats.
"""
import struct,json
import numpy as np,pandas as pd
def decode(b):
 def u(p):return struct.unpack_from('<I',b,p)[0]
 def field(p,n):
  v=p-struct.unpack_from('<i',b,p)[0];s=struct.unpack_from('<H',b,v)[0];off=struct.unpack_from('<H',b,v+4+n*2)[0] if 4+n*2<s else 0
  return p+off if off else None
 def ref(p):return p+u(p)
 def vector(typ,pos):
  p=ref(field(pos,0));n=u(p)
  if typ==5:return json.loads(b[p+4:p+4+n].decode())
  return np.frombuffer(b,dtype={1:'<f4',2:'<i4',3:'<u4',4:'<f8'}[typ],count=n,offset=p+4)
 root=u(0);nr=u(field(root,0));nc=u(field(root,1));c=ref(field(root,2));assert u(c)==nc
 fi=field(root,3);names=vector(b[fi],ref(field(root,4))) if fi else list(range(nc))
 cols={}
 for j in range(nc):
  col=ref(c+4+4*j);vals=vector(b[field(col,0)],ref(field(col,1)));assert len(vals)==nr;cols[names[j]]=vals
 return pd.DataFrame(cols)
