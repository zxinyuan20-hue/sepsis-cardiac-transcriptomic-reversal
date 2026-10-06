"""Verify publisher SHA512, decompress Level5 and extract a local landmark subset.
Packages: h5py, numpy, pandas. Input: frozen GSE70138 release and metadata QC.
Output: outputs/preparation/lincs/*.gctx and eligible_landmarks.h5, audit manifest.
No disease score, compound efficacy ranking, or dose-response inference is computed.
"""
import gzip,hashlib,json,logging,shutil
import h5py,numpy as np,pandas as pd
from pipeline_utils import ROOT,path,sha256,write_json,run_standard_module,now

def digest(p,algorithm="sha512"):
    h=hashlib.new(algorithm)
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
    return h.hexdigest()

def decode(array):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in array])

def main():
    raw=path("raw")/"lincs/GSE70138"
    checksum={}
    with gzip.open(raw/"GSE70138_SHA512SUMS.txt.gz","rt") as f:
        for l in f:
            if l.strip():h,n=l.strip().split(maxsplit=1);checksum[n.lstrip("*")]=h
    checks=[]
    for name,expected in checksum.items():
        p=raw/name
        if p.exists():
            actual=digest(p);checks.append({"file":name,"expected_sha512":expected,"actual_sha512":actual,"pass":actual==expected})
    write_json(path("audit")/"lincs_publisher_checksums.json",checks)
    if not checks or not all(x["pass"] for x in checks):raise ValueError("Publisher checksum mismatch or no matched files")
    src=raw/"GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
    if src.name not in checksum:raise ValueError("No publisher checksum for Level5")
    out=path("preparation")/"lincs";out.mkdir(exist_ok=True)
    dest=out/src.name.removesuffix(".gz")
    manifest=out/"decompression.json"
    if dest.exists():
        info=json.loads(manifest.read_text())
        if sha256(dest)!=info["sha256"]:raise ValueError("Local GCTX changed")
    else:
        if shutil.disk_usage(out).free<15*1024**3:raise ValueError("Insufficient space for decompression")
        temp=dest.with_suffix(".gctx.part")
        if temp.exists():raise ValueError("Incomplete decompression exists; preserve and inspect before retry")
        logging.info("Decompressing verified LINCS Level5")
        with gzip.open(src,"rb") as f,temp.open("xb") as g:shutil.copyfileobj(f,g,8*1024*1024)
        temp.rename(dest)
        write_json(manifest,{"source":str(src.relative_to(ROOT)),"sha256":sha256(dest),"bytes":dest.stat().st_size,"utc":now()})
    genes=pd.read_csv(path("preparation")/"lincs_gene_dictionary.tsv",sep="\t")
    eligible=pd.read_csv(path("preparation")/"lincs_eligible_signatures.tsv.gz",sep="\t")
    if eligible.empty:raise ValueError("Empty eligible signature list; rerun metadata audit after resolving inputs")
    sigmeta=pd.read_csv(raw/"GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",sep="\t")
    landmark=set(genes.loc[genes.pr_is_lm.eq(1),"pr_gene_id"].astype(str))
    with h5py.File(dest,"r") as f:
        rid=decode(f["0/META/ROW/id"][:]);cid=decode(f["0/META/COL/id"][:]);matrix=f["0/DATA/0/matrix"]
        if len(set(rid))!=len(rid) or len(set(cid))!=len(cid):raise ValueError("Duplicate HDF5 axis IDs")
        if set(cid)!=set(sigmeta.sig_id):raise ValueError("Level5 signature IDs disagree with metadata")
        if set(rid)!=set(genes.pr_gene_id.astype(str)):raise ValueError("Level5 gene IDs disagree with metadata")
        gidx=np.flatnonzero(np.isin(rid,list(landmark)));sidx=np.flatnonzero(np.isin(cid,eligible.sig_id))
        if len(gidx)!=len(landmark) or len(sidx)!=len(eligible):raise ValueError("Subset ID coverage is incomplete")
        orientation="signatures_by_genes" if matrix.shape==(len(cid),len(rid)) else "genes_by_signatures" if matrix.shape==(len(rid),len(cid)) else None
        if orientation is None:raise ValueError(f"Unexpected GCTX shape: {matrix.shape}")
        target=out/"eligible_landmarks.h5"
        if target.exists():
            archive=ROOT/"outputs/archive";archive.mkdir(exist_ok=True)
            import datetime
            stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            shutil.copy2(target,archive/f"eligible_landmarks_before_rebuild_{stamp}.h5")
        with h5py.File(target,"w") as h:
            h.create_dataset("gene_id",data=rid[gidx].astype("S"));h.create_dataset("sig_id",data=cid[sidx].astype("S"))
            arr=h.create_dataset("zscore",shape=(len(sidx),len(gidx)),dtype="float32",compression="gzip",chunks=(min(128,len(sidx)),len(gidx)))
            for begin in range(0,len(sidx),128):
                ix=sidx[begin:begin+128]
                block=matrix[ix,:][:,gidx] if orientation=="signatures_by_genes" else matrix[:,ix][gidx,:].T
                if not np.isfinite(block).all():raise ValueError("Nonfinite LINCS values")
                arr[begin:begin+len(ix)]=block
            h.attrs["scope"]="Input preparation only; filtered chemical signatures; no disease score"
            h.attrs["source_accession"]="GSE70138"
            h.attrs["orientation"]="signatures_by_genes"
        result={"full_shape":list(matrix.shape),"orientation":orientation,"gene_ids":len(rid),"signature_ids":len(cid),"subset_shape":[len(sidx),len(gidx)],"subset_sha256":sha256(target),"publisher_checksum_pass":True,"all_subset_values_finite":True}
    write_json(path("audit")/"lincs_matrix_qc.json",result)
    return result
if __name__=="__main__":run_standard_module("04_prepare_lincs_matrix",main)
