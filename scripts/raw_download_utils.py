"""Bounded HTTP ranges with durable resume and publisher MD5 verification.
Dependencies: requests. Inputs: URL, exact byte count, MD5. Outputs stay in project.
Never promotes partial or unverified data. Prior attempt files are preserved.
"""
import hashlib,json,logging,shutil,time,concurrent.futures as cf
from pathlib import Path
import requests
from pipeline_utils import ROOT,sha256,write_json,now

def md5sum(p):
    h=hashlib.md5()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def ranged_acquire(url,dest,expected_bytes,expected_md5,range_workers=1):
    dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
    side=dest.with_name(dest.name+'.provenance.json')
    if dest.exists():
        assert dest.stat().st_size==expected_bytes and md5sum(dest)==expected_md5, str(dest)
        if side.exists():assert json.loads(side.read_text())['sha256']==sha256(dest)
        return {'sha256':sha256(dest),'bytes':expected_bytes,'md5':expected_md5}
    part=dest.with_name(dest.name+'.resume.part')
    if not part.exists():
        candidates=list(dest.parent.glob(dest.name+'.attempt*.part'))
        if candidates:
            source=max(candidates,key=lambda p:p.stat().st_size)
            if 0<source.stat().st_size<=expected_bytes:shutil.copyfile(source,part)
    if not part.exists():part.touch()
    assert part.stat().st_size<=expected_bytes
    if shutil.disk_usage(dest.parent).free < expected_bytes-part.stat().st_size+8*1024**3:
        raise RuntimeError('Disk reserve below 8 GiB')
    last=time.monotonic()
    if range_workers>1:
        def fetch(bounds):
            start,end=bounds
            for attempt in range(6):
                try:
                    with requests.get(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},timeout=(20,60)) as r:
                        r.raise_for_status()
                        assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{expected_bytes}'
                        data=r.content;assert len(data)==end-start+1
                    return data
                except Exception:
                    if attempt==5:raise
                    time.sleep(attempt+1)
        with cf.ThreadPoolExecutor(range_workers) as pool:
            while part.stat().st_size<expected_bytes:
                offset=part.stat().st_size
                ranges=[(a,min(a+4*1024**2-1,expected_bytes-1)) for a in range(offset,min(offset+range_workers*4*1024**2,expected_bytes),4*1024**2)]
                batch=list(pool.map(fetch,ranges))
                with part.open('ab') as f:
                    for data in batch:f.write(data)
                    f.flush()
                if time.monotonic()-last>45:
                    logging.info('PARALLEL_RANGE %s %.1f%%',dest.name,100*part.stat().st_size/expected_bytes);last=time.monotonic()
    with requests.Session() as session:
        while part.stat().st_size<expected_bytes:
            start=part.stat().st_size;end=min(expected_bytes-1,start+4*1024**2-1)
            for attempt in range(6):
                try:
                    with session.get(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},timeout=(20,60)) as r:
                        r.raise_for_status()
                        if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {start}-{end}/{expected_bytes}':
                            raise RuntimeError('Exact bounded Range not honoured')
                        content=r.content
                        assert len(content)==end-start+1
                    with part.open('ab') as f:f.write(content);f.flush()
                    break
                except Exception:
                    if attempt==5:raise
                    logging.warning('RANGE_RETRY %s offset=%s attempt=%s',dest.name,start,attempt+1)
                    time.sleep(min(10,attempt+1))
            if time.monotonic()-last>45:
                logging.info('RANGE %s %.1f%%',dest.name,100*part.stat().st_size/expected_bytes);last=time.monotonic()
    if md5sum(part)!=expected_md5:raise RuntimeError('Publisher MD5 mismatch '+str(part))
    meta={'url':url,'retrieved_utc':now(),'bytes':expected_bytes,'sha256':sha256(part),'publisher_md5':expected_md5,'relative_path':dest.relative_to(ROOT).as_posix(),'source_kind':'public','transfer':'bounded_range_resume'}
    part.rename(dest);write_json(side,meta)
    logging.info('PUBLISHER_CHECKSUM_PASS %s',dest.name)
    return meta
