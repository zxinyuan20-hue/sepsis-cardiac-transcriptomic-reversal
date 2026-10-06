"""Packages: requests,pandas. Inputs/outputs: config/raw_heart_reconstruction_v1.json.
Reproducibility: official ENA MD5 and byte counts, local SHA256, fixed sample manifest.
Acquire only; no changes to the frozen disease/drug analyses. Each run is logged.
"""
import concurrent.futures as cf
import hashlib,json,logging
from pathlib import Path
import pandas as pd
from pipeline_utils import ROOT,acquire,write_json,run_standard_module,sha256,now
from raw_download_utils import ranged_acquire

C=json.loads((ROOT/'config/raw_heart_reconstruction_v1.json').read_text(encoding='utf-8'))
OUT=ROOT/C['out']; SRC=OUT/'sources'
def digest(p,kind='md5'):
    h=hashlib.new(kind)
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def main():
    SRC.mkdir(parents=True,exist_ok=True)
    acquire(C['ena_url'],SRC/'ena_runs_with_md5.tsv')
    old=pd.read_csv(ROOT/C['manifest'],sep='\t')
    ena=pd.read_csv(SRC/'ena_runs_with_md5.tsv',sep='\t')
    m=old.merge(ena,on='run_accession',validate='one_to_one',suffixes=('_prior',''))
    assert len(m)==14 and m.run_accession.nunique()==14
    for col in ['experiment_accession','sample_accession','fastq_ftp','fastq_bytes']:
        assert (m[col].astype(str)==m[col+'_prior'].astype(str)).all(),col
    assert m.title.str.startswith('sepsis').sum()==7
    m[['title','run_accession','experiment_accession','sample_accession','fastq_ftp','fastq_bytes','fastq_md5']].to_csv(OUT/'frozen_fastq_manifest.tsv',sep='\t',index=False)
    def one(row):
        p=SRC/'fastq'/f'{row.run_accession}.fastq.gz'
        meta=ranged_acquire('https://'+row.fastq_ftp,p,int(row.fastq_bytes),row.fastq_md5)
        md5=digest(p)
        assert p.stat().st_size==int(row.fastq_bytes) and md5==row.fastq_md5, row.run_accession
        result=dict(title=row.title,run=row.run_accession,md5=md5,sha256=meta['sha256'],bytes=p.stat().st_size,verified_utc=now())
        write_json(OUT/'checksums'/f'{row.run_accession}.json',result)
        logging.info('ENA_CHECKSUM_PASS %s',row.run_accession)
        return result
    results=[]; failures=[]
    with cf.ThreadPoolExecutor(C['download_workers']) as pool:
        tasks={pool.submit(one,r):r.run_accession for r in m.itertuples()}
        for fut in cf.as_completed(tasks):
            try:results.append(fut.result())
            except Exception as e:failures.append({'run':tasks[fut],'error':str(e)})
    report={'verified':len(results),'expected':14,'bytes_verified':sum(r['bytes'] for r in results),'failures':failures,'analysis_performed':False}
    write_json(OUT/'raw_acquisition_audit.json',report)
    if failures:raise RuntimeError(str(failures))
    return report
if __name__=='__main__':run_standard_module('46_acquire_raw_heart',main)
