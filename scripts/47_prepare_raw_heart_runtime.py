"""Packages: requests,pandas and standard library. Versioned inputs from config.
Outputs: raw_heart_v1 sources, reference, environment and metadata; module log.
Reproducibility: official MD5, exact Subread version, immutable archives.
"""
import concurrent.futures as cf
import gzip,io,json,logging,re,shutil,subprocess,zipfile
import xml.etree.ElementTree as ET
import pandas as pd
import requests
from pipeline_utils import ROOT,acquire,write_json,run_standard_module,sha256,now
from raw_download_utils import ranged_acquire,md5sum
C=json.loads((ROOT/'config/raw_heart_reconstruction_v1.json').read_text(encoding='utf8'))
OUT=ROOT/C['out'];SRC=OUT/'sources'
def main():
    acquire(C['gencode_checksums_url'],SRC/'GENCODE_v47_MD5SUMS')
    sums={line.split()[1]:line.split()[0] for line in (SRC/'GENCODE_v47_MD5SUMS').read_text().splitlines()}
    def get(k,name,md5):
        u=C[k]
        with requests.get(u,headers={'Range':'bytes=0-0'},stream=True,timeout=30) as r:
            r.raise_for_status(); assert r.status_code==206
            size=int(r.headers['Content-Range'].split('/')[-1])
        return ranged_acquire(u,SRC/name,size,md5,range_workers=4 if k=='fasta_url' else 1)
    jobs=[('subread_url','subread-2.1.1-Windows-x86_64.zip',C['subread_md5']),('gtf_url','gencode.v47.primary_assembly.annotation.gtf.gz',sums['gencode.v47.primary_assembly.annotation.gtf.gz']),('fasta_url','GRCh38.primary_assembly.genome.fa.gz',sums['GRCh38.primary_assembly.genome.fa.gz'])]
    with cf.ThreadPoolExecutor(3) as pool: results=list(pool.map(lambda x:get(*x),jobs))
    env=OUT/'environment';env.mkdir(exist_ok=True)
    with zipfile.ZipFile(SRC/'subread-2.1.1-Windows-x86_64.zip') as z:
        for n in z.namelist():
            assert (env/n).resolve().is_relative_to(env.resolve())
        z.extractall(env)
    versions={}
    for name in ['subjunc.exe','featureCounts.exe','subread-buildindex.exe']:
        p=next(env.rglob(name));r=subprocess.run([str(p),'-v'],capture_output=True,text=True)
        versions[name]={'path':p.relative_to(ROOT).as_posix(),'returncode':r.returncode,'version_text':r.stdout+r.stderr,'sha256':sha256(p)}
        assert '2.1.1' in r.stdout+r.stderr,(name,r.stdout,r.stderr)
    ref=OUT/'reference';ref.mkdir(exist_ok=True)
    for name in ['GRCh38.primary_assembly.genome.fa.gz','gencode.v47.primary_assembly.annotation.gtf.gz']:
        dst=ref/name[:-3]
        if not dst.exists():
            partial=dst.with_suffix(dst.suffix+'.part')
            with gzip.open(SRC/name,'rb') as inp,partial.open('wb') as f:shutil.copyfileobj(inp,f)
            partial.rename(dst)
    acquire('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.metadata.EntrezGene.gz',SRC/'gencode.v47.metadata.EntrezGene.gz')
    assert md5sum(SRC/'gencode.v47.metadata.EntrezGene.gz')==sums['gencode.v47.metadata.EntrezGene.gz']
    # Source Table S1 is parsed as row cells rather than flattening layout.
    archive=SRC/'PMC10568675_supplements.zip'
    acquire('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10568675/supplementaryFiles',archive,max_gb=.1)
    with zipfile.ZipFile(archive) as z:
        b=z.read('JCMM-27-3157-s002.docx')
        dst=SRC/'JCMM-27-3157-s002.docx'
        if dst.exists():assert dst.read_bytes()==b
        else:dst.write_bytes(b)
    with zipfile.ZipFile(io.BytesIO(b)) as z:r=ET.fromstring(z.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    rows=[[''.join(t.text or '' for t in c.findall('.//w:t',ns)) for c in row.findall('w:tc',ns)] for row in r.findall('.//w:tbl/w:tr',ns)]
    df=pd.DataFrame(rows[1:],columns=['source_id','age','sex','cause1','cause2','cause3','cause4'])
    df['title']=df.source_id.str.replace(' ','')+'_heart'
    manifest=pd.read_csv(OUT/'frozen_fastq_manifest.tsv',sep='\t')
    df=manifest[['title','run_accession']].merge(df,on='title',validate='one_to_one')
    assert len(df)==14 and df.sex.isin(['M','F']).all()
    df['group']=df.title.str.replace(r'\d+_heart','',regex=True)
    df.to_csv(OUT/'patient_metadata.tsv',sep='\t',index=False)
    report={'version':C['version'],'tools':versions,'reference_sources':results,'metadata_source':'PMC10568675 Table S1 JCMM-27-3157-s002.docx','metadata_source_sha256':sha256(dst),'sex_counts':df.groupby(['group','sex']).size().to_dict().__str__(),'reference_sha256':{p.name:sha256(p) for p in ref.iterdir() if p.is_file()},'ready_utc':now()}
    write_json(OUT/'runtime_reference_audit.json',report)
    return report
if __name__=='__main__':run_standard_module('47_prepare_raw_heart_runtime',main)
