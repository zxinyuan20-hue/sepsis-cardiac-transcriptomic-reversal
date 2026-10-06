"""Packages: Python stdlib/pandas, R oligo/AnnotationDbi.
Inputs: configured raw CEL archive and sample map; outputs: outputs/preprocessing/human.
Seed: centralized config. Safe extraction -> RMA -> annotation, no disease testing.
"""
import json,re,subprocess,tarfile,shutil
from pathlib import Path
import pandas as pd
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,r_environment,run_standard_module

def main():
    out=ROOT/'outputs/preprocessing/human';out.mkdir(parents=True,exist_ok=True)
    folder=out/'cel';folder.mkdir(exist_ok=True)
    source=path('raw')/'geo/GSE79962/GSE79962_RAW.tar'
    meta=json.loads(source.with_name(source.name+'.provenance.json').read_text(encoding='utf-8'))
    if sha256(source)!=meta['sha256']:raise ValueError('CEL archive integrity mismatch')
    samples=pd.read_csv(path('preparation')/'GSE79962_sample_map.tsv',sep='\t')
    entries=[]
    with tarfile.open(source) as t:
        for m in t:
            if not m.isfile() or not m.name.lower().endswith(('.cel','.cel.gz')):continue
            name=Path(m.name).name;match=re.search(r'GSM\d+',name)
            if not match:raise ValueError('No GSM identifier in CEL filename')
            dest=folder/name
            if not dest.exists():
                with t.extractfile(m) as src,dest.open('xb') as dst:shutil.copyfileobj(src,dst)
            if dest.stat().st_size!=m.size:raise ValueError('Truncated CEL extraction')
            entries.append({'gsm':match.group(),'cel_path':dest.relative_to(ROOT).as_posix(),'sha256':sha256(dest)})
    manifest=pd.DataFrame(entries)
    if len(manifest)!=51 or manifest.gsm.duplicated().any() or set(manifest.gsm)!=set(samples.gsm):raise ValueError('CEL/sample mismatch')
    samples.merge(manifest,on='gsm',validate='one_to_one').to_csv(out/'cel_sample_manifest.tsv',sep='\t',index=False)
    with (path('logs')/'09_R_rma_details.log').open('a',encoding='utf-8') as log:
        p=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/09_human_rma.R')],cwd=ROOT,env=r_environment(),stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError('RMA failed; inspect detailed log')
    files={p.name:sha256(p) for p in out.glob('*.tsv.gz')}
    write_json(out/'matrix_hashes.json',files)
    return {'files':files,'scope':'RMA and outcome-independent annotation; no disease contrasts'}
if __name__=='__main__':run_standard_module('09_human_rma',main)
