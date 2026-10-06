"""Packages: stdlib; R limma. Input config/frozen matrix, output screening/bootstrap files.
Fixed seed stratified bootstrap; log each 50 replicates, no drug rankings or rat effects.
"""
import subprocess
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,r_environment,run_standard_module
def main():
    with (path('logs')/'15_R_bootstrap_details.log').open('a',encoding='utf-8') as log:
        p=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/15_human_bootstrap.R')],cwd=ROOT,env=r_environment(),stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError('Bootstrap failed')
    out=path('analysis')/'screening'
    result={p.name:sha256(p) for p in out.glob('bootstrap_*.tsv*')}
    write_json(out/'bootstrap_file_hashes.json',result)
    return result
if __name__=='__main__':run_standard_module('15_human_bootstrap',main)
