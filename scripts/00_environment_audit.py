"""Smoke-check required packages and R execution; save machine-readable evidence."""
import importlib,importlib.metadata,platform,shutil,subprocess,sys
from pipeline_utils import ROOT,RUNTIME,path,r_environment,write_json,run_standard_module,now
def main():
    packages=["requests","numpy","pandas","scipy","sklearn","h5py","matplotlib","seaborn","networkx","yaml","openpyxl","xlrd","nbformat","ipykernel","rdkit"]
    versions={}
    for name in packages:
        m=importlib.import_module(name);versions[name]=getattr(m,"__version__","import_passed")
    import numpy as np,h5py
    from rdkit import Chem
    assert Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O") is not None
    test=path("audit")/"environment_smoke.h5"
    with h5py.File(test,"w") as f:f.create_dataset("test",data=np.eye(3))
    with h5py.File(test,"r") as f:assert float(f["test"][:].sum())==3
    r=subprocess.run([RUNTIME["rscript"],"--vanilla","-e","cat(R.version.string)"],env=r_environment(),capture_output=True,text=True)
    if r.returncode:raise RuntimeError(r.stderr)
    result={"checked_utc":now(),"python":sys.version,"executable":sys.executable,"platform":platform.platform(),"packages":versions,
      "r":r.stdout.strip(),"r_warnings":r.stderr,"disk_free_gb":round(shutil.disk_usage(ROOT).free/1024**3,2),"project_library":RUNTIME["r_library"]}
    write_json(path("audit")/"environment_audit.json",result)
    return result
if __name__=="__main__":run_standard_module("00_environment_audit",main)
