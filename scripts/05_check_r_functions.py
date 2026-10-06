"""Check real CEL readability and R statistical/annotation capability without disease analysis."""
import subprocess,tarfile,shutil
from pathlib import Path
from pipeline_utils import ROOT,RUNTIME,path,r_environment,run_standard_module
def main():
    folder=path("preparation")/"cel_smoke";folder.mkdir(exist_ok=True)
    with tarfile.open(path("raw")/"geo/GSE79962/GSE79962_RAW.tar") as tar:
        members=[m for m in tar.getmembers() if m.isfile() and (m.name.lower().endswith(".cel.gz") or m.name.lower().endswith(".cel"))]
        if len(members)!=51:raise ValueError(f"Expected 51 CEL files; found {len(members)}")
        m=members[0];target=folder/Path(m.name).name
        if not target.exists():
            with tar.extractfile(m) as src,target.open("xb") as dst:shutil.copyfileobj(src,dst)
    with (path("logs")/"05_R_functional_details.log").open("a",encoding="utf-8") as out:
        p=subprocess.run([RUNTIME["rscript"],"--vanilla",str(ROOT/"scripts/05_r_functional_check.R")],env=r_environment(),cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError("R functional check failed; inspect 05_R_functional_details.log")
    return {"core_functionality":"PASS","scope":"real input/annotation and synthetic statistical smoke checks only"}
if __name__=="__main__":run_standard_module("05_check_r_functions",main)
