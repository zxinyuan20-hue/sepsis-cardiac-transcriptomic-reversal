"""Run the R package installer with a project library; no global R changes."""
import subprocess
from pipeline_utils import ROOT,RUNTIME,path,r_environment,run_standard_module
def main():
    with (path("logs")/"03_R_install_details.log").open("a",encoding="utf-8") as out:
        process=subprocess.run([RUNTIME["rscript"],"--vanilla",str(ROOT/"scripts/03_install_r_packages.R")],env=r_environment(),cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
    if process.returncode:raise RuntimeError("R installation failed; inspect outputs/logs/03_R_install_details.log")
    return {"library":RUNTIME["r_library"]}
if __name__=="__main__":run_standard_module("03_prepare_r_environment",main)
