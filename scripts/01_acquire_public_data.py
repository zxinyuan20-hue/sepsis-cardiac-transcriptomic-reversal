"""Acquire study GEO metadata/expression and LINCS Phase II metadata.
Packages: requests. Inputs: study config and public GEO URLs.
Outputs: immutable data/raw files plus SHA256 manifests in outputs/audit.
Seed: config seed. --large fetches CEL and LINCS Level5; no statistical analysis.
"""
import argparse, concurrent.futures, logging
from pipeline_utils import ROOT,CONFIG,path,acquire,write_json,run_standard_module

def resources(large=False):
    import json
    manifest=json.loads((ROOT/"config/data_sources.json").read_text(encoding="utf-8"))
    items=[]
    for row in manifest:
        if row["stage"]=="core" or large:
            target=(ROOT/row["path"]).resolve()
            if not target.is_relative_to((ROOT/"data/raw").resolve()):
                raise ValueError("Download target must remain under data/raw")
            items.append((row["url"],target))
    return items

def main():
    args=argparse.ArgumentParser();args.add_argument("--large",action="store_true");opts=args.parse_args()
    found=[];errors=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=CONFIG["downloads"]["max_workers"]) as pool:
        futures={pool.submit(acquire,u,p,"NCBI GEO public"):u for u,p in resources(opts.large)}
        for future in concurrent.futures.as_completed(futures):
            try:found.append(future.result())
            except Exception as e:errors.append({"url":futures[future],"error":str(e)})
    write_json(path("audit")/("acquisition_large.json" if opts.large else "acquisition_core.json"),{"files":found,"errors":errors})
    if errors:raise RuntimeError(f"{len(errors)} resources failed; details in acquisition manifest")
    return {"files":len(found),"bytes":sum(x["bytes"] for x in found)}

if __name__=="__main__":run_standard_module("01_acquire_public_data",main)
