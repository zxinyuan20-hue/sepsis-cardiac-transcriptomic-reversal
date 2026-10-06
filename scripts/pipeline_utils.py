"""Shared configuration, immutable source acquisition, provenance and module logging.

Packages: requests; other modules declare additional requirements.
Input: config/study_config.json. Output: outputs/logs and module-specific outputs.
Reproducibility: fixed seed from config; source hashes and UTC acquisition records.
"""
from __future__ import annotations
import contextlib, datetime as dt, hashlib, json, logging, os, random, shutil, sys, time
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "config/study_config.json").read_text(encoding="utf-8"))
_runtime_file = ROOT / "config/runtime.local.json"
RUNTIME = json.loads((_runtime_file if _runtime_file.exists() else ROOT / "config/runtime.example.json").read_text(encoding="utf-8"))
RUNTIME["rscript"] = os.environ.get("RSCRIPT", RUNTIME["rscript"])
RUNTIME["python"] = sys.executable
if os.environ.get("H9C2_SOURCE_DIR"):
    RUNTIME["local_rat_source"] = os.environ["H9C2_SOURCE_DIR"]

def path(name):
    p = ROOT / CONFIG["paths"][name]
    p.mkdir(parents=True, exist_ok=True)
    return p

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha256(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024), b""): h.update(b)
    return h.hexdigest()

def write_json(p, value):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding="utf-8")

def run_standard_module(name, function):
    random.seed(CONFIG["seed"])
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(path("logs")/f"{name}.log",encoding="utf-8"),logging.StreamHandler(sys.stdout)])
    logging.info("START %s protocol=%s", name, CONFIG["protocol_version"])
    status={"module":name,"started_utc":now(),"protocol_version":CONFIG["protocol_version"]}
    try:
        result=function()
        status.update(status="PASS",result=result)
        logging.info("FINISH %s",name)
    except Exception as exc:
        status.update(status="FAIL",error=str(exc))
        logging.exception("FAILED %s",name)
        raise
    finally:
        status["finished_utc"]=now()
        write_json(path("audit")/f"{name}_status.json",status)

def acquire(url, dest, source_kind="public", max_gb=None):
    """Download to a separate .part; final files are never overwritten.
    A verified existing source is reused even if the remote source later changes.
    Retry partials only when the server honours the exact Range; otherwise retain
    the partial and use a new attempt file. SHA256 is local integrity, not publisher proof.
    """
    dest=Path(dest); dest.parent.mkdir(parents=True,exist_ok=True)
    side=dest.with_name(dest.name+".provenance.json")
    if dest.exists():
        if not side.exists(): raise RuntimeError(f"Unmanifested existing raw file: {dest}")
        meta=json.loads(side.read_text(encoding="utf-8"))
        if meta["sha256"] != sha256(dest): raise RuntimeError(f"Hash mismatch: {dest}")
        logging.info("REUSE %s", dest.name)
        return meta
    if shutil.disk_usage(dest.parent).free < CONFIG["downloads"]["min_free_gb"]*1024**3:
        raise RuntimeError("Insufficient free space under configured reserve")
    limit=(max_gb or CONFIG["downloads"]["max_automatic_file_gb"])*1024**3
    for attempt in range(CONFIG["downloads"]["max_attempts"]):
        partial=dest.with_name(dest.name+f".attempt{attempt+1}.part")
        offset=partial.stat().st_size if partial.exists() else 0
        headers={"User-Agent":"SepsisCardiacResearch/0.1 public academic data acquisition","Accept-Encoding":"identity"}
        if offset: headers["Range"]=f"bytes={offset}-"
        try:
            with requests.get(url,stream=True,headers=headers,timeout=CONFIG["downloads"]["timeout_seconds"]) as r:
                r.raise_for_status()
                if offset and (r.status_code!=206 or not r.headers.get("Content-Range","").startswith(f"bytes {offset}-")):
                    raise RuntimeError("Server did not honour resume range; preserved partial")
                expected=int(r.headers.get("Content-Length",0))+offset
                if expected>limit: raise RuntimeError(f"File exceeds configured limit: {expected}")
                if expected and shutil.disk_usage(dest.parent).free < expected-offset+5*1024**3:
                    raise RuntimeError("Insufficient space for download and reserve")
                total=offset; last=time.monotonic()
                with partial.open("ab" if offset else "xb") as out:
                    for chunk in r.iter_content(4*1024*1024):
                        if chunk:
                            out.write(chunk); total+=len(chunk)
                            if total>limit: raise RuntimeError("Streaming size exceeds configured limit")
                            if time.monotonic()-last>25:
                                logging.info("DOWNLOAD %s %.1f MiB / %.1f MiB",dest.name,total/2**20,expected/2**20);last=time.monotonic()
                if expected and total!=expected: raise RuntimeError(f"Truncated transfer {total}/{expected}")
                if total==0: raise RuntimeError("Empty transfer")
                with partial.open("rb") as f:prefix=f.read(16)
                if dest.name.endswith(".gz") and prefix[:2]!=b"\x1f\x8b":raise RuntimeError("Expected gzip bytes, received another format")
                if dest.name.endswith(".gctx") and prefix[:8]!=b"\x89HDF\r\n\x1a\n":raise RuntimeError("Expected HDF5 signature")
                meta={"url":url,"resolved_url":r.url,"retrieved_utc":now(),"bytes":total,
                    "sha256":sha256(partial),"etag":r.headers.get("ETag"),"last_modified":r.headers.get("Last-Modified"),
                    "relative_path":dest.relative_to(ROOT).as_posix(),"source_kind":source_kind}
                partial.rename(dest)
                write_json(side,meta)
                logging.info("ACQUIRED %s %.2f MiB",dest.name,total/2**20)
                return meta
        except Exception:
            logging.exception("Attempt %s failed for %s",attempt+1,dest.name)
            if attempt+1==CONFIG["downloads"]["max_attempts"]:raise
            time.sleep(2*(attempt+1))

def r_environment():
    env=os.environ.copy()
    lib=(ROOT/RUNTIME["r_library"]).resolve();lib.mkdir(parents=True,exist_ok=True)
    env.update(R_LIBS_USER=str(lib),LANG="English",LC_ALL="C",STUDY_SEED=str(CONFIG["seed"]),OMP_NUM_THREADS=str(RUNTIME["threads"]),OPENBLAS_NUM_THREADS=str(RUNTIME["threads"]))
    return env
