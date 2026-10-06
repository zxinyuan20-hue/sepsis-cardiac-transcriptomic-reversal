"""Shared ROOT-relative configuration and input-integrity checks for amendment 1."""
import json
import numpy as np
from pipeline_utils import ROOT,sha256
CFG=json.loads((ROOT/'config/robustness_amendment_v1.json').read_text(encoding='utf8'))
OUT=ROOT/CFG['output']
def verify_inputs():
 plan=json.loads((OUT/'plan_freeze.json').read_text(encoding='utf8'))
 for x in plan['files']:
  assert sha256(ROOT/x['path'])==x['sha256'],x['path']
 return plan
def layout(meta,key):
 pairs=sorted(set(zip(meta[key],meta.cell_id)))
 groups=[np.flatnonzero((meta[key].values==e)&(meta.cell_id.values==c)) for e,c in pairs]
 entities=sorted(set(meta[key]));starts=np.r_[0,np.cumsum([len(g) for g in groups])].astype(np.int32)
 compound=np.r_[0,np.cumsum([sum(e==entity for e,c in pairs) for entity in entities])].astype(np.int32)
 return entities,np.concatenate(groups).astype(np.int32),starts,compound
