"""Packages: numpy, pandas, scipy; fixed config and no random draws.
Input/output paths: config/contribution_review_v1.json.
Pipeline: verify frozen packages -> exact sign-flip floors -> BH resolution bounds
-> reconcile original q values -> export source tables and audit summary.
No new biological fitting; results describe mathematical limits, not power.
"""
import json,itertools,math
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
def multipletests(p,method='fdr_bh'):
 return None,false_discovery_control(p,method='bh')
from pipeline_utils import ROOT,run_standard_module,write_json,sha256
CFG=json.loads((ROOT/'config/contribution_review_v1.json').read_text(encoding='utf8'));OUT=ROOT/CFG['output']
def main():
 OUT.mkdir(parents=True,exist_ok=True);(OUT/'tables').mkdir(exist_ok=True);(OUT/'qa').mkdir(exist_ok=True)
 frozen={}
 for rel in ['outputs/publication/workflow_ppt_v3/workflow_manifest.json','outputs/publication/project_review_v1/review_manifest.json','outputs/audit/external_compound_release_manifest.json']:
  data=json.loads((ROOT/rel).read_text(encoding='utf8'))
  if isinstance(data,dict): data=data.get('files',data.get('manifest',[]))
  for r in data: assert sha256(ROOT/r['path'])==r['sha256'],r['path']
  frozen[rel]=len(data)
 d=pd.read_csv(ROOT/CFG['input'],sep='\t');d=d[d.cohort.eq(CFG['cohort'])].copy();M=len(d);assert M==CFG['family_size']
 p=d.iut_p.to_numpy();q=multipletests(p,method='fdr_bh')[1];assert np.allclose(q,d.iut_BH_q,atol=1e-12)
 floor=2.**(-d.donors.to_numpy());assert np.all(p>=floor-1e-12)
 d['resolution_p_floor']=floor;d['floor_BH_q_all_drugs_best_case']=multipletests(floor,method='fdr_bh')[1]
 d[['drug_name','donors','iut_p','iut_BH_q','resolution_p_floor','floor_BH_q_all_drugs_best_case']].to_csv(OUT/'tables/observed_and_best_case.tsv',sep='\t',index=False)
 rows=[]
 for n in CFG['scenario_donors']:
  a=np.arange(1,n+1,dtype=float);s=np.array(list(itertools.product([-1,1],repeat=n)));emp=float(np.mean((s@a)>=a.sum()-1e-12));assert emp==2.**(-n)
  k=math.floor(M*emp/CFG['alpha'])+1 # strict q < alpha as in original gate
  # At k tied floor P values and all others P=1, their BH q = M*p/k.
  if k<=M:assert multipletests([emp]*k+[1.]*(M-k),method='fdr_bh')[1][0]<CFG['alpha']
  rows.append({'donors':n,'assignments':2**n,'minimum_one_sided_p':emp,'best_q_if_only_one_low_p':min(1,M*emp),'minimum_equal_floor_drugs_for_strict_BH':k if k<=M else 'not_attainable','interpretation':'Hypothetical common-n design; all other P=1; not power or observed evidence'})
 pd.DataFrame(rows).to_csv(OUT/'tables/hypothetical_resolution.tsv',sep='\t',index=False)
 f=np.sort(floor);rank=np.arange(1,M+1);ok=f<CFG['alpha']*rank/M
 pd.DataFrame({'rank':rank,'ordered_floor':f,'BH_critical_value':CFG['alpha']*rank/M,'strict_crossing':ok}).to_csv(OUT/'tables/actual_design_floor_order.tsv',sep='\t',index=False)
 first=int(rank[ok][0]) if ok.any() else None
 # Exhaust all floor-or-one allocations of the actual donor count groups (3*16*5=240).
 allocations=[]
 for n6 in range(3):
  for n5 in range(16):
   for n4 in range(5):
    pp=[1/64]*n6+[1/32]*n5+[1/16]*n4+[1.]*(M-n6-n5-n4)
    num=int((multipletests(pp,method='fdr_bh')[1]<CFG['alpha']).sum())
    if num:allocations.append({'n6':n6,'n5':n5,'n4':n4,'floor_drugs':n6+n5+n4,'rejections':num})
 assert min(x['floor_drugs'] for x in allocations)==first==14
 summary={'original_results_unchanged':True,'input_sha256':sha256(ROOT/CFG['input']),'frozen_files_verified':frozen,'donor_distribution':{str(k):int(v) for k,v in d.groupby('donors').size().items()},'observed_min_q':float(q.min()),'original_q_reproduced':True,'best_case_min_q':float(d.floor_BH_q_all_drugs_best_case.min()),'best_case_max_rejections':int((d.floor_BH_q_all_drugs_best_case<CFG['alpha']).sum()),'minimum_number_floor_drugs_actual_design':first,'n4_cannot_reach_unadjusted_005':True,'isolated_floor_hit_min_donors':9,'interpretation':'Original design can support broad coordinated signals but has limited resolution for sparse hits. These are mathematical bounds; cannot attribute all non-rejection to small n, estimate power, or recommend n=9 as adequate sample size.'}
 write_json(OUT/'resolution_summary.json',summary)
 return summary
if __name__=='__main__':run_standard_module('69_exact_resolution_audit',main)
