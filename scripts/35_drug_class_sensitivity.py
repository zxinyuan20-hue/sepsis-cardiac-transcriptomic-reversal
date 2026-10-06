"""Packages pandas/scipy. Post-result descriptive class composition sensitivity.
Inputs fixed 21-drug scores and historical MoA. Outputs all class omissions and
within-class summaries. No re-ranking, new P, or altered primary results.
Paths from ROOT; deterministic; run log via standard module.
"""
import pandas as pd
from scipy.stats import spearmanr
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module

def main():
 out=path('analysis')/'expansion';p=ROOT/'docs/12_扩展结果后药物类别敏感性.md'
 write_json(out/'post_result_class_plan.json',{'recorded_utc':now(),'plan_sha256':sha256(p),'prior_overall_rho_known':True,'status':'POST_RESULT_DESCRIPTIVE_ONLY'})
 x=pd.read_csv(out/'cross_context_comparison.tsv',sep='\t')
 special={'BOR':'proteasome_inhibitors','CAR':'proteasome_inhibitors','DOX':'anthracycline_topoisomerase','AZA':'DNA_methyltransferase_inhibitor','ROS':'PPAR_related_agonist'}
 x['broad_class']=x.state.map(special).fillna('kinase_inhibitors')
 annotation=pd.read_csv(ROOT/'outputs/analysis/screening/annotated_complete_results.tsv',sep='\t')
 m=pd.read_csv(out/'compound_overlap_identity.tsv',sep='\t');m=m[m.structure_supported&m.state.isin(x.state)].merge(annotation[['pert_id','moa','target']],on='pert_id')
 x=x.merge(m[['state','moa','target']],on='state',validate='one_to_one');x.to_csv(out/'drug_class_context.tsv',sep='\t',index=False)
 rows=[]
 def rec(scope,label,d):
  r=float(spearmanr(d.harmonized_lincs_score,d.harmonized_cardiac_score).statistic) if len(d)>2 and d.harmonized_lincs_score.nunique()>1 and d.harmonized_cardiac_score.nunique()>1 else None
  rows.append({'scope':scope,'class':label,'drugs':len(d),'harmonized_spearman':r,'status':'post-result descriptive; no inferential P'})
 for c,g in x.groupby('broad_class'):
  rec('exclude_class',c,x[x.broad_class!=c])
  if len(g)>=4:rec('within_class',c,g)
 pd.DataFrame(rows).to_csv(out/'drug_class_sensitivity.tsv',sep='\t',index=False)
 write_json(out/'drug_class_sensitivity_summary.json',rows)
 return {'classes':x.broad_class.nunique(),'comparisons':len(rows),'results':rows,'scope':'post-result robustness; original scores unchanged'}

if __name__=='__main__':run_standard_module('35_drug_class_sensitivity',main)
