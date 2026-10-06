"""Packages pandas; input config/external_compound_v1.json; output frozen inputs.
Reproducibility seed in shared wrapper, hash checks against prior release manifests.
No drug score computation; identity-based gene intersection only.
"""
import json
import pandas as pd
from pipeline_utils import ROOT,sha256,write_json,run_standard_module,now
C=ROOT/'config/external_compound_v1.json'
def main():
 cfg=json.loads(C.read_text(encoding='utf8'));out=ROOT/cfg['out'];out.mkdir(parents=True,exist_ok=True)
 files=[C,ROOT/'docs/17_独立疾病背景下化合物复核方案.md']
 rel=['outputs/analysis/specificity/specificity_common_genes.tsv','outputs/analysis/specificity/donor_rank_change_vectors.tsv.gz','outputs/analysis/specificity/specificity_results.tsv','outputs/analysis/expansion/expanded_counts.tsv.gz','outputs/analysis/expansion/expanded_analysis_design.tsv','outputs/analysis/raw_heart_v1/GSE237861_all_gene_results.tsv','outputs/analysis/raw_heart_v1/entrez_counts_14.tsv.gz','outputs/analysis/raw_heart_v1/expression_filter.tsv','outputs/analysis/raw_heart_v1/patient_metadata.tsv','outputs/analysis/raw_heart_v1/patient_signature_scores.tsv','outputs/analysis/raw_heart_v1/leave_one_patient_out.tsv','outputs/analysis/external_human/GSE141864_all_gene_results.tsv','outputs/analysis/external_human/GSE141864_gene_expression.tsv.gz','outputs/analysis/external_human/GSE141864_heart_design.tsv','outputs/analysis/external_human/GSE141864_patient_signature_scores.tsv','outputs/analysis/external_human/GSE141864_leave_one_patient_out.tsv','outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv','outputs/audit/raw_heart_final_qa.json']
 files += [ROOT/p for p in rel]
 old={}
 for name in ['specificity_release_manifest.json','external_human_release_manifest.json','raw_heart_release_manifest.json']:
  for r in json.loads((ROOT/'outputs/audit'/name).read_text(encoding='utf8')):old[r['path']]=r['sha256']
 rows=[]
 for p in files:
  k=p.relative_to(ROOT).as_posix();h=sha256(p)
  if k in old:assert h==old[k],k
  rows.append({'path':k,'sha256':h})
 freeze=out/'plan_freeze.json'
 if freeze.exists():assert json.loads(freeze.read_text(encoding='utf8'))['files']==rows
 else:write_json(freeze,{'frozen_utc':now(),'client_date':'2026-10-05','files':rows,'config':cfg})
 common=pd.read_csv(ROOT/rel[0],sep='\t',dtype=str)
 for p in ['outputs/analysis/raw_heart_v1/GSE237861_all_gene_results.tsv','outputs/analysis/external_human/GSE141864_all_gene_results.tsv','outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv']:
  t=pd.read_csv(ROOT/p,sep='\t',dtype={'entrez_id':str});assert t.entrez_id.is_unique
  common=common[common.ENTREZID.isin(t.entrez_id)]
 assert common.ENTREZID.is_unique and common.SYMBOL.is_unique and len(common)>300
 common.to_csv(out/'common_genes.tsv',sep='\t',index=False)
 prior=pd.read_csv(ROOT/'outputs/analysis/specificity/specificity_results.tsv',sep='\t');assert len(prior)==21
 prior[['drug_state','drug_name','donors']].sort_values('drug_name').to_csv(out/'frozen_drugs.tsv',sep='\t',index=False)
 return {'version':cfg['version'],'common_genes':len(common),'drugs':21,'frozen_inputs':len(rows)}
if __name__=='__main__':run_standard_module('52_prepare_external_compound',main)
