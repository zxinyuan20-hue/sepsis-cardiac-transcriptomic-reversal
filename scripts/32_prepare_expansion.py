"""Packages pandas/numpy. Pre-effect cohort and common-landmark preparation.
Inputs frozen metadata/identity coverage and public counts. Outputs prepared counts,
pairing design, gene filter/common background, plan hashes. No drug effect summaries.
Config/paths centralized, deterministic seed and module log via pipeline_utils.
"""
import gzip,json
import numpy as np,pandas as pd
from pipeline_utils import ROOT,path,sha256,write_json,now,run_standard_module

def readj(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 out=path('analysis')/'expansion';cfg=readj(ROOT/'config/expansion_v1.json')
 if (out/'expanded_endpoint_results.tsv').exists():raise RuntimeError('Effects already exist; cannot overwrite pre-effect freeze')
 for item in readj(out/'coverage_freeze.json')['files']:assert sha256(ROOT/item['path'])==item['sha256']
 summary=readj(out/'feasibility_summary.json');assert summary['proceed_comparative_gate']
 cov=pd.read_csv(out/'all_drug_feasibility.tsv',sep='\t');eligible=cov[cov.eligible].sort_values('state')
 d=pd.read_csv(out/'all_samples_metadata_audit.tsv',sep='\t',dtype={'Experiment':str,'Dish':str,'Plate':str})
 design=[]
 for row in eligible.itertuples():
  for donor in row.donors.split('|'):
   drug=d[(d.Cell==donor)&(d.State==row.state)]
   ctrl=d[(d.Cell==donor)&(d.State=='CTRL')]
   cultures=sorted(set(drug.Experiment)&set(ctrl.Experiment));assert len(cultures)>=cfg['minimum_matched_cultures_per_donor']
   selected=pd.concat([drug[drug.Experiment.isin(cultures)],ctrl[ctrl.Experiment.isin(cultures)]]).copy()
   selected['drug_state']=row.state;selected['condition']=np.where(selected.State=='CTRL','CTRL','DRUG');selected['block']=selected.Cell+'.'+selected.Experiment
   design.append(selected)
 design=pd.concat(design,ignore_index=True).sort_values(['drug_state','Cell','Experiment','condition','column_id'])
 assert not design[['drug_state','column_id']].duplicated().any()
 design.to_csv(out/'expanded_analysis_design.tsv',sep='\t',index=False)
 cp=ROOT/'data/raw/geo/GSE217421/GSE217421_Conv-RNAseq-Read-Counts.tsv.gz'
 assert sha256(cp)==readj(cp.with_name(cp.name+'.provenance.json'))['sha256']
 with gzip.open(cp,'rt',encoding='utf-8') as f:header=f.readline().rstrip('\r\n').split('\t')
 keepcols=sorted(design.column_id.unique());x=pd.read_csv(cp,sep='\t',names=['gene_symbol']+header,skiprows=1,usecols=['gene_symbol']+keepcols,index_col=0)
 assert x.index.is_unique and np.isfinite(x.to_numpy()).all() and (x>=0).all().all() and (x.sum()>0).all()
 x.to_csv(out/'expanded_counts.tsv.gz',sep='\t')
 filters=pd.DataFrame(index=x.index);th=cfg['gene_filter']
 for row in eligible.itertuples():
  cols=design.loc[design.drug_state==row.state,'column_id'];v=x[cols]
  filters[row.state]=((v.div(v.sum(axis=0),axis=1)*1e6>=th['minimum_cpm']).sum(axis=1)>=th['minimum_samples']) & (v.sum(axis=1)>=th['minimum_total_count'])
 filters.index.name='gene_symbol';filters.to_csv(out/'expanded_gene_filters.tsv',sep='\t')
 mapping=pd.read_csv(ROOT/'outputs/analysis/cardiac/symbol_entrez_one_to_one.tsv',sep='\t',dtype=str)
 common=set(filters.index[filters.all(axis=1)]);a=np.load(ROOT/'outputs/analysis/screening/screening_inputs.npz');landmark=set(a['genes'])
 commonmap=mapping[mapping.SYMBOL.isin(common)&mapping.ENTREZID.isin(landmark)].sort_values('ENTREZID')
 commonmap.to_csv(out/'common_landmark_background.tsv',sep='\t',index=False)
 covers=[]
 sig=ROOT/'outputs/analysis/human/frozen_signatures'
 for size in [100,50,150]:
  for direction in ['up','down']:
   q=pd.read_csv(sig/f'landmarks_{size}_{direction}.tsv',sep='\t',dtype={'entrez_id':str});n=int(q.entrez_id.isin(commonmap.ENTREZID).sum());ratio=n/len(q)
   assert n>=cfg['minimum_query_genes'] and ratio>=cfg['minimum_query_coverage']
   covers.append({'size':size,'direction':direction,'source_genes':len(q),'common_landmark_genes':n,'coverage':ratio})
 pd.DataFrame(covers).to_csv(out/'common_query_coverage.tsv',sep='\t',index=False)
 results={'drug_count':len(eligible),'unique_records':len(keepcols),'drug_record_uses':int(design.condition.eq('DRUG').sum()),'control_record_uses':int(design.condition.eq('CTRL').sum()),'unique_drug_records':design[design.condition=='DRUG'].column_id.nunique(),'unique_control_records':design[design.condition=='CTRL'].column_id.nunique(),'donor_counts':design.groupby('drug_state').Cell.nunique().to_dict(),'retained_common_all_genes':len(common),'common_measured_landmarks':len(commonmap),'drug_effects_not_inspected_except_prior_VEM':True}
 write_json(out/'expansion_preparation_summary.json',results)
 files=[ROOT/'config/expansion_v1.json',ROOT/'docs/11_多药物心肌扩展方案_效应前.md',out/'all_drug_feasibility.tsv',out/'expanded_analysis_design.tsv',out/'expanded_counts.tsv.gz',out/'expanded_gene_filters.tsv',out/'common_landmark_background.tsv',ROOT/'outputs/analysis/cardiac/symbol_entrez_one_to_one.tsv',ROOT/'outputs/analysis/screening/screening_inputs.npz',ROOT/'outputs/analysis/screening/signature_metadata.tsv.gz',ROOT/'outputs/analysis/screening/complete_screening_results.tsv']
 files += [sig/f'landmarks_{sz}_{di}.tsv' for sz in [50,100,150] for di in ['up','down']]
 write_json(out/'expansion_plan_freeze.json',{'version':cfg['version'],'frozen_utc':now(),'prior_VEM_known':True,'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in files]})
 return results

if __name__=='__main__':run_standard_module('32_prepare_expansion',main)
