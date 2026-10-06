"""Packages numpy/pandas/scipy, standard library. Configured inputs and outputs.
Deterministic independent numeric QA and evidence-source consolidation.
Preserves historical data; records inherited non-input code/document hash drift.
Saves numeric audit, 21-drug evidence table and timestamped extension result manifest.
"""
import json
import numpy as np,pandas as pd
from scipy.spatial.distance import cosine
from pipeline_utils import ROOT,run_standard_module,write_json,sha256,now
from robustness_utils import CFG,OUT,verify_inputs

def main():
 plan=verify_inputs();e=pd.read_csv(OUT/'entities/entity_results.tsv',sep='\t');dist=np.load(OUT/'entities/entity_distributions.npz')
 e=e.set_index('entity_id').loc[dist['entity_ids']];p=(1+(dist['null']>=e.score.to_numpy()).sum(axis=0))/(len(dist['null'])+1);p[e.score<=0]=1
 order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);assert np.allclose(p,e.empirical_p,atol=1e-14);assert np.allclose(q,e.BH_q,atol=1e-14)
 sig=pd.read_csv(ROOT/'outputs/analysis/screening/signature_metadata.tsv.gz',sep='\t');a=np.load(ROOT/'outputs/analysis/screening/screening_inputs.npz');query=a['queries'][list(a['query_names']).index('size_100')]
 from importlib.util import spec_from_file_location,module_from_spec
 spec=spec_from_file_location('multiscore',ROOT/'scripts/79_multiscore_sensitivity.py');mod=module_from_spec(spec);spec.loader.exec_module(mod)
 c,s=mod.alternatives(query,a['z']);coserr=max(abs(c[i]-(cosine(query.astype(float),a['z'][i])-1)) for i in [0,1,10,100,6768]);assert coserr<1e-12
 # Check fixed gene-universe RMS by an explicit squared-sum implementation.
 magn=pd.read_csv(OUT/'toxicity/drug_annotations_and_magnitude.tsv',sep='\t');genes=pd.read_csv(OUT/'toxicity/common_gene_universe.tsv',sep='\t').gene_symbol
 errors=[]
 for row in magn.itertuples():
  z=pd.read_csv(ROOT/f'outputs/analysis/expansion/{row.drug_state}_donor_logFC.tsv.gz',sep='\t',index_col=0).loc[genes].to_numpy();means=z.sum(axis=1)/z.shape[1]
  rms=np.sqrt(np.dot(means,means)/len(means));errors.append(abs(rms-row.RMS_mean_logFC))
 assert max(errors)<1e-12
 allres=pd.read_csv(ROOT/'outputs/analysis/external_compound_v1/all_drug_results.tsv',sep='\t')
 primary=allres[allres.cohort=='GSE237861'].copy().set_index('drug_state');secondary=allres[allres.cohort=='GSE141864'].copy().set_index('drug_state')
 card=pd.read_csv(OUT/'multiscore/cardiac_method_scores.tsv',sep='\t').set_index('drug_state');old=pd.read_csv(ROOT/'outputs/analysis/expansion/cross_context_comparison.tsv',sep='\t').set_index('state')
 evidence=primary[['drug_name','UP_suppression','DOWN_restoration','net_reversal','iut_BH_q']].copy()
 evidence['LINCS_positive']=old.frozen_lincs_score>0
 evidence['cardiac_three_scores_positive']=(card[['weighted_ES','negative_cosine','negative_Spearman']]>0).all(axis=1)
 evidence['primary_mean_bidirectional']=primary.bidirectional
 evidence['secondary_mean_bidirectional']=secondary.bidirectional
 evidence['primary_all_donor_omissions']=primary.all_donor_loo_bidirectional
 evidence['both_contexts_all_patient_omissions']=primary.all_patient_loo_bidirectional.astype(bool)&secondary.all_patient_loo_bidirectional.astype(bool)
 evidence['primary_formal_IUT_BH']=primary.formal_bidirectional_gate
 evidence=evidence.join(magn.set_index('drug_state')[['source_cardiotoxicity','RMS_mean_logFC']]).sort_values('drug_name')
 evidence.to_csv(OUT/'evidence_matrix.tsv',sep='\t')
 checks={'entity_P_BH_independent_reproduction':True,'cosine_scipy_max_error':coserr,'RMS_independent_max_error':max(errors),'frozen_input_hashes_verified':len(plan['files']),'historical_manifest_checks':plan['historical_checks'],'all_21_evidence_rows':len(evidence)==21,'formal_hits_unchanged':int(evidence.primary_formal_IUT_BH.sum())==0,'biological_scope':'post-result exploratory sensitivity; no efficacy or safety validation','reviewed_utc':now()}
 write_json(OUT/'numeric_audit.json',checks)
 return {k:v for k,v in checks.items() if k!='historical_manifest_checks'}

if __name__=='__main__':run_standard_module('81_audit_robustness_extension',main)
