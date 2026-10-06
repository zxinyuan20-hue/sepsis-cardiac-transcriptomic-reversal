"""Packages numpy/pandas/scipy and standard library.
Inputs/outputs: amendment config, original published Supplementary Data 3, frozen donor logFC.
Deterministic post-result descriptive analysis; no independent-drug P/CI or efficacy gate.
Pipeline: source annotation -> common-gene perturbation magnitude -> correlations/OLS.
Outputs source table, full joined metrics, omissions, class summaries and source hashes.
"""
import io,json
from zipfile import ZipFile
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from pipeline_utils import ROOT,run_standard_module,write_json,sha256
from robustness_utils import CFG,OUT,verify_inputs

def main():
 verify_inputs();out=OUT/'toxicity';out.mkdir(exist_ok=True);archive=OUT/'sources/Hansen2024_publisher_supplement.zip'
 provenance=json.loads(archive.with_suffix('.provenance.json').read_text());assert sha256(archive)==provenance['sha256']
 member='Supplementary_data/Supplementary Data 03 - Cardiotoxicity of used drugs.txt'
 with ZipFile(archive) as z:content=z.read(member)
 (out/'Hansen2024_Supplementary_Data03_original.txt').write_bytes(content)
 raw=pd.read_csv(io.StringIO(content.decode('utf-8-sig').replace('\r\r\n','\n')),sep='\t')
 raw.to_csv(out/'source_classification.tsv',sep='\t',index=False)
 annot=raw.rename(columns={'Drug abbreviation':'drug_state','Drug name':'source_drug_name','Drug class':'source_drug_class','Is_cardiotoxic':'source_cardiotoxicity','References for cardiotoxicity':'source_reference','Cardiotoxicity frequency':'source_frequency_code'})
 drugs=pd.read_csv(ROOT/'outputs/analysis/external_compound_v1/frozen_drugs.tsv',sep='\t')
 exp=ROOT/'outputs/analysis/expansion';filters=pd.read_csv(exp/'expanded_gene_filters.tsv',sep='\t',index_col=0)
 common=filters.index[filters[drugs.drug_state].all(axis=1)];assert len(common)>1000
 pd.Series(common,name='gene_symbol').to_csv(out/'common_gene_universe.tsv',sep='\t',index=False)
 magnitude=[]
 for r in drugs.itertuples():
  delta=pd.read_csv(exp/f'{r.drug_state}_donor_logFC.tsv.gz',sep='\t',index_col=0).loc[common]
  assert delta.shape[1]==r.donors and np.isfinite(delta.to_numpy()).all()
  means=delta.mean(axis=1).to_numpy();rms=float(np.sqrt(np.mean(means**2)))
  magnitude.append({'drug_state':r.drug_state,'RMS_mean_logFC':rms,'median_abs_mean_logFC':float(np.median(abs(means))),'mean_donor_RMS':float(np.sqrt((delta**2).mean(axis=0)).mean()),'common_genes':len(common)})
 allres=pd.read_csv(ROOT/'outputs/analysis/external_compound_v1/all_drug_results.tsv',sep='\t')
 metrics=drugs.merge(pd.DataFrame(magnitude),on='drug_state',validate='one_to_one').merge(annot,on='drug_state',how='left',validate='one_to_one')
 assert metrics.source_drug_name.notna().all()
 aliases={'azacytidine':'azacitidine'}
 assert all(aliases.get(a.lower(),a.lower())==b.lower() for a,b in zip(metrics.drug_name,metrics.source_drug_name))
 classes=pd.read_csv(exp/'drug_class_context.tsv',sep='\t')[['state','broad_class']].rename(columns={'state':'drug_state'})
 metrics=metrics.merge(classes,on='drug_state',validate='one_to_one')
 metrics.to_csv(out/'drug_annotations_and_magnitude.tsv',sep='\t',index=False)
 joined=allres.merge(metrics.drop(columns=['drug_name','donors']),on='drug_state',validate='many_to_one');joined.to_csv(out/'all_context_scores.tsv',sep='\t',index=False)
 association=[];omissions=[];classrows=[];models=[]
 for cohort,g in joined.groupby('cohort'):
  for x in ['RMS_mean_logFC','median_abs_mean_logFC']:
   for y in ['net_reversal','UP_suppression','DOWN_restoration']:
    rho=float(spearmanr(g[x],g[y]).statistic)
    loo=[]
    for state in g.drug_state:
     sub=g[g.drug_state!=state];v=float(spearmanr(sub[x],sub[y]).statistic);loo.append(v)
     omissions.append({'cohort':cohort,'magnitude':x,'endpoint':y,'omitted_drug':state,'rho':v})
    association.append({'cohort':cohort,'magnitude':x,'endpoint':y,'n':len(g),'rho':rho,'leave_one_drug_rho_min':min(loo),'leave_one_drug_rho_max':max(loo),'inference':'descriptive only; drugs share donors and controls'})
  for cl,h in g.groupby('source_cardiotoxicity',dropna=False):
   for y in ['net_reversal','UP_suppression','DOWN_restoration','RMS_mean_logFC']:
    classrows.append({'cohort':cohort,'source_class':cl,'n':len(h),'endpoint':y,'median':h[y].median(),'minimum':h[y].min(),'maximum':h[y].max()})
  for cl in sorted(g.broad_class.unique()):
   sub=g[g.broad_class!=cl]
   omissions.append({'cohort':cohort,'magnitude':'RMS_mean_logFC','endpoint':'net_reversal','omitted_drug':'class:'+cl,'rho':float(spearmanr(sub.RMS_mean_logFC,sub.net_reversal).statistic)})
  eligible=g[g.source_cardiotoxicity.isin(['Yes','No'])].copy();counts=eligible.source_cardiotoxicity.value_counts()
  if len(counts)==2 and counts.min()>=4:
   x=np.log(eligible.RMS_mean_logFC.to_numpy());y=eligible.net_reversal.to_numpy();x=(x-x.mean())/x.std(ddof=1);y=(y-y.mean())/y.std(ddof=1)
   design=np.column_stack([np.ones(len(y)),eligible.source_cardiotoxicity.eq('Yes').astype(int),x]);assert np.linalg.matrix_rank(design)==3
   beta=np.linalg.lstsq(design,y,rcond=None)[0];res=y-design@beta
   models.append({'cohort':cohort,'n':len(y),'source_yes':int(counts['Yes']),'source_no':int(counts['No']),'intercept':float(beta[0]),'source_yes_vs_no_coefficient_in_y_SD':float(beta[1]),'log_RMS_coefficient_per_SD':float(beta[2]),'descriptive_R2':float(1-np.sum(res**2)/np.sum((y-y.mean())**2)),'inference':'OLS description only; no causal adjustment, no independent-drug P or CI'})
 pd.DataFrame(association).to_csv(out/'magnitude_correlations.tsv',sep='\t',index=False)
 pd.DataFrame(omissions).to_csv(out/'drug_and_class_omissions.tsv',sep='\t',index=False)
 pd.DataFrame(classrows).to_csv(out/'source_class_summaries.tsv',sep='\t',index=False)
 pd.DataFrame(models).to_csv(out/'descriptive_OLS.tsv',sep='\t',index=False)
 primary=next(x for x in association if x['cohort']=='GSE237861' and x['magnitude']=='RMS_mean_logFC' and x['endpoint']=='net_reversal')
 summary={'source':'Hansen 2024 Supplementary Data 3; source categories preserved verbatim','archive_sha256':sha256(archive),'member':member,'source_table_sha256':sha256(out/'Hansen2024_Supplementary_Data03_original.txt'),'mapped_drugs':len(metrics),'source_class_counts':metrics.source_cardiotoxicity.value_counts().to_dict(),'common_genes':len(common),'primary_association':primary,'source_No_interpretation':'Not classified cardiotoxic in this source; not proof of clinical safety or absence of toxicity at assay exposure','ND_interpretation':'Not determined by source; not a negative label','scope':'Associations cannot distinguish protective from harmful perturbations; no candidate support upgraded'}
 write_json(out/'summary.json',summary);return summary

if __name__=='__main__':run_standard_module('80_toxicity_magnitude_context',main)
