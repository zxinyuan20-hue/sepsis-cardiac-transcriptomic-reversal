"""Packages numpy/pandas/scipy; R limma. Fixed config/input freeze; outputs/external_human.
Patient-level exact tests preserve gene correlation; descriptive cross-gene correlations
have no gene-independent P. No new drug tests. All runs logged and seed fixed.
"""
import json,itertools,subprocess,os
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from pipeline_utils import ROOT,RUNTIME,write_json,sha256,run_standard_module
OUT=ROOT/'outputs/analysis/external_human'
def bh(p):
 p=np.array(p,float);o=np.argsort(p,kind='stable');q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];r=np.empty_like(q);r[o]=np.minimum(q,1);return r
def main():
 cfg=json.loads((ROOT/'config/external_human_v1.json').read_text(encoding='utf-8'))
 assert json.loads((ROOT/'outputs/audit/43_prepare_external_human_status.json').read_text(encoding='utf-8'))['status']=='PASS'
 for r in json.loads((OUT/'input_freeze.json').read_text(encoding='utf-8'))['files']:assert sha256(ROOT/r['path'])==r['sha256']
 h=pd.read_csv(ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
 x=pd.read_csv(OUT/'GSE141864_gene_expression.tsv.gz',sep='\t',index_col=0);x.index=x.index.astype(str)
 alt=pd.read_csv(OUT/'GSE141864_median_gene_expression.tsv.gz',sep='\t',index_col=0);alt.index=alt.index.astype(str)
 s=pd.read_csv(OUT/'GSE141864_heart_design.tsv',sep='\t');s=s[s.primary].reset_index(drop=True);ids=s.gsm.tolist();cases=s.group.eq('sepsis').to_numpy();bg=sorted(set(x.index)&set(h.index));assert len(bg)>1000
 pd.DataFrame({'entrez_id':bg}).to_csv(OUT/'shared_gene_background.tsv',sep='\t',index=False)
 queries={}
 for size in [50,100]:
  queries[size]={d:pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_{size}_{d}.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id.tolist() for d in ['up','down']}
 def scores(mat,size):
  rank=(mat.loc[bg,ids].rank(axis=0,method='average')-.5)/len(bg)
  return rank.loc[rank.index.isin(queries[size]['up'])].mean()-rank.loc[rank.index.isin(queries[size]['down'])].mean()
 v=scores(x,100).to_numpy();effect=float(v[cases].mean()-v[~cases].mean());null=[]
 for sel in itertools.combinations(range(len(v)),int(cases.sum())):
  use=np.zeros(len(v),bool);use[list(sel)]=True;null.append(float(v[use].mean()-v[~use].mean()))
 p=float(np.mean(np.array(null)>=effect-1e-14));q=bh([p,1.])[0]
 pd.DataFrame({'assignment':range(len(null)),'case_minus_control':null}).to_csv(OUT/'GSE141864_exact_permutation_distribution.tsv',sep='\t',index=False)
 tab=s[['gsm','patient_id','group']].copy();tab['primary_disease_rank_score']=v;tab['query50_score']=scores(x,50).values;tab['median_collapse_score']=scores(alt,100).values;tab.to_csv(OUT/'GSE141864_patient_signature_scores.tsv',sep='\t',index=False)
 lod=[]
 for i,r in s.iterrows():
  keep=np.arange(len(s))!=i;z=v[keep];g=cases[keep];lod.append({'omitted_gsm':r.gsm,'omitted_group':r.group,'effect':float(z[g].mean()-z[~g].mean()),'preprocessing':'held fixed; direction diagnostic only'})
 pd.DataFrame(lod).to_csv(OUT/'GSE141864_leave_one_patient_out.tsv',sep='\t',index=False)
 env=os.environ.copy();env['STUDY_SEED']=str(cfg['seed'])
 r=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/44_external_human_analysis.R')],cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace');(OUT/'R_analysis_console.txt').write_text(r.stdout+'\n'+r.stderr,encoding='utf-8');assert r.returncode==0,r.stderr[-2000:]
 genes=pd.read_csv(OUT/'GSE141864_all_gene_results.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id');joint=h[['logFC']].join(genes[['logFC','adj.P.Val']],how='inner',lsuffix='_discovery',rsuffix='_external');joint.to_csv(OUT/'GSE141864_discovery_gene_comparison.tsv',sep='\t')
 rho=float(spearmanr(joint.logFC_discovery,joint.logFC_external).statistic)
 pathway=pd.read_csv(OUT/'GSE141864_fixed_pathway_tests.tsv',sep='\t')
 for mode,idx in pathway.groupby('mode').groups.items():
  vals=pathway.loc[idx,'PValue'].fillna(1).to_numpy();pathway.loc[idx,'planned_10_test_BH_q']=bh(np.r_[vals,np.ones(5)])[:len(vals)]
 pathway.to_csv(OUT/'GSE141864_fixed_pathway_results.tsv',sep='\t',index=False)
 # Author individual contrasts: only separate descriptive direction, no p pooling.
 ann=pd.read_csv(ROOT/'outputs/preparation/human_probe_gene_mapping.tsv',sep='\t',dtype=str).dropna(subset=['SYMBOL','ENTREZID']).drop_duplicates(['SYMBOL','ENTREZID']);ann=ann[ann.groupby('SYMBOL').ENTREZID.transform('nunique').eq(1)].drop_duplicates('SYMBOL');lookup=ann.set_index('SYMBOL').ENTREZID
 rec=[];allrows=[];qc=[]
 for i in range(1,8):
  src=next((OUT/'sources').glob(f'*sepsis{i}.heart.tsv.gz'));t=pd.read_csv(src,sep='\t');assert t.gene_ID.is_unique
  qc.append({'patient':i,'rows':len(t),'FDR_lt_01':int((t.FDR<.1).sum()),'FDR_ge_01':int((t.FDR>=.1).sum()),'positive_logFC':int((t.logFC>0).sum()),'negative_logFC':int((t.logFC<0).sum()),'treated_is_group_label':t.treated.nunique()==1,'not_sample_expression_matrix':True})
  t['entrez_id']=t.gene_ID.map(lookup);t=t.dropna(subset=['entrez_id']);t=t[~t.entrez_id.duplicated(keep=False)].set_index('entrez_id');t['patient']='sepsis'+str(i);allrows.append(t.reset_index())
  j=h[['logFC']].join(t[['logFC']],how='inner',lsuffix='_discovery',rsuffix='_external');u=t[t.index.isin(queries[100]['up'])];d=t[t.index.isin(queries[100]['down'])]
  rec.append({'patient':'sepsis'+str(i),'shared_genes':len(j),'logFC_spearman_descriptive':float(spearmanr(j.logFC_discovery,j.logFC_external).statistic),'up_available':len(u),'up_same_direction_fraction':float((u.logFC>0).mean()),'down_available':len(d),'down_same_direction_fraction':float((d.logFC<0).mean()),'both_query_sets_majority_agree':bool((u.logFC>0).mean()>.5 and (d.logFC<0).mean()>.5),'inference':'none; shared sex-matched controls and single-case fits'})
 pd.DataFrame(rec).to_csv(OUT/'GSE237861_individual_descriptive_concordance.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(OUT/'GSE237861_author_table_schema_audit.tsv',sep='\t',index=False);pd.concat(allrows).to_csv(OUT/'GSE237861_author_mapped_gene_contrasts.tsv.gz',sep='\t',index=False)
 summary={'version':cfg['version'],'formal_cohort_test_completed':'GSE141864 FFPE5+2','primary_rank_difference':effect,'rank_percentage_points':effect*100,'exact_one_sided_p':p,'exact_assignments':len(null),'planned_two_cohort_BH_q':float(q),'GSE237861_placeholder_p':1,'GSE237861_placeholder_is_observed_test':False,'query50_difference':float(tab.query50_score[cases].mean()-tab.query50_score[~cases].mean()),'median_collapse_difference':float(tab.median_collapse_score[cases].mean()-tab.median_collapse_score[~cases].mean()),'all_leave_one_patient_differences_positive':all(z['effect']>0 for z in lod),'all_gene_logFC_spearman':rho,'shared_genes':len(joint),'external_gene_FDR05':int((genes['adj.P.Val']<.05).sum()),'GSE237861_formal_test':'NOT_RUN_INCOMPLETE_MATRIX_2CASES_1CONTROL_AVAILABLE','GSE237861_descriptive_patients_both_query_majorities':sum(z['both_query_sets_majority_agree'] for z in rec),'formal_replication_family_positive':bool(q<.05 and effect>0),'drug_validation_performed':False,'therapeutic_candidates_supported':0}
 write_json(OUT/'external_human_summary.json',summary);return summary
if __name__=='__main__':run_standard_module('44_external_human_analysis',main)
