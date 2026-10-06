"""Packages pandas/numpy. Inputs fixed literature-derived gene list and existing results.
Outputs gene context and descriptive summaries under outputs/analysis/independent_evidence.
No new inference or ranking; source-wide FDR retained. Seed/log via pipeline_utils.
"""
import json
import pandas as pd,numpy as np
from pipeline_utils import ROOT,write_json,sha256,now,run_standard_module
def main():
 out=ROOT/'outputs/analysis/independent_evidence';cp=ROOT/'config/independent_mechanism_context_v1.json';cfg=json.loads(cp.read_text(encoding='utf-8'))
 paths=[cp,ROOT/'outputs/preparation/human_probe_gene_mapping.tsv',ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv']
 paths += [ROOT/f'outputs/analysis/mechanism/{c}_all_gene_results.tsv' for c in cfg['cohorts'] if c not in ['human','DAS','PON','REG']]
 paths += [ROOT/f'outputs/analysis/expansion/{c}_gene_results.tsv' for c in ['DAS','PON','REG']]
 fp=out/'mechanism_context_input_freeze.json'
 if fp.exists():
  for f in json.loads(fp.read_text(encoding='utf-8'))['files']:assert sha256(ROOT/f['path'])==f['sha256']
 else:write_json(fp,{'fixed_utc':now(),'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in paths]})
 mapping=pd.read_csv(paths[1],sep='\t',dtype=str).dropna(subset=['SYMBOL','ENTREZID'])
 tables={'human':pd.read_csv(paths[2],sep='\t',dtype={'entrez_id':str}).rename(columns={'entrez_id':'human_entrez'})}
 for c in cfg['cohorts']:
  if c=='human':continue
  if c in ['DAS','PON','REG']:t=pd.read_csv(ROOT/f'outputs/analysis/expansion/{c}_gene_results.tsv',sep='\t',dtype={'entrez_id':str}).rename(columns={'entrez_id':'human_entrez'})
  else:t=pd.read_csv(ROOT/f'outputs/analysis/mechanism/{c}_all_gene_results.tsv',sep='\t',dtype={'human_entrez':str})
  tables[c]=t
 rows=[]
 for family in cfg['families']:
  for symbol in family['genes']:
   ids=mapping.loc[mapping.SYMBOL.eq(symbol),'ENTREZID'].unique()
   assert len(ids)<=1,(symbol,ids)
   eid=ids[0] if len(ids)==1 else None
   for cohort,t in tables.items():
    z=t[t.human_entrez.eq(eid)] if eid else t.iloc[:0];assert len(z)<=1
    rec={'family':family['id'],'pmid':family['pmid'],'symbol':symbol,'human_entrez':eid,'cohort':cohort,'available':len(z)==1,'logFC':np.nan,'source_genome_FDR':np.nan,'interpretation_limit':family['limitation']}
    if len(z):rec.update(logFC=float(z.logFC.iloc[0]),source_genome_FDR=float(z['adj.P.Val'].iloc[0]))
    rows.append(rec)
 result=pd.DataFrame(rows);result.to_csv(out/'literature_fixed_gene_context.tsv',sep='\t',index=False)
 pivot=result.pivot(index=['family','symbol'],columns='cohort',values='logFC');q=result.pivot(index=['family','symbol'],columns='cohort',values='source_genome_FDR')
 checks=[]
 for (fam,symbol),r in pivot.iterrows():
  shared=pd.notna(r['human']) and pd.notna(r['LOCAL_RAT_LPS_2025']);agree=bool(r['human']*r['LOCAL_RAT_LPS_2025']>0) if shared else None
  rec={'family':fam,'symbol':symbol,'human_logFC':r['human'],'rat_logFC':r['LOCAL_RAT_LPS_2025'],'human_rat_same_direction':agree,'both_disease_genome_FDR_lt_005':bool(q.loc[(fam,symbol),'human']<.05 and q.loc[(fam,symbol),'LOCAL_RAT_LPS_2025']<.05) if shared else None}
  for drug in ['DAS','PON','REG']:
   rec[drug+'_logFC']=r[drug];rec[drug+'_opposes_both_disease']=bool(agree and r[drug]*r['human']<0) if shared and pd.notna(r[drug]) else None
  checks.append(rec)
 check=pd.DataFrame(checks);check.to_csv(out/'literature_fixed_gene_direction_check.tsv',sep='\t',index=False)
 summ=[]
 for fam,g in check.groupby('family'):
  summ.append({'family':fam,'listed':len(g),'human_rat_available':int(g.human_rat_same_direction.notna().sum()),'same_direction':int(g.human_rat_same_direction.fillna(False).sum()),'same_direction_both_genome_FDR':int((g.human_rat_same_direction.fillna(False)&g.both_disease_genome_FDR_lt_005.fillna(False)).sum()),**{d+'_opposes_both':int(g[d+'_opposes_both_disease'].fillna(False).sum()) for d in ['DAS','PON','REG']}})
 pd.DataFrame(summ).to_csv(out/'mechanism_context_descriptive_summary.tsv',sep='\t',index=False)
 return {'listed_genes':len(check),'families':len(summ),'new_significance_tests':0,'mechanism_activity_validated':False,'rows':len(result)}
if __name__=='__main__':run_standard_module('41_independent_mechanism_context',main)
