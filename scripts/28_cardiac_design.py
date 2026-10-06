"""Packages pandas; R AnnotationDbi/org.Hs.eg.db. Input public metadata and counts IDs.
Output design audit, one-to-one annotation, protocol freeze under outputs/analysis/cardiac.
No drug effect inspection. Paths from ROOT, seed/config recorded and module logged.
"""
import gzip,json,subprocess,xml.etree.ElementTree as ET
import pandas as pd
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,r_environment,run_standard_module

def main():
 out=path('analysis')/'cardiac';out.mkdir(exist_ok=True)
 if (out/'donor_reversal_scores.tsv').exists():
  raise RuntimeError('Effects already exist; do not overwrite pre-effect freeze')
 cfg=ROOT/'config/cardiac_v1.json';prep=path('preparation')/'GSE217421'
 raw=path('raw')/'geo/GSE217421'
 for p in raw.glob('*.gz'):
  provenance=json.loads(p.with_name(p.name+'.provenance.json').read_text())
  assert sha256(p)==provenance['sha256'],p.name
 d=pd.read_csv(prep/'VEM_CTRL_design.tsv',sep='\t',dtype=str)
 d['block']=d.Cell+'.'+d.Experiment
 tab=d.groupby(['Cell','Experiment','block','State']).size().unstack(fill_value=0).reset_index()
 tab['matched']=(tab.CTRL>0)&(tab.VEM>0);tab.to_csv(out/'culture_coverage.tsv',sep='\t',index=False)
 donors=sorted(d.loc[d.State=='VEM','Cell'].unique());matched=set(tab.loc[tab.matched,'block'])
 d['primary_included']=d.block.isin(matched)
 d['same_donor_sensitivity']=d.Cell.isin(donors)
 d['exclusion_reason']=d.apply(lambda r:'' if r.primary_included else ('no VEM in donor' if r.Cell not in donors else 'no opposite condition in same culture'),axis=1)
 d.to_csv(out/'analysis_design.tsv',sep='\t',index=False)
 assert len(donors)==5 and int((d.primary_included&(d.State=='VEM')).sum())==16
 assert d.Plate.nunique()==1 and set(d.Time)=={'48'}
 for donor in donors:assert tab[(tab.Cell==donor)&tab.matched].shape[0]>=2
 with gzip.open(raw/'GSE217421_Conv-RNAseq-Configs.tsv.gz','rt',encoding='utf-8') as f:config_text=f.read()
 assert 'field_names_positions_value\t1, 3, 4, 5, 6, 7, 8, 2' in config_text
 evidence={'field_mapping':{'Cell':'Subject','Experiment':'Culture','Dish':'Replicate','Plate':'Measure','Well':'Sample'},'source':'GEO Configs field_names_positions_key/value','cautions':['Culture is not independently verified differentiation batch','All Plate values are 0, no measured plate adjustment possible','No published randomization assumption','Healthy-donor cardiac cultures have documented cell heterogeneity']}
 root=ET.parse(ROOT/'data/reference/mechanism/PMC11390749_fulltext.xml').getroot()
 evidence['title']=root.findtext('.//article-title');evidence['doi']='10.1038/s41467-024-52145-4'
 sections=[]
 for sec in root.findall('.//sec'):
  if sec.findtext('title','') in ['Treatment of cardiomyocyte cell lines','Healthy human subject iPSC lines','Bulk transcriptomics','Identification of differentially expressed genes']:
   sections.append({'title':sec.findtext('title'),'paragraphs':[''.join(el.itertext()) for el in sec.findall('./p')]})
 evidence['methods_excerpts']=sections;write_json(out/'design_source_evidence.json',evidence)
 subprocess.run([RUNTIME['rscript'],str(ROOT/'scripts/28_cardiac_design.R')],cwd=ROOT,env=r_environment(),check=True)
 primary=d[d.primary_included];summary={'donors':donors,'input_records':len(d),'primary_records':len(primary),'primary_VEM':int((primary.State=='VEM').sum()),'primary_CTRL':int((primary.State=='CTRL').sum()),'matched_cultures':len(matched),'cultures_per_donor':tab[tab.matched].groupby('Cell').size().to_dict(),'excluded_records':len(d)-len(primary),'unmatched_VEM':int(((d.State=='VEM')&~d.primary_included).sum()),'same_donor_sensitivity_records':int(d.same_donor_sensitivity.sum()),'independent_units':5,'scope':'pre-effect metadata/ID audit'}
 write_json(out/'design_summary.json',summary)
 files=[cfg,ROOT/'docs/10_心肌药物扰动方案_效应分析前.md',prep/'VEM_CTRL_counts.tsv.gz',prep/'VEM_CTRL_design.tsv',out/'analysis_design.tsv',out/'symbol_entrez_one_to_one.tsv',out/'symbol_entrez_all.tsv',out/'design_source_evidence.json',ROOT/'data/reference/mechanism/PMC11390749_fulltext.xml',ROOT/'outputs/analysis/mechanism/fixed_pathways.tsv']
 signature=json.loads((ROOT/'outputs/analysis/human/frozen_signatures/signature_freeze_manifest.json').read_text())
 for f in signature['files']:
  p=ROOT/f['path'];assert sha256(p)==f['sha256'];files.append(p)
 freeze={'frozen_utc':__import__('pipeline_utils').now(),'version':'cardiac-1.0','drug_effects_not_inspected':True,'prior_LINCS_rank_known':True,'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in files]}
 write_json(out/'cardiac_plan_freeze.json',freeze)
 return summary

if __name__=='__main__':run_standard_module('28_cardiac_design',main)
