"""Packages pandas/numpy/requests/RDKit. Inputs all public cardiac design metadata
and frozen LINCS structures/ranks. Outputs identity/dose/coverage audit; no drug effects.
Paths and repeat thresholds centralized in config. Public API responses immutable with
provenance. Seed via pipeline_utils; log every module run. No private uploads.
"""
import json,re,concurrent.futures
from urllib.parse import quote
import pandas as pd,numpy as np
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from pipeline_utils import ROOT,path,acquire,sha256,write_json,now,run_standard_module

def norm_name(s):return re.sub('[^a-z0-9]','',str(s).lower())
def main():
 cfg=json.loads((ROOT/'config/expansion_feasibility_v1.json').read_text(encoding='utf-8'))
 out=path('analysis')/'expansion';out.mkdir(exist_ok=True);raw=path('reference')/'cardiac_expansion';raw.mkdir(exist_ok=True)
 d=pd.read_csv(ROOT/'outputs/preparation/GSE217421/full_design.tsv',sep='\t',dtype=str)
 m=pd.read_csv(ROOT/'outputs/analysis/mechanism/cardiac_resource_sample_metadata.tsv',sep='\t').fillna('');m=m[m.accession==cfg['dataset']]
 rows=[]
 for r in m.to_dict('records'):
  props=dict(x.strip().split(': ',1) for x in r['!Sample_characteristics_ch1'].split(' | ') if ': ' in x)
  match=re.search(r'MSN\d{2}-\d{2}R-CM',r['!Sample_title']);cell=match.group(0) if match else ''
  parts=[props.get('treatment condition',''),cell,props.get('investigative unit',''),props.get('sample id',''),'0',props.get('well','')]
  rows.append({'column_id':'.'.join(parts),'gsm':r['gsm'],'metadata_state':props.get('treatment condition',''),'drug_name':props.get('drug name',''),'concentration':props.get('concentration',''),'dose_normalized':props.get('concentration','').replace(' ','').lower(),'hours':props.get('treatment time (hrs)','')})
 meta=pd.DataFrame(rows);assert meta.column_id.is_unique,'GEO metadata keys duplicate'
 merged=d.merge(meta,on='column_id',how='left',validate='one_to_one',indicator=True)
 merged['metadata_matched']=merged['_merge'].eq('both');merged.drop(columns='_merge',inplace=True)
 merged.to_csv(out/'all_samples_metadata_audit.tsv',sep='\t',index=False)
 compounds=[]
 for state,g in merged.groupby('State'):
  names=sorted(g.drug_name.dropna().unique());doses=sorted(g.dose_normalized.dropna().unique())
  compounds.append({'state':state,'drug_name':names[0] if len(names)==1 else '|'.join(names),'names_unique':len(names),'doses':'|'.join(doses),'dose_unique':len(doses),'records':len(g),'metadata_complete':bool(g.metadata_matched.all()),'48h_all':bool(g.hours.eq('48').all())})
 drug=pd.DataFrame(compounds);drug.to_csv(out/'all_drug_metadata.tsv',sep='\t',index=False)
 query=drug[(~drug.state.isin(cfg['excluded_biologics']+['CTRL'])) & drug.names_unique.eq(1)]
 def fetch(r):
  dest=raw/f'{r.state}_pubchem_name_properties.json'
  try:
   acquire(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(r.drug_name,safe="")}/property/InChI,InChIKey,Title/JSON',dest,'PubChem named compound identity for public cardiac overlap')
   props=json.loads(dest.read_text(encoding='utf-8'))['PropertyTable']['Properties']
   return {'state':r.state,'query_name':r.drug_name,'properties':props,'status':'acquired'}
  except Exception as e:return {'state':r.state,'query_name':r.drug_name,'status':'failed','error':str(e)}
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(fetch,list(query.itertuples())))
 write_json(out/'pubchem_identity_records.json',records)
 rank=pd.read_csv(ROOT/'outputs/analysis/screening/complete_screening_results.tsv',sep='\t')
 dictionary=pd.read_csv(ROOT/'outputs/preparation/lincs_compound_dictionary.tsv',sep='\t',dtype=str)
 dictionary=dictionary[dictionary.pert_id.isin(rank.pert_id)].drop_duplicates('pert_id')
 te=rdMolStandardize.TautomerEnumerator();identity=[]
 for row in dictionary.itertuples():
  mol=Chem.MolFromSmiles(row.canonical_smiles)
  if mol is None:continue
  key=inchi.MolToInchiKey(mol);taut=inchi.MolToInchiKey(te.Canonicalize(mol))
  identity.append({'pert_id':row.pert_id,'inchi_key':row.inchi_key,'computed_key':key,'tautomer_key':taut,'source_key_agrees':key==row.inchi_key})
 ident=pd.DataFrame(identity);ident.to_csv(out/'screened_structure_keys.tsv',sep='\t',index=False)
 rtab=rank.merge(ident,on='pert_id',validate='one_to_one');rtab['name_key']=rtab.pert_iname.map(norm_name)
 matches=[]
 for rec in records:
  if rec['status']!='acquired':continue
  for prop in rec['properties']:
   pmol=inchi.MolFromInchi(prop['InChI']);pk=prop['InChIKey'];tk=inchi.MolToInchiKey(te.Canonicalize(pmol))
   for row in rtab.itertuples():
    exact=row.inchi_key==pk;taut=row.tautomer_key==tk;name=row.name_key==norm_name(rec['query_name'])
    if exact or taut or name:
     status='full_key' if exact else ('canonical_tautomer_full_key' if taut else 'name_only_structure_unresolved')
     matches.append({'state':rec['state'],'cardiac_name':rec['query_name'],'pubchem_cid':prop['CID'],'pubchem_key':pk,'pert_id':row.pert_id,'lincs_name':row.pert_iname,'lincs_key':row.inchi_key,'rank':row.rank,'lincs_score':row.reversal_score,'identity_status':status,'source_key_agrees':row.source_key_agrees,'structure_supported':(exact or taut) and row.source_key_agrees})
 match=pd.DataFrame(matches).drop_duplicates(['state','pert_id','identity_status']);match.to_csv(out/'compound_overlap_identity.tsv',sep='\t',index=False)
 controls=merged[(merged.State=='CTRL')&merged.metadata_matched];blocks=[];coverage=[]
 for row in drug[drug.state!='CTRL'].itertuples():
  samples=merged[(merged.State==row.state)&merged.metadata_matched]
  for (cell,exp),s in samples.groupby(['Cell','Experiment']):
   ctrl=controls[(controls.Cell==cell)&(controls.Experiment==exp)]
   blocks.append({'state':row.state,'donor':cell,'experiment':exp,'drug_records':len(s),'control_records':len(ctrl),'matched':len(ctrl)>0})
  bb=pd.DataFrame([b for b in blocks if b['state']==row.state]);n=bb[bb.matched].groupby('donor').size()
  qualified=n[n>=cfg['minimum_matched_cultures_per_donor']]
  ids=match[(match.state==row.state)&match.structure_supported]
  reasons=[]
  if len(ids)==0:reasons.append('no structure-supported original-screen overlap')
  if row.dose_unique!=1:reasons.append('inconsistent dose metadata')
  if not row.metadata_complete:reasons.append('incomplete sample metadata link')
  if len(qualified)<cfg['minimum_donors_per_drug']:reasons.append('insufficient donors with >=2 matched cultures')
  if not samples.hours.eq('48').all():reasons.append('not uniformly 48h')
  coverage.append({'state':row.state,'drug_name':row.drug_name,'records':row.records,'doses':row.doses,'structure_supported_entries':len(ids),'pert_ids':'|'.join(sorted(ids.pert_id)),'lincs_ranks':'|'.join(str(int(x)) for x in ids['rank']),'frozen_lincs_score':float(ids.lincs_score.median()) if len(ids) else np.nan,'matched_cultures':int(n.sum()),'qualified_donors':len(qualified),'donors':'|'.join(sorted(qualified.index)),'eligible':not reasons,'exclusion_reasons':'; '.join(reasons)})
 cov=pd.DataFrame(coverage);cov.to_csv(out/'all_drug_feasibility.tsv',sep='\t',index=False);pd.DataFrame(blocks).to_csv(out/'all_culture_coverage.tsv',sep='\t',index=False)
 eligible=cov[cov.eligible];top=rank.head(10);topcov=top[['pert_id','rank','pert_iname']].merge(match[['pert_id','state','identity_status','structure_supported']],on='pert_id',how='left');topcov.to_csv(out/'frozen_top10_cardiac_coverage.tsv',sep='\t',index=False)
 result={'version':cfg['version'],'sample_records':len(d),'metadata_matched':int(merged.metadata_matched.sum()),'drug_states':len(drug)-1,'pubchem_queries':len(records),'pubchem_failed':sum(r['status']!='acquired' for r in records),'structure_supported_shared_drugs':match[match.structure_supported].state.nunique(),'eligible_drugs':len(eligible),'eligible_states':eligible.state.tolist(),'proceed_comparative_gate':len(eligible)>=cfg['minimum_drugs_for_comparative_extension'],'minimum_drugs_gate':cfg['minimum_drugs_for_comparative_extension'],'effects_inspected_this_module':False,'dose_inconsistent_states':drug[(drug.dose_unique!=1)&(drug.state!='CTRL')].state.tolist(),'top10_supported_overlap':topcov[topcov.structure_supported.eq(True)].pert_iname.unique().tolist()}
 write_json(out/'feasibility_summary.json',result)
 manifest_paths=[ROOT/'config/expansion_feasibility_v1.json',ROOT/'outputs/preparation/GSE217421/full_design.tsv',ROOT/'outputs/analysis/mechanism/cardiac_resource_sample_metadata.tsv',out/'all_drug_feasibility.tsv',out/'compound_overlap_identity.tsv',out/'all_samples_metadata_audit.tsv']
 write_json(out/'coverage_freeze.json',{'frozen_utc':now(),'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in manifest_paths]})
 return result

if __name__=='__main__':run_standard_module('31_expansion_coverage',main)
