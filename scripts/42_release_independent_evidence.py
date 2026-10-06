"""Packages pandas/numpy; fixed curated review plus local source responses and stage41.
Outputs evidence/data/exposure matrices, citation QA, source excerpts, Chinese report.
No new drug inference; deterministic release and hashes; logs under outputs/logs.
"""
import json,re,hashlib,xml.etree.ElementTree as ET
from difflib import SequenceMatcher
import pandas as pd,numpy as np
from pipeline_utils import ROOT,write_json,sha256,now,run_standard_module
OUT=ROOT/'outputs/analysis/independent_evidence';AUDIT=ROOT/'outputs/audit'
def j(p):return json.loads(p.read_text(encoding='utf-8'))
def clean(s):return re.sub(r'[^a-z0-9]','',s.lower())
def parse_soft(p):
 rows=[]
 for b in p.read_text(encoding='utf-8').split('^SAMPLE = ')[1:]:
  r={'gsm':b.splitlines()[0].strip()};props=[]
  for l in b.splitlines():
   if not l.startswith('!Sample_') or ' = ' not in l:continue
   k,v=l.split(' = ',1)
   if k=='!Sample_title':r['title']=v
   if k=='!Sample_characteristics_ch1':
    props.append(v)
    if ': ' in v:
     a,z=v.split(': ',1);r[a]=z
  r['characteristics']=' | '.join(props);rows.append(r)
 return rows
def main():
 curated=j(ROOT/'config/independent_evidence_curated_v1.json');records={r['pmid']:r for r in j(OUT/'literature_records.json')}
 for r in j(OUT/'additional_literature_records.json'):records.setdefault(r['pmid'],r)
 # Verify every immutable source and source-provenance pair.
 provenance=[]
 for m in (OUT/'sources').glob('*.source.json'):
  r=j(m);p=m.with_name(m.name.removesuffix('.source.json'));assert sha256(p)==r['sha256'];provenance.append({'path':p.relative_to(ROOT).as_posix(),**r})
 write_json(OUT/'source_manifest.json',provenance)
 checks=[]
 for p in (OUT/'sources').glob('*_crossref.json'):
  pmid=p.name.split('_')[0];r=records[pmid];x=j(p)['message'];ratio=SequenceMatcher(None,clean(r['title']),clean(' '.join(x.get('title',[])))).ratio()
  assert r['doi'].lower()==x['DOI'].lower();assert ratio>.85,(pmid,ratio)
  checks.append({'pmid':pmid,'doi':r['doi'],'title_similarity':ratio,'status':'PASS'})
 pd.DataFrame(checks).to_csv(OUT/'citation_crossref_checks.tsv',sep='\t',index=False)
 # Preserve exact textual excerpts, with machine-readable identifiers and provenance.
 full=[]
 for p in sorted((OUT/'sources').glob('PMC*.xml')):
  t=ET.parse(p);paragraphs=[' '.join(''.join(e.itertext()).split()) for e in t.findall('.//body//p')]
  body='\n\n'.join(paragraphs);full.append({'file':p.name,'body_paragraphs':len(paragraphs),'scope':'full_text_body' if paragraphs else 'metadata_or_abstract_only'})
  if paragraphs:(OUT/(p.stem+'_body.txt')).write_text(body,encoding='utf-8')
 pd.DataFrame(full).to_csv(OUT/'fulltext_availability_audit.tsv',sep='\t',index=False)
 evidence=[];ris=[]
 for e in curated['evidence']:
  r=records[e['pmid']];e={**e,**{k:r[k] for k in ['title','doi','pmcid','year']},'pubmed_url':'https://pubmed.ncbi.nlm.nih.gov/'+e['pmid']+'/'};evidence.append(e)
  ris.extend(['TY  - JOUR','TI  - '+r['title'],'PY  - '+str(r['year']),'DO  - '+r['doi'],'AN  - '+r['pmid'],'UR  - '+e['pubmed_url'],'ER  - ',''])
 pd.DataFrame(evidence).to_csv(OUT/'independent_drug_evidence_matrix.tsv',sep='\t',index=False)
 (OUT/'reviewed_references.ris').write_text('\n'.join(ris),encoding='utf-8')
 pd.DataFrame(curated['datasets']).to_csv(OUT/'independent_dataset_qualification.tsv',sep='\t',index=False)
 # Sample-level audits: do not count FF/FFPE specimens as different people.
 samples=[]
 for acc,filename in [('GSE237861','GSE237861_samples.soft'),('GSE141864','GSE141864_samples_brief.soft'),('GSE314561','GSE314561_samples.soft'),('GSE289264','GSE289264_samples.soft')]:
  samples.extend({'accession':acc,**r} for r in parse_soft(OUT/'sources'/filename))
 s=pd.DataFrame(samples);s.to_csv(OUT/'independent_sample_metadata.tsv',sep='\t',index=False)
 heart=s[s.accession.eq('GSE237861')&s.get('tissue').eq('heart')];assert len(heart)==14;assert heart.title.str.startswith('sepsis').sum()==7;assert heart.title.str.startswith('control').sum()==7
 h2=s[s.accession.eq('GSE141864')&s.get('tissue').eq('Heart')];assert len(h2)==10 and h2['patients id nr'].nunique()==7
 ffpe=h2[h2['type of storage tissue'].eq('FFPE')];assert len(ffpe)==7 and ffpe['patients id nr'].nunique()==7
 assert ffpe.diagnosis.eq('Meningococcal septic shock (MSS)').sum()==5
 h3=s[s.accession.eq('GSE314561')];assert len(h3)==7 and h3.treatment.eq('high-risk serum').sum()==3
 d=s[s.accession.eq('GSE289264')];assert d.treatment.str.contains('Ponatinib',case=False,na=False).sum()==3;assert d.treatment.str.contains('Dasatinib',case=False,na=False).sum()==3;assert not d.treatment.str.contains('Regorafenib',case=False,na=False).any()
 old_gsm=set(pd.read_csv(ROOT/'outputs/preprocessing/human/cel_sample_manifest.tsv',sep='\t').gsm);assert not old_gsm.intersection(set(heart.gsm)|set(h2.gsm))
 # Generic-label fallback resolved brand queries; label version explicitly retained.
 drugmeta=pd.read_csv(ROOT/'outputs/analysis/expansion/all_drug_metadata.tsv',sep='\t');exposure=[]
 for drug,state in [('dasatinib','DAS'),('ponatinib','PON'),('regorafenib','REG')]:
  p=OUT/'sources'/f'{drug}_fda_generic.json';lab=j(p)['results'][0]
  assert any(drug in v.lower() for v in lab.get('openfda',{}).get('generic_name',[]));assert lab['effective_time']<='20261002'
  row=drugmeta[drugmeta.state.eq(state)].iloc[0]
  rec={'drug':drug,'current_culture_concentration':row.doses,'current_culture_hours':48,'label_brand':';'.join(lab['openfda'].get('brand_name',[])),'label_effective_date':lab.get('effective_time'),'label_set_id':lab['set_id'],'label_url':'https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid='+lab['set_id'],'boxed_warning':' '.join(lab.get('boxed_warning',[])),'mechanism_of_action':' '.join(lab.get('mechanism_of_action',[])),'pharmacokinetics':' '.join(lab.get('pharmacokinetics',[])),'warnings_and_cautions':' '.join(lab.get('warnings_and_cautions',[])),'dose_bridge':'Nominal culture concentration is not equivalent to free myocardial clinical exposure; protein binding, metabolites and time differ'}
  exposure.append(rec)
  (OUT/f'{drug}_official_label_sections.txt').write_text('\n\n'.join(k+'\n'+str(v) for k,v in rec.items()),encoding='utf-8')
 pd.DataFrame(exposure).to_csv(OUT/'exposure_safety_matrix.tsv',sep='\t',index=False)
 # Recompute source rows for all selected gene-context entries.
 context=pd.read_csv(OUT/'literature_fixed_gene_context.tsv',sep='\t');assert len(context)==224
 for cohort,g in context.groupby('cohort'):
  if cohort=='human':f=ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv';key='entrez_id'
  elif cohort in ['DAS','PON','REG']:f=ROOT/f'outputs/analysis/expansion/{cohort}_gene_results.tsv';key='entrez_id'
  else:f=ROOT/f'outputs/analysis/mechanism/{cohort}_all_gene_results.tsv';key='human_entrez'
  source=pd.read_csv(f,sep='\t',dtype={key:str})
  for r in g.itertuples():
   if not r.available:continue
   z=source[source[key].eq(str(int(r.human_entrez)))];assert len(z)==1 and np.isclose(z.logFC.iloc[0],r.logFC) and np.isclose(z['adj.P.Val'].iloc[0],r.source_genome_FDR)
 sf=j(ROOT/'outputs/analysis/screening/screening_result_freeze.json')
 for filename,key in [('complete_screening_results.tsv','complete_results_sha256'),('supported_shortlist.tsv','supported_shortlist_sha256'),('exploratory_top10.tsv','exploratory_shortlist_sha256')]:assert sha256(ROOT/'outputs/analysis/screening'/filename)==sf[key]
 for item in j(AUDIT/'specificity_release_manifest.json'):
  if item['path'].startswith(('outputs/analysis/specificity/','outputs/figures/specificity/')):assert sha256(ROOT/item['path'])==item['sha256']
 summary={'version':'independent-evidence-1.0','pubmed_unique_initial':len(j(OUT/'literature_records.json')),'reviewed_evidence_rows':len(evidence),'dataset_qualifications':len(curated['datasets']),'source_responses_hash_verified':len(provenance),'crossref_verified':len(checks),'fixed_gene_context_rows':len(context),'new_external_validation_performed':False,'supported_therapeutic_candidates':0,'decision':'PRIORITIZE_INDEPENDENT_HUMAN_DISEASE_REPLICATION_AND_TOXICITY_DISCRIMINATION','drug_decisions':{'Dasatinib':'Retain conditional immune/cardiac compatibility hypothesis; prior SCM prediction limits novelty; dose-dependent harm','Ponatinib':'Retain as infection-context benefit versus cardiac-toxicity comparison, not preferred therapeutic lead','Regorafenib':'Deprioritize standalone cardioprotective narrative; indirect lung benefit and direct cardiac toxicity conflict'},'new_disease_data':'GSE237861 priority; GSE141864 storage-matched sensitivity','drug_resource':'GSE289264 exposure/context comparison only, one RNA-seq line','limits':['Capped targeted searches, not systematic review','Ponatinib infection-benefit paper full-text body unavailable','No new disease-cohort statistics run','RNA context cannot prove kinase activity or therapeutic benefit']}
 write_json(OUT/'independent_evidence_summary.json',summary)
 report=AUDIT/'三药物独立证据与后续研究决策.md'
 # Main report is maintained as a readable source artifact alongside this release.
 assert report.exists()
 current=ROOT/'CURRENT_RELEASE.json';snapshot=AUDIT/'pre_independent_evidence_release_snapshot.json'
 if not snapshot.exists():snapshot.write_bytes(current.read_bytes())
 status=j(current);preserve_external=bool(status.get('external_human_version'));status.update(status='INDEPENDENT_EVIDENCE_REVIEW_COMPLETE_DISEASE_REPLICATION_PRIORITY_NO_THERAPEUTIC_ADVANCEMENT',updated_utc=now(),independent_evidence_version=summary['version'],continuation_decision=summary['decision'],independent_evidence_report=report.relative_to(ROOT).as_posix(),report=report.relative_to(ROOT).as_posix(),supported_candidates=0)
 for name in ['39_independent_evidence_search','40_independent_source_details','41_independent_mechanism_context']:
  assert j(AUDIT/(name+'_status.json'))['status']=='PASS';status['module_checks'][name]='PASS'
 if not preserve_external:write_json(current,status)
 qa={'source_hashes':'PASS','citation_DOI_and_title':'PASS','human_heart_sample_counts':'PASS','duplicate_storage_patient_detection':'PASS','GSE289264_drug_coverage':'PASS','selected_gene_values_reproduced':'PASS','original_screening_and_specificity_numerical_results':'UNCHANGED','retracted_study':'EXCLUDED','claim_limit':'NO_NEW_EXTERNAL_VALIDATION_OR_EFFICACY','search_errors':'3brandFDA404 resolved via generic labels; PMC bodies selectively recovered; PMID32492702 body remains unavailable'}
 write_json(AUDIT/'independent_evidence_qa.json',qa)
 paths=[p for p in OUT.rglob('*') if p.is_file()]+[report,current,ROOT/'README.md',AUDIT/'independent_evidence_qa.json']
 for folder in ['scripts','docs','config']:paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
 write_json(AUDIT/'independent_evidence_release_manifest.json',[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p),'bytes':p.stat().st_size} for p in paths])
 return summary
if __name__=='__main__':run_standard_module('42_release_independent_evidence',main)
