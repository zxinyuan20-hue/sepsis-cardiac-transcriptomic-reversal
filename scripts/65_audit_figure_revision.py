"""Packages pandas,numpy,Pillow; SVG XML parsing with standard library.
Input/output: config/figure_revision_v2.json; frozen figures and rat result tables.
Reproducibility: shared seed/logger and SHA256 manifests; no inferential testing.
Pipeline: prior integrity -> source reconciliation -> export QA -> revision release.
Outputs: visual/export audit, user-confirmed metadata addendum and new figure pointer.
"""
import json,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
from PIL import Image
from pipeline_utils import ROOT,run_standard_module,write_json,sha256

def main():
 cfg=json.loads((ROOT/'config/figure_revision_v2.json').read_text(encoding='utf8'));out=ROOT/cfg['out'];prev=ROOT/cfg['previous']
 checks=[]
 def check(name,ok):
  assert ok,name
  checks.append(name)
 for manifest in [prev/'draft_manifest.json',ROOT/'outputs/audit/external_compound_release_manifest.json']:
  for row in json.loads(manifest.read_text(encoding='utf8')):check('unchanged:'+row['path'],sha256(ROOT/row['path'])==row['sha256'])
 for path,digest in json.loads((out/'input_sources.json').read_text(encoding='utf8')).items():check('input:'+path,sha256(ROOT/path)==digest)
 exports=[]
 for f in sorted((out/'figures').glob('*.png')):
  record={'figure':f.stem}
  for ext in ['pdf','svg','png','tiff']:check('exists:'+f.stem+ext,f.with_suffix('.'+ext).stat().st_size>1000)
  for ext in ['png','tiff']:
   with Image.open(f.with_suffix('.'+ext)) as im:
    check('600dpi:'+f.stem+ext,all(abs(float(x)-600)<1 for x in im.info['dpi']))
    check('183mm width:'+f.stem+ext,im.width==4320)
    record[ext]={'pixels':im.size,'dpi':tuple(float(x) for x in im.info['dpi'])}
  svg=ET.parse(f.with_suffix('.svg'));ns={'s':'http://www.w3.org/2000/svg'}
  check('editable labels:'+f.stem,len(svg.findall('.//s:text',ns))>5)
  check('PDF header:'+f.stem,f.with_suffix('.pdf').read_bytes().startswith(b'%PDF'))
  exports.append(record)
 check('six figures',len(exports)==6)
 src=out/'source_data';read=lambda p:pd.read_csv(p,sep='\t')
 for name,original in [('Figure2_LINCS.tsv','screening/complete_screening_results.tsv'),('Figure2_cross_context.tsv','expansion/drug_class_context.tsv'),('figure2_source.tsv','external_compound_v1/all_drug_results.tsv')]:
  pd.testing.assert_frame_equal(read(src/name),read(ROOT/'outputs/analysis'/original),check_exact=False,rtol=1e-12,atol=1e-14)
  check('all values retained:'+name,True)
 a=read(src/'Figure2_cardiac.tsv');b=read(ROOT/'outputs/analysis/expansion/expanded_endpoint_results.tsv');b=b[b['size'].eq(100)]
 pd.testing.assert_frame_equal(a.sort_values('drug_state').reset_index(drop=True),b.sort_values('drug_state').reset_index(drop=True),check_exact=False,rtol=1e-12,atol=1e-14)
 check('all 21 cardiac intervals preserved',len(a)==21)
 rat=read(src/'FigureS2_rat_all_genes.tsv');sig=rat['adj.P.Val']<.05
 check('rat retained genes',len(rat)==14894);check('rat up/down counts',(int((sig&(rat.logFC>0)).sum()),int((sig&(rat.logFC<0)).sum()))==(1694,1383))
 orth=read(src/'FigureS2_human_rat_orthologs.tsv');check('shared ortholog genes',len(orth)==10822)
 check('saved rho agrees',np.isclose(orth.logFC_human.corr(orth.logFC_rat,method='spearman'),.06194364825493577))
 for co in ['GSE237861','GSE141864']:
  check('external gene source identical:'+co,sha256(src/(co+'_gene_concordance_source.tsv'))==sha256(prev/'source_data'/(co+'_gene_concordance_source.tsv')))
 write_json(out/'rat_metadata_addendum.json',{'date':cfg['date'],'source':'Direct user statement in this conversation','confirmed':{'species':'Rattus norvegicus','material':'rat cardiomyocytes','LPS_protocol':'standardized, per investigator','count_workflow':'scientifically standardized and consistent with similar studies, per investigator','independent_samples_per_group':6},'specific_parameters_not_yet_supplied':['primary-cell or cell-line identity and source','LPS concentration and duration','actual counting software/parameters and reference','experimental and library batches'],'interpretation':'Method details pending is not a finding of nonstandard experimental practice. No biological results changed.'})
 audit={'status':'PASS','checks':checks,'exports':exports,'visual_review':{'reviewed':['Figure1','Figure2','Figure3','Figure4','FigureS1','FigureS2'],'method':'Viewed complete exported PNGs; re-viewed final Figure1 and S2 after layout repairs','repairs':['Removed crossing flow connector','Shortened S2 y-axis label to clear panel label','Kept all 21 drugs and all rat genes visible'],'scope':'Visual and numerical-source consistency; no journal-specific certification'},'new_inferential_tests':False,'prior_manuscript_status':'Methods/Results v1 unchanged; Figure S2 citation to integrate in next manuscript revision'}
 write_json(out/'qa/revision_checks.json',audit)
 write_json(ROOT/'FIGURE_REVISION.json',{'version':cfg['version'],'date':cfg['date'],'scientific_release':'external-compound-1.0','status':'INTERNAL_FIGURE_QA_PASS','entry':cfg['out']+'/阅读入口.md','manifest':cfg['out']+'/revision_manifest.json','manuscript_base':cfg['previous']})
 paths=[p for p in out.rglob('*') if p.is_file() and p.name!='revision_manifest.json']+[ROOT/'FIGURE_REVISION.json',ROOT/'config/figure_revision_v2.json',ROOT/'scripts/64_revise_publication_figures.py',ROOT/'scripts/65_audit_figure_revision.py']
 rows=[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in sorted(paths)];write_json(out/'revision_manifest.json',rows)
 return {'figures':6,'checks':len(checks),'manifest_files':len(rows),'new_inferential_tests':False,'status':'PASS'}
if __name__=='__main__':run_standard_module('65_audit_figure_revision',main)
