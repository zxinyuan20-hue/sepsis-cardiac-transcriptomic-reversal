"""Packages pandas numpy scipy matplotlib h5py; frozen protocol from module 97.
Inputs public-reference GET responses and query map, outputs donor-level tables,
descriptive localization, expression-scale check, manuscript text and Figure S4.
No cell-level significance test, no cell composition or mechanistic inference.
"""
import json,logging,importlib.util
import numpy as np,pandas as pd,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT,now,sha256,write_json,run_standard_module
from cellxgene_read_utils import decode
BIO=ROOT/'outputs/analysis/biology_extension_v1';CELL=BIO/'cell_localization';SRC=BIO/'sources';PUB=ROOT/'outputs/publication/evidence_hierarchy_revision_v1'
def main():
 proto=json.loads((CELL/'protocol_freeze.json').read_text());acq=json.loads((CELL/'acquisition.json').read_text());obs=pd.read_csv(CELL/'reference_annotations.tsv',sep='\t');var=pd.read_csv(CELL/'reference_genes.tsv',sep='\t').symbol;q=pd.read_csv(CELL/'query_mapping.tsv',sep='\t')
 assert sha256(CELL/'query_mapping.tsv')==proto['query_hash'] and sha256(CELL/'reference_annotations.tsv')==proto['metadata_hash']
 keep=(obs.source=='Nuclei')&~obs.cell_type.isin(['NotAssigned','doublets']);counts=obs[keep].groupby(['cell_type','donor']).size().rename('nuclei').reset_index();counts['eligible_stratum']=counts.nuclei>=20
 n=counts[counts.eligible_stratum].groupby('cell_type').size();types=n[n>=7].index.tolist();counts['eligible_type']=counts.cell_type.isin(types);counts.to_csv(CELL/'donor_cell_counts.tsv',sep='\t',index=False)
 groups={}
 for r in counts[counts.eligible_stratum&counts.eligible_type].itertuples():groups[(r.cell_type,r.donor)]=np.flatnonzero((keep&(obs.cell_type==r.cell_type)&(obs.donor==r.donor)).to_numpy())
 rows=[];inputhash={};sample_api={}
 for name in acq['files']:
  p=SRC/name;inputhash[name]=sha256(p);d=decode(p.read_bytes());assert len(d)==len(obs)
  for col in d:
   sample_api[int(col)]=d[col].to_numpy()[:10]
   values=d[col].to_numpy();assert np.isfinite(values).all() and values.min()>=0
   for (ct,donor),ids in groups.items():rows.append({'gene':var.iloc[col],'cell_type':ct,'donor':donor,'nuclei':len(ids),'mean_expression':float(np.mean(values[ids],dtype=np.float64)),'detection_fraction':float(np.mean(values[ids]>0))})
 donor=pd.DataFrame(rows);donor.to_csv(CELL/'donor_gene_expression.tsv',sep='\t',index=False)
 mean=donor.groupby(['gene','cell_type'])[['mean_expression','detection_fraction']].mean().reset_index();mean.to_csv(CELL/'equal_donor_gene_expression.tsv',sep='\t',index=False)
 profile=[];pergene=[]
 for metric in ['mean_expression','detection_fraction']:
  m=mean.pivot(index='gene',columns='cell_type',values=metric);assert not m.isna().any().any();total=m.sum(axis=1);fr=m.loc[total>0].div(total[total>0],axis=0)
  assert np.allclose(fr.sum(axis=1),1)
  for dr in ['up','down']:
   genes=q.loc[(q.direction==dr)&q.reference_present,'SYMBOL'];genes=genes[genes.isin(fr.index)];assert len(genes)==len(set(genes));v=fr.loc[genes].mean()
   for ct,x in v.items():profile.append({'metric':metric,'direction':dr,'cell_type':ct,'localization_fraction':x,'genes':len(genes),'donors':int(n[ct])})
   pg=fr.loc[genes].reset_index().melt(id_vars='gene',var_name='cell_type',value_name='fraction');pg['direction']=dr;pg['metric']=metric;pergene.append(pg)
 result=pd.DataFrame(profile);result.to_csv(CELL/'localization_profiles.tsv',sep='\t',index=False);pd.concat(pergene).to_csv(CELL/'gene_localization.tsv',sep='\t',index=False)
 # Independently recompute equal-donor means with explicit arithmetic for selected genes.
 err=0
 for (gene,ct),g in donor.groupby(['gene','cell_type']):
  expected=sum(g.mean_expression.tolist())/len(g);actual=mean[(mean.gene==gene)&(mean.cell_type==ct)].mean_expression.iloc[0];err=max(err,abs(expected-actual))
 # Expression scale audit against raw counts for ACTB in first ten cells.
 spec=importlib.util.spec_from_file_location('atlas',ROOT/'scripts/95_inspect_atlas_h5.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 api=decode((SRC/'cellxgene_ACTB.bin').read_bytes());actb=int(api.columns[0]);assert var.iloc[actb]=='ACTB'
 with mod.RangeFile() as f:
  with h5py.File(f,'r') as h:
   ptr=h['X/indptr'][:11];ind=h['X/indices'][ptr[0]:ptr[-1]];val=h['X/data'][ptr[0]:ptr[-1]]
 ratios=[];native=[];query_scale_errors=[]
 for i in range(10):
  ii=ind[ptr[i]:ptr[i+1]];vv=val[ptr[i]:ptr[i+1]];x=float(vv[ii==actb].sum());tot=float(vv.sum());y=float(api.iloc[i,0]);native.append({'cell_row':i,'ACTB_count':x,'total_count':tot,'viewer_expression':y})
  if x>0:ratios.append(float(np.expm1(y)*tot/x))
  for col,values in sample_api.items():
   expected=np.log1p(float(vv[ii==col].sum())*1959/tot);query_scale_errors.append(abs(expected-float(values[i])))
 assert max(query_scale_errors)<1e-5
 scale={'sampled_entries':native,'implied_log1p_normalization_targets':ratios,'matches_log1p_CPTT':bool(np.allclose(ratios,10000,rtol=1e-4)),'matches_log1p_library_target_1959':bool(np.allclose(ratios,1959,rtol=1e-5)),'query_gene_cell_entries_checked':len(query_scale_errors),'query_scale_max_error':max(query_scale_errors),'note':'sampled diagnostic supports log1p library normalization to 1959; official viewer expression is the analyzed scale'};write_json(CELL/'expression_scale_check.json',scale)
 def pretty(ct):return ct.replace('_',' ').replace('Cardiomyocyte','cardiomyocytes').replace('Smooth muscle cells','smooth muscle cells').lower()
 mainres=result[result.metric=='mean_expression'];tops={dr:mainres[mainres.direction==dr].sort_values('localization_fraction',ascending=False).head(3) for dr in ['up','down']}
 def label(dr):return ', '.join(f'{pretty(r.cell_type)} ({r.localization_fraction*100:.1f}%)' for r in tops[dr].itertuples())
 coverage=q.groupby('direction').reference_present.sum().astype(int).to_dict();ngenes={dr:int(mainres[mainres.direction==dr].genes.iloc[0]) for dr in ['up','down']};N=int(keep.sum());included=int(counts.loc[counts.eligible_stratum&counts.eligible_type,'nuclei'].sum())
 resulttext=f'Adult cardiac reference mapping covered {coverage["up"]}/100 UP and {coverage["down"]}/90 DOWN genes. After donor and nuclei coverage requirements, {len(types)} cardiac cell types were eligible, with {int(n.loc[types].min())}–{int(n.loc[types].max())} donors per type. The largest mean relative localization fractions for the UP query were in {label("up")}; for the DOWN query they were in {label("down")}. These fractions summarize relative reference expression, not cell abundance or statistically significant enrichment (Figure S4; Supplementary Note 4).'
 resulttext=resulttext.replace('with 14–14 donors per type','with 14 donors per type')
 resulttext+=f' Positive reference expression was available for {ngenes["up"]} UP and {ngenes["down"]} DOWN genes. Relative to UP, the DOWN query had larger atrial and ventricular cardiomyocyte fractions, whereas UP had larger endothelial and myeloid fractions; these contrasts retained their direction in the detection-frequency sensitivity. Adipocytes had the largest DOWN fraction, and the leading UP type changed to adipocytes under the detection summary.'
 discussion='Adult cardiac reference localization placed the UP query relatively more in endothelial and myeloid expression profiles and the DOWN query relatively more in cardiomyocyte profiles when the two query components were compared. These contrasts persisted in the detection-frequency sensitivity. The distributions were broad: adipocytes had the highest DOWN fraction, and the leading UP type depended on the expression summary. Thus, the query components are not exclusive markers of immune cells or cardiomyocytes. Equal donor weighting avoids treating large cell collections as independent biological replication. Healthy-reference localization cannot distinguish cell-intrinsic sepsis responses from disease-associated composition changes, and broad metabolic genes need not be cardiomyocyte-specific. The analysis supplies reference-based biological annotation rather than single-cell validation of disease or protection.'
 note=f'''The Litviňuková et al. 2020 atlas supplied original global_raw.h5ad annotations and expression from its official global viewer. All 486,134 reference cell identities were matched one-to-one after combining each raw identifier with its documented cell-source suffix. This was an identifier-format difference, not evidence of missing or duplicated cells. Analyses selected nuclei and excluded NotAssigned and doublets, leaving {N:,} eligible nuclei before donor-stratum requirements. Donor–cell-type strata required at least 20 nuclei; cell types required at least seven such donors. The analyzed strata contained {included:,} nuclei in {len(types)} types, with {int(n.loc[types].min())}–{int(n.loc[types].max())} donors per type. Counts for all strata, including excluded ones, are reported.

Reference gene coverage was {coverage['up']}/100 UP and {coverage['down']}/90 DOWN. Ambiguous gene-symbol mappings or absent genes were excluded without replacement. There were {ngenes['up']} UP and {ngenes['down']} DOWN genes with positive total reference expression in the eligible types. The original query is restricted to measured LINCS landmarks, so these profiles do not characterize the entire human myocardial transcriptome. Optional sepsis-dominant and shared-injury modules were not introduced because no separately frozen module definition was used.

For each gene and cell type, we first averaged the official viewer expression within each eligible donor, then averaged donor means equally. Each gene's nonnegative cell-type means were divided by their sum across eligible cell types, giving relative localization fractions summing to one. UP and DOWN profiles average these fractions with equal gene weights. Detection frequency (>0 expression) was summarized identically as a sensitivity. No whole-genome enrichment background, random-gene null, P value, cell-level test or cell-proportion estimate was used. The analyzed universe is the frozen query intersected with reference genes, and the denominator is the eligible cardiac cell-type expression profile. Fractions are not estimated transcript fractions or cell fractions. Unequal availability of donor strata and differences in nuclei capture remain limitations.

The official viewer uses nonnegative transformed expression. A raw-count diagnostic sampled ACTB in the first ten cells and recorded implied library-normalization factors; the complete check is supplied as expression_scale_check.json. Results are described on the official viewer scale rather than assuming a raw-count interpretation. An independent explicit arithmetic reconstruction of equal-donor mean expression agreed within {err:.2g}. Source hashes and exact row-identity checks accompany the result. These healthy-reference profiles cannot identify which cell types changed in sepsis or establish a causal cellular mechanism.'''
 legend=f'Figure S4. Reference-based localization of the frozen human query in the adult heart atlas. (a) Equal-gene averages of relative cell-type expression fractions for the UP and DOWN queries. Expression is averaged within donor and then equally across donors before within-gene normalization across cell types. (b) Detection-frequency sensitivity using the same donor and gene weighting. Labels give eligible donor counts; each donor–cell-type stratum contains at least 20 nuclei and each displayed cell type has at least seven donors. Reference coverage is {coverage["up"]}/100 UP and {coverage["down"]}/90 DOWN genes. Fractions are relative localization summaries, not cell proportions, enrichment probabilities or sepsis effects. The dashed line is 1/{len(types)}, the equal-distribution reference across displayed types; it is not a significance threshold. No inferential P values or cell-level tests are shown.'
 note+='\n\nThe processed reference contained 486,134 rows, whereas the article reports 487,106 cells/nuclei; analyses therefore refer to the specifically archived downloadable/interactive reference version, not a reconstructed count identical to the publication. Exact symbol matching left GFUS (UP) and DIPK1A (DOWN) unmatched; this is a mapping limitation, not evidence of absent expression. CHERP (UP) and ALDOC (DOWN) had zero expression across the retained reference nuclei and were omitted from within-gene fractions. Mesothelial cells had only six qualifying donors and were excluded by the seven-donor gate. All ten retained types had 14 donors. Source-scale checking across 1,880 gene-by-cell entries was consistent with log1p library normalization to 1,959 counts; maximum absolute discrepancy was '+f'{max(query_scale_errors):.2g}. This sampled check does not independently reconstruct the entire reference preprocessing.\n\n'+'Relative localization was distributed across multiple types and depended partly on whether mean expression or detection frequency was summarized. Adipocytes led the DOWN profile under both summaries and the UP detection profile. The main text therefore reports comparative UP-versus-DOWN tendencies rather than exclusive cell-type assignment or statistically significant enrichment.'
 legend+=f' Positive-expression denominators were {ngenes["up"]} UP and {ngenes["down"]} DOWN genes after excluding one zero-expression gene in each query.'
 summary={'recorded_utc':now(),'reference_rows':len(obs),'nuclei_before_stratum_filter':N,'nuclei_included':included,'eligible_types':types,'query_coverage':coverage,'nonzero_genes':ngenes,'donor_range':[int(n.loc[types].min()),int(n.loc[types].max())],'equal_donor_arithmetic_max_error':err,'manuscript_results':resulttext,'manuscript_discussion':discussion,'supplement_methods':note,'figure_legend':legend,'input_hashes':inputhash,'formal_inference':False};write_json(CELL/'summary.json',summary)
 logging.getLogger('fontTools').setLevel(logging.WARNING);plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(1,2,figsize=(7.2,4.6),sharey=True);fig.subplots_adjust(left=.30,right=.98,bottom=.22,top=.85,wspace=.27);order=types
 for ax,metric,title in zip(axes,['mean_expression','detection_fraction'],['a  Expression localization','b  Detection sensitivity']):
  for dr,c,mk,off in [('up','#BE7552','o',-.10),('down','#326C87','D',.10)]:
   z=result[(result.metric==metric)&(result.direction==dr)].set_index('cell_type').loc[order];ax.scatter(z.localization_fraction*100,np.arange(len(order))+off,c=c,marker=mk,s=22,label=dr.upper())
  ax.axvline(100/len(order),c='#9AA3AA',ls='--',lw=.7);ax.set_xlabel('Mean relative localization (%)');ax.set_title(title,loc='left',fontsize=8,fontweight='bold');ax.grid(axis='x',color='#E8EBED',lw=.5);ax.set_xlim(left=0)
 axes[0].set_yticks(range(len(order)),[f'{pretty(ct).capitalize()} (n={int(n[ct])})' for ct in order],fontsize=6.7);axes[0].invert_yaxis();axes[1].legend(frameon=False,loc='best',fontsize=6.5)
 fig.text(.05,.055,f'Healthy adult heart reference; {included:,} nuclei summarized at donor level.\nMapped: UP {coverage["up"]}/100, DOWN {coverage["down"]}/90; positive-expression profiles: {ngenes["up"]} UP, {ngenes["down"]} DOWN.\nDashed line: equal distribution across {len(types)} types. Descriptive localization, not cell proportions.',fontsize=6.7,linespacing=1.4)
 for ext in ['pdf','svg','png','tiff']:
  kw={'dpi':600,'facecolor':'white'}
  if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(PUB/'figures'/f'FigureS4_cell_localization.{ext}',**kw)
 plt.close(fig)
 for name in ['localization_profiles.tsv','donor_cell_counts.tsv','query_mapping.tsv','equal_donor_gene_expression.tsv']:
  import shutil;shutil.copy2(CELL/name,PUB/'source_data'/name)
 return {'eligible_types':len(types),'coverage':coverage,'nuclei':included,'normalization_check':scale['matches_log1p_CPTT']}
if __name__=='__main__':run_standard_module('100_cell_localization',main)
