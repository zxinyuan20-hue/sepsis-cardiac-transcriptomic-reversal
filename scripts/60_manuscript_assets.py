"""Packages numpy/pandas/scipy/matplotlib/Pillow. Inputs frozen result tables.
Central output/export settings: config/manuscript_assets_v1.json. Shared seeded log.
Outputs Figure1-4/S1, source tables, editable table data, captions, integrity checks.
No new biological tests. Existing Figure3/4/S1 exports reused without alteration.
"""
import json,shutil,logging
import numpy as np,pandas as pd
from scipy.stats import spearmanr,false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pipeline_utils import ROOT,run_standard_module,write_json,sha256
CFG=json.loads((ROOT/'config/manuscript_assets_v1.json').read_text(encoding='utf8'))
OUT=ROOT/CFG['out']; FIG=OUT/'figures'; SRC=OUT/'source_data'; INPUTS={}
def read(rel):
 p=ROOT/rel;INPUTS[rel]=sha256(p);return pd.read_csv(p,sep='\t')
def save(fig,name):
 for ext in ['pdf','svg','png','tiff']:fig.savefig(FIG/(name+'.'+ext),dpi=CFG['dpi'],facecolor='white')
 plt.close(fig)
def table_md(df):
 return '| '+' | '.join(df.columns)+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))

def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 for folder in [OUT,FIG,SRC,OUT/'tables',OUT/'qa']:folder.mkdir(parents=True,exist_ok=True)
 for manifest in ['outputs/audit/external_compound_release_manifest.json','outputs/publication/project_review_v1/review_manifest.json']:
  for r in json.loads((ROOT/manifest).read_text(encoding='utf8')):assert sha256(ROOT/r['path'])==r['sha256'],r['path']
 plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'axes.linewidth':.65})
 # Workflow separates independent disease data from one reused drug experiment.
 fig,ax=plt.subplots(figsize=(7.2,6.5));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 def box(x,y,w,h,title,body,color='#eef3f5'):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.008',facecolor=color,edgecolor='#81949e',lw=.7))
  ax.text(x+w/2,y+h-.028,title,ha='center',va='top',fontsize=8,fontweight='bold')
  ax.text(x+w/2,y+h-.070,body,ha='center',va='top',fontsize=7,linespacing=1.45)
 def arrow(x1,y1,x2,y2):ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'->','color':'#697982','lw':1})
 box(.04,.76,.42,.19,'Human myocardial discovery','GSE79962: 51 individuals\nSepsis 20; nonfailing 11; IHD 11; DCM 9\nFrozen query: 100 UP / 90 DOWN')
 box(.55,.76,.41,.19,'Independent disease contexts','GSE237861: 7 sepsis + 7 controls\nGSE141864: 5 sepsis + 2 controls\nPatients remain the analysis unit')
 box(.04,.51,.42,.17,'LINCS screening','482 eligible compound entries\n6,769 perturbation signatures\nMatched random-set null; BH across 482')
 box(.55,.51,.41,.17,'Cross-model disease support','Two mouse studies: 5 + 5 each\nRat cardiomyocytes: 6 control + 6 LPS\nDisease association; no drug rescue','#f4f1ed')
 box(.04,.255,.42,.175,'Cardiac drug reassessment','GSE217421: 21 eligible drugs\n4–6 donors per drug; one dose; 48 h\nMatched culture blocks, then donor means')
 box(.55,.255,.41,.175,'External disease-context reassessment','783 shared genes; all 21 drugs\nUP suppression + DOWN restoration\nDonor IUT; BH across 21')
 box(.14,.04,.72,.12,'Evidence interpretation','Separate direction, statistical support and mechanism\nPatient / donor omission analyses; no efficacy claim','#f4f1ed')
 arrow(.25,.75,.25,.69);arrow(.25,.50,.25,.44);arrow(.47,.34,.535,.34)
 arrow(.76,.75,.98,.72);arrow(.98,.72,.98,.45);arrow(.98,.45,.78,.44)
 arrow(.48,.84,.525,.84);ax.text(.49,.87,'query',fontsize=6,ha='center')
 arrow(.25,.245,.32,.17);arrow(.76,.245,.67,.17)
 ax.text(.5,.205,'Same drug experiment and donors reused across disease contexts',ha='center',fontsize=7,color='#465861')
 fig.subplots_adjust(left=.025,right=.975,bottom=.03,top=.98);save(fig,'Figure1_workflow')
 screen=read('outputs/analysis/screening/complete_screening_results.tsv');ctx=read('outputs/analysis/expansion/drug_class_context.tsv')
 ep=read('outputs/analysis/expansion/expanded_endpoint_results.tsv'); ep=ep[ep['size'].eq(100)]
 assert len(screen)==482 and len(ctx)==21 and len(ep)==21
 assert (screen.empirical_bh_q>=.05).all() and (ep.family_BH_q>=.05).all()
 assert np.allclose(false_discovery_control(screen.empirical_p),screen.empirical_bh_q)
 assert np.allclose(false_discovery_control(ep.exact_one_sided_p),ep.family_BH_q)
 screen.to_csv(SRC/'Figure2_LINCS.tsv',sep='\t',index=False);ctx.to_csv(SRC/'Figure2_cross_context.tsv',sep='\t',index=False);ep.to_csv(SRC/'Figure2_cardiac.tsv',sep='\t',index=False)
 rho=spearmanr(ctx.harmonized_lincs_score,ctx.harmonized_cardiac_score).statistic
 ki=ctx[ctx.broad_class.eq('kinase_inhibitors')]; krho=spearmanr(ki.harmonized_lincs_score,ki.harmonized_cardiac_score).statistic
 assert abs(rho-.4941634888814772)<1e-10 and abs(krho-.2008979233214592)<1e-10
 fig=plt.figure(figsize=(7.2,6.2));gs=fig.add_gridspec(2,2,hspace=.66,wspace=.47);axes=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]
 a,b,c,d=axes
 a.scatter(screen['rank'],screen.reversal_score,s=5,color='#697982',edgecolors='none');a.scatter(screen.loc[screen.reversal_score>0,'rank'],screen.loc[screen.reversal_score>0,'reversal_score'],s=8,color='#276b8c',edgecolors='none')
 a.axhline(0,color='#999999',ls='--',lw=.7);a.set(xlabel='Screening rank (482 compound entries)',ylabel='LINCS reversal score',title='LINCS screening')
 a.text(.97,.96,'23 positive scores\n0 BH-supported entries',ha='right',va='top',transform=a.transAxes,fontsize=7)
 x=ep.sort_values('mean_reversal').reset_index(drop=True)
 b.errorbar(np.arange(1,22),x.mean_reversal*100,yerr=np.c_[(x.mean_reversal-x.t_ci_low)*100,(x.t_ci_high-x.mean_reversal)*100].T,fmt='o',ms=3,lw=.7,color='#697982',capsize=1)
 b.scatter(np.flatnonzero(x.both_directions_opposing)+1,x.loc[x.both_directions_opposing,'mean_reversal']*100,s=14,color='#276b8c',zorder=4)
 b.axhline(0,color='#999999',ls='--',lw=.7);b.set(xlabel='Drug ordered by mean reversal (21 drugs)',ylabel='Cardiac reversal (rank percentage points)',title='Cardiac net reversal; donor 95% t intervals')
 b.set_xticks([1,5,10,15,21])
 b.text(.32,.15,'4–6 donors per drug\nBlue: both mean directions\noppose disease\n0 BH-supported drugs',transform=b.transAxes,fontsize=6)
 palette={'kinase_inhibitors':'#276b8c','proteasome_inhibitors':'#bd6940','anthracycline_topoisomerase':'#bd6940','DNA_methyltransferase_inhibitor':'#697982','PPAR_related_agonist':'#697982'}
 for group,ss in ctx.groupby('broad_class'):
  c.scatter(ss.harmonized_lincs_score,ss.harmonized_cardiac_score,s=22,color=palette.get(group,'#697982'),marker='o' if group=='kinase_inhibitors' else 's',edgecolor='white',lw=.35)
 for state in ['CAR','BOR','DOX','REG','VEM']:
  r=ctx[ctx.state.eq(state)].iloc[0];c.annotate(state,(r.harmonized_lincs_score,r.harmonized_cardiac_score),xytext=(3,3),textcoords='offset points',fontsize=6)
 c.axhline(0,color='#cccccc',lw=.6);c.axvline(0,color='#cccccc',lw=.6);c.set(xlabel='Harmonized LINCS score',ylabel='Harmonized cardiac score',title=f'All 21 drugs: descriptive Spearman ρ = {rho:.3f}')
 c.text(.97,.22,'Circles: kinase inhibitors\nSquares: other drug classes',ha='right',transform=c.transAxes,fontsize=6)
 d.scatter(ki.harmonized_lincs_score,ki.harmonized_cardiac_score,s=24,color='#276b8c',edgecolor='white',lw=.4)
 d.axhline(0,color='#cccccc',lw=.6);d.axvline(0,color='#cccccc',lw=.6);d.set(xlabel='Harmonized LINCS score',ylabel='Harmonized cardiac score',title=f'16 kinase inhibitors: descriptive ρ = {krho:.3f}')
 d.text(.04,.94,'Post-result class sensitivity',va='top',transform=d.transAxes,fontsize=6.5)
 for label,axis in zip('abcd',axes):axis.text(-.17,1.12,label,transform=axis.transAxes,fontweight='bold',fontsize=10)
 fig.subplots_adjust(left=.105,right=.98,bottom=.095,top=.92);save(fig,'Figure2_screening_context')
 # Preserve old exports byte-for-byte and record the manuscript numbering map.
 numbering={'Figure3_disease_replication':'figure1_disease_replication','Figure4_external_compounds':'figure2_all_drug_contexts','FigureS1_gene_pathway':'figureS1_gene_and_pathway_evidence'}
 old=ROOT/'outputs/analysis/external_compound_v1/figures'; reuse=[]
 for target,source in numbering.items():
  for ext in ['pdf','svg','png','tiff']:
   p=old/(source+'.'+ext);dest=FIG/(target+'.'+ext);shutil.copyfile(p,dest);assert sha256(p)==sha256(dest)
   reuse.append({'source':p.relative_to(ROOT).as_posix(),'target':dest.relative_to(ROOT).as_posix(),'sha256':sha256(p)})
 for p in old.glob('*.tsv'):shutil.copyfile(p,SRC/p.name)
 write_json(OUT/'figure_reuse_map.json',reuse)
 # Tables contain analytical datasets only; deferred resources are excluded.
 rows=[
 ['GSE79962','Human myocardium','20 sepsis; 11 nonfailing; 11 IHD; 9 DCM','Individual','Discovery; within-cohort disease comparisons'],
 ['GSE70138','LINCS Phase II','482 compound entries; 6,769 signatures','Compound entry; cell-line aggregation','Transcriptomic screening'],
 ['GSE217421','Healthy donor iPSC cardiac cultures','21 drugs; 4–6 donors/drug; 421 distinct records','Donor','Cardiac perturbation; reused across disease contexts'],
 ['GSE237861','Human myocardium','7 sepsis; 7 noninfectious critically ill controls','Patient','Primary external disease context'],
 ['GSE141864','Human FFPE myocardium','5 sepsis; 2 controls','Patient','Secondary external disease context'],
 ['GSE185754','Mouse myocardium','5 saline; 5 LPS, 24 h','Animal sample','Cross-model disease support'],
 ['GSE267388','Mouse myocardium','5 PBS; 5 LPS, 12 h','Animal sample','Held-out cross-model support'],
 ['Local rat dataset','Rat cardiomyocytes','6 control; 6 LPS','Independent biological sample','Exploratory cross-model support']]
 t1=pd.DataFrame(rows,columns=['Dataset','Material','Design','Unit','Role']);t1.to_csv(OUT/'tables/Table1_datasets.tsv',sep='\t',index=False)
 allr=read('outputs/analysis/external_compound_v1/all_drug_results.tsv');cross=read('outputs/analysis/external_compound_v1/cross_cohort_decision.tsv')
 dl=read('outputs/analysis/external_compound_v1/donor_loo.tsv');pl=read('outputs/analysis/external_compound_v1/patient_loo.tsv');exposure=read('outputs/analysis/external_compound_v1/drug_exposure_audit.tsv')
 def frac(df,drug,cohort):
  x=df[df.drug_state.eq(drug)&df.cohort.eq(cohort)];return f'{int(x.both_positive.sum())}/{len(x)}'
 rr=[]
 for row in cross.sort_values('drug_name').itertuples():
  a=allr[allr.drug_state.eq(row.drug_state)&allr.cohort.eq('GSE237861')].iloc[0]
  b=allr[allr.drug_state.eq(row.drug_state)&allr.cohort.eq('GSE141864')].iloc[0]
  rr.append([row.drug_name,str(row.donors),f'{a.iut_BH_q:.4f}',f'{b.iut_BH_q:.4f}',f'{row.secondary_BH_q:.4f}',frac(dl,row.drug_state,'GSE237861')+'; '+frac(dl,row.drug_state,'GSE141864'),frac(pl,row.drug_state,'GSE237861')+'; '+frac(pl,row.drug_state,'GSE141864')])
 t2=pd.DataFrame(rr,columns=['Drug','Donors','Primary q','Secondary q','Four-component q','Donor omissions A; B','Patient omissions A; B']);assert len(t2)==21
 t2.to_csv(OUT/'tables/Table2_compound_evidence.tsv',sep='\t',index=False)
 tnotes={
 'Table1':'IHD, ischemic heart disease; DCM, dilated cardiomyopathy; FFPE, formalin-fixed paraffin-embedded; LPS, lipopolysaccharide; PBS, phosphate-buffered saline. Counts are study-specific and must not be summed as a single homogeneous cohort. GSE141864 includes three additional frozen-tissue records from already represented patients; these were not added to the FFPE patient count. GSE217421 controls were reused across drugs; records and culture wells are not independent donors. Local rat data were exploratory and experimental metadata remain incomplete.',
 'Table2':'A = GSE237861; B = GSE141864. Primary and secondary q values are BH-adjusted across 21 drugs separately in A and B. Each drug test is the maximum of one-sided donor sign-flip P values for UP suppression and DOWN restoration. The four-component q uses the maximum across both components in both cohorts, with a separate secondary 21-drug BH family. Omission entries are the number of analyses retaining both positive mean components divided by all omissions, in A; B order. Donor omissions use cardiac donors; patient omissions refit disease profiles. All primary q values exceeded 0.05, and no drug met the escalation gate. Stability is descriptive and does not establish efficacy.'}
 for num,frame,title in [('Table1',t1,'Data resources, analytical units and study roles'),('Table2',t2,'Compound evidence and directional stability across external disease contexts')]:
  (OUT/'tables'/f'{num}.md').write_text(f'# {num}. {title}\n\n'+table_md(frame)+'\n\n'+tnotes[num]+'\n',encoding='utf8')
 write_json(OUT/'tables/table_definitions.json',{'Table1':{'title':'Data resources, analytical units and study roles','columns':list(t1.columns),'rows':rows,'note':tnotes['Table1']},'Table2':{'title':'Compound evidence and directional stability across external disease contexts','columns':list(t2.columns),'rows':rr,'note':tnotes['Table2']}})
 exposure.to_csv(SRC/'drug_exposure.tsv',sep='\t',index=False)
 oldcap=(old/'captions.md').read_text(encoding='utf8');INPUTS[(old/'captions.md').relative_to(ROOT).as_posix()]=sha256(old/'captions.md')
 oldcap=oldcap.replace('Figure 1 |','Figure 3 |').replace('Figure 2 |','Figure 4 |')
 captions='''# Figure legends

Figure 1 | Study design and evidence roles. The discovery cohort defined the 100-UP/90-DOWN query used in LINCS screening. Human external cohorts supplied independent disease contexts; cross-model disease data provided separate support. GSE217421 provided one healthy-donor cardiac drug experiment, reused throughout the cardiac and external disease-context analyses. Drug wells were averaged within matched culture blocks and then within donor. Arrows indicate analytical information flow, not a fully prospective sequence: cardiac, disease-specificity and external-context analyses were extensions after the original screening results were known. IUT, intersection–union test; BH, Benjamini–Hochberg; IHD, ischemic heart disease; DCM, dilated cardiomyopathy. LPS, lipopolysaccharide. No drug rescue in diseased cardiac cells was performed.

Figure 2 | Screening evidence and transfer across cell backgrounds. a, All 482 eligible LINCS compound entries ordered by their frozen reversal score. Blue points denote positive scores; the matched random-gene-set test with 10,000 draws and BH correction across 482 entries yielded no supported entry. Scores were medians within compound–cell combinations and then across cells; compound entries are not biological replicates. b, Original cardiac net-reversal endpoint for all 21 drugs, ordered by the mean; dots show equal-donor means and bars show two-sided Student-t 95% confidence intervals across 4–6 donors per drug, conditional on the frozen disease query. Blue dots indicate opposing mean changes in both query directions. The one-sided exhaustive donor sign-flip test for this net endpoint, BH-adjusted across 21 drugs, yielded no supported drug. This endpoint differs from the later two-component IUT in Figure 4. c, Descriptive Spearman agreement of harmonized scores for all 21 drugs; circles indicate kinase inhibitors and squares other classes. Orange squares denote proteasome inhibitors or an anthracycline; gray squares denote other non-kinase classes. Labels identify carfilzomib (CAR), bortezomib (BOR), doxorubicin (DOX), regorafenib (REG) and vemurafenib (VEM). d, The same comparison restricted to 16 kinase inhibitors, a post-result descriptive sensitivity analysis. Correlations have no independence-based P values or confidence intervals; shared controls, compound classes and donor structure preclude interpreting these plots as an independent drug-response validation. Harmonized enrichment scores in c,d are not on the rank-percentage-point scale used in b. Source tables retain all drugs.

'''+oldcap.replace('# Figure captions and interpretation\n\n','')
 (OUT/'Figure_legends.md').write_text(captions,encoding='utf8')
 write_json(OUT/'input_sources.json',INPUTS)
 checks={'status':'PASS','screening_entries':len(screen),'cardiac_drugs':len(ep),'table2_drugs':len(t2),'all_drug_BH_recalculation':'PASS','all_drug_rho':float(rho),'kinase_rho':float(krho),'old_exports_identical':len(reuse),'visual_review':'PENDING','new_biological_tests':False}
 write_json(OUT/'qa/numeric_checks.json',checks);return checks

if __name__=='__main__':run_standard_module('60_manuscript_assets',main)
