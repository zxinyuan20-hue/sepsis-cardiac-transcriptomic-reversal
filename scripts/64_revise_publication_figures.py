"""Packages: numpy,pandas,scipy,matplotlib,Pillow; existing project environment.
Input/output paths and export parameters: config/figure_revision_v2.json.
Fixed seed; shared run_standard_module log. Raw and frozen outputs are read-only.
Pipeline: verify previous package -> redesign saved data -> export figures/source data.
No new inferential analysis. Saved rat statistics are reconciled for visualization.
"""
import importlib.util,json,logging,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Ellipse,Circle,PathPatch
from matplotlib.path import Path as MPath
from scipy.stats import spearmanr
from pipeline_utils import ROOT,run_standard_module,write_json,sha256
CFG=json.loads((ROOT/'config/figure_revision_v2.json').read_text(encoding='utf8'))
OUT=ROOT/CFG['out'];FIG=OUT/'figures';SRC=OUT/'source_data';INPUTS={}
BLUE='#326C87';ORANGE='#BE7552';PURPLE='#887099';GRAY='#76838B';INK='#263B46';PALE='#F1F5F7'
def read(rel,**kwargs):
 p=ROOT/rel;INPUTS[rel]=sha256(p);return pd.read_csv(p,sep='\t',**kwargs)
def save(fig,name):
 for ext in ['pdf','svg','png','tiff']:fig.savefig(FIG/f'{name}.{ext}',dpi=CFG['dpi'],facecolor='white')
 plt.close(fig)
def label(ax,s):ax.text(-.17,1.07,s,transform=ax.transAxes,fontweight='bold',fontsize=10,color=INK)
def style():
 plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.3,'ytick.labelsize':6.3,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
def heart(ax,x,y,s,color):
 verts=np.array([[0,.25],[-.6,.85],[-1,.1],[0,-.65],[1,.1],[.6,.85],[0,.25]])*s+[x,y]
 ax.add_patch(PathPatch(MPath(verts,[MPath.MOVETO]+[MPath.CURVE4]*6),facecolor=color,edgecolor='none',alpha=.8))
 ax.plot(np.array([-.6,-.28,-.14,.02,.16,.3,.55])*s+x,np.array([.02,.02,.23,-.16,.11,.02,.02])*s+y,color='white',lw=.8)
def dish(ax,x,y,s,color):
 ax.add_patch(Ellipse((x,y),s*1.7,s*.8,facecolor='white',edgecolor=color,lw=1.1))
 for dx,dy in [(-.4,.06),(0,.13),(.4,.02),(-.1,-.18)]:
  ax.add_patch(Ellipse((x+dx*s,y+dy*s),s*.43,s*.19,angle=25,facecolor=color,alpha=.5,edgecolor=color,lw=.5))
  ax.add_patch(Circle((x+dx*s,y+dy*s),s*.035,color=color))
def molecules(ax,x,y,s,color):
 pts=np.array([[-.6,0],[-.25,.4],[.15,.15],[.55,.45],[.6,-.25],[.1,-.45]])*s+[x,y]
 for i,j in [(0,1),(1,2),(2,3),(2,4),(2,5)]:ax.plot(pts[[i,j],0],pts[[i,j],1],color=color,lw=1.3)
 for a,b in pts:ax.add_patch(Circle((a,b),s*.10,facecolor=color,edgecolor='white',lw=.6))
def workflow():
 fig,ax=plt.subplots(figsize=(7.2,7.5));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 def panel(x,y,w,h,letter,title,lines,color,icon):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.008,rounding_size=.01',facecolor=PALE if color==BLUE else '#F6F2F7' if color==PURPLE else '#FAF5F1',edgecolor='none'))
  ax.text(x,y+h+.022,letter,fontweight='bold',fontsize=10)
  ax.text(x+.030,y+h-.030,title,va='top',fontsize=9,fontweight='bold',color=color)
  icon(ax,x+.063,y+h*.43,.047,color)
  ax.text(x+.12,y+h-.079,'\n'.join(lines),va='top',ha='left',fontsize=7.2,linespacing=1.55)
 def arrow(a,b,color=GRAY,dashed=False):
  ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'-|>','mutation_scale':10,'lw':1.15,'color':color,'linestyle':'--' if dashed else '-'})
 # The left branch carries screening; the right branch supplies independent contexts.
 ax.text(.04,.983,'DISEASE-DRIVEN SCREENING',fontsize=7.5,fontweight='bold',color=BLUE)
 ax.text(.55,.983,'INDEPENDENT DISEASE CONTEXTS',fontsize=7.5,fontweight='bold',color=PURPLE)
 panel(.04,.735,.41,.205,'a','Human myocardial discovery',['GSE79962','Sepsis 20 · nonfailing 11','IHD 11 · DCM 9','Frozen query: 100 UP / 90 DOWN'],BLUE,heart)
 panel(.55,.735,.41,.205,'b','Human myocardial cohorts',['GSE237861: 7 sepsis / 7 controls','GSE141864: 5 sepsis / 2 controls','Independent patient profiles','Analyzed separately'],PURPLE,heart)
 panel(.04,.49,.41,.19,'c','Transcriptomic screening',['LINCS Phase II · GSE70138','482 compound entries','6,769 perturbation signatures','Cell-line aggregation'],BLUE,molecules)
 panel(.55,.49,.41,.19,'d','Cross-cohort drug reassessment',['All 21 drugs · 783 shared genes','UP suppression','DOWN restoration','Patient and donor omissions'],PURPLE,molecules)
 panel(.04,.245,.41,.19,'e','Cardiac drug perturbations',['GSE217421 · 21 drugs','Healthy donor iPSC cardiac cells','4–6 donors per drug · 48 h','Matched wells → donor means'],BLUE,dish)
 panel(.55,.245,.41,.19,'f','Cross-model disease evidence',['Rat cardiomyocytes: 6 / 6','Two mouse datasets: 5 / 5 each','LPS versus control','Orthologs and fixed pathways'],ORANGE,dish)
 arrow((.245,.725),(.245,.69));arrow((.245,.48),(.245,.445))
 ax.plot([.46,.50,.50],[.34,.34,.58],color=GRAY,lw=1.1);arrow((.50,.58),(.54,.58))
 ax.text(.502,.456,'Same donors',fontsize=6,ha='center',color=GRAY,rotation=90,bbox={'facecolor':'white','edgecolor':'none','pad':2})
 arrow((.755,.725),(.755,.69),PURPLE)
 # Route drug reassessment outside the independent model panel, with no crossing.
 ax.plot([.97,.986,.986],[.58,.58,.195],color=PURPLE,lw=1.1);arrow((.986,.195),(.91,.153),PURPLE)
 arrow((.755,.235),(.755,.153),ORANGE,True)
 ax.add_patch(FancyBboxPatch((.04,.063),.92,.087,boxstyle='round,pad=.008,rounding_size=.008',facecolor='#EDF2F3',edgecolor='none'))
 ax.text(.50,.125,'Integrated evidence interpretation',fontsize=9,fontweight='bold',ha='center')
 ax.text(.50,.088,'Direction and stability  ·  Multiplicity-adjusted support  ·  Pathway context',fontsize=7.4,ha='center')
 fig.subplots_adjust(left=.02,right=.98,bottom=.01,top=.99);save(fig,'Figure1_workflow')
def screening():
 s=read('outputs/analysis/screening/complete_screening_results.tsv');ctx=read('outputs/analysis/expansion/drug_class_context.tsv');e=read('outputs/analysis/expansion/expanded_endpoint_results.tsv');e=e[e['size'].eq(100)].sort_values('mean_reversal')
 fig=plt.figure(figsize=(7.2,7.2));gs=fig.add_gridspec(2,2,width_ratios=[1.05,1],height_ratios=[1,1.55],hspace=.58,wspace=.60)
 a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[0,1]);d=fig.add_subplot(gs[1,1])
 a.scatter(s['rank'],s.reversal_score,s=5,color='#ACB6BD',lw=0);a.scatter(s.loc[s.reversal_score>0,'rank'],s.loc[s.reversal_score>0,'reversal_score'],s=8,color=BLUE,lw=0)
 a.axhline(0,color=GRAY,ls='--',lw=.6);a.set(xlabel='Compound rank',ylabel='LINCS reversal score',title='LINCS compound screening')
 a.text(.98,.92,'482 entries\n23 positive scores',transform=a.transAxes,ha='right',va='top',fontsize=7)
 y=np.arange(len(e));x=e.mean_reversal.to_numpy()*100
 b.errorbar(x,y,xerr=np.vstack([(e.mean_reversal-e.t_ci_low)*100,(e.t_ci_high-e.mean_reversal)*100]),fmt='o',color=GRAY,ms=2.7,lw=.65,capsize=1.1)
 sel=e.both_directions_opposing.to_numpy();b.scatter(x[sel],y[sel],s=15,color=BLUE,zorder=4)
 b.set_yticks(y,[f'{r.drug_name} ({r.donors})' for r in e.itertuples()]);b.tick_params(axis='y',length=0,labelsize=6)
 b.set_ylim(-.7,20.7);b.axvline(0,color=GRAY,ls='--',lw=.6);b.set(xlabel='Net reversal (rank percentage points)',title='Cardiac donor effects')
 b.grid(axis='x',color='#E9EDEF',lw=.5);b.text(0,-.19,'Dots: means; bars: 95% donor t intervals\nBlue: both mean directions oppose disease',transform=b.transAxes,fontsize=6.2)
 ki=ctx[ctx.broad_class.eq('kinase_inhibitors')]
 for ax,data,title in [(c,ctx,'All cardiac drugs'),(d,ki,'Kinase inhibitors')]:
  for kinase,sub in data.groupby(data.broad_class.eq('kinase_inhibitors')):
   ax.scatter(sub.harmonized_lincs_score,sub.harmonized_cardiac_score,s=23,marker='o' if kinase else 's',color=BLUE if kinase else ORANGE,edgecolor='white',lw=.5)
  rho=spearmanr(data.harmonized_lincs_score,data.harmonized_cardiac_score).statistic
  ax.axhline(0,color='#B7C0C4',lw=.6);ax.axvline(0,color='#B7C0C4',lw=.6);ax.set(xlabel='Harmonized LINCS score',ylabel='Harmonized cardiac score',title=f'{title} (n = {len(data)})')
  ax.text(.03,.96,f'Spearman ρ = {rho:.3f}',transform=ax.transAxes,va='top',fontsize=7)
 d.set_box_aspect(1);d.text(0,-.28,'Descriptive, post-result class sensitivity',transform=d.transAxes,fontsize=6.2)
 c.text(0,-.32,'● Kinase inhibitors    ■ Other classes',transform=c.transAxes,fontsize=6.1)
 for ax,l in [(a,'a'),(b,'b'),(c,'c'),(d,'d')]:label(ax,l)
 fig.subplots_adjust(left=.18,right=.975,top=.94,bottom=.125);save(fig,'Figure2_screening_context')
 for name,data in [('Figure2_LINCS',s),('Figure2_cardiac',e),('Figure2_cross_context',ctx)]:data.to_csv(SRC/f'{name}.tsv',sep='\t',index=False)
def refresh_existing():
 # Reuse the established plotting module's exact data mappings; apply only visual edits.
 spec=importlib.util.spec_from_file_location('frozen_figures',ROOT/'scripts/55_external_compound_figures.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 mod.FIG=SRC
 mapping={'figure1_disease_replication':'Figure3_disease_replication','figure2_all_drug_contexts':'Figure4_external_compounds','figureS1_gene_and_pathway_evidence':'FigureS1_gene_pathway'}
 def export(fig,name):
  if name=='figure2_all_drug_contexts':
   for text in list(fig.texts):
    if 'All 21 fixed drugs' in text.get_text():text.remove()
   for tick in fig.axes[0].get_yticklabels():tick.set_fontweight('normal')
   for axis in fig.axes[:2]:
    for i in range(0,21,2):axis.axhspan(i-.5,i+.5,color='#87949B',alpha=.035,zorder=-2)
  if name=='figureS1_gene_and_pathway_evidence':
   for axis in fig.axes:
    for table in axis.tables:
     for (ri,ci),cell in table.get_celld().items():
      cell.set_facecolor('white');cell.set_edgecolor('#BAC5CA');cell.set_linewidth(.55);cell.visible_edges='TB' if ri==0 else 'B' if ri==5 else ''
      if ri==0:cell.get_text().set_fontweight('bold')
  save(fig,mapping[name])
 mod.save=export;mod.main();style()
def rat_figure():
 rat=read('outputs/analysis/mechanism/LOCAL_RAT_LPS_2025_all_gene_results.tsv',dtype={'human_entrez':str})
 human=read('outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv',dtype={'entrez_id':str})
 concord=read('outputs/analysis/mechanism/cross_model_concordance.tsv')
 shared=rat.dropna(subset=['human_entrez']).merge(human[['entrez_id','logFC']],left_on='human_entrez',right_on='entrez_id',suffixes=('_rat','_human'),validate='one_to_one')
 sig=rat['adj.P.Val']<.05;up=int((sig&(rat.logFC>0)).sum());down=int((sig&(rat.logFC<0)).sum())
 rho=spearmanr(shared.logFC_human,shared.logFC_rat).statistic
 assert (up,down,len(shared))==(1694,1383,10822) and abs(rho-.06194364825493577)<1e-12
 fig=plt.figure(figsize=(7.2,6.15));gs=fig.add_gridspec(2,2,hspace=.63,wspace=.48,height_ratios=[1,1.05]);axs=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)];a,b,c,d=axs
 a.set(xlim=(0,1),ylim=(0,1));a.axis('off');a.set_title('Rat cardiomyocyte model',loc='left',pad=12)
 dish(a,.23,.75,.16,GRAY);dish(a,.75,.75,.16,ORANGE)
 a.text(.23,.51,'Control\nn = 6',ha='center',fontsize=8);a.text(.75,.51,'LPS\nn = 6',ha='center',fontsize=8)
 a.plot([.23,.23,.75,.75],[.45,.36,.36,.45],color=GRAY,lw=.8);a.annotate('',xy=(.49,.25),xytext=(.49,.36),arrowprops={'arrowstyle':'-|>','color':GRAY,'lw':.8})
 a.text(.49,.16,'Bulk transcriptome comparison',ha='center',fontsize=8,fontweight='bold');a.text(.49,.025,'Independent biological samples',ha='center',fontsize=7)
 colors=np.where(sig,np.where(rat.logFC>0,ORANGE,BLUE),'#C4CDD1');y=-np.log10(rat['adj.P.Val'].clip(lower=np.nextafter(0.,1.)))
 b.scatter(rat.logFC,y,c=colors,s=2,alpha=.6,lw=0,rasterized=True);b.axhline(-np.log10(.05),ls='--',lw=.6,color=GRAY)
 b.set(xlabel='LPS versus control log2 fold change',ylabel='−log10(BH-adjusted P)',title='Differential expression')
 b.text(.02,.98,f'{down:,} down',color=BLUE,transform=b.transAxes,va='top',fontsize=7);b.text(.98,.98,f'{up:,} up',color=ORANGE,transform=b.transAxes,va='top',ha='right',fontsize=7)
 h=c.hexbin(shared.logFC_human,shared.logFC_rat,gridsize=40,mincnt=1,cmap='Blues',bins='log',linewidths=0,rasterized=True)
 c.axhline(0,color=GRAY,lw=.55);c.axvline(0,color=GRAY,lw=.55);c.set(xlabel='Human myocardium log2 fold change',ylabel='Rat log2 fold change',title='Human–rat ortholog agreement')
 c.text(.03,.96,f'{len(shared):,} genes; Spearman ρ = {rho:.3f}',transform=c.transAxes,va='top',fontsize=6.7)
 cb=fig.colorbar(h,ax=c,location='bottom',fraction=.055,pad=.33);cb.set_label('Genes per hexagon (log scale)',fontsize=6)
 r=concord[concord.cohort.str.startswith('LOCAL_RAT')].copy();labels=['All samples\n6 control / 6 LPS','Omit C1/C2\n4 control / 6 LPS'];yy=[1,0]
 d.scatter(r.logfc_spearman,yy,color=[BLUE,GRAY],s=30);d.axvline(0,color=GRAY,lw=.6,ls='--')
 d.set_yticks(yy,labels);d.tick_params(axis='y',length=0,labelsize=6.5);d.set(xlim=(-.015,.12),ylim=(-.7,1.7),xlabel='Descriptive Spearman ρ',title='Sample-omission sensitivity')
 for val,y0 in zip(r.logfc_spearman,yy):d.text(val+.006,y0,f'{val:.3f}',va='center',fontsize=7)
 d.grid(axis='x',color='#E9EDEF',lw=.5)
 for ax,l in zip(axs,'abcd'):label(ax,l)
 fig.subplots_adjust(left=.11,right=.975,top=.91,bottom=.14);save(fig,'FigureS2_rat_cardiomyocytes')
 rat.to_csv(SRC/'FigureS2_rat_all_genes.tsv',sep='\t',index=False);shared.to_csv(SRC/'FigureS2_human_rat_orthologs.tsv',sep='\t',index=False);r.to_csv(SRC/'FigureS2_omission_summary.tsv',sep='\t',index=False)
 return {'up':up,'down':down,'shared_genes':len(shared),'spearman':rho}
def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 for folder in [OUT,FIG,SRC,OUT/'qa']:folder.mkdir(parents=True,exist_ok=True)
 previous=json.loads((ROOT/CFG['previous']/'draft_manifest.json').read_text(encoding='utf8'))
 for row in previous:assert sha256(ROOT/row['path'])==row['sha256'],row['path']
 style();workflow();screening();refresh_existing();rat=rat_figure()
 cap=(ROOT/CFG['previous']/'Figure_legends.md').read_text(encoding='utf8')
 start=cap.index('Figure 2 |');end=cap.index('Figure 3 |')
 cap=cap[:start]+'''Figure 2 | Initial compound screening and cardiac-context comparison. a, All 482 eligible LINCS compound entries, ordered by the frozen screening rank; blue points indicate the 23 positive reversal scores. No entry passed the 482-entry BH family. b, All 21 cardiac drugs, ordered by mean net reversal; numbers in parentheses indicate independent donors. Dots show equal-donor means and bars show two-sided 95% Student-t intervals conditional on the frozen disease query. Blue dots indicate both mean components opposing disease; no drug passed the original net-endpoint 21-drug BH family. c, Harmonized enrichment scores for all 21 drugs; circles indicate kinase inhibitors and squares indicate other classes. d, The corresponding post-result descriptive comparison restricted to 16 kinase inhibitors. Spearman correlations are descriptive; shared controls and donor structure preclude interpreting this as independent drug-response validation. Scores in a,c,d differ from the rank-percentage-point endpoint in b. Source tables retain all drugs.

'''+cap[end:]
 cap=cap.replace('Bold names mark the three previously selected hypotheses and are not a selected testing family.','All drug names receive the same typographic weight; the full 21-drug testing family is retained.')
 cap+='''

Figure S2 | Existing rat cardiomyocyte transcriptomes and their agreement with human myocardium. a, Schematic of the existing LPS and control groups, each with six independent biological samples as confirmed by the investigator. Icons are conceptual and are not microscopy images. b, All 14,894 retained rat genes from the saved TMM/voom/limma analysis; blue and orange denote BH-adjusted P < 0.05 with negative and positive log2 fold change, respectively (1,383 downregulated and 1,694 upregulated genes). The dashed line denotes BH q = 0.05; no additional fold-change cutoff was applied. c, All 10,822 human-mapped genes shared with the discovery contrast; hexagon color indicates gene density on a logarithmic scale. Spearman rho is descriptive and genes are not biological replicates. d, Saved human–rat agreement in the main analysis and the sensitivity omitting control samples C1/C2. The full 12-sample analysis remains primary. These results describe an LPS-associated cellular response and cross-model agreement, not a drug-treatment rescue experiment. All values are from existing frozen analyses; no new differential-expression model or candidate selection was run for this figure.
'''
 (OUT/'Figure_legends_v2.md').write_text(cap,encoding='utf8')
 write_json(OUT/'input_sources.json',INPUTS)
 write_json(OUT/'qa/authoring_checks.json',{'previous_manifest_files_unchanged':len(previous),'rat_reconciliation':rat,'figures':6,'visual_review':'PENDING','new_inferential_tests':False})
 return {'figures':6,'rat_reconciliation':rat,'visual_review':'PENDING'}
if __name__=='__main__':run_standard_module('64_revise_publication_figures',main)
