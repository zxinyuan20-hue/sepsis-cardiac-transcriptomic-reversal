"""Packages numpy pandas matplotlib Pillow. Input/output paths from biology JSON.
Fixed placement, no biological data simulation; saved source tables, 600dpi figures.
Reuses frozen patient scores and evidence matrix; Hallmark from module 92.
"""
import json,shutil,logging
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from pipeline_utils import ROOT,write_json,run_standard_module
CFG=json.loads((ROOT/'config/biology_extension_v1.json').read_text(encoding='utf8'))
OUT=ROOT/CFG['publication'];FIG=OUT/'figures';SRC=OUT/'source_data';BIO=ROOT/CFG['analysis']
def save(fig,name):
 for ext in ['pdf','svg','png','tiff']:
  kw={'dpi':600,'facecolor':'white'}
  if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(FIG/f'{name}.{ext}',**kw)
 plt.close(fig)
def main():
 FIG.mkdir(parents=True,exist_ok=True);SRC.mkdir(exist_ok=True)
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'legend.frameon':False})
 # Copy unchanged assets without replacing any assets built in this revision.
 for pattern in ['figures/*','source_data/*','Supplementary_Note_[12].*','Tables.docx']:
  for p in (ROOT/CFG['source_publication']).glob(pattern):
   if p.is_file():
    dest=OUT/p.relative_to(ROOT/CFG['source_publication']);dest.parent.mkdir(parents=True,exist_ok=True)
    if not dest.exists():shutil.copy2(p,dest)
 # Graded replication on the exact frozen score scale.
 f=plt.figure(figsize=(7.2,5.7));gs=f.add_gridspec(2,2,height_ratios=[1.1,1],hspace=.70,wspace=.38);axes=[f.add_subplot(gs[0,0]),f.add_subplot(gs[0,1]),f.add_subplot(gs[1,:])]
 allscores=[];allloo=[]
 for i,(cohort,folder,prefix) in enumerate([('GSE237861','raw_heart_v1',''),('GSE141864','external_human','GSE141864_')]):
  path=ROOT/'outputs/analysis'/folder;s=pd.read_csv(path/f'{prefix}patient_signature_scores.tsv',sep='\t');s['cohort']=cohort;allscores.append(s)
  for j,group in enumerate(['control','sepsis']):
   v=s.loc[s.group==group,'primary_disease_rank_score'].to_numpy()*100;color=['#657984','#BE7552'][j]
   axes[i].scatter(j+np.linspace(-.10,.10,len(v)),v,s=24,c=color,edgecolors='white',linewidths=.5);axes[i].plot([j-.2,j+.2],[v.mean()]*2,c=color,lw=1.5)
  axes[i].set(xticks=[0,1],xticklabels=[f'Control (n={(s.group=="control").sum()})',f'Sepsis (n={(s.group=="sepsis").sum()})'],ylabel='Disease rank score (percentage points)',xlim=(-.5,1.5))
  axes[i].set_title(f'{chr(97+i)}  {cohort}\n'+['Stable directional replication','Weaker, omission-sensitive support'][i],loc='left',pad=10)
  axes[i].text(.02,.98,['P = 0.025; BH q = 0.050*','P = 0.286; BH q = 0.286'][i],transform=axes[i].transAxes,va='top',fontsize=6.7);axes[i].margins(y=.23)
  l=pd.read_csv(path/f'{prefix}leave_one_patient_out.tsv',sep='\t');v=l['case_minus_control' if i==0 else 'effect'].to_numpy()*100
  full=(s.loc[s.group=='sepsis','primary_disease_rank_score'].mean()-s.loc[s.group=='control','primary_disease_rank_score'].mean())*100
  axes[2].scatter(v,i+np.linspace(-.09,.09,len(v)),s=19,c=['#326C87','#887099'][i]);axes[2].scatter(full,i,c='black',marker='D',s=24)
  allloo.extend({'cohort':cohort,'omission':j,'effect':float(x)} for j,x in enumerate(v))
 axes[2].axvline(0,c='#999999',ls='--',lw=.7);axes[2].set(yticks=[0,1],yticklabels=['GSE237861','GSE141864'],ylim=(1.35,-.35),xlabel='Sepsis minus control score (percentage points)');axes[2].set_title('c  Every patient omission; diamonds denote full-cohort effects',loc='left')
 f.subplots_adjust(left=.13,right=.98,bottom=.17,top=.87);f.text(.13,.04,'*Unrounded q = 0.0501166, above the prespecified q < 0.05 threshold.\nPoints are patients or omission diagnostics; horizontal lines in a–b are group means.',fontsize=7)
 save(f,'Figure2_disease_replication');pd.concat(allscores).to_csv(SRC/'Figure2_patient_scores.tsv',sep='\t',index=False);pd.DataFrame(allloo).to_csv(SRC/'Figure2_omissions.tsv',sep='\t',index=False)
 # Evidence hierarchy: group headers and strong formal-inference separator.
 e=pd.read_csv(ROOT/'outputs/analysis/robustness_amendment_v1/evidence_matrix.tsv',sep='\t')
 fields=['LINCS_positive','cardiac_three_scores_positive','primary_mean_bidirectional','secondary_mean_bidirectional','primary_all_donor_omissions','both_contexts_all_patient_omissions','primary_formal_IUT_BH']
 labels=['LINCS\npositive','Cardiac\n3 scores > 0','Primary\nUP & DOWN','Secondary\nUP & DOWN','Primary\ndonor LOO','Both contexts\npatient LOO','Primary\nIUT / BH','Historical\ncardiotoxicity']
 f,ax=plt.subplots(figsize=(7.2,6.7));f.subplots_adjust(left=.20,right=.98,bottom=.17,top=.80)
 for i,r in e.iterrows():
  ax.axhspan(i-.5,i+.5,color='#F5F6F7' if i%2==0 else 'white',zorder=0)
  for j,col in enumerate(fields):
   yes=bool(r[col]);ax.scatter(j,i,s=34,facecolor='#287E8E' if yes else 'white',edgecolor='#287E8E' if yes else '#AEB7BD',lw=.7)
  ax.text(7.35,i,r.source_cardiotoxicity,ha='center',va='center',fontsize=7,color='#A05E2A' if r.source_cardiotoxicity=='Yes' else '#505B64')
 for x,lw in [(3.5,.7),(5.5,1.5),(6.5,.6)]:ax.axvline(x,color='#64737D',lw=lw)
 ax.set(xlim=(-.5,7.95),ylim=(20.65,-.65),yticks=range(21),yticklabels=e.drug_name,xticks=list(range(7))+[7.35],xticklabels=labels);ax.tick_params(axis='x',labeltop=True,labelbottom=False,length=0,pad=8);ax.tick_params(axis='y',length=0)
 for sp in ax.spines.values():sp.set_visible(False)
 for x,txt in [(1.5,'Directional evidence'),(4.5,'Robustness\nevidence'),(6.75,'Formal inference /\nsafety context')]:
  ax.text(x,1.14,txt,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontweight='bold',fontsize=7.4)
 f.text(.025,.98,'Evidence dimensions across 21 cardiac drugs',fontsize=10,weight='bold',va='top')
 f.text(.20,.09,'Filled: criterion met   Open: criterion not met   LOO: leave one out\nPrimary: GSE237861; secondary: GSE141864. IUT/BH: 21-drug family, q < 0.05.\nNo drug passed formal inference. Historical No does not establish safety.\nND: not determined. Columns are not independent or an additive evidence score.',fontsize=6.8,va='center',linespacing=1.5)
 save(f,'Figure5_evidence_matrix');e.to_csv(SRC/'Figure5_evidence_matrix.tsv',sep='\t',index=False)
 # Full predefined collection; no significance-selected pathways.
 m=pd.read_csv(BIO/'common_NES_matrix.tsv',sep='\t',index_col=0);names=[r['name'] for r in CFG['datasets']];m=m[names]
 f=plt.figure(figsize=(7.2,9.6));gs=f.add_gridspec(2,1,height_ratios=[6.2,1.7],hspace=.36);ax=f.add_subplot(gs[0]);bx=f.add_subplot(gs[1]);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#E4E7E9');limit=np.nanmax(abs(m.to_numpy()));im=ax.imshow(m,aspect='auto',cmap=cmap,vmin=-limit,vmax=limit)
 display=[]
 for x in m.index:
  s=x.removeprefix('HALLMARK_').replace('_',' ').capitalize()
  for a,b in [('Dna','DNA'),('Il2 stat5','IL2–STAT5'),('Il6 jak stat3','IL6–JAK–STAT3'),('Kras','KRAS'),('Mtorc1','mTORC1'),('Myc','MYC'),('P53','p53'),('Pi3k akt mtor','PI3K–AKT–mTOR'),('Tgf beta','TGFβ'),('Tnfa signaling via nfkb','TNFα signaling via NF-κB'),('Uv','UV')]:s=s.replace(a,b)
  display.append(s)
 ax.set_yticks(range(50),display,fontsize=6.0);ax.set_xticks(range(6),['Human\n79962','Human\n237861','Human\n141864','Rat\nH9c2','Mouse\n185754','Mouse\n267388'],fontsize=6.5);ax.xaxis.tick_top();ax.tick_params(length=0);ax.set_title('a  Hallmark program correspondence',loc='left',fontweight='bold',pad=32)
 f.subplots_adjust(left=.34,right=.91,top=.90,bottom=.10);cax=f.add_axes([.935,.38,.018,.45]);cb=f.colorbar(im,cax=cax);cb.set_label('Normalized enrichment score (NES)',fontsize=7);cb.ax.tick_params(labelsize=6)
 c=pd.read_csv(BIO/'cross_model_concordance.tsv',sep='\t');c=c[(c.scope=='common')&(c.dataset_a=='GSE79962')]
 y=np.arange(5);bx.hlines(y,c.gene_signed_z_rho,c.pathway_NES_rho,color='#9DA6AB',lw=1);bx.scatter(c.gene_signed_z_rho,y,c='#72818B',s=23,label='Gene signed z');bx.scatter(c.pathway_NES_rho,y,c='#326C87',marker='D',s=23,label='Pathway NES')
 bx.set(yticks=y,yticklabels=c.dataset_b,xlim=(-.05,1),xlabel='Descriptive Spearman correlation with GSE79962');bx.invert_yaxis();bx.axvline(0,c='#BBBBBB',lw=.6);bx.set_title('b  Matched gene background',loc='left',fontweight='bold');bx.legend(loc='lower right',fontsize=6.5)
 f.text(.05,.025,'All 50 Hallmark sets displayed; gray = below 15-gene coverage (pancreas beta cells: 7 genes).\nNES concordance uses 49 eligible sets over 8,258 shared genes. Positive NES: higher in sepsis/LPS.\nGene and pathway correlations compare different units; their difference is descriptive, not a formal test.',fontsize=6.7,linespacing=1.4)
 save(f,'Figure6_Hallmark_concordance')
 for name in ['hallmark_enrichment.tsv','hallmark_coverage.tsv','cross_model_concordance.tsv','common_NES_matrix.tsv','full_NES_matrix.tsv']:shutil.copy2(BIO/name,SRC/name)
 p=FIG/'Figure1_workflow.png';im=Image.open(p).convert('RGB');im.save(p,dpi=(600,600));im.save(FIG/'Figure1_workflow.tiff',dpi=(600,600),compression='tiff_lzw')
 write_json(OUT/'figure_revision_checks.json',{'backend':'Python plots; existing editable PowerPoint workflow retained','all_50_pathways_displayed':True,'NA_not_zero':True,'formal_drug_hits_unchanged':0,'visual_review':'pending'})
 return {'figures_revised':[1,2,5],'figures_added':[6]}
if __name__=='__main__':run_standard_module('96_evidence_hierarchy_figures',main)
