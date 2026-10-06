"""Python-only matplotlib figures; packages pandas,numpy,matplotlib.
Inputs audited54 tables; outputs external_compound_v1/figures, source data, captions.
Seed20261005; 183mm width, white background, editable PDF/SVG,600dpi raster.
"""
import json,logging,textwrap
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline_utils import ROOT,write_json,run_standard_module
OUT=ROOT/'outputs/analysis/external_compound_v1';FIG=OUT/'figures'
def save(fig,name):
 for suffix in ['pdf','svg','png','tiff']:fig.savefig(FIG/(name+'.'+suffix),dpi=600,facecolor='white')
 plt.close(fig)
def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 assert json.loads((OUT/'numeric_qa.json').read_text())['status']=='PASS'
 FIG.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.65,'legend.frameon':False})
 colors={'control':'#697982','sepsis':'#BD6940'}
 fig=plt.figure(figsize=(7.2,5.35));gs=fig.add_gridspec(2,2,height_ratios=[1.25,1],hspace=.62,wspace=.42);axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
 cohorts=['GSE237861','GSE141864'];source=[];loo=[]
 for i,cohort in enumerate(cohorts):
  folder=ROOT/'outputs/analysis'/('raw_heart_v1' if i==0 else 'external_human');sf='patient_signature_scores.tsv' if i==0 else 'GSE141864_patient_signature_scores.tsv';lf='leave_one_patient_out.tsv' if i==0 else 'GSE141864_leave_one_patient_out.tsv'
  s=pd.read_csv(folder/sf,sep='\t');s['cohort']=cohort;source.append(s);ax=axs[i]
  for j,group in enumerate(['control','sepsis']):
   values=s.loc[s.group.eq(group),'primary_disease_rank_score'].to_numpy()*100;offset=np.linspace(-.11,.11,len(values))
   ax.scatter(j+offset,values,s=23,color=colors[group],edgecolor='white',linewidth=.4,zorder=3);ax.plot([j-.2,j+.2],[values.mean()]*2,color=colors[group],lw=1.5)
  ax.set_xticks([0,1],[f'Control (n={sum(s.group.eq("control"))})',f'Sepsis (n={sum(s.group.eq("sepsis"))})']);ax.set_xlim(-.5,1.5);ax.set_ylabel('Disease rank score (percentage points)');ax.set_title(cohort+'\n'+('P = 0.025058; BH q = 0.050117' if i==0 else 'P = 0.285714; BH q = 0.285714'),loc='left');ax.grid(axis='y',color='#eeeeee',lw=.5);ax.text(-.19,1.17,'ab'[i],transform=ax.transAxes,fontweight='bold',fontsize=10)
  l=pd.read_csv(folder/lf,sep='\t');vals=l['case_minus_control' if i==0 else 'effect'].to_numpy()*100
  for j,value in enumerate(vals):loo.append({'cohort':cohort,'omission_index':j,'effect_percentage_points':value})
  axs[2].scatter(vals,i+np.linspace(-.10,.10,len(vals)),s=18,color=['#276B8C','#9873A1'][i],zorder=3)
  full=(s.loc[s.group.eq('sepsis'),'primary_disease_rank_score'].mean()-s.loc[s.group.eq('control'),'primary_disease_rank_score'].mean())*100
  axs[2].scatter([full],[i],marker='D',s=25,color='black',zorder=4)
 axs[2].axvline(0,color='#999999',ls='--',lw=.8);axs[2].set_yticks([0,1],cohorts);axs[2].set_ylim(1.4,-.4);axs[2].set_xlabel('Case minus control score (percentage points)');axs[2].set_title('Leave-one-patient-out direction; diamonds show full-cohort effects',loc='left');axs[2].grid(axis='x',color='#eeeeee',lw=.5);axs[2].text(-.08,1.12,'c',transform=axs[2].transAxes,fontweight='bold',fontsize=10)
 fig.subplots_adjust(left=.13,right=.97,top=.86,bottom=.12);save(fig,'figure1_disease_replication')
 pd.concat(source,ignore_index=True).to_csv(FIG/'figure1_patient_source.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(FIG/'figure1_omission_source.tsv',sep='\t',index=False)
 r=pd.read_csv(OUT/'all_drug_results.tsv',sep='\t');names=r[['drug_name','drug_state','donors']].drop_duplicates().sort_values('drug_name');order=names.drug_state.tolist()
 net=r.pivot(index='drug_state',columns='cohort',values='net_reversal').loc[order,['GSE79962','GSE237861','GSE141864']]*100
 main=r[r.cohort.eq('GSE237861')].set_index('drug_state').loc[order]
 fig=plt.figure(figsize=(7.2,8.2));gs=fig.add_gridspec(1,2,width_ratios=[1,1.5],wspace=.2);ax=fig.add_subplot(gs[0]);bx=fig.add_subplot(gs[1]);lim=float(abs(net.to_numpy()).max())
 im=ax.imshow(net,aspect='auto',cmap='RdBu',vmin=-lim,vmax=lim,interpolation='nearest')
 ax.set_yticks(range(21),[f'{r.drug_name} (n={r.donors})' for r in names.itertuples()]);ax.set_xticks(range(3),['Discovery\n79962','External\n237861','External\n141864']);ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=5)
 for label,state in zip(ax.get_yticklabels(),order):
  if state in ['DAS','PON','REG']:label.set_fontweight('bold')
 for i in range(21):
  for j in range(3):ax.text(j,i,f'{net.iloc[i,j]:.2f}',ha='center',va='center',fontsize=6,color='white' if abs(net.iloc[i,j])>.6*lim else '#222222')
 yy=np.arange(21)
 for k,name,color,label in [(0,'UP_suppression','#276B8C','UP suppression'),(1,'DOWN_restoration','#BD6940','DOWN restoration')]:
  point=main[name].to_numpy()*100;lo=main[name+'_ci_low'].to_numpy()*100;hi=main[name+'_ci_high'].to_numpy()*100
  bx.errorbar(point,yy+(k-.5)*.18,xerr=np.c_[point-lo,hi-point].T,fmt='o',ms=3,elinewidth=.65,capsize=1.5,color=color,label=label)
 bx.set_ylim(20.5,-.5);bx.set_yticks(yy,[]);bx.tick_params(axis='y',length=0);bx.axvline(0,color='#999999',lw=.8,ls='--');bx.grid(axis='x',color='#eeeeee',lw=.5);bx.set_xlabel('Component score (rank percentage points)');bx.set_title('GSE237861 components\nConditional donor t intervals (95%)',fontsize=8,pad=20)
 bx.legend(loc='upper center',bbox_to_anchor=(.5,-.065),ncol=1,fontsize=7)
 fig.text(.04,.96,'a',fontweight='bold',fontsize=10);fig.text(.55,.96,'b',fontweight='bold',fontsize=10);fig.text(.5,.985,'All 21 fixed drugs: primary bidirectional BH hits = 0; minimum q = 0.328125',ha='center',fontsize=8)
 fig.subplots_adjust(left=.22,right=.985,top=.89,bottom=.15)
 box=ax.get_position();cax=fig.add_axes([box.x0,.075,box.width,.014]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Net reversal (rank percentage points)',fontsize=6.5)
 assert abs(ax.get_position().y0-bx.get_position().y0)<1e-12 and abs(ax.get_position().y1-bx.get_position().y1)<1e-12
 save(fig,'figure2_all_drug_contexts')
 r.to_csv(FIG/'figure2_source.tsv',sep='\t',index=False)
 # Supplementary evidence contract: partial all-gene concordance does not establish fixed pathways.
 fig=plt.figure(figsize=(7.2,6.0));gs=fig.add_gridspec(2,2,height_ratios=[1,1.25],hspace=.75,wspace=.4)
 for i,(cohort,rel) in enumerate([('GSE237861','raw_heart_v1/gene_effect_comparison.tsv'),('GSE141864','external_human/GSE141864_discovery_gene_comparison.tsv')]):
  t=pd.read_csv(ROOT/'outputs/analysis'/rel,sep='\t');ax=fig.add_subplot(gs[0,i]);rho=t.logFC_discovery.corr(t.logFC_external,method='spearman')
  ax.scatter(t.logFC_discovery,t.logFC_external,s=1.3,alpha=.25,color='#276B8C',rasterized=True,edgecolors='none');ax.axhline(0,color='#bbbbbb',lw=.5);ax.axvline(0,color='#bbbbbb',lw=.5)
  ax.set_xlabel('GSE79962 log2 fold change');ax.set_ylabel(cohort+' log2 fold change');ax.set_title(f'{cohort}: {len(t):,} shared genes\nDescriptive Spearman rho = {rho:.3f}',loc='left');ax.text(-.17,1.16,'ab'[i],transform=ax.transAxes,fontweight='bold',fontsize=10)
  t.to_csv(FIG/(cohort+'_gene_concordance_source.tsv'),sep='\t',index=False)
 paths=pd.read_csv(ROOT/'outputs/analysis/raw_heart_v1/two_cohort_pathway_family.tsv',sep='\t');paths=paths[paths['mode'].eq('estimated_residual_correlation')]
 labels=pd.read_csv(ROOT/'outputs/analysis/mechanism/fixed_pathways.tsv',sep='\t').drop_duplicates('pathway_id').set_index('pathway_id').pathway_name
 ax=fig.add_subplot(gs[1,:]);ax.axis('off');cells=[]
 for pid in sorted(paths.pathway_id.unique()):
  values=[textwrap.fill(labels[pid],36)]
  for cohort in ['GSE237861','GSE141864']:
   row=paths[paths.pathway_id.eq(pid)&paths.cohort.eq(cohort)].iloc[0]
   values.append(f'{row.Direction}; q={row.planned_10_test_BH_q:.3f}' if pd.notna(row.PValue) else f'Not tested\n{int(row.NGenes)} genes; coverage {row.coverage:.1%}')
  cells.append(values)
 table=ax.table(cellText=cells,colLabels=['Fixed pathway','GSE237861','GSE141864'],colWidths=[.50,.25,.25],loc='center',cellLoc='left',bbox=[0,0,1,1]);table.auto_set_font_size(False);table.set_fontsize(6.5)
 for (ri,ci),cell in table.get_celld().items():cell.set_edgecolor('#dddddd');cell.set_linewidth(.5);cell.set_facecolor('#edf0f2' if ri==0 else 'white')
 ax.set_title('Fixed pathway tests with estimated gene correlation\nBH family: two cohorts x five sets',loc='left',pad=10);ax.text(-.08,1.16,'c',transform=ax.transAxes,fontweight='bold',fontsize=10)
 fig.subplots_adjust(left=.13,right=.97,bottom=.055,top=.88);save(fig,'figureS1_gene_and_pathway_evidence');paths.to_csv(FIG/'figureS1_pathway_source.tsv',sep='\t',index=False)
 captions='''# Figure captions and interpretation

Figure 1 | Independent human myocardial disease-pattern assessment. Panels a and b show all individual patients, with group means as horizontal lines; 7 sepsis and 7 noninfectious critically ill controls in GSE237861, and 5 sepsis and 2 controls in the FFPE subset of GSE141864. Scores use the frozen 100UP/90DOWN query, with cohort-specific measured/shared backgrounds. Their absolute levels should not be compared across panels. One-sided patient-label permutation tests use 3432 and 21 assignments, respectively; BH adjustment includes the two planned cohorts. GSE237861 q=0.050117 does not meet q<0.05. Panel c shows every patient omission with fixed original score preprocessing; diamonds mark full-cohort effects. These omission results are direction diagnostics, not refitted disease models. No patient was removed from the main analysis.

Figure 2 | Cardiac compound reversal across disease contexts. All21 previously frozen drugs are displayed alphabetically. Panel a gives equal-donor mean net reversal on the same783-gene intersection; drug ranks remain defined over the original876 measured genes. Net reversal is the sum of UP suppression and DOWN restoration and cannot alone satisfy the bidirectional endpoint. Positive blue values indicate opposition to the disease expression profile, not clinical benefit. Panel b shows the two components for the primary GSE237861 context, with means and two-sided95% Student-t intervals across4–6 independent donors. Intervals condition on the disease weights, are not multiplicity-adjusted and do not incorporate patient-profile uncertainty. Primary significance requires the maximum of two one-sided exhaustive donor sign-flip P values (IUT), adjusted across21 drugs by BH; no drug passes. Bold names mark the three previously selected hypotheses and are not a selected testing family. Drug and matched control wells are first averaged within matched culture blocks and then equally within donor. The same healthy iPSC-CM drug experiment is reused; only the disease context is new. Drug-state doses/times are preserved from the frozen input design. No dose-response or diseased-cell rescue claim is made.

Figure S1 | Gene-level agreement and fixed pathway evidence. Panels a,b show every shared fitted gene, with descriptive Spearman correlations and no gene-independent significance test. The GSE237861 effects are sex-adjusted edgeR coefficients; the GSE141864 and GSE79962 effects are limma coefficients. Panel c retains the five previously fixed Reactome sets and the primary camera estimates of residual gene correlation. The planned two-cohort ten-test BH family is unchanged. An ineligible low-coverage set is labelled not tested; its correction placeholder must not be read as an observedP. No eligible set passes the primary BH threshold. Cohorts have different controls/platforms and are not pooled.

All figures use actual saved results, without simulated biological data. Source tables accompany each figure. VectorPDF/SVG retain text, rasterPNG/TIFF are600dpi; target-journal formatting review has not been performed.
'''
 (FIG/'captions.md').write_text(captions,encoding='utf8');return {'figures':3,'formats':['pdf','svg','png','tiff'],'width_mm':182.88,'backend':'Python matplotlib','visual_review':'PENDING'}
if __name__=='__main__':run_standard_module('55_external_compound_figures',main)
