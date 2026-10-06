"""Packages matplotlib/numpy/pandas/Pillow; Python-only existing figure workflow.
Input/output: amendment config and audited tables. Deterministic visual placement.
Contract: distinguish directional/stability criteria, formal support and source toxicity;
quantitative grid, white background, 183 mm width, Arial, PDF/PNG/TIFF and source tables.
No inferential regression bands or P values for non-independent drug associations.
"""
import json,shutil,logging
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from pipeline_utils import ROOT,run_standard_module,write_json
from robustness_utils import CFG,OUT

def main():
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 assert json.loads((OUT/'numeric_audit.json').read_text())['formal_hits_unchanged']
 pub=ROOT/CFG['publication'];figdir=pub/'figures';figdir.mkdir(parents=True,exist_ok=True);(pub/'source_data').mkdir(exist_ok=True)
 plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False})
 contract={'conclusion':'High overall scoring concordance coexists with algorithm-sensitive shortlists and no formal cardiac drug support; source toxicity labels do not establish effects at assay exposure.','archetype':'quantitative grid and evidence matrix','backend':'Python matplotlib','width_mm':183,'exports':['PDF','600dpi PNG','600dpi TIFF','SVG'],'inference':'descriptive drug associations; no independent-drug P/CI; no causal toxicity interpretation'}
 write_json(pub/'figure_contract.json',contract)
 def save(fig,name):
  for ext in ['pdf','svg','png','tiff']:
   kw={'dpi':600,'facecolor':'white'}
   if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
   fig.savefig(figdir/f'{name}.{ext}',**kw)
  plt.close(fig)
 evidence=pd.read_csv(OUT/'evidence_matrix.tsv',sep='\t');shutil.copy2(OUT/'evidence_matrix.tsv',pub/'source_data/evidence_matrix.tsv')
 fields=['LINCS_positive','cardiac_three_scores_positive','primary_mean_bidirectional','secondary_mean_bidirectional','primary_all_donor_omissions','both_contexts_all_patient_omissions','primary_formal_IUT_BH']
 labels=['LINCS\npositive','Cardiac\n3 scores > 0','Primary\nUP & DOWN','Secondary\nUP & DOWN','Primary\ndonor LOO','Both contexts\npatient LOO','Primary\nIUT / BH']
 f,ax=plt.subplots(figsize=(183/25.4,166/25.4));f.subplots_adjust(left=.205,right=.98,bottom=.17,top=.86)
 for i,r in evidence.iterrows():
  ax.axhspan(i-.5,i+.5,color='#F5F6F7' if i%2==0 else 'white',zorder=0)
  for j,col in enumerate(fields):
   yes=bool(r[col]);ax.scatter(j,i,s=37,facecolor='#287E8E' if yes else 'white',edgecolor='#287E8E' if yes else '#B4BAC0',linewidth=.7,zorder=2)
  ax.text(7.3,i,r.source_cardiotoxicity,ha='center',va='center',fontsize=7,color='#A05E2A' if r.source_cardiotoxicity=='Yes' else '#505B64')
 ax.axvline(5.5,color='#ADB4BA',lw=.7);ax.axvline(6.5,color='#ADB4BA',lw=.7)
 ax.set_xlim(-.5,7.85);ax.set_ylim(20.65,-.65);ax.set_yticks(range(21),evidence.drug_name);ax.set_xticks(range(7),labels);ax.tick_params(axis='both',length=0);ax.tick_params(axis='x',labeltop=True,labelbottom=False,pad=8);ax.text(7.3,-1.05,'Source\ncardiotoxicity',ha='center',va='bottom',fontsize=7)
 for sp in ax.spines.values():sp.set_visible(False)
 f.text(.025,.975,'Evidence dimensions across 21 cardiac drugs',fontsize=10,weight='bold',va='top')
 f.text(.205,.09,'Filled: criterion met     Open: criterion not met\nLOO: leave one out. Primary: GSE237861; secondary: GSE141864.\nIUT/BH: both components, 21-drug family, q < 0.05; no drug passed.\nToxicity: Hansen et al. source category; No does not establish safety; ND = not determined.',fontsize=7,va='center',linespacing=1.5)
 save(f,'Figure5_evidence_matrix')
 entry=pd.read_csv(OUT/'multiscore/entry_method_scores.tsv',sep='\t');joined=pd.read_csv(OUT/'toxicity/all_context_scores.tsv',sep='\t');x=joined[joined.cohort=='GSE237861'].copy();corr=pd.read_csv(OUT/'multiscore/method_concordance.tsv',sep='\t')
 shutil.copy2(OUT/'multiscore/entry_method_scores.tsv',pub/'source_data/LINCS_method_scores.tsv');x.to_csv(pub/'source_data/primary_toxicity_context.tsv',sep='\t',index=False)
 f,axes=plt.subplots(2,2,figsize=(183/25.4,157/25.4));f.subplots_adjust(left=.10,right=.98,top=.94,bottom=.15,hspace=.51,wspace=.36)
 for i,(method,label) in enumerate([('negative_cosine','Negative cosine'),('negative_Spearman','Negative Spearman')]):
  ax=axes[0,i];ax.scatter(entry.weighted_ES,entry[method],s=5,c='#427E93',alpha=.48,rasterized=True);ax.axhline(0,color='#ADB4BA',lw=.6);ax.axvline(0,color='#ADB4BA',lw=.6);ax.set_xlabel('Original weighted ES reversal');ax.set_ylabel(label+' reversal')
  rho=corr[(corr.scope=='482 LINCS entries')&(corr.method_a=='weighted_ES')&(corr.method_b==method)].rank_spearman.iloc[0]
  ax.set_title(f'{chr(97+i)}  LINCS algorithm comparison',loc='left',fontweight='bold');ax.text(.04,.94,f'n = 482 entries\nDescriptive ρ = {rho:.3f}',transform=ax.transAxes,va='top',fontsize=7)
 colors={'Yes':'#B57542','No':'#397E91','ND':'#8C9298'};ax=axes[1,0]
 for cl,g in x.groupby('source_cardiotoxicity'):ax.scatter(g.RMS_mean_logFC,g.net_reversal,s=21,c=colors[cl],edgecolor='white',linewidth=.4,label=cl)
 ax.set_xscale('log');ax.set_xticks([.2,.4,.8,1.6],['0.2','0.4','0.8','1.6']);ax.xaxis.set_minor_formatter(NullFormatter());ax.set_ylim(-.145,.035);ax.axhline(0,color='#ADB4BA',lw=.6);ax.set_xlabel('Perturbation magnitude (RMS logFC; log scale)');ax.set_ylabel('Net reversal (percentile-rank units)');ax.set_title('c  Perturbation magnitude',loc='left',fontweight='bold')
 for state in ['PON','REG','CAB','CAR','DOX']:
  r=x[x.drug_state==state].iloc[0];offset={'PON':(4,5),'REG':(4,-9),'CAB':(4,5),'CAR':(-18,-10),'DOX':(3,-10)}[state];ax.annotate(state,(r.RMS_mean_logFC,r.net_reversal),xytext=offset,textcoords='offset points',fontsize=6)
 ax.text(.03,.06,'n = 21; descriptive ρ = −0.471',transform=ax.transAxes,fontsize=7)
 ax=axes[1,1]
 for j,cl in enumerate(['No','Yes','ND']):
  g=x[x.source_cardiotoxicity==cl].sort_values('drug_name');xx=j+np.linspace(-.12,.12,len(g));ax.scatter(xx,g.net_reversal,s=20,color=colors[cl],alpha=.85);med=g.net_reversal.median();ax.plot([j-.2,j+.2],[med,med],color='#283239',lw=1.4)
 ax.axhline(0,color='#ADB4BA',lw=.6);ax.set_xticks([0,1,2],['No\n(n = 11)','Yes\n(n = 8)','ND\n(n = 2)']);ax.set_xlabel('Source cardiotoxicity category');ax.set_ylabel('Net reversal (percentile-rank units)');ax.set_title('d  Historical toxicity categories',loc='left',fontweight='bold')
 handles=[Line2D([0],[0],marker='o',ls='',color=v,label=k,markersize=4) for k,v in colors.items()];f.legend(handles=handles,title='Source cardiotoxicity',loc='lower center',bbox_to_anchor=(.5,.042),ncol=3,fontsize=7,title_fontsize=7)
 f.text(.1,.012,'Panels c–d: GSE237861 profile; shared donors/controls. Lines in d show medians. No does not establish safety.',fontsize=6.5)
 save(f,'FigureS3_scoring_and_perturbation')
 legend='''# New figure legends

Figure 5. Directional, stability and inferential evidence across the fixed 21-drug cardiac panel. Filled circles denote satisfaction of the column-specific criterion; open circles denote its absence. LINCS positivity refers to the original frozen cross-cell median score. Cardiac three-score positivity requires positive median donor scores under the harmonized weighted enrichment, negative cosine and negative Spearman methods on 897 common measured landmarks. Primary and secondary mean bidirectionality refer to positive UP suppression and DOWN restoration under GSE237861 and GSE141864, respectively. Donor omission requires bidirectionality after every donor omission in the primary context. Patient omission requires bidirectionality after every patient omission in both external contexts. Formal support requires the primary intersection–union test with BH q < 0.05 across all 21 drugs; none met this criterion. Columns reuse data and are not independent replications or an additive evidence score. Toxicity categories reproduce Hansen et al. Supplementary Data 3: Yes, No and ND (not determined). No does not establish clinical safety or absence of toxicity at the studied exposure.

Figure S3. Fixed-query scoring sensitivity and perturbation context. (a,b) Original weighted enrichment reversal versus negative cosine and negative Spearman scores for all 482 eligible LINCS entries. All methods use the same signed 100-UP/90-DOWN query over 978 measured landmarks and median aggregation within cell followed by median aggregation across cells. Spearman rho is descriptive; no independent-entry test is shown. (c) RMS of equal-donor mean drug-minus-control log2-expression changes over 14,940 genes retained for all 21 drugs versus net reversal under the GSE237861 disease profile. The x-axis is logarithmic. PON, ponatinib; REG, regorafenib; CAB, cabozantinib; CAR, carfilzomib; DOX, doxorubicin. (d) The same endpoint by historical source cardiotoxicity category; horizontal segments mark medians. No indicates absence of classification as cardiotoxic in this source, not established safety. Shared donors and controls, heterogeneous drugs and exposures preclude treating these drug summaries as independent causal observations. No drug-efficacy or safety conclusion follows from these descriptive relationships.
'''
 (pub/'Figure_legends_addendum.md').write_text(legend,encoding='utf8')
 return {'new_figures':2,'formats':4,'visual_review':'pending'}

if __name__=='__main__':run_standard_module('82_robustness_figures',main)
