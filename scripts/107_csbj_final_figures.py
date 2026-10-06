"""Required packages: matplotlib, numpy, pandas, Pillow; Python figure backend.
Input/output paths: config/csbj_delivery_v1.json, audited figure sources.
Deterministic geometry; no simulated data or statistical changes.
Main pipeline: retain quantitative originals; redraw workflow and evidence matrix;
save native vector PDF/SVG and RGB 600-dpi TIFF/PNG plus export metadata.
"""
import json,shutil,logging,io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Ellipse,Circle,RegularPolygon,PathPatch
from matplotlib.path import Path as MP
from PIL import Image
from pipeline_utils import ROOT,write_json,run_standard_module
CFG=json.loads((ROOT/'config/csbj_delivery_v1.json').read_text())
OUT=ROOT/CFG['output'];FIG=OUT/'figure_sources'
B='#376F8A';R='#BB6856';P='#807195';I='#263C48';G='#79878F'

def save(f,name):
 for ext in ['pdf','svg','png']:
  f.savefig(FIG/f'{name}.{ext}',dpi=CFG['dpi'],facecolor='white')
 plt.close(f)
 with Image.open(FIG/f'{name}.png') as im:im.convert('RGB').save(FIG/f'{name}.tiff',dpi=(600,600),compression='tiff_lzw')

def main():
 FIG.mkdir(parents=True,exist_ok=True)
 logging.getLogger('fontTools').setLevel(logging.WARNING)
 plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
 src=ROOT/CFG['scientific_source']/'figures'
 for p in src.iterdir():
  if p.suffix in ['.pdf','.svg','.png','.tiff'] and not p.name.startswith(('Figure1_','Figure5_')):shutil.copy2(p,FIG/p.name)
 # Schematic-led composite: discovery and external context feed separate evidence levels.
 # Native primitives only; no raster generative illustrations and no PPT conversion.
 f,ax=plt.subplots(figsize=(7.2,4.5));f.subplots_adjust(0,0,1,1);ax.set(xlim=(0,960),ylim=(570,0));ax.axis('off')
 def txt(x,y,s,size=7,c=I,bold=False,ha='center'):
  ax.text(x,y,s,fontsize=size,color=c,weight='bold' if bold else 'normal',ha=ha,va='center',linespacing=1.45)
 def line(x,y,xx,yy,c=G,dash=False,arrow=True):
  ax.annotate('',(xx,yy),(x,y),arrowprops=dict(arrowstyle='->' if arrow else '-',lw=.8,color=c,linestyle='--' if dash else '-'))
 def heart(x,y):
  v=np.array([[0,.2],[-.6,.85],[-1,.1],[0,-.65],[1,.1],[.6,.85],[0,.2]])*np.array([45,-45])+[x,y]
  ax.add_patch(PathPatch(MP(v,[MP.MOVETO]+[MP.CURVE4]*6),fc=B,ec=B,lw=.7))
  ax.plot(np.array([-.60,-.28,-.14,.02,.16,.30,.55])*45+x,np.array([.02,.02,.23,-.16,.11,.02,.02])*45+y,c='white',lw=.7)
 txt(25,25,'(a)  Human myocardial signature and compound screening',9,B,True,'left')
 txt(25,49,'Discovery / prespecified screening',6.5,G,False,'left')
 heart(105,111)
 txt(105,176,'Human myocardium',8,B,True);txt(105,217,'GSE79962\n20 sepsis / 11 nonfailing\nIHD 11 / DCM 9')
 ax.arrow(273,140,0,-49,width=10,head_width=27,head_length=16,length_includes_head=True,fc=R,ec=R)
 ax.arrow(326,89,0,49,width=10,head_width=27,head_length=16,length_includes_head=True,fc=B,ec=B)
 txt(271,158,'100 UP',7,R,True);txt(331,158,'90 DOWN',7,B,True)
 txt(301,185,'Disease signature',8,B,True);txt(301,219,'Differential expression\nMeasured landmark genes')
 for x,y in [(490,106),(531,126),(572,106)]:ax.add_patch(RegularPolygon((x,y),6,radius=19,orientation=0,fc='white',ec=B,lw=1))
 line(505,113,515,119,B,arrow=False);line(547,119,558,113,B,arrow=False)
 txt(532,176,'Compound screening',8,B,True);txt(532,217,'LINCS Phase II\n482 compound entries\n6,769 signatures')
 ax.add_patch(FancyBboxPatch((734,80),132,66,boxstyle='round,pad=0,rounding_size=6',fc='#F3F7F8',ec=B,lw=.8))
 for y in [93,113,133]:
  for x in [747,768,789,810,831,852]:ax.add_patch(Circle((x,y),5,fc='#D5E5EA',ec=B,lw=.5))
 txt(799,176,'Cardiac perturbations',8,B,True);txt(799,217,'GSE217421\n21 drugs / 4–6 donors per drug\nHealthy iPSC cardiac cultures')
 txt(799,61,'External evaluation',6.5,G)
 for x,xx in [(160,241),(367,464),(615,719)]:line(x,113,xx,113)
 line(25,263,935,263,'#D9E1E5',arrow=False)
 txt(25,291,'(b)  H9c2 cellular response',8.4,R,True,'left');txt(25,312,'Cross-model context',6.3,G,False,'left')
 ax.add_patch(Ellipse((91,375),104,42,fc='#FAF0EC',ec=R,lw=.8))
 for x,y in [(63,372),(84,379),(109,367),(116,381),(83,364)]:ax.add_patch(Ellipse((x,y),15,6,angle=25,fc='#D7A191',ec=R,lw=.4))
 txt(94,438,'H9c2 cardiomyoblasts\n6 control / 6 LPS\n100 μg/mL; 48 h',6.8)
 line(155,375,199,375,R)
 txt(290,367,'Transcriptomic response\nand one-to-one orthologs',7.4,R,True)
 txt(290,429,'Gene-level correspondence\nMouse contexts:\nGSE185754 / GSE267388',6.8)
 txt(426,291,'(c)  Human cohorts',8.4,P,True,'left');txt(426,312,'External evaluation',6.3,G,False,'left')
 txt(544,352,'GSE237861\n7 sepsis / 7 controls',7.4);txt(544,407,'GSE141864\n5 sepsis / 2 controls',7.4)
 line(544,434,544,460,P);txt(544,481,'Independent disease profiles',7,P,True)
 txt(724,291,'(d)  Evidence integration',8.4,B,True,'left')
 txt(812,347,'UP suppression\nDOWN restoration',8,B,True)
 txt(812,409,'Donor and patient stability\nMultiplicity-adjusted support',6.9)
 txt(812,479,'Hallmark correspondence\nAdult cardiac cell localization',6.8,B,True)
 txt(812,514,'Post-result sensitivity',6.3,G)
 line(799,248,799,272,B);line(665,480,706,480,P)
 line(290,472,290,545,R,True,False);line(290,545,812,545,R,True,False);line(812,545,812,528,R,True)
 save(f,'Figure1_workflow')
 # Frozen evidence matrix, identical rows/criteria; only title removal and native re-export.
 e=pd.read_csv(ROOT/CFG['scientific_source']/'source_data/Figure5_evidence_matrix.tsv',sep='\t')
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
 for x,txt_ in [(1.5,'Directional evidence'),(4.5,'Robustness\nevidence'),(6.75,'Formal inference /\nsafety context')]:ax.text(x,1.14,txt_,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontweight='bold',fontsize=7.4)
 f.text(.20,.09,'Filled: criterion met   Open: criterion not met   LOO: leave one out\nPrimary: GSE237861; secondary: GSE141864. IUT/BH: 21-drug family, q < 0.05.\nNo drug passed formal inference. Historical No does not establish safety.\nND: not determined. Columns are not independent or an additive evidence score.',fontsize=6.8,va='center',linespacing=1.5)
 save(f,'Figure5_evidence_matrix')
 # Normalize raster production copies; no resampling, no metadata-only upscaling.
 audit=[]
 for p in FIG.glob('*.tiff'):
  with Image.open(io.BytesIO(p.read_bytes())) as im:
   dpi=tuple(float(v) for v in im.info.get('dpi'));size=im.size;rgb=im.convert('RGB')
  rgb.save(p,dpi=dpi,compression='tiff_lzw')
  audit.append({'name':p.name,'pixels':size,'dpi':dpi,'mode':'RGB','compression':'LZW','bytes':p.stat().st_size})
 write_json(OUT/'figure_export_audit.json',{'native_redraw':[1,5],'raster_AI_assets_included':False,'PPT_conversion_used':False,'statistics_changed':False,'figures':audit,'visual_review':'pending'})
 return {'figures':len(audit),'formats':['pdf','svg','png','tiff']}

if __name__=='__main__':run_standard_module('107_csbj_final_figures',main)
