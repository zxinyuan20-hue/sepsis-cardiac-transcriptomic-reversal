"""Packages: numpy/pandas/scipy/sklearn/matplotlib; R edgeR for rat input QC.
Inputs: outputs/preprocessing/{human,rat}; outputs: same and outputs/figures/preprocessing.
Seed: study config. PCA/correlation/distribution, technical flags, no automatic exclusions.
"""
import json,subprocess,logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from pipeline_utils import ROOT,CONFIG,RUNTIME,path,sha256,write_json,r_environment,run_standard_module

COLORS={'control':'#4477AA','lps':'#CC6677','nonfailing':'#4477AA','sepsis':'#CC6677','IHD':'#228833','DCM':'#AA3377'}

def dataset_qc(mat,samples,folder,label):
    samples=samples.set_index('sample_id').loc[list(mat.columns)].rename_axis('sample_id').reset_index()
    a=mat.to_numpy(float)
    if not np.isfinite(a).all():raise ValueError('Nonfinite normalized values')
    corr=np.corrcoef(a.T)
    if not np.isfinite(corr).all():raise ValueError('Invalid correlation')
    off=corr.copy();np.fill_diagonal(off,np.nan)
    score=PCA(n_components=5,svd_solver='full').fit(a.T)
    pcs=score.transform(a.T)
    metrics=samples[['sample_id','group']].copy()
    metrics['median_expression']=np.median(a,axis=0)
    metrics['iqr_expression']=np.quantile(a,.75,axis=0)-np.quantile(a,.25,axis=0)
    metrics['median_sample_correlation']=np.nanmedian(off,axis=1)
    center=np.median(metrics.median_sample_correlation)
    mad=np.median(np.abs(metrics.median_sample_correlation-center))*1.4826
    cutoff=center-3*mad
    metrics['low_correlation_review_flag']=metrics.median_sample_correlation<cutoff
    for k in range(5):metrics[f'PC{k+1}']=pcs[:,k]
    metrics.to_csv(folder/'sample_qc_metrics.tsv',sep='\t',index=False)
    pd.DataFrame(corr,index=mat.columns,columns=mat.columns).to_csv(folder/'sample_correlations.tsv',sep='\t')
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),layout='constrained')
    for group in samples.group.unique():
        mask=(samples.group==group).values
        axes[0].scatter(pcs[mask,0],pcs[mask,1],s=42,color=COLORS[group],label=f'{group} (n={mask.sum()})',edgecolor='white',linewidth=.4)
    if len(samples)<=15:
        for i,sid in enumerate(samples.sample_id):axes[0].annotate(sid,(pcs[i,0],pcs[i,1]),xytext=(4,4),textcoords='offset points',fontsize=7)
    axes[0].set(xlabel=f'PC1 ({score.explained_variance_ratio_[0]:.1%})',ylabel=f'PC2 ({score.explained_variance_ratio_[1]:.1%})',title='A  Sample structure')
    axes[0].legend(fontsize=8,frameon=False)
    boxes=axes[1].boxplot(a,showfliers=False,patch_artist=True,widths=.65)
    for patch,group in zip(boxes['boxes'],samples.group):patch.set_facecolor(COLORS[group]);patch.set_alpha(.7)
    axes[1].set(xlabel='Samples in manifest order',ylabel='Log2 expression' if label=='Human' else 'TMM log2 CPM',title='B  Expression distributions')
    axes[1].set_xticks(range(1,len(samples)+1));axes[1].set_xticklabels(samples.sample_id if len(samples)<=15 else range(1,len(samples)+1),rotation=90,fontsize=6)
    im=axes[2].imshow(corr,cmap='viridis',vmin=float(np.nanmin(off)),vmax=1,interpolation='nearest')
    axes[2].set(title='C  Sample correlations',xlabel='Samples in manifest order',ylabel='Samples in manifest order')
    fig.colorbar(im,ax=axes[2],shrink=.8,label='Pearson r')
    fig.suptitle(label+' expression quality overview',fontsize=13)
    figures=ROOT/'outputs/figures/preprocessing';figures.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(figures/(label.lower()+'_qc.'+ext),dpi=300,facecolor='white')
    plt.close(fig)
    result={'samples':len(samples),'genes':len(mat),'pca_variance_ratio':score.explained_variance_ratio_.tolist(),'minimum_pair_correlation':float(np.nanmin(off)),'median_pair_correlation':float(np.nanmedian(off)),'review_flags':metrics.loc[metrics.low_correlation_review_flag,'sample_id'].tolist(),'correlation_review_cutoff':float(cutoff),'flag_rule':'median sample correlation below global median minus 3 scaled MAD; review only, disease groups can differ','excluded_samples':[]}
    write_json(folder/'sample_qc_summary.json',result)
    return result

def main():
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    with (path('logs')/'10_R_rat_qc_details.log').open('a',encoding='utf-8') as log:
        p=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/10_rat_qc.R')],cwd=ROOT,env=r_environment(),stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError('Rat QC R failed')
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'font.size':9})
    h=ROOT/'outputs/preprocessing/human';r=ROOT/'outputs/preprocessing/rat'
    human=pd.read_csv(h/'human_rma_genes.tsv.gz',sep='\t',index_col=0)
    hs=pd.read_csv(h/'cel_sample_manifest.tsv',sep='\t').rename(columns={'gsm':'sample_id'})
    rat=pd.read_csv(r/'rat_tmm_logcpm_qc.tsv.gz',sep='\t',index_col=0)
    rs=pd.read_csv(r/'sample_manifest.tsv',sep='\t')
    results={'human':dataset_qc(human,hs,h,'Human'),'rat':dataset_qc(rat,rs,r,'Rat')}
    # Same-sample comparison with author's submitted matrix, common probes only.
    probes=pd.read_csv(h/'human_rma_probes.tsv.gz',sep='\t',index_col=0)
    submitted=pd.read_csv(path('preparation')/'GSE79962_submitted_expression.tsv.gz',sep='\t',index_col=0)
    probes.index=probes.index.astype(str);submitted.index=submitted.index.astype(str)
    shared=probes.index.intersection(submitted.index)
    if len(shared)<1000 or set(probes.columns)!=set(submitted.columns):raise ValueError('Submitted/RMA identity mismatch')
    pairs=[{'sample_id':sid,'common_probes':len(shared),'spearman_r':float(spearmanr(probes.loc[shared,sid],submitted.loc[shared,sid]).statistic)} for sid in probes.columns]
    pd.DataFrame(pairs).to_csv(h/'submitted_rma_concordance.tsv',sep='\t',index=False)
    results['submitted_comparison']={'common_probes':len(shared),'min_spearman_r':min(x['spearman_r'] for x in pairs),'median_spearman_r':float(np.median([x['spearman_r'] for x in pairs])),'note':'preprocessing pipelines differ; rank agreement is not scientific validation'}
    write_json(path('audit')/'preprocessing_qc.json',results)
    return results
if __name__=='__main__':run_standard_module('10_preprocessing_qc',main)
