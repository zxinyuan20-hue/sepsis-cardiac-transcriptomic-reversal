"""Packages pandas,numpy,scipy; R edgeR/limma. Fixed seed20261002.
Inputs: config/raw_heart_quantification_v1.json and audited14counts.
Outputs under raw_heart_v1; original results immutable. No drug re-ranking.
"""
import gzip,itertools,json,re,subprocess
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from pipeline_utils import ROOT,RUNTIME,sha256,write_json,run_standard_module,r_environment,now
OUT=ROOT/'outputs/analysis/raw_heart_v1'
def bh(values):
    v=np.asarray(values,float);ix=np.argsort(v,kind='stable');r=np.empty_like(v)
    r[ix]=np.minimum(1,np.minimum.accumulate((v[ix]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);return r
def main():
    qa=json.loads((OUT/'quantification_audit.json').read_text());assert qa['disease_inference_ready'],qa
    assert sha256(OUT/'gene_counts_14.tsv.gz')==qa['counts_sha256']
    freeze=json.loads((OUT/'quantification_plan_freeze.json').read_text());assert sha256(ROOT/'config/raw_heart_quantification_v1.json')==freeze['signature']['config_sha256']
    for p,h in freeze['signature']['inherited_analysis_inputs'].items():assert sha256(ROOT/p)==h,p
    for f in ['GSE141864_patient_signature_scores.tsv','external_human_summary.json']:
        assert (ROOT/'outputs/analysis/external_human'/f).exists()
    # Transcript IDs match exact annotation versions; no stripping/rescue mapping.
    tx=pd.read_csv(OUT/'sources/gencode.v47.metadata.EntrezGene.gz',sep='\t',header=None,names=['transcript_id','entrez_id'],dtype=str)
    lookup=tx.groupby('transcript_id').entrez_id.agg(set).to_dict();mapping={}
    with (OUT/'reference/gencode.v47.primary_assembly.annotation.gtf').open() as f:
        for line in f:
            if line.startswith('#'):continue
            a=line.rstrip().split('\t')
            if a[2]!='transcript':continue
            att=dict(re.findall(r'(\w+) "([^"]+)"',a[8]));gid=att['gene_id'];tid=att['transcript_id']
            mapping.setdefault(gid,set()).update(lookup.get(tid,set()))
    rows=[{'gene_id':g,'entrez_id':next(iter(v)) if len(v)==1 else '', 'n_entrez':len(v)} for g,v in mapping.items()]
    mp=pd.DataFrame(rows);mp.to_csv(OUT/'gencode_gene_entrez_mapping.tsv',sep='\t',index=False)
    counts=pd.read_csv(OUT/'gene_counts_14.tsv.gz',sep='\t',index_col=0)
    unique=mp.loc[mp.n_entrez.eq(1)].set_index('gene_id').entrez_id
    counts['entrez_id']=counts.index.map(unique);entrez=counts.dropna(subset=['entrez_id']).groupby('entrez_id').sum()
    assert len(entrez)>10000 and entrez.shape[1]==14
    entrez.to_csv(OUT/'entrez_counts_14.tsv.gz',sep='\t',compression='gzip')
    env=r_environment();r=subprocess.run([RUNTIME['rscript'],'--vanilla',str(ROOT/'scripts/49_analyze_reconstructed_heart.R')],capture_output=True,text=True,encoding='utf8',errors='replace',env=env,cwd=ROOT)
    (OUT/'R_analysis_console.txt').write_text(r.stdout+'\n'+r.stderr,encoding='utf8');assert r.returncode==0,r.stderr[-2500:]
    x=pd.read_csv(OUT/'logCPM_14.tsv.gz',sep='\t',index_col=0);x.index=x.index.astype(str)
    h=pd.read_csv(ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
    s=pd.read_csv(OUT/'patient_metadata.tsv',sep='\t');ids=s.title.tolist();case=s.group.eq('sepsis').to_numpy();bg=sorted(set(h.index)&set(x.index))
    queries={size:{d:pd.read_csv(ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_{size}_{d}.tsv',sep='\t',dtype={'entrez_id':str}).entrez_id.tolist() for d in ['up','down']} for size in [50,100]}
    coverage={d:len(set(queries[100][d])&set(bg))/len(queries[100][d]) for d in ['up','down']}
    write_json(OUT/'query_coverage.json',coverage);assert min(coverage.values())>=.70,coverage
    ranks=(x.loc[bg,ids].rank(axis=0,method='average')-.5)/len(bg)
    score=lambda size:(ranks.loc[ranks.index.isin(queries[size]['up'])].mean()-ranks.loc[ranks.index.isin(queries[size]['down'])].mean()).to_numpy()
    v=score(100);effect=float(v[case].mean()-v[~case].mean());null=[]
    for sel in itertools.combinations(range(14),7):
        g=np.zeros(14,bool);g[list(sel)]=True;null.append(float(v[g].mean()-v[~g].mean()))
    assert len(null)==3432
    p=float(np.mean(np.asarray(null)>=effect-1e-14))
    pd.DataFrame({'case_minus_control':null}).to_csv(OUT/'unadjusted_permutation_distribution.tsv',sep='\t',index=False)
    female=s.sex.eq('F').to_numpy();coef=lambda g:float(np.linalg.lstsq(np.c_[np.ones(14),female.astype(int),g.astype(int)],v,rcond=None)[0][-1])
    adjusted=coef(case);strat=[];males=np.where(~female)[0];females=np.where(female)[0]
    for men in itertools.combinations(males,int(case[~female].sum())):
        for women in itertools.combinations(females,int(case[female].sum())):
            g=np.zeros(14,bool);g[list(men)+list(women)]=True;strat.append(coef(g))
    assert len(strat)==840
    ps=float(np.mean(np.asarray(strat)>=adjusted-1e-14))
    pd.DataFrame({'sex_adjusted_group_coefficient':strat}).to_csv(OUT/'sex_stratified_permutation_distribution.tsv',sep='\t',index=False)
    s['primary_disease_rank_score']=v;s['query50_score']=score(50);s.to_csv(OUT/'patient_signature_scores.tsv',sep='\t',index=False)
    loo=[]
    for i in range(14):
        k=np.arange(14)!=i;z=v[k];g=case[k];loo.append({'omitted_title':ids[i],'case_minus_control':float(z[g].mean()-z[~g].mean())})
    pd.DataFrame(loo).to_csv(OUT/'leave_one_patient_out.tsv',sep='\t',index=False)
    old=json.loads((ROOT/'outputs/analysis/external_human/external_human_summary.json').read_text())
    family=pd.DataFrame({'cohort':['GSE141864','GSE237861'],'p':[old['exact_one_sided_p'],p],'effect':[old['primary_rank_difference'],effect]});family['BH_q']=bh(family.p);family.to_csv(OUT/'two_cohort_primary_family.tsv',sep='\t',index=False)
    genes=pd.read_csv(OUT/'GSE237861_all_gene_results.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id')
    joint=h[['logFC']].join(genes[['logFC','FDR']],how='inner',lsuffix='_discovery',rsuffix='_external');joint.to_csv(OUT/'gene_effect_comparison.tsv',sep='\t')
    newpath=pd.read_csv(OUT/'GSE237861_fixed_pathway_tests.tsv',sep='\t');newpath['cohort']='GSE237861'
    oldpath=pd.read_csv(ROOT/'outputs/analysis/external_human/GSE141864_fixed_pathway_tests.tsv',sep='\t');oldpath['cohort']='GSE141864'
    paths=pd.concat([oldpath,newpath],ignore_index=True)
    for mode,idx in paths.groupby('mode').groups.items():
        assert len(idx)==10;paths.loc[idx,'planned_10_test_BH_q']=bh(paths.loc[idx,'PValue'].fillna(1))
    paths.to_csv(OUT/'two_cohort_pathway_family.tsv',sep='\t',index=False)
    result={'version':'raw-heart-1.0','completed_utc':now(),'samples':14,'entrez_mapped_genes':len(entrez),'filtered_genes':len(x),'shared_genes':len(bg),'query_coverage':coverage,'primary_effect':effect,'rank_percentage_points':effect*100,'primary_one_sided_exact_p':p,'two_cohort_BH_q':float(family.loc[family.cohort.eq('GSE237861'),'BH_q'].iloc[0]),'primary_permutations':3432,'sex_adjusted_effect':adjusted,'sex_stratified_sensitivity_p':ps,'sex_stratified_permutations':840,'all_loo_effects_positive':all(r['case_minus_control']>0 for r in loo),'query50_effect':float(score(50)[case].mean()-score(50)[~case].mean()),'gene_logFC_spearman_descriptive':float(spearmanr(joint.logFC_discovery,joint.logFC_external).statistic),'genes_FDR05':int(genes.FDR.lt(.05).sum()),'therapeutic_candidates_supported':0,'drug_tests_performed':False,'release_status':'AWAITING_INDEPENDENT_NUMERICAL_AUDIT'}
    write_json(OUT/'disease_replication_summary.json',result);return result
if __name__=='__main__':run_standard_module('49_analyze_reconstructed_heart',main)
