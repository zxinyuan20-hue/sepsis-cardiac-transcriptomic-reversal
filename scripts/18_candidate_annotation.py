"""Packages: pandas. Inputs frozen rankings, dated Broad annotation, human gene effects.
Outputs screening annotation and target-context tables; exact-name lookup only.
No rat effects, no efficacy inference, no remote queries containing private data.
"""
import json
import pandas as pd
from pipeline_utils import ROOT,path,acquire,write_json,sha256,run_standard_module

def main():
    cfg=json.loads((ROOT/'config/compound_annotation.json').read_text())
    out=path('analysis')/'screening'
    freeze=json.loads((out/'screening_result_freeze.json').read_text())
    if sha256(out/'complete_screening_results.tsv')!=freeze['complete_results_sha256']:raise ValueError('Frozen ranking changed')
    source=acquire(cfg['source_url'],path('reference')/cfg['local_file'],'Broad Repurposing Hub dated mechanism annotation')
    annotation=pd.read_csv(path('reference')/cfg['local_file'],sep='\t',comment='!',dtype=str).fillna('')
    annotation['match_name']=annotation.pert_iname.str.strip().str.casefold()
    # Preserve any multiple records as combined database fields, never multiply compounds.
    fields=['moa','target','clinical_phase','disease_area','indication']
    unique=annotation.groupby('match_name')[fields].agg(lambda s:' | '.join(sorted(set(v for v in s if v))))
    rank=pd.read_csv(out/'complete_screening_results.tsv',sep='\t');rank['match_name']=rank.pert_iname.str.strip().str.casefold()
    rank=rank.merge(unique,on='match_name',how='left',validate='many_to_one')
    rank['annotation_match']=rank.match_name.isin(unique.index)
    rank['annotation_date']=cfg['source_date'];rank['identity_verification']='exact_name_only_not_structure_confirmed'
    rank.to_csv(out/'annotated_complete_results.tsv',sep='\t',index=False)
    rank[rank['rank']<=10].to_csv(out/'annotated_top10.tsv',sep='\t',index=False)
    rank[rank.supported_candidate].head(10).to_csv(out/'annotated_supported_shortlist.tsv',sep='\t',index=False)
    mapping=pd.read_csv(ROOT/'outputs/preprocessing/human/probe_mapping_all.tsv',sep='\t',dtype=str).dropna(subset=['ENTREZID','SYMBOL'])
    mapping=mapping[['SYMBOL','ENTREZID']].drop_duplicates();valid=mapping.groupby('SYMBOL').ENTREZID.nunique();mapping=mapping[mapping.SYMBOL.isin(valid[valid==1].index)].drop_duplicates('SYMBOL').set_index('SYMBOL')
    effects={name:pd.read_csv(path('analysis')/f'human/primary_{name}.tsv',sep='\t',dtype={'entrez_id':str}).set_index('entrez_id') for name in ['sepsis_vs_nonfailing','sepsis_vs_IHD','sepsis_vs_DCM']}
    candidates=rank[(rank['rank']<=10)|rank.supported_candidate]
    records=[]
    for row in candidates.itertuples():
        targets=[] if pd.isna(row.target) else sorted(set(t.strip() for t in row.target.split('|') if t.strip()))
        for target in targets:
            entrez=mapping.loc[target,'ENTREZID'] if target in mapping.index else None
            rec={'rank':row.rank,'pert_id':row.pert_id,'pert_iname':row.pert_iname,'target_symbol':target,'entrez_id':entrez,'database_moa':row.moa,'target_mapping':'unambiguous_symbol' if entrez else 'unmapped_or_ambiguous'}
            for name,df in effects.items():
                rec[name+'_logFC']=df.loc[entrez,'logFC'] if entrez in df.index else None
                rec[name+'_FDR']=df.loc[entrez,'adj.P.Val'] if entrez in df.index else None
            records.append(rec)
    pd.DataFrame(records).to_csv(out/'candidate_target_expression_context.tsv',sep='\t',index=False)
    summary={'annotation_snapshot':cfg['source_date'],'source_sha256':source['sha256'],'source_url':cfg['source_url'],'matched_compounds':int(rank.annotation_match.sum()),'total_compounds':len(rank),'top10_matched':int(rank.loc[rank['rank']<=10,'annotation_match'].sum()),'target_rows':len(records),'scope':'Database MoA/target annotation plus human transcript levels; target activity, compound identity, cardiac toxicity and current development status not yet independently validated; original rankings unchanged'}
    write_json(out/'annotation_audit.json',summary)
    return summary
if __name__=='__main__':run_standard_module('18_candidate_annotation',main)
