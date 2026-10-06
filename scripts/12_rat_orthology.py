"""Packages: pandas, requests. Inputs: official NCBI references and private rat gene IDs only.
Outputs: outputs/preprocessing/rat crosswalks; no expression values sent over network.
Seed from shared config. Snapshot -> ambiguity filtering -> coverage audit; no effects.
"""
import json
import pandas as pd
from pipeline_utils import ROOT,path,acquire,write_json,sha256,run_standard_module

def main():
    cfg=json.loads((ROOT/'config/local_rat.json').read_text(encoding='utf-8'))
    dest=path('reference')/'Rattus_norvegicus.gene_info.gz'
    source=acquire(cfg['ncbi_gene_info_url'],dest,'NCBI rat gene information snapshot')
    info=pd.read_csv(dest,sep='\t',dtype=str)
    mappings=[]
    for row in info.itertuples(index=False):
        d=dict(zip(info.columns,row))
        if d['#tax_id']!='10116':continue
        for ref in d['dbXrefs'].split('|'):
            if ref.startswith('Ensembl:ENSRNOG'):
                mappings.append({'rat_ensembl':ref.split(':',1)[1],'rat_entrez':d['GeneID'],'symbol':d['Symbol']})
    allmap=pd.DataFrame(mappings).drop_duplicates()
    cross=allmap[~allmap.rat_ensembl.duplicated(keep=False)&~allmap.rat_entrez.duplicated(keep=False)]
    pairs=[]
    orth=path('reference')/'gene_orthologs.gz'
    meta=json.loads(orth.with_name(orth.name+'.provenance.json').read_text())
    if sha256(orth)!=meta['sha256']:raise ValueError('Orthology snapshot mismatch')
    for chunk in pd.read_csv(orth,sep='\t',dtype=str,chunksize=200000):
        chunk.columns=[c.lstrip('#') for c in chunk.columns]
        direct=chunk[(chunk.tax_id=='9606')&(chunk.Other_tax_id=='10116')]
        reverse=chunk[(chunk.tax_id=='10116')&(chunk.Other_tax_id=='9606')]
        if len(direct):pairs.append(direct[['GeneID','Other_GeneID']].rename(columns={'GeneID':'human_entrez','Other_GeneID':'rat_entrez'}))
        if len(reverse):pairs.append(reverse[['GeneID','Other_GeneID']].rename(columns={'GeneID':'rat_entrez','Other_GeneID':'human_entrez'}))
    if not pairs:raise ValueError('No human rat orthologs')
    allpairs=pd.concat(pairs,ignore_index=True).drop_duplicates()
    one=allpairs[~allpairs.human_entrez.duplicated(keep=False)&~allpairs.rat_entrez.duplicated(keep=False)]
    out=ROOT/'outputs/preprocessing/rat'
    cross.to_csv(out/'rat_ensembl_entrez_unambiguous.tsv',sep='\t',index=False)
    allmap.to_csv(out/'rat_ensembl_entrez_all.tsv',sep='\t',index=False)
    one.to_csv(out/'human_rat_ortholog_one_to_one.tsv',sep='\t',index=False)
    allpairs.to_csv(out/'human_rat_ortholog_all.tsv',sep='\t',index=False)
    counts=pd.read_csv(out/'rat_counts.tsv.gz',sep='\t',index_col=0)
    filtered=pd.read_csv(out/'low_expression_filter.tsv',sep='\t')
    joined=cross.merge(one,on='rat_entrez',validate='one_to_one')
    joined=joined[joined.rat_ensembl.isin(counts.index)]
    retained=set(filtered.loc[filtered.retained,'feature_id'])
    joined['passes_expression_qc_filter']=joined.rat_ensembl.isin(retained)
    joined.to_csv(out/'rat_to_human_mapping.tsv',sep='\t',index=False)
    human=pd.read_csv(ROOT/'outputs/preprocessing/human/selected_probe_per_gene.tsv',sep='\t',dtype=str)
    shared=joined[joined.passes_expression_qc_filter&joined.human_entrez.isin(human.ENTREZID)]
    result={'human_rat_one_to_one_pairs':len(one),'local_genes_mapped_to_human':len(joined),'filtered_genes_mapped_to_human':int(joined.passes_expression_qc_filter.sum()),'filtered_genes_shared_with_human_array':len(shared),'rat_source':source,'orthology_sha256':meta['sha256'],'limitation':'Current NCBI crosswalk versus vendor Ensembl release 112; unmapped genes remain excluded, no symbol-based rescue','disease_effects_examined':False}
    write_json(out/'orthology_qc.json',result)
    return {k:v for k,v in result.items() if k!='rat_source'}
if __name__=='__main__':run_standard_module('12_rat_orthology',main)
