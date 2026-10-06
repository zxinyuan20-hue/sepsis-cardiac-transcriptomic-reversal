"""Freeze ortholog, pathway and literature sources for later analysis.
Packages: requests,pandas. Inputs: NCBI/Reactome/Europe PMC public endpoints.
Outputs: immutable data/reference files; human/mouse mappings and sources in outputs.
No public redistribution licence is inferred from download accessibility.
"""
import concurrent.futures,gzip,json,logging
import pandas as pd
from pipeline_utils import ROOT,path,acquire,write_json,run_standard_module

def main():
    refs=path("reference")
    resources=[
      ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",refs/"gene_orthologs.gz","NCBI ortholog snapshot"),
      ("https://reactome.org/download/current/NCBI2Reactome.txt",refs/"NCBI2Reactome.txt","Reactome direct mapping snapshot"),
      ("https://reactome.org/download/current/ReactomePathways.txt",refs/"ReactomePathways.txt","Reactome pathway names snapshot")]
    records=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(acquire,*x) for x in resources]
        for f in concurrent.futures.as_completed(futures):records.append(f.result())
    pair=[]
    for chunk in pd.read_csv(refs/"gene_orthologs.gz",sep="\t",dtype=str,chunksize=200000):
        chunk.columns=[c.lstrip("#") for c in chunk.columns]
        direct=chunk[(chunk.tax_id=="9606")&(chunk.Other_tax_id=="10090")].copy()
        reverse=chunk[(chunk.tax_id=="10090")&(chunk.Other_tax_id=="9606")].copy()
        if len(direct):pair.append(direct[["GeneID","Other_GeneID"]].rename(columns={"GeneID":"human_entrez","Other_GeneID":"mouse_entrez"}))
        if len(reverse):pair.append(reverse[["GeneID","Other_GeneID"]].rename(columns={"GeneID":"mouse_entrez","Other_GeneID":"human_entrez"})[["human_entrez","mouse_entrez"]])
    if not pair:raise ValueError("No human-mouse ortholog pairs")
    allpairs=pd.concat(pair,ignore_index=True).drop_duplicates()
    one=allpairs[~allpairs.human_entrez.duplicated(keep=False)&~allpairs.mouse_entrez.duplicated(keep=False)]
    allpairs.to_csv(path("preparation")/"human_mouse_ortholog_all_pairs.tsv",sep="\t",index=False)
    one.to_csv(path("preparation")/"human_mouse_ortholog_one_to_one.tsv",sep="\t",index=False)
    react=pd.read_csv(refs/"NCBI2Reactome.txt",sep="\t",header=None,dtype=str,names=["entrez","pathway_id","url","pathway_name","evidence","species"])
    human=react[react.species=="Homo sapiens"].drop_duplicates(["entrez","pathway_id"])
    if human.empty:raise ValueError("No human Reactome pathways parsed")
    human.to_csv(path("preparation")/"reactome_human_direct_membership.tsv.gz",sep="\t",index=False)
    # Current direct mapping is deliberately distinguished from all-level mapping.
    literature=["EXT_ID:28067713","DOI:10.1038/s41598-023-47699-0","DOI:10.3390/ph18010043","DOI:10.3390/ph18071040","DOI:10.1016/j.cell.2017.10.049"]
    from urllib.parse import urlencode
    bib=[]
    for i,q in enumerate(literature,1):
        url="https://www.ebi.ac.uk/europepmc/webservices/rest/search?"+urlencode({"query":q,"format":"json","resultType":"core"})
        target=refs/f"literature_{i:02d}.json"
        records.append(acquire(url,target,"Europe PMC citation verification"))
        data=json.loads(target.read_text(encoding="utf-8"))
        for r in data.get("resultList",{}).get("result",[]):
            bib.append({k:r.get(k) for k in ["title","authorString","doi","pmid","pmcid","firstPublicationDate","abstractText"]})
    write_json(path("preparation")/"verified_literature.json",bib)
    summary={"human_mouse_pairs":len(allpairs),"one_to_one_pairs":len(one),"human_reactome_pathways":int(human.pathway_id.nunique()),"literature_records":len(bib),"sources":records}
    write_json(path("audit")/"reference_resources_qc.json",summary)
    return {k:v for k,v in summary.items() if k!="sources"}
if __name__=="__main__":run_standard_module("06_acquire_reference_resources",main)
