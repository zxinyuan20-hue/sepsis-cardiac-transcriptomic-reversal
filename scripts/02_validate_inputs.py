"""Audit GEO sample units, expression formats, joins and LINCS metadata.
Packages: pandas/numpy. Inputs: immutable raw GEO and LINCS metadata.
Outputs: outputs/preparation sample maps/count matrices; outputs/audit QC JSON/CSV.
Seed: config. Does NOT perform differential expression or compound ranking.
"""
import csv,gzip,io,json,logging,re
from collections import defaultdict
import numpy as np
import pandas as pd
from pipeline_utils import ROOT,CONFIG,path,write_json,run_standard_module

def soft_samples(filename):
    records=[]; current=None
    with gzip.open(filename,"rt",encoding="utf-8") as f:
        for line in f:
            if line.startswith("^"):
                if current is not None:records.append(dict(current));current=None
                if line.startswith("^SAMPLE = "):
                    current=defaultdict(list);current["gsm"]=line.strip().split(" = ",1)[1]
            elif current is not None and line.startswith("!Sample_") and " = " in line:
                k,v=line.rstrip("\r\n").split(" = ",1);current[k.removeprefix("!Sample_")].append(v)
    if current is not None:records.append(dict(current))
    return records

def meta_frame(records,acc):
    rows=[]
    for rec in records:
        chars={}
        for x in rec.get("characteristics_ch1",[]):
            if ": " in x:
                k,v=x.split(": ",1);chars[k.lower()]=v
        rows.append({"dataset":acc,"gsm":rec["gsm"],"title":" | ".join(rec.get("title",[])),
            "source":" | ".join(rec.get("source_name_ch1",[])),"organism":" | ".join(rec.get("organism_ch1",[])),
            "characteristics_json":json.dumps(chars,ensure_ascii=False),"processing":" | ".join(rec.get("data_processing",[])),**chars})
    out=pd.DataFrame(rows)
    if out.gsm.duplicated().any():raise ValueError(f"Duplicate GSM in {acc}")
    return out

def parse_series_matrix(filename):
    with gzip.open(filename,"rt",encoding="utf-8") as f:
        lines=[];inside=False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):inside=True;continue
            if line.startswith("!series_matrix_table_end"):break
            if inside:lines.append(line)
    if not lines:raise ValueError("Series matrix lacks expression table")
    return pd.read_csv(io.StringIO("".join(lines)),sep="\t",index_col=0)

def matrix_profile(df):
    arr=df.to_numpy(dtype=float)
    return {"features":len(df),"samples":len(df.columns),"duplicate_feature_ids":int(df.index.duplicated().sum()),
      "duplicate_sample_ids":int(df.columns.duplicated().sum()),"missing_values":int(np.isnan(arr).sum()),
      "nonfinite_values":int((~np.isfinite(arr)).sum()),"negative_values":int((arr<0).sum()),
      "noninteger_values":int((np.abs(arr-np.rint(arr))>1e-6).sum()),
      "quantiles":dict(zip(["min","q01","median","q99","max"],map(float,np.nanquantile(arr,[0,.01,.5,.99,1])))),
      "exact_duplicate_sample_pairs":[[str(df.columns[i]),str(df.columns[j])] for i in range(df.shape[1]) for j in range(i) if np.array_equal(arr[:,i],arr[:,j],equal_nan=True)]}

def parse_dose_um(values):
    """Accept only explicitly micromolar dose units; do not reinterpret mg/mL."""
    return pd.to_numeric(values.astype(str).str.extract(r"(?i)^([0-9.eE+-]+)\s+um$")[0],errors="coerce")

def align_mouse(df,meta,acc):
    mapping=[]
    for _,r in meta.iterrows():
        if acc=="GSE185754":
            m=re.search(r"_(\d+)\s*$",r.title)
            if not m:raise ValueError(f"No replicate suffix: {r.title}")
            group="control" if r.get("stress")=="saline" else "LPS" if "lps" in str(r.get("stress")).lower() else None
            col=("Saline_" if group=="control" else "LPS_24h_")+m[1]+"_count"
        else:
            m=re.fullmatch(r"Heart_(PBS|LPS)_rep(\d+)",r.title)
            if not m:raise ValueError(f"Unknown sample title: {r.title}")
            group="control" if m[1]=="PBS" else "LPS";col=f"read_count_WT_{m[1]}_{m[2]}"
        if group is None or col not in df:raise ValueError(f"Unmatched sample {r.gsm}: {col}")
        mapping.append({"dataset":acc,"gsm":r.gsm,"title":r.title,"group":group,"expression_column":col,"mapping_evidence":"GEO sample title/characteristics and exact expression column label"})
    mp=pd.DataFrame(mapping)
    if mp.expression_column.duplicated().any():raise ValueError("Multiple GSMs mapped to same column")
    cols=[c for c in df if c.endswith("_count")] if acc=="GSE185754" else [c for c in df if c.startswith("read_count_WT_")]
    if set(cols)!=set(mp.expression_column):raise ValueError("Unmapped count columns")
    expr=df[mp.expression_column].copy();expr.columns=mp.gsm
    return mp,expr

def read_metrics(p):
    # GEO release has an extra row index field without a matching header field.
    with gzip.open(p,"rt") as f:
        header=f.readline().rstrip().split("\t");first=f.readline().rstrip().split("\t")
    if len(first)==len(header)+1:
        df=pd.read_csv(p,sep="\t",index_col=0)
    elif len(first)==len(header):df=pd.read_csv(p,sep="\t")
    else:raise ValueError("Unexpected LINCS metrics header width")
    if not df.sig_id.astype(str).str.contains(":").all():raise ValueError("LINCS sig_id parsing failed")
    return df

def audit_lincs():
    folder=path("raw")/"lincs/GSE70138"
    def read(s):return pd.read_csv(folder/f"GSE70138_Broad_LINCS_{s}.txt.gz",sep="\t",low_memory=False)
    sig=read("sig_info_2017-03-06");met=read_metrics(folder/"GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz")
    gene=read("gene_info_2017-03-06");pert=read("pert_info_2017-03-06");cells=read("cell_info_2017-04-28")
    for name,df,key in [("sig",sig,"sig_id"),("metrics",met,"sig_id"),("genes",gene,"pr_gene_id"),("pert",pert,"pert_id"),("cells",cells,"cell_id")]:
        if df[key].duplicated().any():raise ValueError(f"Duplicate {key} in {name}")
    quality_cols=["sig_id","distil_cc_q75","distil_ss","distil_nsample","tas"]
    merged=sig.merge(met[quality_cols],on="sig_id",how="left",validate="one_to_one",indicator=True)
    merged["dose_um"]=parse_dose_um(merged.pert_idose)
    merged["time_h"]=pd.to_numeric(merged.pert_itime.str.extract(r"^([0-9.eE+-]+)\s+h$")[0],errors="coerce")
    c=CONFIG["lincs"]
    eligible=(merged.pert_type.eq(c["pert_type"]) & merged.distil_cc_q75.ge(c["min_distil_cc_q75"]) & merged.distil_ss.ge(c["min_distil_ss"]) & merged.distil_nsample.ge(c["min_distil_nsample"]) & merged.dose_um.gt(0) & merged.dose_um.le(c["primary_max_dose_um"]) & merged.time_h.isin(c["primary_times_hours"]))
    filtered=merged.loc[eligible].copy()
    if filtered.empty:raise ValueError("No eligible LINCS signatures: inspect units, metadata and gates; do not continue")
    by=filtered.groupby("pert_id").agg(signatures=("sig_id","size"),cell_lines=("cell_id","nunique"))
    eligible_drugs=by.index[(by.signatures>=c["minimum_high_quality_signatures"])&(by.cell_lines>=c["minimum_distinct_cell_lines"])]
    filtered["compound_replicability_eligible"]=filtered.pert_id.isin(eligible_drugs)
    filtered.to_csv(path("preparation")/"lincs_eligible_signatures.tsv.gz",sep="\t",index=False)
    by.to_csv(path("preparation")/"lincs_compound_coverage.tsv",sep="\t")
    gene.to_csv(path("preparation")/"lincs_gene_dictionary.tsv",sep="\t",index=False)
    pert.to_csv(path("preparation")/"lincs_compound_dictionary.tsv",sep="\t",index=False)
    summary={"release":"GSE70138_2017-03-06","signatures":len(sig),"genes":len(gene),"landmark_genes":int(gene.pr_is_lm.eq(1).sum()),
        "chemical_perturbagens":int(pert.pert_type.eq("trt_cp").sum()),"metrics_missing":int(merged._merge.ne("both").sum()),
        "eligible_signatures":len(filtered),"eligible_compounds_before_multicell_gate":int(filtered.pert_id.nunique()),
        "eligible_compounds_after_multicell_gate":len(eligible_drugs),"eligible_cell_lines":sorted(filtered.cell_id.unique().tolist()),
        "quality_filter":c["signature_quality"],"is_hiq_available":("is_hiq" in met),
        "sentinel_minus666_cc":int(merged.distil_cc_q75.eq(-666).sum()),
        "chemical_signatures_with_unparsed_dose":int((merged.pert_type.eq("trt_cp") & merged.dose_um.isna()).sum()),
        "note":"Coverage only; no disease signatures or drug rankings computed. Tumour-cell activity is not cardiac efficacy."}
    write_json(path("audit")/"lincs_metadata_qc.json",summary)
    return summary

def main():
    summaries=[];metadata={}
    for d in CONFIG["datasets"]:
        acc=d["accession"];meta=meta_frame(soft_samples(path("raw")/"geo"/acc/f"{acc}_family.soft.gz"),acc)
        if len(meta)!=d["expected_samples"]:raise ValueError(f"Unexpected sample count {acc}: {len(meta)}")
        meta.to_csv(path("preparation")/f"{acc}_sample_metadata.tsv",sep="\t",index=False);metadata[acc]=meta
    acc="GSE79962";meta=metadata[acc]
    labels={"nonischemic dilated cardiomyopathy":"DCM","ischemic heart disease":"IHD","nonfailing heart":"nonfailing","septic cardiomyopathy":"sepsis"}
    meta["group"]=meta.condition.map(labels)
    if meta.group.isna().any():raise ValueError(f"Unknown human condition: {meta.condition.unique()}")
    expr=parse_series_matrix(path("raw")/"geo"/acc/f"{acc}_series_matrix.txt.gz")
    if set(expr.columns)!=set(meta.gsm):raise ValueError("Human matrix/GSM mismatch")
    hp=matrix_profile(expr)
    if hp["nonfinite_values"] or hp["duplicate_feature_ids"] or hp["duplicate_sample_ids"] or hp["exact_duplicate_sample_pairs"]:
        raise ValueError("Human matrix failed uniqueness/completeness gate")
    meta=meta.set_index("gsm").loc[expr.columns].reset_index();meta=meta.rename(columns={meta.columns[0]:"gsm"}) if "gsm" not in meta else meta
    meta.to_csv(path("preparation")/f"{acc}_sample_map.tsv",sep="\t",index=False)
    expr.to_csv(path("preparation")/f"{acc}_submitted_expression.tsv.gz",sep="\t",index_label="probe_id")
    summaries.append({"dataset":acc,"kind":"submitted normalized microarray","groups":meta.group.value_counts().to_dict(),**matrix_profile(expr)})
    for acc,name,key in [("GSE185754","GSE185754_Raw_gene_counts_matrix.txt.gz","gene_id"),("GSE267388","GSE267388_gene_expression.xls.gz","gene_id")]:
        data=pd.read_csv(path("raw")/"geo"/acc/name,sep="\t",index_col=key,low_memory=False)
        mp,expr=align_mouse(data,metadata[acc],acc)
        prof=matrix_profile(expr)
        if prof["nonfinite_values"] or prof["negative_values"] or prof["duplicate_feature_ids"] or prof["exact_duplicate_sample_pairs"]:raise ValueError(f"Invalid count data in {acc}")
        mp.to_csv(path("preparation")/f"{acc}_sample_map.tsv",sep="\t",index=False)
        expr.to_csv(path("preparation")/f"{acc}_counts.tsv.gz",sep="\t",index_label="gene_id")
        ann=[x for x in data if x in ["gene_name","gene_symbol","gene_biotype","gene_chr","gene_description"]]
        data[ann].to_csv(path("preparation")/f"{acc}_gene_annotations.tsv.gz",sep="\t")
        summaries.append({"dataset":acc,"kind":"estimated counts" if prof["noninteger_values"] else "integer counts","groups":mp.group.value_counts().to_dict(),**prof})
    with gzip.open(path("raw")/"geo/GSE79962/GSE79962_Clinical_metadata.xls.gz","rb") as f:
        excel=pd.ExcelFile(io.BytesIO(f.read()),engine="xlrd")
        clinical={"sheets":excel.sheet_names}
        for sheet in excel.sheet_names:
            df=pd.read_excel(excel,sheet_name=sheet,header=None)
            clinical[sheet]={"shape":list(df.shape),"gsm_ids":sorted(set(re.findall(r"GSM\d+"," ".join(df.fillna("").astype(str).to_numpy().ravel()))))}
            df.to_csv(path("preparation")/f"GSE79962_clinical_sheet_{len(clinical)-2}.tsv",sep="\t",index=False,header=False)
    write_json(path("audit")/"geo_input_qc.json",{"datasets":summaries,"clinical_workbook":clinical,
        "single_cell":"Metadata only; CD45 selection, pooling and genotype/time design require review before using cell proportions or group tests."})
    lincs=audit_lincs()
    return {"geo_samples":sum(len(m) for m in metadata.values()),"bulk_profiles":summaries,"lincs":lincs}

if __name__=="__main__":run_standard_module("02_validate_inputs",main)
