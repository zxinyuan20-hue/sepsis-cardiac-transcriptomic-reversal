# Functional smoke checks, NOT scientific differential expression.
# Input: prepared matrices, project libraries and one original CEL.
# Output: outputs/audit/r_functional_checks.json, probe mapping, environment/renv.lock.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
lib <- file.path(root,"environment","R-library")
.libPaths(c(lib,.Library))
set.seed(as.integer(Sys.getenv("STUDY_SEED")))
suppressPackageStartupMessages({library(limma);library(edgeR);library(AnnotationDbi);library(hugene10sttranscriptcluster.db);library(oligo)})
audit <- file.path(root,"outputs","audit")
prep <- file.path(root,"outputs","preparation")
# Synthetic contrast test catches a reversed coefficient direction.
group <- factor(rep(c("control","treated"),each=5),levels=c("control","treated"))
y <- matrix(rnorm(200*10,sd=.2),nrow=200)
y[1:10,group=="treated"] <- y[1:10,group=="treated"]+3
design <- model.matrix(~group)
fit <- eBayes(lmFit(y,design))
stopifnot(all(fit$coefficients[1:10,2]>2))
counts <- as.matrix(read.delim(gzfile(file.path(prep,"GSE267388_counts.tsv.gz")),row.names=1,check.names=FALSE))
storage.mode(counts) <- "double"
dg <- DGEList(counts=counts)
dg <- calcNormFactors(dg)
stopifnot(all(is.finite(dg$samples$norm.factors)),all(dg$samples$norm.factors>0))
# Map existing probe identifiers without inspecting disease effects.
human <- read.delim(gzfile(file.path(prep,"GSE79962_submitted_expression.tsv.gz")),row.names=1,check.names=FALSE)
probes <- rownames(human)
valid <- intersect(probes,keys(hugene10sttranscriptcluster.db,keytype="PROBEID"))
mapping <- AnnotationDbi::select(hugene10sttranscriptcluster.db,keys=valid,keytype="PROBEID",columns=c("ENTREZID","SYMBOL"))
write.table(mapping,file.path(prep,"human_probe_gene_mapping.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="")
lm <- read.delim(file.path(prep,"lincs_gene_dictionary.tsv"),check.names=FALSE)
mapped <- unique(na.omit(mapping$ENTREZID))
landmark <- as.character(lm$pr_gene_id[lm$pr_is_lm==1])
cel <- list.files(file.path(prep,"cel_smoke"),pattern="[.][Cc][Ee][Ll]([.]gz)?$",full.names=TRUE)
if(length(cel)!=1)stop("Expected exactly one CEL smoke input")
arr <- oligo::read.celfiles(cel,verbose=FALSE)
result <- list(limma_synthetic_direction_test=TRUE,edgeR_estimated_count_input=TRUE,
   human_probes=length(probes),annotated_probe_ids=length(unique(mapping$PROBEID)),mapped_entrez_genes=length(mapped),
   platform_covered_lincs_landmarks=length(intersect(mapped,landmark)),lincs_landmarks_total=length(landmark),
   cel_read_success=TRUE,cel_probe_features=nrow(arr),cel_annotation=annotation(arr),
   note="Software/input preparation only; no disease contrast or compound rank produced. Mapping includes ambiguous rows for later filtering.")
jsonlite::write_json(result,file.path(audit,"r_functional_checks.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(audit,"R_functional_sessionInfo.txt"))
options(repos=c(CRAN="https://cloud.r-project.org"),renv.config.auto.snapshot=FALSE)
renv::snapshot(project=root,lockfile=file.path(root,"environment","renv.lock"),library=c(lib,.Library),
  packages=c("BiocManager","limma","edgeR","oligo","pd.hugene.1.0.st.v1","AnnotationDbi","org.Hs.eg.db","org.Mm.eg.db","hugene10sttranscriptcluster.db","jsonlite","data.table","ggplot2","readxl","renv"),prompt=FALSE,force=TRUE)
cat("R_FUNCTIONAL_CHECKS_PASS\n")
