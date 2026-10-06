# Packages: oligo, AnnotationDbi, hugene10sttranscriptcluster.db, jsonlite.
# Input/output paths: project root, outputs/preprocessing/human. Seed from config.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
root <- normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library))
set.seed(as.integer(Sys.getenv('STUDY_SEED')))
suppressPackageStartupMessages({library(oligo);library(AnnotationDbi);library(hugene10sttranscriptcluster.db)})
out <- file.path(root,'outputs','preprocessing','human')
samples <- read.delim(file.path(out,'cel_sample_manifest.tsv'),check.names=FALSE)
raw <- read.celfiles(file.path(root,samples$cel_path),verbose=TRUE)
sampleNames(raw) <- samples$gsm
raw_quantiles <- t(apply(exprs(raw),2,quantile,probs=c(.01,.25,.5,.75,.99),na.rm=TRUE))
write.table(data.frame(gsm=samples$gsm,raw_quantiles,check.names=FALSE),file.path(out,'raw_intensity_quantiles.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
eset <- rma(raw,target='core')
x <- exprs(eset)
stopifnot(ncol(x)==51,all(is.finite(x)),identical(colnames(x),samples$gsm))
saveRDS(eset,file.path(out,'human_rma_expression_set.rds'))
write_matrix <- function(mat,name){con<-gzfile(file.path(out,name),'wt');write.table(data.frame(feature_id=rownames(mat),mat,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE);close(con)}
write_matrix(x,'human_rma_probes.tsv.gz')
valid <- intersect(rownames(x),keys(hugene10sttranscriptcluster.db,keytype='PROBEID'))
mapping <- unique(AnnotationDbi::select(hugene10sttranscriptcluster.db,keys=valid,keytype='PROBEID',columns=c('ENTREZID','SYMBOL')))
mapping <- mapping[!is.na(mapping$ENTREZID),]
ngenes <- tapply(mapping$ENTREZID,mapping$PROBEID,function(z)length(unique(z)))
unambiguous <- names(ngenes)[ngenes==1]
eligible <- unique(mapping[mapping$PROBEID %in% unambiguous,c('PROBEID','ENTREZID')])
eligible$mean_expression <- rowMeans(x)[eligible$PROBEID]
eligible <- eligible[order(-eligible$mean_expression,eligible$PROBEID),]
selected <- eligible[!duplicated(eligible$ENTREZID),]
selected <- selected[order(as.numeric(selected$ENTREZID)),]
genes <- x[selected$PROBEID,,drop=FALSE];rownames(genes)<-selected$ENTREZID
write_matrix(genes,'human_rma_genes.tsv.gz')
write.table(selected,file.path(out,'selected_probe_per_gene.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(mapping,file.path(out,'probe_mapping_all.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
lm <- read.delim(file.path(root,'outputs','preparation','lincs_gene_dictionary.tsv'))
covered <- intersect(rownames(genes),as.character(lm$pr_gene_id[lm$pr_is_lm==1]))
jsonlite::write_json(list(samples=ncol(x),core_probes=nrow(x),probes_without_entrez=nrow(x)-length(unique(mapping$PROBEID)),ambiguous_probes=sum(ngenes>1),unambiguous_probes=nrow(eligible),unique_genes=nrow(genes),landmarks_covered=length(covered),normalization='oligo RMA core; log2',annotation_package=as.character(packageVersion('hugene10sttranscriptcluster.db')),excluded_samples=character(0)),file.path(out,'rma_annotation_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('RMA_AND_ANNOTATION_PASS\n')
