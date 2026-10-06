# Packages AnnotationDbi/org.Mm.eg.db; map mouse Ensembl to Entrez without phenotype effects.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library))
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Mm.eg.db)})
out<-file.path(root,'outputs','analysis','mechanism');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-read.delim(gzfile(file.path(root,'outputs','preparation','GSE185754_counts.tsv.gz')),row.names=1)
ids<-rownames(x);valid<-intersect(ids,keys(org.Mm.eg.db,keytype='ENSEMBL'))
m<-unique(AnnotationDbi::select(org.Mm.eg.db,keys=valid,keytype='ENSEMBL',columns='ENTREZID'));m<-m[!is.na(m$ENTREZID),]
write.table(m,file.path(out,'mouse_ensembl_entrez_all.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
m<-m[!duplicated(m$ENSEMBL)&!duplicated(m$ENSEMBL,fromLast=TRUE)&!duplicated(m$ENTREZID)&!duplicated(m$ENTREZID,fromLast=TRUE),]
write.table(m,file.path(out,'mouse_ensembl_entrez_one_to_one.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
