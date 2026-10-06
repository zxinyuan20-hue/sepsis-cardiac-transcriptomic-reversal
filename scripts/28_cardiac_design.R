# Packages AnnotationDbi/org.Hs.eg.db. Only feature IDs; no expression effects.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library))
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db)})
x<-read.delim(gzfile(file.path(root,'outputs/preparation/GSE217421/VEM_CTRL_counts.tsv.gz')),check.names=FALSE)
ids<-intersect(x$gene_symbol,keys(org.Hs.eg.db,keytype='SYMBOL'))
m<-unique(AnnotationDbi::select(org.Hs.eg.db,keys=ids,keytype='SYMBOL',columns='ENTREZID'));m<-m[!is.na(m$ENTREZID),]
out<-file.path(root,'outputs/analysis/cardiac')
write.table(m,file.path(out,'symbol_entrez_all.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
m<-m[!duplicated(m$SYMBOL)&!duplicated(m$SYMBOL,fromLast=TRUE)&!duplicated(m$ENTREZID)&!duplicated(m$ENTREZID,fromLast=TRUE),]
write.table(m,file.path(out,'symbol_entrez_one_to_one.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'mapping_R_sessionInfo.txt'))
cat('EXACT_SYMBOL_MAPPING_PASS',nrow(m),'\n')
