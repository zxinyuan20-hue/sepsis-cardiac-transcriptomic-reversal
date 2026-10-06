# edgeR/jsonlite; fixed retained genes, seed; no plotting. Paths relative to root.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(20261005)
suppressPackageStartupMessages(library(edgeR))
src<-file.path(root,'outputs','analysis','raw_heart_v1');out<-file.path(root,'outputs','analysis','external_compound_v1')
x<-as.matrix(read.delim(gzfile(file.path(src,'entrez_counts_14.tsv.gz')),row.names=1,check.names=FALSE))
f<-read.delim(file.path(src,'expression_filter.tsv'));x<-x[as.character(f$entrez_id[f$retained]),,drop=FALSE]
s<-read.delim(file.path(src,'patient_metadata.tsv'));x<-x[,s$title,drop=FALSE]
common<-read.delim(file.path(out,'common_genes.tsv'),colClasses='character')$ENTREZID
result<-matrix(NA_real_,length(common),15,dimnames=list(common,c('full',s$title)))
for(i in 0:14){
 keep<-if(i==0) seq_len(14) else setdiff(seq_len(14),i)
 z<-s[keep,];z$sex<-factor(z$sex);z$group<-factor(z$group,levels=c('control','sepsis'))
 # Removing the only female case leaves female controls: design still identifiable.
 design<-model.matrix(~sex+group,z);stopifnot(qr(design)$rank==ncol(design))
 y<-calcNormFactors(DGEList(x[,keep,drop=FALSE]));y<-estimateDisp(y,design,robust=FALSE)
 fit<-glmQLFit(y,design,robust=FALSE);t<-topTags(glmQLFTest(fit,coef='groupsepsis'),n=Inf,sort.by='none')$table
 result[,i+1]<-t[common,'logFC']
}
write.table(data.frame(entrez_id=rownames(result),result,check.names=FALSE),file.path(out,'GSE237861_patient_loo_logFC.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
