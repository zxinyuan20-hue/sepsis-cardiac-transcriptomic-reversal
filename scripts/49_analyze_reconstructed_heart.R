# Packages edgeR,limma,jsonlite; project library. Fixed sex+group model, seed.
# Input raw_heart_v1/entrez_counts_14.tsv.gz; outputs all-gene/pathway results.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(20261002)
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
out<-file.path(root,'outputs','analysis','raw_heart_v1')
x<-as.matrix(read.delim(gzfile(file.path(out,'entrez_counts_14.tsv.gz')),row.names=1,check.names=FALSE))
s<-read.delim(file.path(out,'patient_metadata.tsv'));x<-x[,s$title,drop=FALSE]
stopifnot(ncol(x)==14,!anyDuplicated(s$title),sum(s$group=='sepsis')==7,all(x>=0),all(x==round(x)))
s$sex<-factor(s$sex);s$group<-factor(s$group,levels=c('control','sepsis'));design<-model.matrix(~sex+group,s)
stopifnot(qr(design)$rank==ncol(design))
# Label-independent filter, then TMM. Not length normalized for QuantSeq.
keep<-rowSums(cpm(x)>=1)>=7
y<-DGEList(x[keep,,drop=FALSE]);y<-calcNormFactors(y,method='TMM')
y<-estimateDisp(y,design,robust=FALSE);fit<-glmQLFit(y,design,robust=FALSE)
res<-topTags(glmQLFTest(fit,coef='groupsepsis'),n=Inf,sort.by='none')$table;res$entrez_id<-rownames(res)
write.table(res,file.path(out,'GSE237861_all_gene_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
logcpm<-cpm(y,log=TRUE,prior.count=2)
write.table(data.frame(entrez_id=rownames(logcpm),logcpm,check.names=FALSE),gzfile(file.path(out,'logCPM_14.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(title=rownames(y$samples),y$samples),file.path(out,'normalization_factors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(entrez_id=rownames(x),retained=keep),file.path(out,'expression_filter.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
v<-voom(y,design,plot=FALSE)
ps<-read.delim(file.path(root,'outputs','analysis','mechanism','fixed_pathways.tsv'),colClasses='character');sets<-split(ps$entrez,ps$pathway_id)
cfg<-jsonlite::fromJSON(file.path(root,'config','mechanism_v1.json'))
idx<-lapply(sets,function(a)which(rownames(y)%in%a));n<-lengths(idx);coverage<-n/lengths(sets)
eligible<-n>=cfg$minimum_pathway_genes & coverage>=cfg$minimum_pathway_coverage
rows<-list()
for(mode in c('estimated_residual_correlation','fixed_001_sensitivity')){
 z<-data.frame(pathway_id=names(sets),NGenes=n,coverage=coverage,Direction=NA_character_,PValue=NA_real_,Correlation=NA_real_,mode=mode)
 if(any(eligible)){
  a<-camera(v,idx[eligible],design,contrast=which(colnames(design)=='groupsepsis'),inter.gene.cor=if(mode=='estimated_residual_correlation') NA else .01)
  ii<-match(rownames(a),z$pathway_id);z$Direction[ii]<-a$Direction;z$PValue[ii]<-a$PValue
  if('Correlation'%in%names(a))z$Correlation[ii]<-a$Correlation
 }
 rows[[mode]]<-z
}
write.table(do.call(rbind,rows),file.path(out,'GSE237861_fixed_pathway_tests.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('RECONSTRUCTED_HEART_ANALYSIS_COMPLETED\n')
