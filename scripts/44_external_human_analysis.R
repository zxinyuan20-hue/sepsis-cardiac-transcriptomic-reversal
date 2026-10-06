# Packages limma/jsonlite from project library. Inputs submitted log2 RMA matrix/design.
# Outputs all-gene disease contrasts and fixed-pathway tests under outputs/analysis/external_human.
# Fixed patient design; no data-driven exclusion or drug ranking.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(20261002)
suppressPackageStartupMessages(library(limma))
out<-file.path(root,'outputs','analysis','external_human')
x<-as.matrix(read.delim(gzfile(file.path(out,'GSE141864_gene_expression.tsv.gz')),row.names=1,check.names=FALSE))
s<-read.delim(file.path(out,'GSE141864_heart_design.tsv'));s<-s[s$primary %in% c(TRUE,'True','TRUE'),];x<-x[,s$gsm,drop=FALSE]
stopifnot(ncol(x)==7,!anyDuplicated(s$patient_id),sum(s$group=='sepsis')==5)
g<-factor(s$group,levels=c('control','sepsis'));design<-model.matrix(~g)
fit<-eBayes(lmFit(x,design),trend=FALSE,robust=FALSE)
t<-topTable(fit,coef=2,number=Inf,sort.by='none',adjust.method='BH');t$entrez_id<-rownames(t)
t$moderated_SE<-fit$stdev.unscaled[,2]*sqrt(fit$s2.post);t$df_total<-fit$df.total
write.table(t,file.path(out,'GSE141864_all_gene_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
ps<-read.delim(file.path(root,'outputs','analysis','mechanism','fixed_pathways.tsv'),colClasses='character');sets<-split(ps$entrez,ps$pathway_id)
oldcfg<-jsonlite::fromJSON(file.path(root,'config','mechanism_v1.json'))
idx<-lapply(sets,function(v)which(rownames(x)%in%v));n<-lengths(idx);coverage<-n/lengths(sets)
eligible<-n>=oldcfg$minimum_pathway_genes & coverage>=oldcfg$minimum_pathway_coverage
rows<-list()
for(mode in c('estimated_residual_correlation','fixed_001_sensitivity')){
 z<-data.frame(pathway_id=names(sets),NGenes=n,coverage=coverage,Direction=NA_character_,PValue=NA_real_,Correlation=NA_real_,mode=mode)
 if(any(eligible)){
  a<-camera(x,idx[eligible],design,contrast=2,inter.gene.cor=if(mode=='estimated_residual_correlation') NA else .01)
  ii<-match(rownames(a),z$pathway_id);z$Direction[ii]<-a$Direction;z$PValue[ii]<-a$PValue
  if('Correlation'%in%names(a))z$Correlation[ii]<-a$Correlation
 }
 rows[[mode]]<-z
}
write.table(do.call(rbind,rows),file.path(out,'GSE141864_fixed_pathway_tests.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('EXTERNAL_HUMAN_GENE_ANALYSIS_PASS\n')
