# Packages limma/jsonlite; ROOT-relative config, input matrix and header metadata.
# Output sensitivity model coefficients/designs under configured outputs directory.
# Fixed seed; original full-data preprocessing retained; no compound scoring.
args<-commandArgs(FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library))
suppressPackageStartupMessages(library(limma))
cfg<-jsonlite::fromJSON(file.path(root,'config','discovery_scan_sensitivity_v1.json'));set.seed(cfg$seed)
out<-file.path(root,cfg$output);dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-as.matrix(read.delim(gzfile(file.path(root,'outputs/preprocessing/human/human_rma_genes.tsv.gz')),row.names=1,check.names=FALSE))
meta<-read.delim(file.path(root,'outputs/analysis/discovery_exchangeability_v1/cel_scan_metadata.tsv'))
meta<-meta[match(colnames(x),meta$gsm),];stopifnot(identical(meta$gsm,colnames(x)),nrow(x)==18866)
shared<-intersect(meta$scan_day[meta$group=='sepsis'],meta$scan_day[meta$group=='nonfailing'])
summaries<-list()
for(label in c('primary_reproduction','scan_adjusted','shared_days_adjusted')){
 keep<-if(label=='shared_days_adjusted')meta$scan_day%in%shared else rep(TRUE,nrow(meta))
 m<-meta[keep,];g<-factor(m$group,levels=c('nonfailing','sepsis','IHD','DCM'));day<-factor(m$scan_day)
 design<-if(label=='primary_reproduction')model.matrix(~0+g) else model.matrix(~0+g+day)
 colnames(design)[1:4]<-levels(g);stopifnot(qr(design)$rank==ncol(design))
 contrast<-rep(0,ncol(design));contrast[1:2]<-c(-1,1)
 fit<-eBayes(contrasts.fit(lmFit(x[,keep,drop=FALSE],design),contrast),trend=FALSE,robust=FALSE)
 tab<-topTable(fit,coef=1,number=Inf,sort.by='none',adjust.method='BH');tab$entrez_id<-rownames(tab)
 write.table(tab,file.path(out,paste0(label,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(gsm=m$gsm,design,check.names=FALSE),file.path(out,paste0(label,'_design.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 summaries[[label]]<-list(n=nrow(design),design_rank=qr(design)$rank,residual_df=nrow(design)-ncol(design),groups=as.list(table(g)),genes=nrow(tab),up=sum(tab$adj.P.Val<.05&tab$logFC>0),down=sum(tab$adj.P.Val<.05&tab$logFC<0))
}
jsonlite::write_json(summaries,file.path(out,'model_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
