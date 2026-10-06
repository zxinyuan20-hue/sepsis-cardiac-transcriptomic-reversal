# Packages: limma,jsonlite. Input human matrix and fixed queries; output bootstrap masks.
# Resample samples within four groups, full-gene BH, then landmark selection. No rat effects.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library))
suppressPackageStartupMessages(library(limma))
cfg<-jsonlite::fromJSON(file.path(root,'config','screening_v1.json'))
set.seed(cfg$seed)
out<-file.path(root,'outputs','analysis','screening')
x<-as.matrix(read.delim(gzfile(file.path(root,'outputs','preprocessing','human','human_rma_genes.tsv.gz')),row.names=1,check.names=FALSE))
samples<-read.delim(file.path(root,'outputs','preprocessing','human','cel_sample_manifest.tsv'))
samples<-samples[match(colnames(x),samples$gsm),]
group<-factor(samples$group,levels=c('nonfailing','sepsis','IHD','DCM'))
design<-model.matrix(~0+group);colnames(design)<-levels(group)
contrast<-makeContrasts(sepsis-nonfailing,levels=design)
genes<-read.delim(file.path(out,'gene_order.tsv'),colClasses='character')$gene_id
landmark<-which(rownames(x)%in%genes)
queries<-matrix(0L,nrow=cfg$bootstrap_iterations,ncol=length(genes),dimnames=list(NULL,genes))
sampled<-matrix(0L,nrow=cfg$bootstrap_iterations,ncol=ncol(x))
quality<-data.frame(replicate=seq_len(cfg$bootstrap_iterations),up=0L,down=0L,valid=FALSE)
for(b in seq_len(cfg$bootstrap_iterations)){
 idx<-unlist(lapply(levels(group),function(g)sample(which(group==g),sum(group==g),replace=TRUE)))
 sampled[b,]<-idx
 fit<-eBayes(contrasts.fit(lmFit(x[,idx,drop=FALSE],design[idx,,drop=FALSE]),contrast),trend=FALSE,robust=FALSE)
 p<-p.adjust(fit$p.value[,1],method='BH')
 for(sign in c(1L,-1L)){
  valid<-landmark[p[landmark]<.05 & fit$coefficients[landmark,1]*sign>0]
  chosen<-head(valid[order(-abs(fit$t[valid,1]),rownames(x)[valid])],100)
  queries[b,match(rownames(x)[chosen],genes)]<-sign
 }
 quality$up[b]<-sum(queries[b,]==1);quality$down[b]<-sum(queries[b,]==-1)
 quality$valid[b]<-min(quality$up[b],quality$down[b])>=15
 if(b%%50==0){cat('BOOTSTRAP',b,'/',cfg$bootstrap_iterations,'\n');flush.console()}
}
con<-gzfile(file.path(out,'bootstrap_queries.tsv.gz'),'wt');write.table(queries,con,sep='\t',quote=FALSE,row.names=FALSE);close(con)
write.table(quality,file.path(out,'bootstrap_query_quality.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(sampled,file.path(out,'bootstrap_sample_indices.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(samples$gsm,file.path(out,'bootstrap_sample_index_reference.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'bootstrap_R_sessionInfo.txt'))
cat('BOOTSTRAP_QUERIES_PASS\n')
