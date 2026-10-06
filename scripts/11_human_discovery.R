# Packages: limma/jsonlite. Inputs: normalized human data/config. Outputs: outputs/analysis/human.
# Reproducibility: fixed seed; contrasts and BH per contrast; no private rat disease effects.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(as.integer(Sys.getenv('STUDY_SEED')))
suppressPackageStartupMessages(library(limma))
cfg<-jsonlite::fromJSON(file.path(root,'config','study_config.json'))
stopifnot(isTRUE(cfg$analysis$formal_analysis_enabled))
input<-file.path(root,'outputs','preprocessing','human');out<-file.path(root,'outputs','analysis','human')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
readmat<-function(p)as.matrix(read.delim(gzfile(p),row.names=1,check.names=FALSE))
x<-readmat(file.path(input,'human_rma_genes.tsv.gz'))
samples<-read.delim(file.path(input,'cel_sample_manifest.tsv'))
samples<-samples[match(colnames(x),samples$gsm),]
group<-factor(samples$group,levels=c('nonfailing','sepsis','IHD','DCM'))
design<-model.matrix(~0+group);colnames(design)<-levels(group)
contrast<-makeContrasts(sepsis_vs_nonfailing=sepsis-nonfailing,sepsis_vs_IHD=sepsis-IHD,sepsis_vs_DCM=sepsis-DCM,levels=design)
stopifnot(qr(design)$rank==4,ncol(x)==51,all(is.finite(x)))
write.table(data.frame(sample_id=samples$gsm,design),file.path(out,'design_matrix.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(contrast,file.path(out,'contrast_matrix.tsv'),sep='\t',quote=FALSE,col.names=NA)
run_fit<-function(mat,design,label){
 fit<-eBayes(contrasts.fit(lmFit(mat,design),contrast),trend=FALSE,robust=FALSE)
 summary<-list()
 for(k in seq_len(ncol(contrast))){
  result<-topTable(fit,coef=k,number=Inf,sort.by='none',adjust.method='BH')
  result$entrez_id<-rownames(result)
  result$moderated_SE<-fit$stdev.unscaled[,k]*sqrt(fit$s2.post)
  result$df_total<-fit$df.total
  result$ci_low<-result$logFC-qt(.975,fit$df.total)*result$moderated_SE
  result$ci_high<-result$logFC+qt(.975,fit$df.total)*result$moderated_SE
  result<-result[,c('entrez_id','logFC','moderated_SE','ci_low','ci_high','df_total','AveExpr','t','P.Value','adj.P.Val','B')]
  write.table(result,file.path(out,paste0(label,'_',colnames(contrast)[k],'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  summary[[colnames(contrast)[k]]]<-list(tested=nrow(result),fdr_up=sum(result$adj.P.Val<cfg$analysis$fdr&result$logFC>0),fdr_down=sum(result$adj.P.Val<cfg$analysis$fdr&result$logFC<0))
 }
 summary
}
stats<-list(primary=run_fit(x,design,'primary'))
flagged<-cfg$analysis$qc_flagged_samples
keep<-!(samples$gsm %in% flagged)
if(length(flagged))stats$qc_sensitivity<-run_fit(x[,keep,drop=FALSE],design[keep,,drop=FALSE],'without_qc_flag')
submitted<-readmat(file.path(root,'outputs','preparation','GSE79962_submitted_expression.tsv.gz'))
selected<-read.delim(file.path(input,'selected_probe_per_gene.tsv'),colClasses='character')
sub<-submitted[selected$PROBEID,samples$gsm,drop=FALSE];rownames(sub)<-selected$ENTREZID
stopifnot(identical(rownames(x),rownames(sub)))
stats$submitted_sensitivity<-run_fit(sub,design,'submitted')
jsonlite::write_json(stats,file.path(out,'discovery_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('HUMAN_DISCOVERY_PASS\n')
