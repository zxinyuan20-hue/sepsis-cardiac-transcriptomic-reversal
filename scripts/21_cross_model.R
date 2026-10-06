# Packages edgeR/limma/jsonlite. Fixed cohort inputs/config, outputs/analysis/mechanism.
# Original counts preserved; exploratory estimated-count inference, not drug treatment effects.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(as.integer(Sys.getenv('STUDY_SEED')))
suppressPackageStartupMessages({library(edgeR);library(limma)})
cfg<-jsonlite::fromJSON(file.path(root,'config','mechanism_v1.json'));out<-file.path(root,'outputs','analysis','mechanism')
cohorts<-jsonlite::fromJSON(file.path(out,'cohort_inputs.json'))
sets<-read.delim(file.path(out,'fixed_pathways.tsv'),colClasses='character')
paths<-split(sets$entrez,sets$pathway_id)
allcamera<-list();estimated_camera<-list();qc<-list()
write_matrix<-function(x,file){con<-gzfile(file,'wt');write.table(data.frame(feature_id=rownames(x),x,check.names=FALSE),con,sep='\t',row.names=FALSE,quote=FALSE);close(con)}
run_camera<-function(expr,design,contrast,ids,dataset,primary=TRUE){
 indices<-lapply(paths,function(g)which(ids%in%g));n<-lengths(indices);coverage<-n/lengths(paths)
 eligible<-n>=cfg$minimum_pathway_genes & coverage>=cfg$minimum_pathway_coverage
 result<-data.frame(pathway_id=names(paths),NGenes=n,coverage=coverage,Direction=NA_character_,PValue=NA_real_,cohort=dataset,primary=primary)
 if(any(eligible)){
  fit<-camera(expr,index=indices[eligible],design=design,contrast=contrast,inter.gene.cor=cfg$camera_inter_gene_correlation)
  idx<-match(rownames(fit),result$pathway_id);result$Direction[idx]<-fit$Direction;result$PValue[idx]<-fit$PValue
 }
 allcamera[[dataset]]<<-result
 # Added post-result sensitivity: estimate each set's residual correlation.
 alternate<-result;alternate$Direction<-NA_character_;alternate$PValue<-NA_real_;alternate$Correlation<-NA_real_
 if(any(eligible)){
  fit2<-camera(expr,index=indices[eligible],design=design,contrast=contrast,inter.gene.cor=NA)
  idx<-match(rownames(fit2),alternate$pathway_id);alternate$Direction[idx]<-fit2$Direction;alternate$PValue[idx]<-fit2$PValue
  if('Correlation'%in%names(fit2))alternate$Correlation[idx]<-fit2$Correlation
 }
 estimated_camera[[dataset]]<<-alternate
}
# Human pathway context using full four-group primary design, but annotated background.
human<-as.matrix(read.delim(gzfile(file.path(root,'outputs','preprocessing','human','human_rma_genes.tsv.gz')),row.names=1,check.names=FALSE))
hs<-read.delim(file.path(root,'outputs','preprocessing','human','cel_sample_manifest.tsv'));hs<-hs[match(colnames(human),hs$gsm),]
g<-factor(hs$group,levels=c('nonfailing','sepsis','IHD','DCM'));design<-model.matrix(~0+g);colnames(design)<-levels(g)
contrast<-makeContrasts(sepsis-nonfailing,levels=design)[,1]
run_camera(human,design,contrast,rownames(human),'human')
for(i in seq_len(nrow(cohorts))){
 row<-cohorts[i,];dataset<-row$dataset
 x<-as.matrix(read.delim(gzfile(file.path(root,row$counts)),row.names=1,check.names=FALSE));storage.mode(x)<-'double'
 s<-read.delim(file.path(root,row$samples));s<-s[match(colnames(x),s$sample_id),]
 mapping<-read.delim(file.path(root,row$mapping),colClasses='character');stopifnot(!anyDuplicated(mapping$feature_id),!anyDuplicated(mapping$human_entrez))
 for(sensitivity in if(dataset=='LOCAL_RAT_LPS_2025')c(FALSE,TRUE) else FALSE){
  label<-if(sensitivity)paste0(dataset,'_without_C1_C2') else dataset
  use<-if(sensitivity)!(s$sample_id%in%cfg$rat_sensitivity_excluded_samples) else rep(TRUE,nrow(s))
  g<-factor(s$group[use],levels=c('control','lps'));des<-model.matrix(~g)
  y<-DGEList(counts=x[,use,drop=FALSE],group=g);keep<-filterByExpr(y,design=des);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
  pdf(file.path(out,paste0(label,'_voom.pdf')),width=6,height=5);v<-voom(y,des,plot=TRUE);dev.off()
  fit<-eBayes(lmFit(v,des),trend=FALSE,robust=FALSE)
  tab<-topTable(fit,coef=2,number=Inf,sort.by='none',adjust.method='BH');tab$feature_id<-rownames(tab)
  tab$moderated_SE<-fit$stdev.unscaled[,2]*sqrt(fit$s2.post);tab$df_total<-fit$df.total
  tab$human_entrez<-mapping$human_entrez[match(tab$feature_id,mapping$feature_id)]
  write.table(tab,file.path(out,paste0(label,'_all_gene_results.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='')
  mapped<-which(!is.na(tab$human_entrez));vv<-v[mapped,];ids<-tab$human_entrez[mapped]
  write_matrix(v$E,file.path(out,paste0(label,'_voom_expression.tsv.gz')))
  write.table(data.frame(sample_id=colnames(y),y$samples),file.path(out,paste0(label,'_library_metrics.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  run_camera(vv,des,2,ids,label,!sensitivity)
  qc[[label]]<-list(input_genes=nrow(x),retained_genes=sum(keep),mapped_human_genes=length(ids),samples=ncol(y),fdr_up=sum(tab$adj.P.Val<.05&tab$logFC>0),fdr_down=sum(tab$adj.P.Val<.05&tab$logFC<0),fractional_input_elements=sum(x!=floor(x)),role=if(sensitivity)'QC sensitivity only' else 'cross-model disease association; no drug exposure',estimated_count_caution=dataset!='GSE185754')
 }
}
write.table(do.call(rbind,allcamera),file.path(out,'camera_all_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
write.table(do.call(rbind,estimated_camera),file.path(out,'camera_estimated_correlation_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
jsonlite::write_json(qc,file.path(out,'cross_model_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('CROSS_MODEL_ANALYSIS_PASS\n')
