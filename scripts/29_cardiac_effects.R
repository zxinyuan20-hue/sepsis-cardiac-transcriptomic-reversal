# Packages edgeR, limma, jsonlite. Inputs frozen cardiac-1.0 design and count matrix.
# Outputs donor-level expression, exploratory differential genes/pathways, library QC.
# Fixed seed from central runtime. No outcome-driven sample removal.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(as.integer(Sys.getenv('STUDY_SEED')))
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
cfg<-fromJSON(file.path(root,'config/cardiac_v1.json'));out<-file.path(root,'outputs/analysis/cardiac')
s<-read.delim(file.path(out,'analysis_design.tsv'),colClasses='character',check.names=FALSE)
s$primary_included<-tolower(s$primary_included)=='true';s$same_donor_sensitivity<-tolower(s$same_donor_sensitivity)=='true'
x<-as.matrix(read.delim(gzfile(file.path(root,'outputs/preparation/GSE217421/VEM_CTRL_counts.tsv.gz')),row.names=1,check.names=FALSE));storage.mode(x)<-'double'
stopifnot(identical(colnames(x),s$column_id),!anyDuplicated(rownames(x)),all(is.finite(x)),all(x>=0),all(colSums(x)>0))
mapping<-read.delim(file.path(out,'symbol_entrez_one_to_one.tsv'),colClasses='character')
sets<-read.delim(file.path(root,'outputs/analysis/mechanism/fixed_pathways.tsv'),colClasses='character');paths<-split(sets$entrez,sets$pathway_id)
thresholds<-cfg$expression_filter
xp<-x[,s$primary_included,drop=FALSE]
keep<-rowSums(cpm(xp)>=thresholds$minimum_cpm)>=thresholds$minimum_samples & rowSums(xp)>=thresholds$minimum_total_count
write.table(data.frame(gene_symbol=rownames(x),retained=keep),file.path(out,'gene_filter.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_matrix<-function(m,name){con<-gzfile(file.path(out,name),'wt');write.table(data.frame(gene_symbol=rownames(m),m,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE);close(con)}
metadata<-list()
for(scenario in c('primary','all_same_donors')){
 use<-if(scenario=='primary')s$primary_included else s$same_donor_sensitivity
 ss<-s[use,];y<-calcNormFactors(DGEList(counts=x[keep,use,drop=FALSE]));e<-cpm(y,log=TRUE,prior.count=0.5)
 write_matrix(e,paste0(scenario,'_sample_logCPM.tsv.gz'))
 write.table(data.frame(column_id=ss$column_id,donor=ss$Cell,condition=ss$State,culture=ss$Experiment,y$samples),file.path(out,paste0(scenario,'_library_metrics.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 # Average wells in each culture/condition, then cultures within each donor equally.
 donors<-sort(unique(ss$Cell)); profiles<-list();pd<-list()
 for(donor in donors)for(condition in c('CTRL','VEM')){
  sub<-ss[ss$Cell==donor & ss$State==condition,]
  cultures<-sort(unique(sub$block))
  vals<-sapply(cultures,function(b)rowMeans(e[,sub$column_id[sub$block==b],drop=FALSE]))
  id<-paste(donor,condition,sep=':');profiles[[id]]<-rowMeans(as.matrix(vals));pd[[id]]<-data.frame(profile_id=id,donor=donor,condition=condition,n_cultures=length(cultures),n_wells=nrow(sub))
 }
 z<-do.call(cbind,profiles);rownames(z)<-rownames(e);designmeta<-do.call(rbind,pd);rownames(designmeta)<-NULL
 designmeta$donor<-factor(designmeta$donor);designmeta$condition<-factor(designmeta$condition,levels=c('CTRL','VEM'))
 des<-model.matrix(~donor+condition,designmeta);stopifnot(qr(des)$rank==ncol(des),nrow(des)-ncol(des)==4)
 fit<-eBayes(lmFit(z,des),trend=TRUE,robust=FALSE)
 tab<-topTable(fit,coef='conditionVEM',number=Inf,sort.by='none',adjust.method='BH');tab$gene_symbol<-rownames(tab)
 tab$moderated_SE<-fit$stdev.unscaled[,'conditionVEM']*sqrt(fit$s2.post);tab$df_total<-fit$df.total;tab$df_residual<-fit$df.residual
 tab$entrez_id<-mapping$ENTREZID[match(tab$gene_symbol,mapping$SYMBOL)]
 tab$ci_low<-tab$logFC-qt(.975,tab$df_total)*tab$moderated_SE;tab$ci_high<-tab$logFC+qt(.975,tab$df_total)*tab$moderated_SE
 delta<-sapply(donors,function(d)z[,paste(d,'VEM',sep=':')]-z[,paste(d,'CTRL',sep=':')]);rownames(delta)<-rownames(z)
 stopifnot(max(abs(tab$logFC-rowMeans(delta)))<1e-8)
 tab$n_donors_up<-rowSums(delta>0);tab$n_donors_down<-rowSums(delta<0)
 write.table(tab,file.path(out,paste0(scenario,'_gene_results.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='')
 write_matrix(z,paste0(scenario,'_donor_condition_logCPM.tsv.gz'));write_matrix(delta,paste0(scenario,'_donor_logFC.tsv.gz'))
 write.table(designmeta,file.path(out,paste0(scenario,'_donor_profile_design.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(data.frame(profile_id=designmeta$profile_id,des,check.names=FALSE),file.path(out,paste0(scenario,'_model_matrix.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 if(scenario=='primary'){
  mapped<-which(!is.na(tab$entrez_id));ids<-tab$entrez_id[mapped];zm<-z[mapped,,drop=FALSE]
  index<-lapply(paths,function(g)which(ids%in%g));n<-lengths(index);coverage<-n/lengths(paths);eligible<-n>=10 & coverage>=.25
  cam<-data.frame(pathway_id=names(paths),NGenes=n,coverage=coverage,Direction=NA_character_,PValue=NA_real_,Correlation=NA_real_,FDR=NA_real_)
  if(any(eligible)){
   fitc<-camera(zm,index[eligible],design=des,contrast=ncol(des),inter.gene.cor=NA)
   ix<-match(rownames(fitc),cam$pathway_id);cam[ix,c('Direction','PValue','Correlation','FDR')]<-fitc[,c('Direction','PValue','Correlation','FDR')]
  }
  write.table(cam,file.path(out,'primary_pathways.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
 }
 metadata[[scenario]]<-list(wells=nrow(ss),donors=length(donors),profiles=ncol(z),retained_genes=nrow(z),mapped_genes=sum(!is.na(tab$entrez_id)),df_residual=4,fdr_up=sum(tab$adj.P.Val<.05 & tab$logFC>0),fdr_down=sum(tab$adj.P.Val<.05 & tab$logFC<0),minimum_library=min(colSums(x[,use,drop=FALSE])),normalization='TMM; logCPM prior.count=0.5; equal cultures then donors',gene_model='paired donor fixed effects limma-trend; exploratory')
}
write_json(metadata,file.path(out,'gene_analysis_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'effects_R_sessionInfo.txt'))
cat('DONOR_LEVEL_CARDIAC_ANALYSIS_PASS\n')
