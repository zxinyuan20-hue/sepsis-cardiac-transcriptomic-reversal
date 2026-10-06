# Required packages edgeR/limma/jsonlite. Inputs frozen expansion config/design/counts.
# Outputs per-drug donor logFC, gene audit tables and sample rank endpoints under outputs.
# Same normalization and paired-donor model as cardiac-1.0; seeded, no outcome exclusions.
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
.libPaths(c(file.path(root,'environment','R-library'),.Library));set.seed(as.integer(Sys.getenv('STUDY_SEED')))
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
cfg<-fromJSON(file.path(root,'config/expansion_v1.json'));out<-file.path(root,'outputs/analysis/expansion')
s<-read.delim(file.path(out,'expanded_analysis_design.tsv'),colClasses='character',check.names=FALSE)
x<-as.matrix(read.delim(gzfile(file.path(out,'expanded_counts.tsv.gz')),row.names=1,check.names=FALSE));storage.mode(x)<-'double'
filters<-read.delim(file.path(out,'expanded_gene_filters.tsv'),row.names=1,check.names=FALSE,colClasses='character')
stopifnot(identical(rownames(x),rownames(filters)),all(colSums(x)>0))
mapping<-read.delim(file.path(root,'outputs/analysis/cardiac/symbol_entrez_one_to_one.tsv'),colClasses='character')
write_matrix<-function(m,name){con<-gzfile(file.path(out,name),'wt');write.table(data.frame(gene_symbol=rownames(m),m,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE);close(con)}
sample_scores<-list();coverage<-list();qc<-list()
for(drug in sort(unique(s$drug_state))){
 cat('DRUG_START',drug,'\n');ss<-s[s$drug_state==drug,];counts<-x[,ss$column_id,drop=FALSE]
 keep<-tolower(filters[[drug]])=='true';th<-cfg$gene_filter
 recomputed<-rowSums(cpm(counts)>=th$minimum_cpm)>=th$minimum_samples & rowSums(counts)>=th$minimum_total_count
 stopifnot(identical(keep,as.logical(recomputed)))
 y<-calcNormFactors(DGEList(counts=counts[keep,,drop=FALSE]));e<-cpm(y,log=TRUE,prior.count=.5)
 ent<-mapping$ENTREZID[match(rownames(e),mapping$SYMBOL)];ranks<-(apply(e,2,rank,ties.method='average')-.5)/nrow(e)
 for(size in c(cfg$primary_query_size,cfg$sensitivity_query_sizes)){
  sets<-list()
  for(direction in c('up','down')){
   query<-read.delim(file.path(root,paste0('outputs/analysis/human/frozen_signatures/landmarks_',size,'_',direction,'.tsv')),colClasses='character')
   mask<-ent%in%query$entrez_id;n<-sum(mask);ratio<-n/nrow(query)
   stopifnot(n>=cfg$minimum_query_genes,ratio>=cfg$minimum_query_coverage)
   sets[[direction]]<-mask;coverage[[paste(drug,size,direction)]]<-data.frame(drug_state=drug,size=size,direction=direction,source_genes=nrow(query),covered_genes=n,coverage=ratio)
  }
  z<-data.frame(drug_state=drug,size=size,column_id=ss$column_id,donor=ss$Cell,block=ss$block,condition=ss$condition,up_rank=colMeans(ranks[sets$up,,drop=FALSE]),down_rank=colMeans(ranks[sets$down,,drop=FALSE]))
  z$disease_rank<-z$up_rank-z$down_rank;sample_scores[[paste(drug,size)]]<-z
 }
 profiles<-list();pd<-list();donors<-sort(unique(ss$Cell))
 for(donor in donors)for(condition in c('CTRL','DRUG')){
  sub<-ss[ss$Cell==donor & ss$condition==condition,];cultures<-sort(unique(sub$block));id<-paste(donor,condition,sep=':')
  profiles[[id]]<-rowMeans(sapply(cultures,function(b)rowMeans(e[,sub$column_id[sub$block==b],drop=FALSE])))
  pd[[id]]<-data.frame(profile_id=id,donor=donor,condition=condition)
 }
 z<-do.call(cbind,profiles);rownames(z)<-rownames(e);dm<-do.call(rbind,pd);dm$donor<-factor(dm$donor);dm$condition<-factor(dm$condition,levels=c('CTRL','DRUG'))
 des<-model.matrix(~donor+condition,dm);stopifnot(qr(des)$rank==ncol(des),nrow(des)-ncol(des)==length(donors)-1)
 fit<-eBayes(lmFit(z,des),trend=TRUE,robust=FALSE)
 tab<-topTable(fit,coef='conditionDRUG',number=Inf,sort.by='none',adjust.method='BH');tab$gene_symbol<-rownames(tab);tab$entrez_id<-ent
 tab$moderated_SE<-fit$stdev.unscaled[,'conditionDRUG']*sqrt(fit$s2.post);tab$df_total<-fit$df.total;tab$df_residual<-fit$df.residual
 delta<-sapply(donors,function(d)z[,paste(d,'DRUG',sep=':')]-z[,paste(d,'CTRL',sep=':')]);rownames(delta)<-rownames(z)
 stopifnot(max(abs(tab$logFC-rowMeans(delta)))<1e-8)
 write.table(tab,file.path(out,paste0(drug,'_gene_results.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='')
 write_matrix(delta,paste0(drug,'_donor_logFC.tsv.gz'));write_matrix(z,paste0(drug,'_donor_condition_logCPM.tsv.gz'))
 write.table(data.frame(column_id=ss$column_id,y$samples),file.path(out,paste0(drug,'_library_metrics.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 qc[[drug]]<-list(drug=drug,donors=length(donors),cultures=length(unique(ss$block)),wells=nrow(ss),genes=nrow(e),df_residual=length(donors)-1,per_drug_gene_fdr_count=sum(tab$adj.P.Val<.05),warning='gene BH is within-drug only, not global across drug-gene pairs')
 cat('DRUG_DONE',drug,'donors',length(donors),'genes',nrow(e),'\n')
}
write.table(do.call(rbind,sample_scores),file.path(out,'expanded_sample_rank_scores.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,coverage),file.path(out,'expanded_query_coverage.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write_json(qc,file.path(out,'expanded_gene_qc.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'expanded_R_sessionInfo.txt'))
cat('EXPANDED_DONOR_ANALYSIS_PASS\n')
