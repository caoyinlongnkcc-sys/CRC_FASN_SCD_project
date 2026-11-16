#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."; raw <- file.path(root,"data_raw"); out <- file.path(root,"results_tables"); log <- file.path(root,"logs"); audit <- file.path("..","feasibility_audit")
anchors <- c("FASN","SREBF1","ACACA","SCD","ELOVL6","FDFT1","TECR")
mf <- list.files(file.path(raw,"project_score_extract"),pattern="fold_change_values",full.names=TRUE)[1]
if(is.na(mf)||!file.exists(mf))stop("Missing Project Score fold-change file")
z <- read.delim(mf,header=FALSE,check.names=FALSE,quote="",fill=TRUE,stringsAsFactors=FALSE)
rn <- as.character(z[1,]); ri <- as.character(z[2,]); src <- as.character(z[3,]); qc <- as.character(z[4,])
idx <- which(seq_along(rn)>3 & !is.na(rn)&nzchar(rn)&!is.na(ri)&nzchar(ri)); sidx <- idx[src[idx]=="Sanger" & qc[idx]=="TRUE"]
sids <- ri[sidx]; unique_sids <- unique(sids)
sym <- as.character(z[-c(1,2,3,4,5),2]); gi <- which(sym%in%anchors)
if(!length(gi))stop("No frozen anchor genes detected in Project Score matrix")
vals <- matrix(NA_real_,nrow=length(unique_sids),ncol=length(anchors),dimnames=list(unique_sids,anchors))
for(j in seq_along(anchors)){rows<-which(sym==anchors[j]); if(!length(rows))next; for(i in seq_along(unique_sids)){cc<-sidx[sids==unique_sids[i]]; vals[i,j]<-median(as.numeric(z[rows[1]+5,cc]),na.rm=TRUE)}}
ml <- read.csv(file.path(raw,"Sanger_CellModelPassports_model_list_latest.csv"),check.names=FALSE,stringsAsFactors=FALSE)
dep_ids <- read.delim(file.path(out,"depmap_crc_model_ids.tsv"),sep="\t",stringsAsFactors=FALSE)$ModelID
dm <- read.csv(file.path(audit,"DepMap_24Q2_Model.csv"),check.names=FALSE,stringsAsFactors=FALSE)
maprow <- match(unique_sids,ml$model_id); crc_ids <- ml$model_id[ml$cancer_type=="Colorectal Carcinoma"]
is_crc <- unique_sids%in%crc_ids; is_model_list <- !is.na(maprow)
sanger_to_dep <- unique(dm$ModelID[dm$SangerModelID%in%unique_sids & dm$ModelID%in%dep_ids]); overlap <- intersect(unique_sids,dm$SangerModelID[dm$ModelID%in%dep_ids])
qc_tab <- data.frame(fitness_model_id=unique_sids,model_list_match=is_model_list,cancer_type=ifelse(is.na(maprow),NA,ml$cancer_type[maprow]),cancer_type_detail=ifelse(is.na(maprow),NA,ml$cancer_type_detail[maprow]),tissue=ifelse(is.na(maprow),NA,ml$tissue[maprow]),depmap_crc_overlap=unique_sids%in%overlap,stringsAsFactors=FALSE)
write.table(qc_tab,file.path(out,"project_score_sanger_model_qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
stats<-lapply(seq_along(anchors),function(j){a<-vals[is_crc,j];b<-vals[!is_crc,j];a<-a[is.finite(a)];b<-b[is.finite(b)]; p<-if(length(a)>=10&&length(b)>=10)suppressWarnings(wilcox.test(a,b,exact=FALSE)$p.value)else NA_real_; data.frame(gene=anchors[j],n_crc=length(a),n_noncrc=length(b),median_crc=if(length(a))median(a)else NA_real_,median_noncrc=if(length(b))median(b)else NA_real_,median_difference=if(length(a)&&length(b))median(a)-median(b)else NA_real_,prop_crc_below_neg0_5=if(length(a))mean(a <= -0.5)else NA_real_,p=p)})
sr<-do.call(rbind,stats); sr$BH<-p.adjust(sr$p,"BH")
dep<-read.delim(file.path(out,"depmap_crc_lipid_dependency_results.tsv"),sep="\t",check.names=FALSE,stringsAsFactors=FALSE); dep<-dep[match(anchors,dep$gene),c("gene","median_difference","fdr_bh")]; names(dep)<-c("gene","depmap_median_difference","depmap_fdr")
sr<-merge(sr,dep,by="gene",sort=FALSE); sr$direction_concordant<-sign(sr$median_difference)==sign(sr$depmap_median_difference); write.table(sr,file.path(out,"project_score_sanger_anchor_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
sumtab<-data.frame(raw_sanger_columns=sum(src[idx]=="Sanger"),qc_pass_sanger_columns=sum(src[idx]=="Sanger"&qc[idx]=="TRUE"),unique_sanger_model_ids=length(unique_sids),model_list_matches=sum(is_model_list),sanger_crc_models=sum(is_crc),sanger_noncrc_models=sum(!is_crc),depmap_crc_model_overlap=length(overlap),anchors_detected=sum(anchors%in%sym),stringsAsFactors=FALSE); write.table(sumtab,file.path(out,"project_score_sanger_audit_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Project Score Sanger audit completed.",paste("unique Sanger models",length(unique_sids)),paste("model-list matches",sum(is_model_list)),paste("CRC models",sum(is_crc)),paste("non-CRC models",sum(!is_crc)),paste("DepMap CRC model overlap",length(overlap)),paste("anchors detected",sum(anchors%in%sym))),file.path(log,"project_score_sanger_audit.log"))
