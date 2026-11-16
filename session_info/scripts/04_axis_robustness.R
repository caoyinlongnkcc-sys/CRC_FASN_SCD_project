#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."; audit <- file.path("..","feasibility_audit"); raw <- file.path(root,"data_raw"); meta <- file.path(root,"data_metadata"); out <- file.path(root,"results_tables"); log <- file.path(root,"logs")
anchors <- c("FASN","SREBF1","ACACA","SCD","ELOVL6","FDFT1","TECR")
depfile <- file.path(audit,"depmap24q2","CRISPRGeneEffect.csv")
dh <- names(read.csv(depfile,nrows=0,check.names=FALSE)); dc <- sub(" \\(.*$","",dh); di <- which(dc%in%anchors)
dclass <- rep("NULL",length(dh)); dclass[1]<-"character"; dclass[di]<-"numeric"; dep <- read.csv(depfile,check.names=FALSE,colClasses=dclass); ids<-as.character(dep[[1]])
md <- read.csv(file.path(audit,"DepMap_24Q2_Model.csv"),check.names=FALSE,stringsAsFactors=FALSE); crc_ids<-md$ModelID[md$OncotreePrimaryDisease=="Colorectal Adenocarcinoma" & md$ModelID%in%ids];
if(length(crc_ids)!=61)stop("QC stop: CRC ID count changed")
dep<-dep[match(crc_ids,ids),,drop=FALSE]; ids<-crc_ids

read_selected <- function(file, genes, sep=","){
 line<-readLines(file,n=1,warn=FALSE); h<-strsplit(line,sep,fixed=TRUE)[[1]]; h[1]<-"ModelID"; c<-sub(" \\(.*$","",h); ii<-which(c%in%genes); cl<-rep("NULL",length(h)); cl[1]<-"character"; cl[ii]<-"numeric"; z<-read.csv(file,header=FALSE,skip=1,check.names=FALSE,colClasses=cl,sep=sep); names(z)<-h[cl!="NULL"]; list(data=z,genes=sub(" \\(.*$","",names(z)[-1]))
}
ex<-read_selected(file.path(audit,"depmap24q2","OmicsExpressionProteinCodingGenesTPMLogp1.csv"),anchors); cn<-read_selected(file.path(raw,"DepMap_24Q2_OmicsAbsoluteCNGene.csv"),anchors)
em<-ex$data; cm<-cn$data; names(em)[1]<-"ModelID"; names(cm)[1]<-"ModelID"; me<-match(ids,em$ModelID); mc<-match(ids,cm$ModelID); keep<-!is.na(me)&!is.na(mc); ids<-ids[keep]; dep<-dep[keep,,drop=FALSE]; em<-em[me[keep],,drop=FALSE]; cm<-cm[mc[keep],,drop=FALSE]
if(!all(keep)) warning("Some CRC models lack expression or CNV data; complete matched cases are used and denominators are reported.")
getcol<-function(z,g){j<-which(sub(" \\(.*$","",names(z))==g); if(!length(j))return(rep(NA_real_,nrow(z))); as.numeric(z[[j[1]]])}
rr<-lapply(anchors,function(g){d<-getcol(dep,g); e<-getcol(em,g); c<-getcol(cm,g); ok<-is.finite(d)&is.finite(e); okc<-ok&is.finite(c); sp<-if(sum(ok)>=10)suppressWarnings(cor.test(d[ok],e[ok],method="spearman",exact=FALSE))else NULL; lm0<-if(sum(okc)>=15)summary(lm(d[okc]~e[okc]+c[okc]))else NULL; data.frame(gene=g,n_crc=sum(ok),spearman_rho=if(!is.null(sp))unname(sp$estimate)else NA_real_,spearman_p=if(!is.null(sp))sp$p.value else NA_real_,n_cnv_model=sum(okc),cnv_median=if(any(is.finite(c)))median(c[is.finite(c)])else NA_real_,expression_median=if(any(is.finite(e)))median(e[is.finite(e)])else NA_real_,lm_expression_beta=if(!is.null(lm0))coef(lm0)[2,1]else NA_real_,lm_expression_p=if(!is.null(lm0))coef(lm0)[2,4]else NA_real_,lm_cnv_beta=if(!is.null(lm0))coef(lm0)[3,1]else NA_real_,lm_cnv_p=if(!is.null(lm0))coef(lm0)[3,4]else NA_real_)})
rr<-do.call(rbind,rr); rr$BH_spearman<-p.adjust(rr$spearman_p,"BH"); write.table(rr,file.path(out,"depmap_crc_lipid_axis_expression_cnv_robustness.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# Background subgroup counts and candidate-gene tests; small groups are descriptive only.
mm<-md[match(ids,md$ModelID),,drop=FALSE]; groups<-list(OncotreeCode=c("COAD","READ"),PrimaryOrMetastasis=c("Primary","Metastatic"),LegacyMolecularSubtype=c("MSI","MSS"))
gc<-do.call(rbind,lapply(names(groups),function(v){z<-mm[[v]]; data.frame(variable=v,group=groups[[v]],n=vapply(groups[[v]],function(g)sum(z==g,na.rm=TRUE),integer(1)))})); write.table(gc,file.path(out,"depmap_crc_background_group_counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
subres<-list(); k<-1L
for(v in names(groups)){z<-mm[[v]]; for(gene in anchors){xv<-getcol(dep,gene); ok<-!is.na(z)&z%in%groups[[v]]&is.finite(xv); if(sum(ok)>=15&&length(unique(z[ok]))==2){wt<-suppressWarnings(wilcox.test(xv[ok]~z[ok],exact=FALSE)); subres[[k]]<-data.frame(variable=v,gene=gene,n=sum(ok),group1=groups[[v]][1],group2=groups[[v]][2],median_group1=median(xv[ok&z==groups[[v]][1]]),median_group2=median(xv[ok&z==groups[[v]][2]]),p=wt$p.value);k<-k+1L}}}
if(length(subres)){sr<-do.call(rbind,subres);sr$BH<-p.adjust(sr$p,"BH")}else sr<-data.frame(variable=character(),gene=character(),n=integer(),p=numeric(),BH=numeric()); write.table(sr,file.path(out,"depmap_crc_axis_background_tests.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Stage 4 robustness analysis completed.",paste("CRC models complete expression+CNV",length(ids)),paste("anchors",paste(anchors,collapse=",")),paste("subgroup tests",nrow(sr))),file.path(log,"stage4_axis_robustness.log"))
