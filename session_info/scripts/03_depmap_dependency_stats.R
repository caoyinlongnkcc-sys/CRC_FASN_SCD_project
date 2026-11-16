#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."; audit <- file.path("..","feasibility_audit"); meta <- file.path(root,"data_metadata"); out <- file.path(root,"results_tables"); log <- file.path(root,"logs")
gs <- read.delim(file.path(meta,"lipid_gene_set.tsv"),sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
gs <- gs[gs$tier=="core",,drop=FALSE]
genes <- unique(gs$gene)
f <- file.path(audit,"depmap24q2","CRISPRGeneEffect.csv")
h <- names(read.csv(f,nrows=0,check.names=FALSE)); clean <- sub(" \\(.*$","",h); gi <- which(clean%in%genes)
classes <- rep("NULL",length(h)); classes[1]<-"character"; classes[gi]<-"numeric"
x <- read.csv(f,check.names=FALSE,colClasses=classes); ids <- as.character(x[[1]])
md <- read.csv(file.path(audit,"DepMap_24Q2_Model.csv"),check.names=FALSE,stringsAsFactors=FALSE)
crc <- ids %in% md$ModelID[md$OncotreePrimaryDisease=="Colorectal Adenocarcinoma"]
if(sum(crc)!=61) stop(paste("QC stop: expected 61 exact CRC models, observed",sum(crc)))
noncrc <- !crc
if(sum(noncrc)<100) stop("QC stop: non-CRC comparator is unexpectedly small")
cols <- names(x)[-1]; g <- sub(" \\(.*$","",cols)
common_dep <- read.delim(file.path(root,"data_raw","DepMap_24Q2_CRISPRInferredCommonEssentials.csv"),header=TRUE,stringsAsFactors=FALSE)$Essentials
common_dep <- sub(" \\(.*$","",common_dep)
res <- lapply(seq_along(g),function(j){
  v <- as.numeric(x[[j+1]]); a <- v[crc]; b <- v[noncrc]; a <- a[is.finite(a)]; b <- b[is.finite(b)]
  p <- if(length(a)>=10&&length(b)>=10) suppressWarnings(wilcox.test(a,b,exact=FALSE)$p.value) else NA_real_
  data.frame(gene=g[j], n_crc=length(a), n_noncrc=length(b), missing_crc=1-length(a)/sum(crc),
    median_crc=if(length(a))median(a) else NA_real_, mean_crc=if(length(a))mean(a) else NA_real_,
    sd_crc=if(length(a)>1)sd(a) else NA_real_, prop_crc_below_neg0_5=if(length(a))mean(a< -0.5) else NA_real_,
    median_noncrc=if(length(b))median(b) else NA_real_, median_difference=if(length(a)&&length(b))median(a)-median(b) else NA_real_,
    p_wilcoxon=p, depmap_common_essential=as.integer(g[j]%in%common_dep), stringsAsFactors=FALSE)
})
r <- do.call(rbind,res); r$fdr_bh <- p.adjust(r$p_wilcoxon,method="BH"); r$crc_selective <- !is.na(r$fdr_bh)&r$fdr_bh<0.05&r$median_difference<0
r <- merge(gs,r,by="gene",all.x=TRUE,sort=FALSE)
write.table(r,file.path(out,"depmap_crc_lipid_dependency_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
summ <- data.frame(total_core_genes=length(genes),detected_genes=sum(!is.na(r$n_crc)),tested_genes=sum(!is.na(r$p_wilcoxon)),fdr_lt_0_05=sum(r$crc_selective,na.rm=TRUE),median_crc_models=sum(crc),median_noncrc_models=sum(noncrc),stringsAsFactors=FALSE)
write.table(summ,file.path(out,"depmap_crc_lipid_dependency_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
top <- r[order(r$fdr_bh, r$median_crc),]; write.table(head(top,100),file.path(out,"depmap_crc_lipid_dependency_top100.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("DepMap Stage 3 completed.",paste("CRC models",sum(crc)),paste("non-CRC comparator",sum(noncrc)),paste("core genes",length(genes)),paste("tested genes",sum(!is.na(r$p_wilcoxon))),paste("CRC-selective BH<0.05",sum(r$crc_selective,na.rm=TRUE))),file.path(log,"stage3_depmap_stats.log"))
