suppressPackageStartupMessages(library(data.table))
root <- "."; audit <- file.path("..","feasibility_audit"); raw <- file.path(root,"data_raw"); out <- file.path(root,"results_tables"); logd <- file.path(root,"logs")
anchors <- c("FASN","SCD","ACACA","SREBF1")
depfile <- file.path(audit,"depmap24q2","CRISPRGeneEffect.csv")
exprfile <- file.path(audit,"depmap24q2","OmicsExpressionProteinCodingGenesTPMLogp1.csv")
modfile <- file.path(audit,"DepMap_24Q2_Model.csv")
if(!all(file.exists(c(depfile,exprfile,modfile)))) stop("Required official DepMap files are missing")
dh <- names(fread(depfile,nrows=0,check.names=FALSE)); didx <- c(1,which(sub(" \\(.*$","",dh)%in%anchors)); dep <- fread(depfile,select=didx,check.names=FALSE); names(dep)[1] <- "ModelID"; names(dep)[-1] <- sub(" \\(.*$","",names(dep)[-1])
md <- fread(modfile); crc_ids <- intersect(md[OncotreePrimaryDisease=="Colorectal Adenocarcinoma",ModelID],dep$ModelID); dep <- dep[ModelID%in%crc_ids]; if(nrow(dep)!=61) stop("CRC model count changed: ",nrow(dep)); for(g in anchors) dep[[g]] <- as.numeric(dep[[g]])
expr <- fread(exprfile,check.names=FALSE); names(expr)[1] <- "ModelID"; expr <- expr[ModelID%in%crc_ids]; if(nrow(expr)<50) stop("Too few CRC expression models")
names(expr) <- make.unique(sub(" \\(.*$","",names(expr))); m <- match(dep$ModelID,expr$ModelID); keep <- !is.na(m); dep <- dep[keep]; expr <- expr[m[keep]]; if(nrow(dep)<50) stop("Too few matched CRC models")
X <- as.matrix(expr[,-1,with=FALSE]); storage.mode(X) <- "numeric"; genes <- colnames(X)
rank_test <- function(groupvar,label){
  med <- median(dep[[groupvar]],na.rm=TRUE); grp <- ifelse(dep[[groupvar]]<=med,"high_dependency","low_dependency"); tab <- table(grp); if(any(tab<15)) stop(label," group below 15")
  hi <- grp=="high_dependency"; lo <- grp=="low_dependency"; st <- pv <- lfc <- rep(NA_real_,ncol(X))
  for(j in seq_len(ncol(X))){a<-X[hi,j]; b<-X[lo,j]; a<-a[is.finite(a)]; b<-b[is.finite(b)]; if(length(a)>=5&&length(b)>=5){tt<-suppressWarnings(t.test(a,b)); st[j]<-unname(tt$statistic); pv[j]<-tt$p.value; lfc[j]<-mean(a)-mean(b)}}
  z <- data.table(gene=genes,grouping=label,high_dependency_n=sum(hi),low_dependency_n=sum(lo),mean_difference=lfc,t_stat=st,p_value=pv); z[,p_adj:=p.adjust(p_value,"BH")]; z[order(p_adj,p_value)]
}
de_scd <- rank_test("SCD","SCD_dependency_median_split"); de_fasn <- rank_test("FASN","FASN_dependency_median_split")
fwrite(de_scd,file.path(out,"depmap_crc_scd_dependency_expression_results.tsv"),sep="\t",na="NA"); fwrite(de_fasn,file.path(out,"depmap_crc_fasn_dependency_expression_results.tsv"),sep="\t",na="NA")
fwrite(rbind(de_scd[gene%in%anchors],de_fasn[gene%in%anchors]),file.path(out,"depmap_crc_anchor_mechanism_checks.tsv"),sep="\t",na="NA")
gmt <- file.path(raw,"ReactomePathways.gmt.zip"); if(!file.exists(gmt)) stop("Reactome GMT missing")
tmpdir <- tempfile("gmt_"); dir.create(tmpdir); unzip(gmt,exdir=tmpdir); gfile <- file.path(tmpdir,"ReactomePathways.gmt"); lines <- readLines(gfile,warn=FALSE); sets <- lapply(lines,function(z){q<-strsplit(z,"\t",fixed=TRUE)[[1]]; unique(q[-c(1,2)])}); names(sets) <- vapply(lines,function(z)strsplit(z,"\t",fixed=TRUE)[[1]][1],character(1))
path_test <- function(de,label){s<-de$t_stat; names(s)<-de$gene; s<-s[is.finite(s)]; ps<-lapply(sets,function(g)intersect(g,names(s))); ps<-ps[vapply(ps,length,integer(1))>=10 & vapply(ps,length,integer(1))<=500]; z<-as.data.table(fgsea::fgsea(pathways=ps,stats=s,minSize=10,maxSize=500,nperm=1000)); setnames(z,c("pathway","size","ES","NES","pval","padj"),c("pathway","n_genes","ES","NES","p_value","p_adj"),skip_absent=TRUE); z[,grouping:=label]; z[order(p_adj,p_value)]}
pt <- rbind(path_test(de_scd,"SCD_dependency_median_split"),path_test(de_fasn,"FASN_dependency_median_split")); fwrite(pt,file.path(out,"depmap_crc_mechanism_reactome_rank_audit.tsv"),sep="\t",na="NA")
avail <- data.table(resource=c("DepMap_CRISPR_dependency","DepMap_RNA_expression","DepMap_model_metadata","DepMap_mutation_table","Reactome_GMT"),available=c(TRUE,TRUE,TRUE,FALSE,TRUE),matched_crc_models=c(nrow(dep),nrow(expr),sum(crc_ids%in%md$ModelID),NA_integer_,NA_integer_),note=c("official 24Q2 CRISPRGeneEffect; exact CRC","official TPMLogp1; matched by ModelID","official model metadata","not present in audited local directory; mutation analysis not run","official downloaded Reactome GMT")); fwrite(avail,file.path(out,"mechanism_feasibility_audit.tsv"),sep="\t",na="NA")
capture.output(sessionInfo(),file=file.path(logd,"mechanism_feasibility_sessionInfo.txt")); cat("Mechanism feasibility audit completed; matched CRC models:",nrow(dep),"\n")
