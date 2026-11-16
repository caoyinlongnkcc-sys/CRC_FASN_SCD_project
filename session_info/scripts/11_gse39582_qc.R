suppressPackageStartupMessages({library(data.table); library(hgu133plus2.db); library(AnnotationDbi)})
root <- "."; f <- file.path(root,"data_raw","GSE39582_series_matrix.txt.gz"); out<-file.path(root,"results_tables");logd<-file.path(root,"logs")
if(!file.exists(f))stop("GSE39582 matrix is missing")
con <- gzfile(f,"rt"); raw <- readLines(con,warn=FALSE); close(con)
begin <- which(grepl("!series_matrix_table_begin",raw,fixed=TRUE)); end <- which(grepl("!series_matrix_table_end",raw,fixed=TRUE)); if(length(begin)!=1||length(end)!=1)stop("Series Matrix boundaries not found")
tab <- fread(text=paste(raw[(begin+1):(end-1)],collapse="\n"),sep="\t",header=TRUE,check.names=FALSE,quote="\"",data.table=TRUE)
sample_line <- raw[grepl("!Sample_geo_accession",raw,fixed=TRUE)][1]; sample_ids <- strsplit(sub("^!Sample_geo_accession\\t","",sample_line),"\\t")[[1]]; sample_ids <- gsub('^"|"$','',sample_ids)
if(ncol(tab)-1 != length(sample_ids))stop("Expression columns do not match GEO sample IDs")
chars <- raw[grepl("!Sample_characteristics_ch1",raw,fixed=TRUE)]; char_mat <- lapply(chars,function(z){q<-strsplit(sub("^!Sample_characteristics_ch1\\t","",z),"\\t")[[1]];gsub('^"|"$','',q)})
char_df <- as.data.table(do.call(rbind,char_mat)); char_df <- transpose(char_df); setnames(char_df,paste0("char",seq_len(ncol(char_df)))); char_df[,sample_id:=sample_ids]
probe <- as.character(tab[[1]]); expr <- as.matrix(tab[,-1,with=FALSE]); storage.mode(expr)<-"numeric"; colnames(expr)<-sample_ids
map <- AnnotationDbi::select(hgu133plus2.db,keys=probe,keytype="PROBEID",columns=c("SYMBOL","ENTREZID")); map <- as.data.table(map)[!is.na(SYMBOL)&SYMBOL!=""]; gene_map <- map[SYMBOL%in%c("FASN","SCD"),.(PROBEID,SYMBOL,ENTREZID)]; fwrite(gene_map,file.path(out,"GSE39582_FASN_SCD_probe_mapping.tsv"),sep="\t",na="NA")
probe_idx <- match(gene_map$PROBEID,probe); detected <- gene_map[!is.na(probe_idx)]
qc <- data.table(accession="GSE39582",file_sha256="5B88FA71C5268FA1A6AC52B062EE6E676DEF6C997AB4A09C36EA7ACD90EF3968",n_probes=nrow(tab),n_samples=ncol(expr),n_sample_ids=length(sample_ids),duplicate_sample_ids=anyDuplicated(sample_ids),FASN_probes=sum(detected$SYMBOL=="FASN"),SCD_probes=sum(detected$SYMBOL=="SCD"),note="Series Matrix and GPL570 Bioconductor annotation; clinical fields require manual parsing")
fwrite(qc,file.path(out,"GSE39582_qc_summary.tsv"),sep="\t",na="NA"); fwrite(char_df,file.path(out,"GSE39582_sample_characteristics_raw.tsv"),sep="\t",na="NA")
capture.output(sessionInfo(),file=file.path(logd,"gse39582_qc_sessionInfo.txt")); cat("GSE39582 matrix QC completed; samples",ncol(expr),"probes",nrow(tab),"\n")
