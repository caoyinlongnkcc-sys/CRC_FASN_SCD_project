suppressPackageStartupMessages({library(httr); library(jsonlite); library(data.table)})
options(stringsAsFactors=FALSE)
root <- "."
raw <- file.path(root,"data_raw","tcga_cptac_validation")
res <- file.path(root,"results_tables")
dir.create(raw,recursive=TRUE,showWarnings=FALSE); dir.create(res,recursive=TRUE,showWarnings=FALSE)
genes <- data.table(entrezGeneId=c(2194,6319,31,6720), gene=c("FASN","SCD","ACACA","SREBF1"))
api <- "https://www.cbioportal.org/api"
get_json <- function(url){r<-GET(url); stop_for_status(r); fromJSON(content(r,"text",encoding="UTF-8"),flatten=TRUE)}
post_mol <- function(profile, sample_list, ids){
  u <- paste0(api,"/molecular-profiles/",profile,"/molecular-data/fetch?projection=SUMMARY")
  body <- toJSON(list(entrezGeneIds=as.integer(ids),sampleListId=sample_list),auto_unbox=TRUE)
  r <- POST(u,content_type_json(),body=body); stop_for_status(r)
  x <- fromJSON(content(r,"text",encoding="UTF-8"),flatten=TRUE)
  if(length(x)==0) return(data.table())
  as.data.table(x)
}
write_tsv <- function(x,f){fwrite(x,f,sep="\t",na="NA")}
study_info <- data.table()
for(st in c("coad_tcga_gdc","read_tcga_gdc","coad_cptac_2019")){
  s <- get_json(paste0(api,"/studies/",st)); study_info <- rbind(study_info, data.table(studyId=s$studyId,name=s$name,description=s$description),fill=TRUE)
}
write_tsv(study_info,file.path(raw,"study_info.tsv"))

fetch_tcga <- function(study, prefix){
  sl <- get_json(paste0(api,"/studies/",study,"/sample-lists")); sl <- as.data.table(sl)
  # Prefer all-sample list; explicitly record selected list.
  all <- sl[grepl("_all$",sampleListId)][1]
  if(is.na(all$sampleListId)) all <- sl[1]
  expr <- post_mol(paste0(study,"_mrna_seq_tpm"),all$sampleListId,genes$entrezGeneId)
  cna <- post_mol(paste0(study,"_cna"),all$sampleListId,genes$entrezGeneId)
  expr <- merge(expr,genes,by="entrezGeneId",all.x=TRUE); cna <- merge(cna,genes,by="entrezGeneId",all.x=TRUE)
  write_tsv(expr,file.path(raw,paste0(prefix,"_tpm.tsv"))); write_tsv(cna,file.path(raw,paste0(prefix,"_cna.tsv")))
  data.table(dataset=prefix,studyId=study,sampleListId=all$sampleListId,profile_expr=paste0(study,"_mrna_seq_tpm"),profile_cna=paste0(study,"_cna"),expr_rows=nrow(expr),cna_rows=nrow(cna),expr_samples=uniqueN(expr$sampleId),cna_samples=uniqueN(cna$sampleId))
}
tcga_manifest <- rbind(fetch_tcga("coad_tcga_gdc","TCGA_COAD"),fetch_tcga("read_tcga_gdc","TCGA_READ"),fill=TRUE)

fetch_cptac <- function(){
  study <- "coad_cptac_2019"; sl <- as.data.table(get_json(paste0(api,"/studies/",study,"/sample-lists"))); all <- sl[grepl("_all$",sampleListId)][1]; if(is.na(all$sampleListId)) all<-sl[1]
  p <- as.data.table(get_json(paste0(api,"/studies/",study,"/molecular-profiles")))
  prof <- p[molecularProfileId=="coad_cptac_2019_protein_quantification",molecularProfileId][1]
  x <- post_mol(prof,all$sampleListId,genes$entrezGeneId); x <- merge(x,genes,by="entrezGeneId",all.x=TRUE)
  write_tsv(x,file.path(raw,"CPTAC_COAD_protein.tsv"))
  data.table(dataset="CPTAC_COAD",studyId=study,sampleListId=all$sampleListId,profile_protein=prof,protein_rows=nrow(x),protein_samples=uniqueN(x$sampleId))
}
cptac_manifest <- fetch_cptac()
manifest <- rbind(tcga_manifest,cptac_manifest,fill=TRUE); write_tsv(manifest,file.path(res,"tcga_cptac_validation_manifest.tsv"))

fetch_clinical <- function(study, prefix, sample_list){
  u <- paste0(api,"/studies/",study,"/clinical-data?sampleListId=",sample_list)
  x <- as.data.table(get_json(u)); write_tsv(x,file.path(raw,paste0(prefix,"_clinical.tsv"))); x
}
fetch_patient_clinical <- function(study, prefix){
  u <- paste0(api,"/studies/",study,"/clinical-data?clinicalDataType=PATIENT")
  x <- as.data.table(get_json(u)); write_tsv(x,file.path(raw,paste0(prefix,"_patient_clinical.tsv"))); x
}
clin_coad <- fetch_clinical("coad_tcga_gdc","TCGA_COAD",tcga_manifest[dataset=="TCGA_COAD"]$sampleListId)
clin_read <- fetch_clinical("read_tcga_gdc","TCGA_READ",tcga_manifest[dataset=="TCGA_READ"]$sampleListId)
clin_cptac <- fetch_clinical("coad_cptac_2019","CPTAC_COAD",cptac_manifest$sampleListId)
patient_coad <- fetch_patient_clinical("coad_tcga_gdc","TCGA_COAD")
patient_read <- fetch_patient_clinical("read_tcga_gdc","TCGA_READ")
patient_cptac <- fetch_patient_clinical("coad_cptac_2019","CPTAC_COAD")

qc_one <- function(f,dataset,layer){
  x <- fread(f); if(!"gene"%in%names(x)) stop("gene field missing in ",f); x[,value_num:=as.numeric(value)]
  out <- x[,.(dataset=dataset,layer=layer,gene=gene,rows=.N,unique_samples=uniqueN(sampleId),missing_values=sum(is.na(value_num)),missing_rate=mean(is.na(value_num)),min= suppressWarnings(min(value_num,na.rm=TRUE)),max=suppressWarnings(max(value_num,na.rm=TRUE))),by=gene]
  out[!is.finite(min),min:=NA_real_]; out[!is.finite(max),max:=NA_real_]; out
}
qc <- rbind(qc_one(file.path(raw,"TCGA_COAD_tpm.tsv"),"TCGA_COAD","TPM"),qc_one(file.path(raw,"TCGA_COAD_cna.tsv"),"TCGA_COAD","CNA"),qc_one(file.path(raw,"TCGA_READ_tpm.tsv"),"TCGA_READ","TPM"),qc_one(file.path(raw,"TCGA_READ_cna.tsv"),"TCGA_READ","CNA"),qc_one(file.path(raw,"CPTAC_COAD_protein.tsv"),"CPTAC_COAD","protein"),fill=TRUE)
write_tsv(qc,file.path(res,"tcga_cptac_validation_qc.tsv"))

cor_stats <- function(expr_file,cna_file,dataset){
  e <- fread(expr_file); c <- fread(cna_file); out <- list()
  for(g in genes$gene){
    ee <- e[gene==g,.(sampleId,value_expr=as.numeric(value))]; cc <- c[gene==g,.(sampleId,value_cna=as.numeric(value))]; z <- merge(ee,cc,by="sampleId");
    if(nrow(z)>=10) { ct <- suppressWarnings(cor.test(z$value_expr,z$value_cna,method="spearman",exact=FALSE)); out[[g]] <- data.table(dataset=dataset,gene=g,n=nrow(z),spearman_rho=unname(ct$estimate),p_value=ct$p.value) }
  }; rbindlist(out,fill=TRUE)
}
cor_tbl <- rbind(cor_stats(file.path(raw,"TCGA_COAD_tpm.tsv"),file.path(raw,"TCGA_COAD_cna.tsv"),"TCGA_COAD"),cor_stats(file.path(raw,"TCGA_READ_tpm.tsv"),file.path(raw,"TCGA_READ_cna.tsv"),"TCGA_READ"),fill=TRUE)
cor_tbl[,p_adj:=p.adjust(p_value,"BH"),by=dataset]; write_tsv(cor_tbl,file.path(res,"tcga_expression_cna_correlations.tsv"))

clinical_summary <- rbindlist(lapply(list(clin_coad,clin_read,clin_cptac,patient_coad,patient_read,patient_cptac),function(x){
  x[clinicalAttributeId %in% c("SAMPLE_TYPE","PATH_STAGE","PATH_T_STAGE","PATH_N_STAGE","AJCC_PATHOLOGIC_STAGE","MSI_STATUS","STAGE","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS"),.(n_records=.N,n_samples=if("sampleId"%in%names(x)) uniqueN(sampleId) else NA_integer_,n_patients=if("patientId"%in%names(x)) uniqueN(patientId) else NA_integer_),by=.(clinicalAttributeId)]
}),fill=TRUE); write_tsv(clinical_summary,file.path(res,"tcga_cptac_clinical_attribute_availability.tsv"))
capture.output(sessionInfo(),file=file.path(root,"logs","tcga_cptac_validation_sessionInfo.txt"))
cat("Completed TCGA/CPTAC validation retrieval and QC\n")
