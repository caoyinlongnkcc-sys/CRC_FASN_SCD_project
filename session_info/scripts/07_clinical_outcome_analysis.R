suppressPackageStartupMessages({library(data.table); library(survival)})
root <- "."; raw <- file.path(root,"data_raw","tcga_cptac_validation"); out <- file.path(root,"results_tables"); logd <- file.path(root,"logs")
genes <- c("FASN","SCD","ACACA","SREBF1")
read_clin <- function(prefix, patient=FALSE){fread(file.path(raw,paste0(prefix,if(patient)"_patient_clinical.tsv" else "_clinical.tsv")))}
wide_clin <- function(x){
  x <- x[clinicalAttributeId %in% c("OS_MONTHS","OS_STATUS","DFS_MONTHS","DFS_STATUS","PATH_STAGE","PATH_T_STAGE","PATH_N_STAGE","STAGE","MSI_STATUS","SAMPLE_TYPE","AGE","SEX")]
  dcast(x, if("sampleId"%in%names(x)) sampleId+patientId ~ clinicalAttributeId else patientId ~ clinicalAttributeId, value.var="value", fun.aggregate=function(z) z[1])
}
select_samples <- function(expr, sampleclin){
  ids <- unique(expr[,.(sampleId,patientId)])
  z <- merge(ids,sampleclin[,intersect(names(sampleclin),c("sampleId","SAMPLE_TYPE")),with=FALSE],by="sampleId",all.x=TRUE)
  if("SAMPLE_TYPE"%in%names(z) && any(z$SAMPLE_TYPE=="Primary Tumor",na.rm=TRUE)) z <- z[SAMPLE_TYPE=="Primary Tumor"]
  z[order(patientId,sampleId),.SD[1],by=patientId]
}
make_expr <- function(file, sampleclin){
  x <- fread(file); x[,value:=as.numeric(value)]; x[,expr:=log2(value+1)]
  ids <- select_samples(x,sampleclin); x <- merge(x,ids[,.(patientId,sampleId)],by=c("patientId","sampleId"));
  w <- dcast(x,patientId~gene,value.var="expr"); w
}
surv_one <- function(d, timevar="OS_MONTHS", statusvar="OS_STATUS", predictors=c(genes,"FASN_SCD_score","Lipid4_score"), dataset=""){
  ans <- list(); ph <- list()
  for(v in predictors){
    if(!v%in%names(d)) next
    z <- d[is.finite(get(v)) & !is.na(get(timevar)) & !is.na(get(statusvar))]
    z[,time:=as.numeric(get(timevar))]; z[,event:=as.integer(grepl("^1",get(statusvar)))]
    z <- z[is.finite(time) & time>=0 & event%in%c(0,1)]
    if(nrow(z)<30 || sum(z$event)<10) next
    z[,x:=as.numeric(scale(get(v)))]; fit <- coxph(Surv(time,event)~x,data=z)
    s <- summary(fit); cz <- tryCatch(cox.zph(fit),error=function(e)NULL); ph_p <- if(!is.null(cz)) cz$table["x","p"] else NA_real_
    ans[[v]] <- data.table(dataset=dataset,endpoint=timevar,feature=v,n=nrow(z),events=sum(z$event),hr_per_sd=unname(s$coefficients[1,"exp(coef)"]),lower95=unname(s$conf.int[1,"lower .95"]),upper95=unname(s$conf.int[1,"upper .95"]),p_value=unname(s$coefficients[1,"Pr(>|z|)"]),ph_test_p=ph_p)
  }; if(!length(ans)) return(data.table(dataset=character(),endpoint=character(),feature=character(),n=integer(),events=integer(),hr_per_sd=numeric(),lower95=numeric(),upper95=numeric(),p_value=numeric(),ph_test_p=numeric())); rbindlist(ans,fill=TRUE)
}
collapse_stage <- function(x, type){
  x <- toupper(trimws(as.character(x))); out <- rep(NA_character_,length(x))
  if(type=="stage") out[grepl("^STAGE I",x)] <- "I"; if(type=="stage") out[grepl("^STAGE II",x)] <- "II"; if(type=="stage") out[grepl("^STAGE III",x)] <- "III"; if(type=="stage") out[grepl("^STAGE IV",x)] <- "IV"
  if(type=="t") out[grepl("^TIS",x)] <- "Tis"; if(type=="t") out[grepl("^T[1]",x)] <- "T1"; if(type=="t") out[grepl("^T[2]",x)] <- "T2"; if(type=="t") out[grepl("^T[3]",x)] <- "T3"; if(type=="t") out[grepl("^T4",x)] <- "T4"
  if(type=="n") out[grepl("^N0",x)] <- "N0"; if(type=="n") out[grepl("^N1",x)] <- "N1"; if(type=="n") out[grepl("^N2",x)] <- "N2"; out
}
stage_one <- function(d, var, features=c(genes,"FASN_SCD_score","Lipid4_score"),dataset=""){
  if(!var%in%names(d)) return(data.table())
  ans<-list();
  for(v in features){if(!v%in%names(d))next; z<-d[!is.na(get(v)) & !is.na(get(var)) & get(var)!="",.(x=get(v),g=get(var))]; if(uniqueN(z$g)<2)next; kt<-kruskal.test(x~g,data=z); ans[[v]]<-data.table(dataset=dataset,clinical_field=var,feature=v,n=nrow(z),groups=uniqueN(z$g),p_value=kt$p.value,group_values=paste(sort(unique(z$g)),collapse=";"))}
  rbindlist(ans,fill=TRUE)
}
analyze_tcga <- function(prefix,exprfile,clinfile,patientfile,dataset){
  e <- fread(exprfile); sc <- wide_clin(fread(clinfile)); pc <- wide_clin(fread(patientfile));
  w <- make_expr(exprfile,sc); d <- merge(w,pc,by="patientId",all.x=FALSE)
  for(v in genes) if(v%in%names(d)) d[[v]] <- as.numeric(d[[v]])
  d[,FASN_SCD_score:=rowMeans(.SD,na.rm=FALSE),.SDcols=c("FASN","SCD")]; d[,Lipid4_score:=rowMeans(.SD,na.rm=FALSE),.SDcols=genes]
  if("PATH_STAGE"%in%names(d)) d[,PATH_STAGE_GROUP:=collapse_stage(PATH_STAGE,"stage")]; if("PATH_T_STAGE"%in%names(d)) d[,PATH_T_STAGE_GROUP:=collapse_stage(PATH_T_STAGE,"t")]; if("PATH_N_STAGE"%in%names(d)) d[,PATH_N_STAGE_GROUP:=collapse_stage(PATH_N_STAGE,"n")]
  fwrite(d,file.path(out,paste0(dataset,"_clinical_analysis_matrix.tsv")),sep="\t",na="NA")
  s <- surv_one(d,dataset=dataset); if(nrow(s)) s[,p_adj:=p.adjust(p_value,"BH")]
  st <- rbind(stage_one(d,"PATH_STAGE_GROUP",dataset=dataset),stage_one(d,"PATH_T_STAGE_GROUP",dataset=dataset),stage_one(d,"PATH_N_STAGE_GROUP",dataset=dataset)); if(nrow(st)) st[,p_adj:=p.adjust(p_value,"BH"),by=clinical_field]
  list(surv=s,stage=st,n=nrow(d))
}
coad <- analyze_tcga("TCGA_COAD",file.path(raw,"TCGA_COAD_tpm.tsv"),file.path(raw,"TCGA_COAD_clinical.tsv"),file.path(raw,"TCGA_COAD_patient_clinical.tsv"),"TCGA_COAD")
read <- analyze_tcga("TCGA_READ",file.path(raw,"TCGA_READ_tpm.tsv"),file.path(raw,"TCGA_READ_clinical.tsv"),file.path(raw,"TCGA_READ_patient_clinical.tsv"),"TCGA_READ")
surv <- rbind(coad$surv,read$surv,fill=TRUE); stage <- rbind(coad$stage,read$stage,fill=TRUE)
fwrite(surv,file.path(out,"tcga_os_cox_results.tsv"),sep="\t",na="NA"); fwrite(stage,file.path(out,"tcga_stage_association_results.tsv"),sep="\t",na="NA")

# CPTAC protein survival: only assayed genes are analyzed; no imputation for absent proteins.
prot <- fread(file.path(raw,"CPTAC_COAD_protein.tsv")); pc <- wide_clin(read_clin("CPTAC_COAD",TRUE));
prot[,expr:=as.numeric(value)]; pw <- dcast(prot,patientId~gene,value.var="expr"); d <- merge(pw,pc,by="patientId");
if(all(c("FASN","ACACA")%in%names(d))) d[,FASN_ACACA_score:=rowMeans(.SD,na.rm=FALSE),.SDcols=c("FASN","ACACA")]
cptac_surv <- surv_one(d,predictors=intersect(c("FASN","ACACA","FASN_ACACA_score"),names(d)),dataset="CPTAC_COAD_protein")
if(nrow(cptac_surv)) cptac_surv[,p_adj:=p.adjust(p_value,"BH")]
cptac_stage <- rbind(stage_one(d,"STAGE",features=intersect(c("FASN","ACACA","FASN_ACACA_score"),names(d)),dataset="CPTAC_COAD_protein"),stage_one(d,"PATHOLOGY_STATUS",features=intersect(c("FASN","ACACA","FASN_ACACA_score"),names(d)),dataset="CPTAC_COAD_protein")); if(nrow(cptac_stage)) cptac_stage[,p_adj:=p.adjust(p_value,"BH"),by=clinical_field]
fwrite(cptac_surv,file.path(out,"cptac_protein_os_cox_results.tsv"),sep="\t",na="NA"); fwrite(cptac_stage,file.path(out,"cptac_protein_stage_association_results.tsv"),sep="\t",na="NA")
capture.output(sessionInfo(),file=file.path(logd,"clinical_outcome_sessionInfo.txt")); cat("Clinical outcome analysis completed\n")
