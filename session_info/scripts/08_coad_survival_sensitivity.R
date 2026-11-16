suppressPackageStartupMessages({library(data.table); library(survival)})
root <- "."; out <- file.path(root,"results_tables"); logd <- file.path(root,"logs")
set.seed(20261112)
d <- fread(file.path(out,"TCGA_COAD_clinical_analysis_matrix.tsv"))
need <- c("patientId","OS_MONTHS","OS_STATUS","SCD","FASN_SCD_score","Lipid4_score","AGE","SEX","PATH_STAGE_GROUP")
missing <- setdiff(need,names(d)); if(length(missing)) stop("Missing required fields: ",paste(missing,collapse=", "))
d[,time:=as.numeric(OS_MONTHS)]; d[,event:=as.integer(grepl("^1",OS_STATUS))]; d[,AGE_num:=as.numeric(AGE)]; d[,SEX:=factor(SEX)]; d[,stage:=factor(PATH_STAGE_GROUP,levels=c("I","II","III","IV"))]
d[,SCD_z:=as.numeric(scale(SCD))]; d[,FASN_SCD_z:=as.numeric(scale(FASN_SCD_score))]; d[,Lipid4_z:=as.numeric(scale(Lipid4_score))]
base <- d[is.finite(time)&time>=0&event%in%c(0,1)&is.finite(AGE_num)&!is.na(SEX)&!is.na(stage)]
fit_adj <- function(xvar,label){
  z <- base[is.finite(get(xvar))]; f <- as.formula(paste0("Surv(time,event) ~ ",xvar," + AGE_num + SEX + stage")); m <- coxph(f,data=z); sm<-summary(m); cz<-cox.zph(m); row<-which(rownames(sm$coefficients)==xvar); data.table(model=label,n=nrow(z),events=sum(z$event),hr=sm$coefficients[row,"exp(coef)"],lower95=sm$conf.int[row,"lower .95"],upper95=sm$conf.int[row,"upper .95"],p_value=sm$coefficients[row,"Pr(>|z|)"],global_ph_p=cz$table["GLOBAL","p"])
}
adj <- rbind(fit_adj("SCD_z","SCD_adjusted"),fit_adj("FASN_SCD_z","FASN_SCD_adjusted"),fit_adj("Lipid4_z","Lipid4_adjusted"),fill=TRUE); adj[,p_adj:=p.adjust(p_value,"BH")]
fwrite(adj,file.path(out,"tcga_coad_adjusted_cox_results.tsv"),sep="\t",na="NA")

z <- base[is.finite(SCD_z)]
tv <- coxph(Surv(time,event) ~ SCD_z + tt(SCD_z) + AGE_num + SEX + stage,data=z,tt=function(x,t,...) x*log(pmax(t,0.1)))
tvsm <- summary(tv); tvrow <- data.table(term=rownames(tvsm$coefficients),coef=tvsm$coefficients[,"coef"],hr=tvsm$coefficients[,"exp(coef)"],lower95=tvsm$conf.int[,"lower .95"],upper95=tvsm$conf.int[,"upper .95"],p_value=tvsm$coefficients[,"Pr(>|z|)"],n=nrow(z),events=sum(z$event)); fwrite(tvrow,file.path(out,"tcga_coad_scd_time_varying_cox.tsv"),sep="\t",na="NA")

tau <- 60; if(max(z$time,na.rm=TRUE)<tau) stop("Prespecified tau=60 months exceeds observed follow-up")
z[,group:=ifelse(SCD<=median(SCD,na.rm=TRUE),"low_or_equal","high")]; z[,group:=factor(group,levels=c("low_or_equal","high"))]
rmst_one <- function(dd,tau){ sf<-survfit(Surv(time,event)~1,data=dd); ss<-summary(sf,times=sort(unique(c(0,sf$time[sf$time<tau],tau))),extend=TRUE); tt<-ss$time; survv<-ss$surv; if(length(tt)<2) return(NA_real_); sum(diff(tt)*survv[-length(survv)]) }
r0<-rmst_one(z[group=="low_or_equal"],tau); r1<-rmst_one(z[group=="high"],tau); diff0<-r1-r0
boot <- numeric(500); for(i in seq_along(boot)){ii<-sample.int(nrow(z),nrow(z),replace=TRUE); b<-z[ii]; boot[i]<-rmst_one(b[group=="high"],tau)-rmst_one(b[group=="low_or_equal"],tau)}
rmst <- data.table(dataset="TCGA_COAD",feature="SCD",tau_months=tau,cutpoint=median(z$SCD,na.rm=TRUE),low_n=sum(z$group=="low_or_equal"),high_n=sum(z$group=="high"),low_events=sum(z[group=="low_or_equal"]$event),high_events=sum(z[group=="high"]$event),rmst_low=r0,rmst_high=r1,rmst_difference_high_minus_low=diff0,bootstrap_lower95=quantile(boot,0.025,na.rm=TRUE),bootstrap_upper95=quantile(boot,0.975,na.rm=TRUE),bootstrap_reps=sum(is.finite(boot)))
fwrite(rmst,file.path(out,"tcga_coad_scd_rmst_60m.tsv"),sep="\t",na="NA")
capture.output(sessionInfo(),file=file.path(logd,"coad_survival_sensitivity_sessionInfo.txt")); cat("COAD adjusted, time-varying and RMST sensitivity analyses completed\n")
