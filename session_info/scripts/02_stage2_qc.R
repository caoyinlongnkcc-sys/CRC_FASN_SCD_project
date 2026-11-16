#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."; audit <- file.path("..","feasibility_audit"); meta <- file.path(root,"data_metadata"); out <- file.path(root,"results_tables"); log <- file.path(root,"logs"); dir.create(out,showWarnings=FALSE)
gs <- read.delim(file.path(meta,"lipid_gene_set.tsv"),sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
genes <- unique(gs$gene); report <- character(); add <- function(...) report <<- c(report,paste0(...))

# ---- DepMap CRISPR structure and joins ----
df <- file.path(audit,"depmap24q2","CRISPRGeneEffect.csv")
dh <- names(read.csv(df,nrows=0,check.names=FALSE))
dclean <- sub(" \\(.*$","",dh)
gidx <- which(dclean %in% genes); keep <- sort(unique(c(1,gidx)))
classes <- rep("NULL",length(dh)); classes[1] <- "character"; classes[gidx] <- "numeric"
dep <- read.csv(df,header=TRUE,check.names=FALSE,colClasses=classes)
dep_ids <- as.character(dep[[1]])
md <- read.csv(file.path(audit,"DepMap_24Q2_Model.csv"),check.names=FALSE,stringsAsFactors=FALSE)
dup_dep <- unique(dep_ids[duplicated(dep_ids)&!is.na(dep_ids)&nzchar(dep_ids)])
miss_meta <- setdiff(dep_ids,md$ModelID)
crc_meta <- md$OncotreePrimaryDisease=="Colorectal Adenocarcinoma"
crc_ids <- intersect(dep_ids,md$ModelID[crc_meta])
add("DepMap raw rows=",length(dep_ids),"; columns retained=",length(keep),"; frozen genes detected=",length(gidx),"/",length(genes))
add("DepMap duplicate model IDs=",length(dup_dep),"; IDs absent from model metadata=",length(miss_meta))
add("DepMap exact colorectal adenocarcinoma models after ID join=",length(crc_ids))
if(length(gidx)>0){ dgene_names <- sub(" \\(.*$","",names(dep)[-1]); dsel <- seq_along(dgene_names)[dgene_names %in% genes] + 1L; miss <- data.frame(gene=dgene_names[dsel-1L], missing_rate=colMeans(is.na(dep[,dsel,drop=FALSE])), n_nonmissing=colSums(!is.na(dep[,dsel,drop=FALSE]))); write.table(miss,file.path(out,"depmap_gene_missingness.tsv"),sep="\t",quote=FALSE,row.names=FALSE) }
crc_dep <- dep[match(crc_ids,dep_ids),,drop=FALSE]
crc_dep <- cbind(ModelID=crc_ids,as.data.frame(crc_dep[,-1,drop=FALSE],check.names=FALSE))
if(nrow(crc_dep)>0){ cm <- data.frame(gene=sub(" \\(.*$","",names(crc_dep)[-1]),missing_rate=colMeans(is.na(crc_dep[,-1,drop=FALSE])),n_nonmissing=colSums(!is.na(crc_dep[,-1,drop=FALSE]))); write.table(cm,file.path(out,"depmap_crc_gene_missingness.tsv"),sep="\t",quote=FALSE,row.names=FALSE); add("DepMap CRC genes with >20% missingness=",sum(cm$missing_rate>0.2),"/",nrow(cm)) }
write.table(data.frame(ModelID=crc_ids),file.path(out,"depmap_crc_model_ids.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# ---- Sanger organoid fitness structure ----
of <- file.path(audit,"organoid2026","fitness","Organoids_fitness_scores_fold_change_values_20260723.tsv")
oraw <- read.delim(of,header=FALSE,check.names=FALSE,quote="",fill=TRUE,stringsAsFactors=FALSE)
first <- as.character(oraw[1,]); second <- as.character(oraw[2,]);
idx <- which(seq_along(first)>3 & !is.na(first) & nzchar(first) & !is.na(second) & nzchar(second))
onames <- first[idx]; oids <- second[idx]
qcpass <- as.character(oraw[4,idx]); add("Organoid QC-pass flags=",sum(qcpass=="TRUE"),"/",length(qcpass))
olab <- as.character(oraw[-c(1,2,3,4),2]); ogene <- sub(" \\(.*$","",olab)
ogidx <- which(ogene %in% genes)
add("Organoid fitness raw dimensions=",nrow(oraw),"x",ncol(oraw),"; model columns with name+ID=",length(idx),"; frozen genes detected=",length(unique(ogene[ogidx])),"/",length(genes))
add("Organoid duplicate model IDs=",sum(duplicated(oids)),"; duplicate model names=",sum(duplicated(onames)),"; missing model IDs=",sum(is.na(oids)|!nzchar(oids)))
if(length(ogidx)>0){ vals <- oraw[ogidx+4,idx,drop=FALSE]; vals <- as.data.frame(lapply(vals,as.numeric)); names(vals)<-oids; rownames(vals)<-ogene[ogidx]; miss <- data.frame(gene=rownames(vals),missing_rate=rowMeans(is.na(vals)),n_nonmissing=rowSums(!is.na(vals))); write.table(miss,file.path(out,"organoid_gene_missingness.tsv"),sep="\t",quote=FALSE,row.names=FALSE) }

# Published organoid model metadata and clinical annotation join audit.
ml <- read.csv(gzfile(file.path(audit,"organoid2026","model_list_latest.csv.gz")),check.names=FALSE,stringsAsFactors=FALSE)
add("Organoid model-list rows=",nrow(ml),"; columns=",paste(names(ml),collapse=","))
add("Organoid fitness IDs present in model-list model_id=",sum(oids%in%ml$model_id),"/",length(oids))
add("Organoid fitness names present in model-list model_name=",sum(onames%in%ml$model_name),"/",length(onames))
add("Organoid fitness IDs present in model-list sample_id=",sum(oids%in%ml$sample_id),"/",length(oids))
short_names <- sub("^WTSI-","",onames)
syn_hits <- vapply(short_names,function(x)sum(grepl(x,ml$synonyms,fixed=TRUE)),integer(1))
add("Organoid fitness names with model-list synonym hit=",sum(syn_hits==1),"/",length(syn_hits),"; ambiguous synonym hits=",sum(syn_hits>1))
mapped_names <- vapply(seq_along(short_names),function(i)if(syn_hits[i]==1)ml$model_name[which(grepl(short_names[i],ml$synonyms,fixed=TRUE))[1]] else NA_character_,character(1))
add("Synonym-mapped names matching CMS sample_ID=",sum(mapped_names%in%read.csv(file.path(audit,"organoid2026","data_CMS-subtypes_organoids.csv"),check.names=FALSE,stringsAsFactors=FALSE)$sample_ID,na.rm=TRUE),"/",length(mapped_names))
map_method <- rep("unmatched",length(oids)); map_row <- rep(NA_integer_,length(oids));
for(i in seq_along(oids)){a<-which(ml$model_id==oids[i]); if(length(a)==1){map_method[i]<-"direct_model_id";map_row[i]<-a} else if(syn_hits[i]==1){map_method[i]<-"unique_synonym";map_row[i]<-which(grepl(short_names[i],ml$synonyms,fixed=TRUE))[1]}}
maptab <- data.frame(fitness_name=onames,fitness_model_id=oids,mapping_method=map_method,model_list_model_id=ifelse(is.na(map_row),NA,ml$model_id[map_row]),model_list_model_name=ifelse(is.na(map_row),NA,ml$model_name[map_row]),sample_id=ifelse(is.na(map_row),NA,ml$sample_id[map_row]),tissue=ifelse(is.na(map_row),NA,ml$tissue[map_row]),cancer_type=ifelse(is.na(map_row),NA,ml$cancer_type[map_row]),cancer_type_detail=ifelse(is.na(map_row),NA,ml$cancer_type_detail[map_row]))
write.table(maptab,file.path(out,"organoid_model_crosswalk_qc.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
add("Crosswalk direct+unique synonym mapped=",sum(map_method!="unmatched"),"/",length(map_method),"; mapped colorectal carcinoma=",sum(maptab$cancer_type=="Colorectal Carcinoma",na.rm=TRUE))
for(f in c("data_CMS-subtypes_organoids.csv","data_CRIS-subtypes_organoids.csv","data_covariates_ALL-organoids.csv")){
 x<-read.csv(file.path(audit,"organoid2026",f),check.names=FALSE,stringsAsFactors=FALSE); add(f," rows=",nrow(x),"; columns=",paste(names(x),collapse=","),"; fitness IDs matching first column=",sum(oids%in%x[[1]]),"/",length(oids))
}
writeLines(report,file.path(log,"stage2_structure_qc.txt"))
write.table(data.frame(dataset=c("DepMap","Organoid"),raw_models=c(length(dep_ids),length(oids)),duplicate_ids=c(length(dup_dep),sum(duplicated(oids))),missing_ids=c(sum(is.na(dep_ids)|!nzchar(dep_ids)),sum(is.na(oids)|!nzchar(oids)))),file.path(out,"stage2_model_qc_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Stage 2 QC completed.",report),file.path(log,"stage2_qc_status.txt"))
