#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- "."; raw <- file.path(root,"data_raw"); meta <- file.path(root,"data_metadata"); dir.create(meta,recursive=TRUE,showWarnings=FALSE)
stamp <- "2026-11-11"
rt <- function(f,...) read.delim(f,quote="",comment.char="",check.names=FALSE,stringsAsFactors=FALSE,...)
h <- rt(file.path(raw,"hgnc_complete_set.txt"),sep="\t"); approved <- unique(h$symbol[h$status=="Approved"])
emap <- list(); umap <- list()
for(i in seq_len(nrow(h))){s<-h$symbol[i]; eid<-as.character(h$entrez_id[i]); uid<-as.character(h$uniprot_ids[i]); if(!is.na(eid)&&nzchar(eid))for(x in strsplit(eid,"[|, ]")[[1]])if(nzchar(x))emap[[x]]<-s; if(!is.na(uid)&&nzchar(uid))for(x in strsplit(uid,"[|, ]")[[1]])if(nzchar(x))umap[[x]]<-s}
rows<-list(); rejected<-list(); add<-function(g,src,id,nm,ev=""){g<-as.character(g); ok<-!is.na(g)&g%in%approved; if(any(ok))for(x in unique(g[ok]))rows[[length(rows)+1L]]<<-data.frame(gene=x,source=src,source_id=id,source_name=nm,evidence=ev,retrieved=stamp); if(any(!ok))rejected[[length(rejected)+1L]]<<-data.frame(source=src,source_id=id,candidate=g[!ok],reason="not_HGNC_approved")}
# Reactome human lipid scope was fixed before dependency review.
gmt<-file.path(raw,"reactome_gmt_tmp","ReactomePathways.gmt"); if(!file.exists(gmt))stop("Missing extracted Reactome GMT")
rx<-"(?i)(fatty acid|glycerolipid|glycerophospholipid|sphingolipid|cholesterol|steroid|eicosanoid|lipoprotein|triglyceride|phospholipid|ether lipid|bile acid|lipid metabolism|lipid homeostasis)"
for(ln in readLines(gmt,warn=FALSE)){p<-strsplit(ln,"\t",fixed=TRUE)[[1]]; if(length(p)>=3&&grepl("^R-HSA-",p[2])&&grepl(rx,p[1]))add(p[-c(1,2)],"Reactome",p[2],p[1],"curated_pathway")}
# KEGG IDs are predeclared in protocol/config.yml; empty responses are retained as QC findings.
kids<-c("hsa00061","hsa00062","hsa00071","hsa00072","hsa00100","hsa00120","hsa00121","hsa00140","hsa00561","hsa00564","hsa00565","hsa00590","hsa00591","hsa00592","hsa00600","hsa01040","hsa01212","hsa04979")
kn<-setNames(rep(NA_character_,length(kids)),kids); kp<-rt(file.path(raw,"kegg_hsa_pathways.tsv"),header=FALSE,sep="\t"); for(i in seq_len(nrow(kp)))kn[as.character(kp[i,1])]<-sub(" - Homo sapiens \\(human\\)$","",as.character(kp[i,2]))
for(id in kids){f<-file.path(raw,paste0("kegg_",id,"_genes.tsv")); if(!file.exists(f)||file.info(f)$size<=1)next; z<-rt(f,header=FALSE,sep="\t"); gn<-unname(unlist(emap[sub("^hsa:","",as.character(z[,2]))])); add(gn,"KEGG",id,kn[id],"curated_pathway")}
# GO root plus descendants in go-basic, excluding ND annotations.
obo<-readLines(file.path(raw,"go-basic.obo"),warn=FALSE); par<-list(); active<-FALSE; cur<-NA_character_
for(ln in obo){if(ln=="[Term]"){active<-TRUE;cur<-NA_character_;next}; if(active&&grepl("^id: GO:",ln)){cur<-sub("^id: ","",ln);par[[cur]]<-character();next}; if(active&&!is.na(cur)&&grepl("^is_a: GO:",ln))par[[cur]]<-c(par[[cur]],sub("^is_a: (GO:[0-9]+).*$","\\1",ln)); if(ln=="")active<-FALSE}
desc<-"GO:0006629"; repeat{new<-names(par)[vapply(par,function(x)any(x%in%desc),logical(1))]; new<-setdiff(new,desc); if(!length(new))break; desc<-c(desc,new)}
gz<-gzfile(file.path(raw,"HUMAN-uniprot.gaf.gz"),"rt"); gaf<-readLines(gz,warn=FALSE);close(gz)
for(ln in gaf[!grepl("^!",gaf)]){p<-strsplit(ln,"\t",fixed=TRUE)[[1]]; if(length(p)<7||!(p[5]%in%desc)||p[7]=="ND")next; g<-umap[[sub("^UniProtKB:","",p[2])]]; if(is.null(g))next; add(g,"GO",p[5],"GO:0006629 descendant",ifelse(p[7]=="IEA","IEA","curated_or_experimental"))}
long<-if(length(rows))unique(do.call(rbind,rows))else data.frame(); write.table(long,file.path(meta,"lipid_gene_set_long.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(length(rejected))write.table(unique(do.call(rbind,rejected)),file.path(meta,"lipid_gene_set_rejected.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(nrow(long)){sp<-split(long,long$gene); sm<-do.call(rbind,lapply(sp,function(d){s<-unique(d$source);data.frame(gene=d$gene[1],sources=paste(sort(s),collapse=";"),n_sources=length(s),reactome=as.integer("Reactome"%in%s),kegg=as.integer("KEGG"%in%s),go=as.integer("GO"%in%s),go_curated=as.integer(any(d$source=="GO"&d$evidence!="IEA")),go_iea=as.integer(any(d$source=="GO"&d$evidence=="IEA")),tier=ifelse(any(s%in%c("Reactome","KEGG"))||any(d$source=="GO"&d$evidence!="IEA"),"core","extended"))})); ce<-function(f)sub(" \\(.*$","",readLines(f)[-1]); sm$depmap_common_essential<-as.integer(sm$gene%in%ce(file.path(raw,"DepMap_24Q2_CRISPRInferredCommonEssentials.csv"))); sm$sanger_common_essential<-as.integer(sm$gene%in%ce(file.path(raw,"Sanger_Organoid_CRISPRInferredCommonEssentials.csv"))); write.table(sm,file.path(meta,"lipid_gene_set.tsv"),sep="\t",quote=FALSE,row.names=FALSE); writeLines(c(paste("retrieved",stamp),paste("genes",nrow(sm)),paste("core",sum(sm$tier=="core")),paste("extended",sum(sm$tier=="extended")),paste("GO_descendants",length(desc)),paste("DepMap_common_essential",sum(sm$depmap_common_essential)),paste("Sanger_common_essential",sum(sm$sanger_common_essential))),file.path(root,"logs","gene_set_qc.txt"))}
write.table(data.frame(resource=c("Reactome GMT","Reactome pathway names","KEGG REST","GO basic + human GAF","HGNC complete set"),retrieved=stamp),file.path(meta,"gene_set_sources.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("KEGG IDs predeclared:",kids,"","Empty response files:",kids[file.info(file.path(raw,paste0("kegg_",kids,"_genes.tsv")))$size<=1]),file.path(meta,"kegg_empty_pathways.txt"))
