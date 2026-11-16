suppressPackageStartupMessages(library(data.table))
root <- "."; audit <- file.path("..","feasibility_audit"); raw <- file.path(root,"data_raw"); out <- file.path(root,"results_tables"); logd <- file.path(root,"logs")
depfile <- file.path(audit,"depmap24q2","CRISPRGeneEffect.csv"); exprfile <- file.path(audit,"depmap24q2","OmicsExpressionProteinCodingGenesTPMLogp1.csv"); modfile <- file.path(audit,"DepMap_24Q2_Model.csv")
anchors <- c("FASN","SCD")
dh <- names(fread(depfile,nrows=0,check.names=FALSE)); didx <- c(1,which(sub(" \\(.*$","",dh)%in%anchors)); dep <- fread(depfile,select=didx,check.names=FALSE); names(dep)[1]<-"ModelID"; names(dep)[-1]<-anchors
md <- fread(modfile); crc_ids <- intersect(md[OncotreePrimaryDisease=="Colorectal Adenocarcinoma",ModelID],dep$ModelID); dep <- dep[ModelID%in%crc_ids]; if(nrow(dep)!=61)stop("CRC dependency count changed")
expr <- fread(exprfile,check.names=FALSE); names(expr)[1]<-"ModelID"; expr<-expr[ModelID%in%crc_ids]; names(expr)<-make.unique(sub(" \\(.*$","",names(expr))); m<-match(dep$ModelID,expr$ModelID); keep<-!is.na(m); dep<-dep[keep]; expr<-expr[m[keep]]; if(nrow(dep)<50)stop("Too few matched CRC models")
X<-as.matrix(expr[,-1,with=FALSE]); storage.mode(X)<-"numeric"; genes<-colnames(X)
gmt <- file.path(raw,"reactome_gmt_tmp","ReactomePathways.gmt"); if(!file.exists(gmt)){td<-tempfile();dir.create(td);unzip(file.path(raw,"ReactomePathways.gmt.zip"),exdir=td);gmt<-file.path(td,"ReactomePathways.gmt")}
lines<-readLines(gmt,warn=FALSE); parts<-lapply(lines,function(z)strsplit(z,"\t",fixed=TRUE)[[1]]); names(parts)<-vapply(parts,`[[`,character(1),1)
targets<-c("Fatty acid metabolism","Regulation of cholesterol biosynthesis by SREBP (SREBF)","Unfolded Protein Response (UPR)","ATF4 activates genes in response to endoplasmic reticulum  stress","Mitochondrial Fatty Acid Beta-Oxidation")
sets<-lapply(targets,function(n){q<-parts[[n]];if(is.null(q))character()else intersect(q[-c(1,2)],genes)});names(sets)<-targets
score<-function(g){z<-X[,g,drop=FALSE];z<-apply(z,2,function(v){as.numeric(scale(v))});rowMeans(z,na.rm=TRUE)}
res<-list();cov<-list();k<-1L
for(nm in names(sets)){g<-sets[[nm]];cov[[k]]<-data.table(pathway=nm,n_genes_requested=if(is.null(parts[[nm]]))0L else length(parts[[nm]])-2,n_genes_expression=length(g),n_models=nrow(dep));if(length(g)>=5){sc<-score(g);for(a in anchors){ok<-is.finite(sc)&is.finite(dep[[a]]);ct<-suppressWarnings(cor.test(sc[ok],dep[[a]][ok],method="spearman",exact=FALSE));res[[length(res)+1L]]<-data.table(pathway=nm,gene=a,n_genes_expression=length(g),n_models=sum(ok),spearman_rho=unname(ct$estimate),p_value=ct$p.value)}};k<-k+1L}
r<-rbindlist(res,fill=TRUE);r[,p_adj:=p.adjust(p_value,"BH")];r[,spearman_rho_dependency_strength:=-spearman_rho];fwrite(r,file.path(out,"depmap_crc_targeted_pathway_dependency_correlations.tsv"),sep="\t",na="NA");fwrite(rbindlist(cov,fill=TRUE),file.path(out,"depmap_crc_targeted_pathway_coverage.tsv"),sep="\t",na="NA")
limitation<-data.table(resource="KEGG_ferroptosis_hsa04216",available=FALSE,note="KEGG REST SSL request failed in this environment; no ferroptosis gene set was invented or analyzed")
fwrite(limitation,file.path(out,"targeted_pathway_data_limitations.tsv"),sep="\t",na="NA");capture.output(sessionInfo(),file=file.path(logd,"targeted_pathway_audit_sessionInfo.txt"));cat("Targeted pathway audit completed; matched models:",nrow(dep),"\n")
