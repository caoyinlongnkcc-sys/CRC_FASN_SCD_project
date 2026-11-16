suppressPackageStartupMessages(library(data.table))
out <- "results_tables"
checks <- list()
add_check <- function(name, expected, observed) {
  checks[[length(checks) + 1L]] <<- data.table(check = name,
    expected = as.character(expected), observed = as.character(observed),
    passed = identical(as.character(expected), as.character(observed)))
}

dep <- fread(file.path(out, "depmap_crc_lipid_dependency_summary.tsv"))
add_check("DepMap_total_core_genes", 1482, dep$total_core_genes)
add_check("DepMap_selective_genes_FDR_lt_0_05", 77, dep$fdr_lt_0_05)
add_check("DepMap_CRC_models", 61, dep$median_crc_models)
add_check("DepMap_nonCRC_models", 1089, dep$median_noncrc_models)

ps <- fread(file.path(out, "project_score_sanger_audit_summary.tsv"))
add_check("ProjectScore_CRC_models", 29, ps$sanger_crc_models)
add_check("ProjectScore_nonCRC_models", 321, ps$sanger_noncrc_models)
add_check("ProjectScore_DepMap_overlap", 28, ps$depmap_crc_model_overlap)
add_check("ProjectScore_anchors_detected", 7, ps$anchors_detected)

g39582 <- fread(file.path(out, "GSE39582_validation_qc.tsv"))
add_check("GSE39582_tumor_samples", 566, g39582$tumor_samples)
add_check("GSE39582_stage_complete", 562, g39582$stage_complete)
add_check("GSE39582_OS_events", 191, g39582$os_events)

g18105 <- fread(file.path(out, "GSE18105_validation_qc.tsv"))
add_check("GSE18105_cancer_LCM_samples", 77, g18105$cancer_lcm_samples)

psr <- fread(file.path(out, "project_score_sanger_anchor_results.tsv"))
add_check("ProjectScore_FASN_BH", "0.000868559108956435", sprintf("%.15g", psr[gene == "FASN", BH]))
add_check("ProjectScore_SCD_BH", "8.45969256024827e-06", sprintf("%.15g", psr[gene == "SCD", BH]))

tc <- fread(file.path(out, "tcga_expression_cna_correlations.tsv"))
add_check("TCGA_COAD_paired_rows", 4, nrow(tc[dataset == "TCGA_COAD"]))
add_check("TCGA_READ_paired_rows", 4, nrow(tc[dataset == "TCGA_READ"]))
add_check("TCGA_COAD_n", 454, unique(tc[dataset == "TCGA_COAD", n]))
add_check("TCGA_READ_n", 167, unique(tc[dataset == "TCGA_READ", n]))

figs <- file.exists(file.path("figures", c("Figure3_ProjectScore_anchor_replication.png", "Figure4_TCGA_expression_CNA_correlations.png")))
add_check("Annotated_Figure3_and_Figure4_exist", "TRUE", all(figs))
fig_info <- file.info(file.path("figures", c("Figure3_ProjectScore_anchor_replication.png", "Figure4_TCGA_expression_CNA_correlations.png")))
add_check("Annotated_Figure3_and_Figure4_nonzero", "TRUE", all(fig_info$size > 0))

res <- rbindlist(checks)
fwrite(res, file.path(out, "final_numeric_consistency_checks.tsv"), sep = "\t")
writeLines(c(
  "# Final numeric consistency audit",
  "",
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("All checks passed: ", all(res$passed)),
  "",
  "The audit compares manuscript-critical counts and annotated-figure source values with the existing audited result tables."
), "logs/final_numeric_consistency_report.md")
if (!all(res$passed)) stop("Numeric consistency audit failed")
cat("All numeric consistency checks passed\n")
