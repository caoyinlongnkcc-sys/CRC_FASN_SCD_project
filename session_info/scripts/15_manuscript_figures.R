suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})
root <- "."
audit <- file.path("..", "feasibility_audit")
fig <- file.path(root, "figures")
out <- file.path(root, "results_tables")
dir.create(fig, showWarnings = FALSE)
theme_set(theme_classic(base_size = 11))
pal <- c("CRC" = "#0072B2", "non-CRC" = "#D55E00")
save_png <- function(plot, path, width, height) {
  png(path, width = round(width * 300), height = round(height * 300), res = 300)
  print(plot)
  dev.off()
}

depfile <- file.path(audit, "depmap24q2", "CRISPRGeneEffect.csv")
dh <- names(fread(depfile, nrows = 0, check.names = FALSE))
anchors <- c("FASN", "SCD", "ACACA", "SREBF1")
idx <- which(sub(" \\(.*$", "", dh) %in% anchors)
dep <- fread(depfile, select = c(1, idx), check.names = FALSE)
names(dep)[1] <- "ModelID"
names(dep)[-1] <- sub(" \\(.*$", "", names(dep)[-1])
md <- fread(file.path(audit, "DepMap_24Q2_Model.csv"))
crc_ids <- md[OncotreePrimaryDisease == "Colorectal Adenocarcinoma", ModelID]
dep[, group := ifelse(ModelID %in% crc_ids, "CRC", "non-CRC")]
long <- melt(dep, id.vars = c("ModelID", "group"), measure.vars = anchors,
             variable.name = "gene", value.name = "gene_effect")
long[, group := factor(group, levels = c("CRC", "non-CRC"))]
long[, gene := factor(gene, levels = anchors)]
fwrite(long, file.path(out, "figure2_depmap_anchor_plot_data.tsv"), sep = "\t", na = "NA")
p1 <- ggplot(long[is.finite(gene_effect)], aes(group, gene_effect, fill = group)) +
  geom_violin(trim = FALSE, alpha = .35, color = NA) +
  geom_boxplot(width = .22, outlier.shape = NA, alpha = .8) +
  facet_wrap(~gene, nrow = 1) + scale_fill_manual(values = pal) +
  labs(x = NULL, y = "DepMap CRISPR gene effect", fill = NULL) +
  theme(legend.position = "top")
save_png(p1, file.path(fig, "Figure2_DepMap_anchor_distributions.png"), 10, 3.8)

g <- fread(file.path(out, "GSE39582_FASN_SCD_clinical_matrix.tsv"))
g[, stage := factor(stage, levels = c("1", "2", "3", "4"),
                    labels = c("I", "II", "III", "IV"))]
gl <- melt(g[!is.na(stage)], id.vars = c("sample_id", "stage"),
           measure.vars = c("FASN", "SCD"), variable.name = "gene",
           value.name = "expression")
fwrite(gl, file.path(out, "figure5_gse39582_stage_plot_data.tsv"), sep = "\t", na = "NA")
p2 <- ggplot(gl, aes(stage, expression, fill = stage)) +
  geom_boxplot(outlier.shape = NA, alpha = .75) +
  geom_jitter(width = .12, alpha = .18, size = .45) +
  facet_wrap(~gene, nrow = 1) + scale_fill_brewer(palette = "Blues") +
  labs(x = "TNM stage", y = "GSE39582 expression", fill = NULL) +
  theme(legend.position = "none")
save_png(p2, file.path(fig, "Figure5_GSE39582_stage_expression.png"), 7.2, 3.8)

a <- fread(file.path(out, "GSE39582_adjusted_os_cox_results.tsv"))
a[, hr := hr_per_sd]
a[, feature := factor(feature, levels = rev(c("FASN_z", "SCD_z", "axis_z")))]
p3 <- ggplot(a, aes(y = feature, x = hr)) +
  geom_vline(xintercept = 1, lty = 2, color = "grey50") +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y",
                height = .16, color = "#0072B2") +
  geom_point(size = 2.5, color = "#0072B2") + scale_x_log10() +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL) + theme_classic()
save_png(p3, file.path(fig, "Figure5_GSE39582_adjusted_OS_forest.png"), 6.2, 3.6)

manifest <- c("Figure2_DepMap_anchor_distributions.png",
              "Figure5_GSE39582_stage_expression.png",
              "Figure5_GSE39582_adjusted_OS_forest.png")
writeLines(manifest, file.path(out, "manuscript_figure_manifest.txt"))
capture.output(sessionInfo(), file = file.path(root, "logs", "manuscript_figures_sessionInfo.txt"))
cat("Manuscript figures generated\n")
