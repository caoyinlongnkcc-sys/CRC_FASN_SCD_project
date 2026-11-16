suppressPackageStartupMessages({library(data.table); library(ggplot2)})
fig <- "figures"
dir.create(fig, showWarnings = FALSE)
save_png <- function(plot, path, width, height) {
  png(path, width = round(width * 300), height = round(height * 300), res = 300)
  print(plot)
  dev.off()
}

ps <- fread("results_tables/project_score_sanger_anchor_results.tsv")
ps[, gene := factor(gene, levels = rev(c("FASN", "SCD", "ACACA", "SREBF1")))]
ps_plot <- ps[gene %in% c("FASN", "SCD", "ACACA", "SREBF1")]
ps_plot[, label := sprintf("n=%d/%d\\nΔ=%.3f\\nBH=%.3g", n_crc, n_noncrc, median_difference, BH)]
p3 <- ggplot(ps_plot,
             aes(x = median_difference, y = gene, color = direction_concordant)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  geom_point(size = 3) +
  geom_text(aes(label = label), hjust = 1.05, nudge_x = -0.03,
            color = "black", size = 2.5, lineheight = 0.85) +
  scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#999999")) +
  labs(x = "Project Score CRC minus non-CRC median gene effect",
       y = NULL, color = "Direction concordant") + theme_classic()
save_png(p3, file.path(fig, "Figure3_ProjectScore_anchor_replication.png"), 7, 3.8)

tc <- fread("results_tables/tcga_expression_cna_correlations.tsv")
tc[, cohort := fifelse(dataset == "TCGA_COAD", "COAD", "READ")]
tc[, gene := factor(gene, levels = rev(c("FASN", "SCD", "ACACA", "SREBF1")))]
tc[, label := sprintf("n=%d\\nBH=%.3g", n, p_adj)]
p4 <- ggplot(tc, aes(x = spearman_rho, y = gene, color = cohort)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  geom_point(size = 2.8, position = position_dodge(width = .35)) +
  geom_text(aes(label = label), position = position_dodge(width = .35),
            hjust = -0.05, color = "black", size = 2.5, lineheight = 0.85) +
  scale_color_manual(values = c("COAD" = "#0072B2", "READ" = "#D55E00")) +
  labs(x = "Spearman correlation: expression versus CNA", y = NULL, color = NULL) +
  theme_classic()
save_png(p4, file.path(fig, "Figure4_TCGA_expression_CNA_correlations.png"), 7, 3.8)
cat("Additional figures generated\n")
