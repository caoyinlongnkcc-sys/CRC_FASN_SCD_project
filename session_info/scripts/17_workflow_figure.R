suppressPackageStartupMessages(library(ggplot2))
workflow <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = 1,
  label = c("Frozen lipid\ngene universe\n1,539 genes", "DepMap 24Q2\nCRC discovery\n61 models", "Project Score\nplatform check\n29 CRC models", "TCGA/CPTAC\npatient multi-omics", "GEO cohorts\nGSE39582 / GSE18105")
)
p <- ggplot(workflow, aes(x, y)) +
  geom_segment(data = data.frame(x = 1:4, xend = 2:5, y = 1, yend = 1),
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.18, "cm")), color = "grey45") +
  geom_label(aes(label = label), size = 3.2, linewidth = 0.3,
             fill = "#EAF2F8", color = "#1B4F72", lineheight = 0.95) +
  scale_x_continuous(limits = c(.5, 5.5), breaks = NULL) +
  scale_y_continuous(limits = c(.7, 1.3), breaks = NULL) +
  labs(title = "Evidence hierarchy for CRC lipid-dependency analysis") +
  theme_void() + theme(plot.title = element_text(hjust = .5, face = "bold"))
png("figures/Figure1_study_workflow.png", width = 3600, height = 1200, res = 300)
print(p)
dev.off()
cat("Workflow figure generated\n")
