# Preregistered analysis plan (Stage 0–1)

## Question

Identify lipid-metabolism genes with reproducible functional dependency in colorectal-cancer models, then test whether dependencies are conserved across DepMap cell lines and Sanger organoids or are associated with CMS/MSI/lineage context.

The previously considered tegavivint-sensitivity analysis is excluded because a verifiable cell-line-level AUC endpoint with stable model identifiers was not available.

## Gene-set freeze

Membership is created by `scripts/01_build_lipid_gene_set.R` from the dated files listed in `data_metadata/gene_set_sources.tsv`:

- Reactome human GMT: retain `R-HSA` pathways whose names match the prespecified lipid scope in the script.
- KEGG: use the 18 pathway IDs listed in `protocol/config.yml`; empty API responses remain documented and are not replaced post hoc.
- GO: `GO:0006629` and descendants in `go-basic.obo`; human GAF annotations with evidence code `ND` are excluded. IEA and non-IEA annotations are retained and labelled.
- HGNC: current approved-symbol table is used for identifier validation.

The primary screen uses `tier=core` (Reactome/KEGG or non-IEA GO support). `tier=extended` is a sensitivity set containing genes supported only by IEA GO annotations. Dependency values are not read by the gene-set script.

## Primary statistical rules (frozen before new result review)

1. Use model-level dependency values, not a gene score fabricated from the outcome.
2. Define CRC using the model metadata tissue lineage; report colon and rectum separately if sample size permits.
3. Summarise gene dependency by median and proportion below the prespecified -0.5 threshold; do not use the threshold as the only result.
4. Compare CRC against other lineages using a two-sided Wilcoxon test and BH correction across the frozen gene set.
5. For organoids, use the published fold-change/fitness endpoint and repeat the same directionally defined summaries. No cross-platform numeric pooling is allowed.
6. Context analyses (CMS, CRIS/MSI and metastatic status) require at least 15 CRC models per tested subgroup; otherwise they are descriptive only.
7. Essential genes are marked using the downloaded DepMap and Sanger `CRISPRInferredCommonEssentials.csv` files; essentiality is not used to select genes.
8. Any missingness, duplicate model identifier, lineage mismatch, or failed identifier mapping is reported and triggers a QC stop rather than silent imputation.

## Required validation before manuscript claims

- Model IDs and lineage labels are one-to-one after joins.
- At least one independent platform supports any gene described as “reproducible”.
- BH-adjusted results, effect sizes and denominators are reported together.
- All analyses can be rerun from the manifest and scripts without editing the frozen gene list.

