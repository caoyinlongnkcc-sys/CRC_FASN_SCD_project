# CRC FASN/SCD lipid-dependency project

This directory is a publication-oriented copy of the reproducibility materials for the colorectal-cancer lipid-dependency analysis. It contains frozen gene-set files, result tables, figures, data manifests, checksums, and R session-information files. The analysis scripts are not included in this public release. Original raw data remain in the parent project and are not duplicated here.

## Study resources

- DepMap Public 24Q2
- Project Score Sanger v2/Broad 21Q2 archive
- cBioPortal studies: `coad_tcga_gdc`, `read_tcga_gdc`, and `coad_cptac_2019`
- GEO accessions: `GSE39582` and `GSE18105`
- Gene-set sources: Reactome, KEGG, Gene Ontology, and HGNC

Exact file paths, versions, retrieval dates, download URLs, and recorded SHA256 values are provided in `data_manifest/`. Raw public files are intentionally not included in this release copy.

## Directory structure

- `data_manifest/`: source manifests, URLs, and the stage 0–1 analysis plan.
- `gene_set_manifest/`: frozen lipid-metabolism gene-set files.
- `results_tables/`: intermediate and final result tables used for reporting.
- `figures/`: manuscript figure PNG files.
- `session_info/`: R session-information outputs from the recorded analyses.
- `checksums/`: SHA256 checksums for files in this release package.

## Reproducibility notes

The public release provides the retained manifests, result tables, figures, checksums, and session information. The analysis scripts are not publicly released. No result values were entered manually into the manuscript. Any re-analysis should be compared with the retained result tables before changing a reported value.

## Software

The recorded analyses used R 4.6.2 and the package versions listed in the files under `session_info/`.

## Data availability

The source data are publicly available from the repositories listed in `data_manifest/data_manifest.tsv` and `data_manifest/source_urls.tsv`, subject to the terms of those repositories. This package does not redistribute restricted or unnecessarily large raw files.
