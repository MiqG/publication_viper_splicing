## source code directory - project-wide scripts
This directory contains scripts used for data processing, statistical analysis, and gene/protein activity computation. Each script performs a specific function within the analysis pipeline. They are here because they may be used in multiple workflows.

### outline
- `compute_gene_sets_enrichment_score.R`: Calculates enrichment scores for predefined gene sets. This is used to assess the overrepresentation of specific gene sets in the data.

- `compute_protein_activity.R`: computes protein activity using different approaches (GSEA, Pearson correlation, Spearman correlation, VIPER) levels based on exon inclusion signatures.

- `evaluate_activity.R`: evaluates computed protein activity scores.

- `run_deseq2.R`: perfoms differential expression analysis using the DESeq2 R package. This script processes count data to identify significantly differentially expressed genes across conditions.

- `MannWhitneyU.py`: implements the Mann-Whitney U test, a nonparametric test for assessing whether two independent samples come from the same distribution. Helper scripts of this script:
    - `PairwiseTest.py`
    - `process_inputs.py`
    - `util.py`
