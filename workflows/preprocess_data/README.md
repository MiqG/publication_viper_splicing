# Preprocess downloaded data for the project

This workflow (`snakefile`) processes and cleans data in the `data/raw` directory to populate the `data/prep` directory with data ready for downstream analysis.

## Outline (`snakefile`)
1. Generate list of splicing factors
2. Preprocess various datasets:
    - CCLE: generates clean event PSI, gene expression TPM, and mutation data.
    - ENCORE KD and KO: generates clean event PSI and gene expression tables.
    - TCGA: generates clean event PSI and gene expression tables for each cancer cohort.
    - ENA hand-curated datasets: generates clean event PSI and gene expression tables.
    - DepMap Demeter2: cleans up DEMETER2 gene dependency scores.
    - STRINGDB: prepares protein protein interaction data.
    - articles: generates clean event PSI and gene expression tables for specific studies:
        - Nijhuis2020
        - Lu2021
        - CardosoMoreira2020
        - tumorigenesis (Danielsson2013)
        - Riaz2017

3. Compute delta PSI for perturbation experiments
    - ENCORE KD and KO
    - ENA hand-curated dataset
    
4. Benchmark imputation of event PSI with KNN algorithm in CardosoMoreira2020 for posterior computational splicing factor network inference

5. Compute summary statistics for each dataset

6. Scripts to make figures for
    - `scripts/figures_eda.R`: overview of splicing factors
    - `scripts/figures_benchmark_imputation.R`: evaluate effect of choice of k in KNN imputation for PSI values

