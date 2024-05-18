# Preprocess downloaded data for the project

## Outline (`snakefile`)
1. List of splicing factors
2. Preprocess datasets:
    - CCLE
    - ENCORE KD and KO
    - TCGA
    - ENA hand-curated datasets
    - DepMap Demeter2
    - STRINGDB
    - articles
        - Nijhuis2020
        - Lu2021
        - CardosoMoreira2020
        - tumorigenesis (Danielsson2013)
        - Riaz2017

3. Compute delta PSI from perturbation experiments
    - ENCORE KD and KO
    - ENA hand-curated dataset
    
4. Impute event PSI in CardosoMoreira2020 for SF-exon network inference with ARACNe and MLR

5. Compute summary stats

6. EDA list of splicing factors

## Important remarks

Make sure to have all packages required to run the scripts.

## Recommendations
Run the workflows using
```
snakemake --cores 6
```
In case you want to run the workflows on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).