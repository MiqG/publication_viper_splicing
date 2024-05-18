# Validation of SF activity estimation with VIPER and SF-exon networks

## Outline
1. Regulon inference
    - `aracne_regulons.smk`
    - `mlr_regulons.smk`
    - `experimentally_derived_regulons.smk`
    - `inferred_and_experimental_regulons.smk`
        - ARACNe MI + experimental MoR
        - MLR beta + experimental MoR
        - ARACNe MI + ARACNe MoR (negative) + MLR MoR (positive)
2. Regulon evaluation
    - `regulon_robustness.smk`: generate empirical networks considering top N edges
    - `regulon_thresholds.smk`: generate empirical netwroks considering PSI thresholds
    - `regulon_evaluation.smk`: wihtin and across recall scores for all types of networks
3. Regulon EDA
    - `regulon_eda.smk`: EDA of empirical networks
    
## Important remarks

Make sure to have all packages required to run the scripts.

To run `aracne_regulons.smk`, adapted ARACNe-AP was obtained from https://github.com/chaolinzhanglab/ARACNe-AP and `scripts/aracne/filter_arachne_bootstraps.pl` and `scripts/aracne/estimate_mor.R` were obtained from https://github.com/chaolinzhanglab/mras.

## Recommendations
Run the workflows using
```
# regulon inference
snakemake -s aracne_regulons.smk --cores 6
snakemake -s mlr_regulons.smk --cores 6
snakemake -s experimentally_derived_regulons.smk --cores 6
snakemake -s inferred_and_experimental_regulons.smk --cores 6
# regulon evaluation
snakemake -s regulon_robustness.smk --cores 6
snakemake -s regulon_thresholds.smk --cores 6
snakemake -s regulon_evaluation.smk --cores 6
# regulon EDA
snakemake -s regulon_eda.smk --cores 6
```
In case you want to run the workflows on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).
