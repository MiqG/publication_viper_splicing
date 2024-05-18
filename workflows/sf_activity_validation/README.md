# Validation of SF activity estimation with VIPER and SF-exon networks

## Outline
1. Protein depletion: identification of protein-level perturbations on splicing factors from transcriptomic changes
    - `protein_depletion_Nijhuis2020.smk`: transcriptomics and proteomics of Indisulam treatment, which targets splicing factor RBM39.
    - `protein_depletion_Lu2021.smk`: transcriptomics of indisulam treatment.
2. `combinatorial_perturbations.smk`: identification of combinatorially perturbed splicing factors
3. `sf3b_complex.smk`: detection of SF3b complex protein-protein interactions

## Important remarks

Make sure to have all packages required to run the scripts.

## Recommendations
Run the workflows using
```
snakemake -s protein_depletion_Nijhuis2020.smk --cores 6
snakemake -s protein_depletion_Lu2021.smk --cores 6
snakemake -s combinatorial_perturbations.smk --cores 6
snakemake -s sf3b_complex.smk --cores 6
```
In case you want to run the workflows on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).