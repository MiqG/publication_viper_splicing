# Validation of SF activity estimation with VIPER and SF-exon networks

This set of Snakemake workflows validate estimation of splicing factor activity using empirical splicing factor networks and VIPER for different datasets of increasing complexity. These workflows use raw, preprocessed files and splicing factor networks from `results/regulon_inference` to populate `results/sf_activity_activation`.

## Outline
1. Protein depletion: 
    - identification of protein-level perturbations on splicing factors from transcriptomic changes
    - `protein_depletion_Nijhuis2020.smk`: transcriptomics and proteomics of Indisulam treatment, which targets splicing factor RBM39.
    - `protein_depletion_Lu2021.smk`: transcriptomics of indisulam treatment.

2. `combinatorial_perturbations.smk`: identification of combinatorially perturbed splicing factors

3. `sf3b_complex.smk`: detection of SF3b complex protein-protein interactions

## Expected outputs
```{shell}
$ tree results/sf_activity_validation/
results/sf_activity_validation/
├── figures
│   ├── protein_depletion-Lu2021-EX
│   │   ├── activity-indisulam_vs_dmso-box-ENSG00000131051.pdf
│   │   ├── activity-indisulam_vs_dmso-ranking-scatter-ENSG00000131051.pdf
│   │   ├── activity-MS023_vs_dmso-box-ENSG00000131051.pdf
│   │   ├── activity-MS023_vs_dmso-ranking-scatter-ENSG00000131051.pdf
│   │   ├── figdata
│   │   │   └── validation_drug_target_activity
│   │   │       ├── genexpr.tsv.gz
│   │   │       ├── metadata.tsv.gz
│   │   │       └── protein_activity.tsv.gz
│   │   ├── genexpr-indisulam_vs_dmso-box-ENSG00000131051.pdf
│   │   └── genexpr-MS023_vs_dmso-box-ENSG00000131051.pdf
│   ├── protein_depletion-Nijhuis2020-EX
│   │   ├── activity-indisulam_vs_dmso-ranking-scatter-ENSG00000131051.pdf
│   │   ├── diff_proteomics-indisulam_vs_dmso-volcano-RBM39.pdf
│   │   ├── figdata
│   │   │   └── validation_drug_target_activity
│   │   │       ├── diff_proteomics.tsv.gz
│   │   │       ├── genexpr.tsv.gz
│   │   │       ├── protein_activity.tsv.gz
│   │   │       └── proteomics.tsv.gz
│   │   ├── genexpr-indisulam_vs_dmso-bar-ENSG00000131051.pdf
│   │   ├── genexpr-indisulam_vs_dmso-box-ENSG00000131051.pdf
│   │   └── proteomics-indisulam_vs_dmso-box-RBM39.pdf
│   ├── validation_combinatorial_perturbations
│   │   ├── activity-double_perturbation_combined-ranking-scatter.pdf
│   │   ├── activity-double_perturbation_rep-ranking-scatter.pdf
│   │   └── figdata
│   │       └── validation_combinatorial_perturbation
│   │           └── protein_activity.tsv.gz
│   ├── validation_drug_targets-EX
│   │   ├── activity-indisulam_vs_dmso-ranking-scatter-ENSG00000131051.pdf
│   │   ├── figdata
│   │   │   └── validation_drug_target_activity
│   │   │       ├── genexpr.tsv.gz
│   │   │       ├── protein_activity.tsv.gz
│   │   │       └── proteomics.tsv.gz
│   │   ├── genexpr-indisulam_vs_dmso-box-ENSG00000131051.pdf
│   │   └── proteomics-indisulam_vs_dmso-box-RBM39.pdf
│   ├── validation_ptms
│   │   ├── activity_acetylation-double_perturbation_combined-ranking-scatter.pdf
│   │   ├── activity_acetylation-double_perturbation_rep-ranking-scatter.pdf
│   │   ├── activity_methylation-prmt_inhibitors-scatter_line.pdf
│   │   ├── activity_phosphorylation-sr_proteins-activity_vs_phosphoproteomics-scatter.pdf
│   │   ├── activity_phosphorylation-sr_proteins-KH-CB19-heatmap.pdf
│   │   ├── activity_phosphorylation-sr_proteins-KH-CB19-scatter_line.pdf
│   │   ├── activity_phosphorylation-sr_proteins-PALBOCICLIB-heatmap.pdf
│   │   ├── activity_phosphorylation-sr_proteins-PALBOCICLIB-scatter_line.pdf
│   │   ├── activity_phosphorylation-sr_proteins-T025-heatmap.pdf
│   │   ├── activity_phosphorylation-sr_proteins-T025-scatter_line.pdf
│   │   ├── activity_phosphorylation-sr_proteins-T3-heatmap.pdf
│   │   └── activity_phosphorylation-sr_proteins-T3-scatter_line.pdf
│   └── validation_sf3b_complex
│       ├── activity_drugs-sf3b_complex-scatter_line.pdf
│       ├── activity_drugs-sf3b_complex_vs_shortest_paths-box.pdf
│       ├── activity_drugs-sf3b_complex_vs_shortest_paths_vs_u2snrnp-box.pdf
│       ├── activity_drugs-shortest_paths-bar.pdf
│       ├── activity_mutations-sf3b_complex-scatter_line.pdf
│       └── figdata
│           └── validation_sf3b_complex
│               ├── protein_activity.tsv.gz
│               └── shortest_paths.tsv.gz
└── files
    ├── metadata
    │   ├── ptms-EX.tsv.gz
    │   └── sf3b_complex-EX.tsv.gz
    ├── ppi
    │   ├── shortest_path_lengths_to_sf3b_complex.tsv.gz
    │   └── shortest_path_lengths_to_u2snrnp_complex.tsv.gz
    ├── protein_activity
    │   ├── combinatorial_perturbations-EX.tsv.gz
    │   ├── protein_depletion-EX.tsv.gz
    │   ├── protein_depletion-Lu2021-EX.tsv.gz
    │   ├── protein_depletion-Nijhuis2020-EX.tsv.gz
    │   ├── ptms-EX.tsv.gz
    │   └── sf3b_complex-EX.tsv.gz
    └── signatures
        ├── combinatorial_perturbations-EX.tsv.gz
        ├── protein_depletion-EX.tsv.gz
        ├── protein_depletion-Lu2021-EX.tsv.gz
        ├── protein_depletion-Nijhuis2020-EX.tsv.gz
        ├── ptms-EX.tsv.gz
        └── sf3b_complex-EX.tsv.gz
```