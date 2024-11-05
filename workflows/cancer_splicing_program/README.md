# Identifying a recurrent cancer splicing program

This set of Snakemake workflows explore how splicing factor activities change across cancer cohorts systematically. These workflows use raw, preprocessed files and splicing factor networks from `results/regulon_inference` to populate `results/cancer_splicing_program`.

## Outline
1. `program_definition.smk`: define cancer program and survival analysis.
2. `immune_evasion.smk`: study how program recapitulates immune evasion.
3. `proliferation.smk`: study how program recapitulates cell proliferation.
4. `tumorigenesis.smk`: study how program recapitulates tumorigenesis.        
        
## Expected outputs
```{shell}
$ tree results/cancer_splicing_program/
results/cancer_splicing_program/
├── figures
│   ├── cancer_program
│   │   ├── comparison-correlation_by_cancer-violin.pdf
│   │   ├── comparison-correlation_median_pancan-violin.pdf
│   │   ├── comparison-diff_analysis-scatter.pdf
│   │   ├── comparison-survival-scatter.pdf
│   │   ├── driver_selection-drivers_vs_cancer_type-activity_drivers_vs_activity-violin.pdf
│   │   ├── driver_selection-drivers_vs_cancer_type-activity_drivers_vs_genexpr-violin.pdf
│   │   ├── driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_activity-violin.pdf
│   │   ├── driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_genexpr-violin.pdf
│   │   ├── driver_selection-n_signif_vs_driver_type-activity-bar.pdf
│   │   ├── driver_selection-n_signif_vs_driver_type_by_genexpr-genexpr-bar.pdf
│   │   ├── driver_selection-n_signif_vs_driver_type-genexpr-bar.pdf
│   │   ├── driver_selection-sf_vs_cancer_type-heatmap.pdf
│   │   ├── eda_programs-rbpdb_family-bar.pdf
│   │   ├── eda_programs-spliceosome_db_complex-bar.pdf
│   │   ├── figdata
│   │   │   ├── cancer_program
│   │   │   │   └── enrichments_reactome.tsv.gz
│   │   │   └── tcga_tumorigenesis
│   │   │       ├── assocs_gene_dependency.tsv.gz
│   │   │       ├── diff_protein_activity.tsv.gz
│   │   │       ├── sf_cross_regulation.tsv.gz
│   │   │       └── survival_roc_activity.tsv.gz
│   │   ├── n_samples-balloon.pdf
│   │   ├── prolif_driver-driver_type_vs_demeter2-box.pdf
│   │   ├── prolif_driver-driver_type_vs_demeter2-essential-box.pdf
│   │   ├── prolif_driver-driver_type_vs_demeter2-not_essential-box.pdf
│   │   ├── reactome_enrichments-bar.pdf
│   │   ├── reactome_enrichments-immune_screen-scatter.pdf
│   │   ├── sf_cross_regulation-correlations-violin-activity.pdf
│   │   ├── sf_cross_regulation-correlations-violin-genexpr.pdf
│   │   ├── sf_cross_regulation-regulon_similarity-violin-activity.pdf
│   │   ├── sf_cross_regulation-regulon_similarity-violin-genexpr.pdf
│   │   ├── survival_analysis-cancer_all-examples-bar-activity.pdf
│   │   ├── survival_analysis-cancer_all-examples-bar-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis-cancer_all-examples-bar-genexpr.pdf
│   │   ├── survival_analysis-cancer_all-examples-bar-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis-cancers_all-roc_curves-activity.pdf
│   │   ├── survival_analysis-cancers_all-roc_curves-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis-cancers_all-roc_curves-genexpr.pdf
│   │   ├── survival_analysis-cancers_all-roc_curves-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis-cancers_all-violin-activity.pdf
│   │   ├── survival_analysis-cancers_all-violin-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis-cancers_all-violin-genexpr.pdf
│   │   ├── survival_analysis-cancers_all-violin-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis-cancers_all_vs_thresh-violin-activity.pdf
│   │   ├── survival_analysis-cancers_all_vs_thresh-violin-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis-cancers_all_vs_thresh-violin-genexpr.pdf
│   │   ├── survival_analysis-cancers_all_vs_thresh-violin-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis-cancers_differential-roc_curves-activity.pdf
│   │   ├── survival_analysis-cancers_differential-roc_curves-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis-cancers_differential-roc_curves-genexpr.pdf
│   │   ├── survival_analysis-cancers_differential-roc_curves-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis-cancers_differential-violin-activity.pdf
│   │   ├── survival_analysis-cancers_differential-violin-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis-cancers_differential-violin-genexpr.pdf
│   │   ├── survival_analysis-cancers_differential-violin-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis_conf-cancer_all-examples-bar-activity.pdf
│   │   ├── survival_analysis_conf-cancer_all-examples-bar-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis_conf-cancer_all-examples-bar-genexpr.pdf
│   │   ├── survival_analysis_conf-cancer_all-examples-bar-genexpr_w_activity_labs.pdf
│   │   ├── survival_analysis_conf-cancers_all-violin-activity.pdf
│   │   ├── survival_analysis_conf-cancers_all-violin-activity_w_genexpr_labs.pdf
│   │   ├── survival_analysis_conf-cancers_all-violin-genexpr.pdf
│   │   └── survival_analysis_conf-cancers_all-violin-genexpr_w_activity_labs.pdf
│   ├── immune_evasion-EX
│   │   ├── diff_response-immune_screen_score_reglators_oi-bar.pdf
│   │   ├── diff_response-median_diff_vs_immune_screen_score-scatter.pdf
│   │   ├── diff_response-median_diff_vs_pvalue-scatter.pdf
│   │   ├── diff_response-psi_vs_survival-km-HsaEX1036341_SEC22B.pdf
│   │   ├── figdata
│   │   │   └── tumorigenesis
│   │   │       ├── diff_response.tsv.gz
│   │   │       ├── immune_screen.tsv.gz
│   │   │       └── splicing.tsv.gz
│   │   └── reactome_enrichments-targets_regulators_exon_oi-bar.pdf
│   ├── proliferation-EX
│   │   ├── figdata
│   │   │   └── proliferation
│   │   │       └── protein_activity_ccle.tsv.gz
│   │   ├── hallmarks_ccle-mki67_vs_activity_median_correlations-violin.pdf
│   │   └── hallmarks_ccle-mki67_vs_activity_median-scatter.pdf
│   └── tumorigenesis-EX
│       ├── figdata
│       │   └── tumorigenesis
│       │       └── protein_activity.tsv.gz
│       ├── tumorigenesis-cell_line_vs_activity-random-violin.pdf
│       ├── tumorigenesis-cell_line_vs_activity-violin.pdf
│       └── tumorigenesis-cell_line_vs_genexpr_fc-violin.pdf
└── files
    ├── diff_genexpr_counts_deseq2
    │   ├── BLCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── BRCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── COAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── HNSC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KICH-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KIRC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KIRP-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LIHC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LUAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LUSC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── PRAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── THCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   └── UCEC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    ├── diff_genexpr_tpm
    │   ├── BLCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── BRCA-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── BRCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── COAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── HNSC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KICH-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KIRC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KIRP-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LIHC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LUAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LUSC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── PRAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── SKCM-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── STAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── THCA-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── THCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   └── UCEC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    ├── diff_protein_activity
    │   ├── BLCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── BRCA-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── BRCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── COAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── HNSC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KICH-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KIRC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── KIRP-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LIHC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LUAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── LUSC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── PRAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── SKCM-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── STAD-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── THCA-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── THCA-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   └── UCEC-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    ├── metadata
    │   └── driver_mutations-EX.tsv.gz
    ├── PANCAN
    │   ├── cancer_program.tsv.gz
    │   ├── genexpr_counts_deseq2-deseq2-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── genexpr_tpm-mannwhitneyu-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── genexpr_tpm-sf_cross_regulation.tsv.gz
    │   ├── genexpr_tpm-survival_analysis-cat.tsv.gz
    │   ├── genexpr_tpm-survival_analysis-surv.tsv.gz
    │   ├── genexpr_tpm-survival_analysis_with_confounders-cat.tsv.gz
    │   ├── genexpr_tpm-survival_analysis_with_confounders-surv.tsv.gz
    │   ├── genexpr_tpm_vs_activity.tsv.gz
    │   ├── mannwhitneyu-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── protein_activity-mannwhitneyu-Metastatic_vs_PrimaryTumor.tsv.gz
    │   ├── protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz
    │   ├── protein_activity-sf_cross_regulation.tsv.gz
    │   ├── protein_activity-survival_analysis-cat.tsv.gz
    │   ├── protein_activity-survival_analysis-surv.tsv.gz
    │   ├── protein_activity-survival_analysis_with_confounders-cat.tsv.gz
    │   ├── protein_activity-survival_analysis_with_confounders-surv.tsv.gz
    │   ├── sf_cross_regulation.tsv.gz
    │   ├── survival_analysis-cat.tsv.gz
    │   ├── survival_analysis-surv.tsv.gz
    │   └── survival_analysis.tsv.gz
    ├── program_enrichment_scores
    │   ├── CCLE-EX.tsv.gz
    │   └── PANCAN-SolidTissueNormal-EX.tsv.gz
    ├── protein_activity
    │   ├── ACC-PrimaryTumor-EX.tsv.gz
    │   ├── Bian2018-EX.tsv.gz
    │   ├── BLCA-PrimaryTumor-EX.tsv.gz
    │   ├── BLCA-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── BLCA-SolidTissueNormal-EX.tsv.gz
    │   ├── BRCA-Metastatic-EX.tsv.gz
    │   ├── BRCA-Metastatic_vs_PrimaryTumor-EX.tsv.gz
    │   ├── BRCA-PrimaryTumor-EX.tsv.gz
    │   ├── BRCA-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── BRCA-SolidTissueNormal-EX.tsv.gz
    │   ├── CCLE-EX.tsv.gz
    │   ├── CESC-PrimaryTumor-EX.tsv.gz
    │   ├── CESC-SolidTissueNormal-EX.tsv.gz
    │   ├── CHOL-PrimaryTumor-EX.tsv.gz
    │   ├── CHOL-SolidTissueNormal-EX.tsv.gz
    │   ├── COAD-PrimaryTumor-EX.tsv.gz
    │   ├── COAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── COAD-SolidTissueNormal-EX.tsv.gz
    │   ├── DLBC-PrimaryTumor-EX.tsv.gz
    │   ├── driver_mutations-EX.tsv.gz
    │   ├── ESCA-PrimaryTumor-EX.tsv.gz
    │   ├── ESCA-SolidTissueNormal-EX.tsv.gz
    │   ├── GBM-PrimaryTumor-EX.tsv.gz
    │   ├── GBM-RecurrentTumor-EX.tsv.gz
    │   ├── GBM-SolidTissueNormal-EX.tsv.gz
    │   ├── HNSC-PrimaryTumor-EX.tsv.gz
    │   ├── HNSC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── HNSC-SolidTissueNormal-EX.tsv.gz
    │   ├── KICH-PrimaryTumor-EX.tsv.gz
    │   ├── KICH-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── KICH-SolidTissueNormal-EX.tsv.gz
    │   ├── KIRC-PrimaryTumor-EX.tsv.gz
    │   ├── KIRC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── KIRC-SolidTissueNormal-EX.tsv.gz
    │   ├── KIRP-PrimaryTumor-EX.tsv.gz
    │   ├── KIRP-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── KIRP-SolidTissueNormal-EX.tsv.gz
    │   ├── LAML-PrimaryBloodDerivedCancerPeripheralBlood-EX.tsv.gz
    │   ├── LGG-PrimaryTumor-EX.tsv.gz
    │   ├── LGG-RecurrentTumor-EX.tsv.gz
    │   ├── LIHC-PrimaryTumor-EX.tsv.gz
    │   ├── LIHC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── LIHC-RecurrentTumor-EX.tsv.gz
    │   ├── LIHC-SolidTissueNormal-EX.tsv.gz
    │   ├── LUAD-PrimaryTumor-EX.tsv.gz
    │   ├── LUAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── LUAD-SolidTissueNormal-EX.tsv.gz
    │   ├── LUSC-PrimaryTumor-EX.tsv.gz
    │   ├── LUSC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── LUSC-SolidTissueNormal-EX.tsv.gz
    │   ├── MESO-PrimaryTumor-EX.tsv.gz
    │   ├── OV-PrimaryTumor-EX.tsv.gz
    │   ├── OV-RecurrentTumor-EX.tsv.gz
    │   ├── PAAD-PrimaryTumor-EX.tsv.gz
    │   ├── PAAD-SolidTissueNormal-EX.tsv.gz
    │   ├── PANCAN-SolidTissueNormal-EX.tsv.gz
    │   ├── PCPG-PrimaryTumor-EX.tsv.gz
    │   ├── PCPG-SolidTissueNormal-EX.tsv.gz
    │   ├── PRAD-PrimaryTumor-EX.tsv.gz
    │   ├── PRAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── PRAD-SolidTissueNormal-EX.tsv.gz
    │   ├── READ-PrimaryTumor-EX.tsv.gz
    │   ├── READ-SolidTissueNormal-EX.tsv.gz
    │   ├── Riaz2017-NIV3-NAIVE-EX.tsv.gz
    │   ├── Riaz2017-NIV3-PROG-EX.tsv.gz
    │   ├── Riaz2017-ON_CombPD1_CTLA4-EX.tsv.gz
    │   ├── Riaz2017-ON-EX.tsv.gz
    │   ├── Riaz2017-ON_PD1-EX.tsv.gz
    │   ├── Riaz2017-PRE_CombPD1_CTLA4-EX.tsv.gz
    │   ├── Riaz2017-PRE-EX.tsv.gz
    │   ├── Riaz2017-PRE_PD1-EX.tsv.gz
    │   ├── SARC-PrimaryTumor-EX.tsv.gz
    │   ├── SARC-RecurrentTumor-EX.tsv.gz
    │   ├── SKCM-Metastatic-EX.tsv.gz
    │   ├── SKCM-Metastatic_vs_PrimaryTumor-EX.tsv.gz
    │   ├── SKCM-PrimaryTumor-EX.tsv.gz
    │   ├── STAD-PrimaryTumor-EX.tsv.gz
    │   ├── STAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── STAD-SolidTissueNormal-EX.tsv.gz
    │   ├── TGCT-PrimaryTumor-EX.tsv.gz
    │   ├── THCA-Metastatic-EX.tsv.gz
    │   ├── THCA-Metastatic_vs_PrimaryTumor-EX.tsv.gz
    │   ├── THCA-PrimaryTumor-EX.tsv.gz
    │   ├── THCA-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── THCA-SolidTissueNormal-EX.tsv.gz
    │   ├── THYM-PrimaryTumor-EX.tsv.gz
    │   ├── tumorigenesis-EX.tsv.gz
    │   ├── UCEC-PrimaryTumor-EX.tsv.gz
    │   ├── UCEC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── UCEC-SolidTissueNormal-EX.tsv.gz
    │   ├── UCS-PrimaryTumor-EX.tsv.gz
    │   └── UVM-PrimaryTumor-EX.tsv.gz
    ├── sf_activity_regulation
    │   ├── genexpr_tpm_vs_activity-ACC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-BLCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-BRCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-CESC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-CHOL-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-COAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-DLBC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-ESCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-GBM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-HNSC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-KICH-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-KIRC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-KIRP-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-LAML-PrimaryBloodDerivedCancerPeripheralBlood.tsv.gz
    │   ├── genexpr_tpm_vs_activity-LGG-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-LIHC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-LUAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-LUSC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-MESO-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-OV-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-PAAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-PCPG-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-PRAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-READ-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-SARC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-SKCM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-STAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-TGCT-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-THCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-THYM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-UCEC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm_vs_activity-UCS-PrimaryTumor.tsv.gz
    │   └── genexpr_tpm_vs_activity-UVM-PrimaryTumor.tsv.gz
    ├── sf_cross_regulation
    │   ├── ACC-PrimaryTumor.tsv.gz
    │   ├── BLCA-PrimaryTumor.tsv.gz
    │   ├── BRCA-PrimaryTumor.tsv.gz
    │   ├── CESC-PrimaryTumor.tsv.gz
    │   ├── CHOL-PrimaryTumor.tsv.gz
    │   ├── COAD-PrimaryTumor.tsv.gz
    │   ├── DLBC-PrimaryTumor.tsv.gz
    │   ├── ESCA-PrimaryTumor.tsv.gz
    │   ├── GBM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-ACC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-BLCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-BRCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-CESC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-CHOL-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-COAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-DLBC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-ESCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-GBM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-HNSC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-KICH-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-KIRC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-KIRP-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-LAML-PrimaryBloodDerivedCancerPeripheralBlood.tsv.gz
    │   ├── genexpr_tpm-LGG-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-LIHC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-LUAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-LUSC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-MESO-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-OV-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-PAAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-PCPG-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-PRAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-READ-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-SARC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-SKCM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-STAD-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-TGCT-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-THCA-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-THYM-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-UCEC-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-UCS-PrimaryTumor.tsv.gz
    │   ├── genexpr_tpm-UVM-PrimaryTumor.tsv.gz
    │   ├── HNSC-PrimaryTumor.tsv.gz
    │   ├── KICH-PrimaryTumor.tsv.gz
    │   ├── KIRC-PrimaryTumor.tsv.gz
    │   ├── KIRP-PrimaryTumor.tsv.gz
    │   ├── LAML-PrimaryBloodDerivedCancerPeripheralBlood.tsv.gz
    │   ├── LGG-PrimaryTumor.tsv.gz
    │   ├── LIHC-PrimaryTumor.tsv.gz
    │   ├── LUAD-PrimaryTumor.tsv.gz
    │   ├── LUSC-PrimaryTumor.tsv.gz
    │   ├── MESO-PrimaryTumor.tsv.gz
    │   ├── OV-PrimaryTumor.tsv.gz
    │   ├── PAAD-PrimaryTumor.tsv.gz
    │   ├── PCPG-PrimaryTumor.tsv.gz
    │   ├── PRAD-PrimaryTumor.tsv.gz
    │   ├── protein_activity-ACC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-BLCA-PrimaryTumor.tsv.gz
    │   ├── protein_activity-BRCA-PrimaryTumor.tsv.gz
    │   ├── protein_activity-CESC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-CHOL-PrimaryTumor.tsv.gz
    │   ├── protein_activity-COAD-PrimaryTumor.tsv.gz
    │   ├── protein_activity-DLBC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-ESCA-PrimaryTumor.tsv.gz
    │   ├── protein_activity-GBM-PrimaryTumor.tsv.gz
    │   ├── protein_activity-HNSC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-KICH-PrimaryTumor.tsv.gz
    │   ├── protein_activity-KIRC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-KIRP-PrimaryTumor.tsv.gz
    │   ├── protein_activity-LAML-PrimaryBloodDerivedCancerPeripheralBlood.tsv.gz
    │   ├── protein_activity-LGG-PrimaryTumor.tsv.gz
    │   ├── protein_activity-LIHC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-LUAD-PrimaryTumor.tsv.gz
    │   ├── protein_activity-LUSC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-MESO-PrimaryTumor.tsv.gz
    │   ├── protein_activity-OV-PrimaryTumor.tsv.gz
    │   ├── protein_activity-PAAD-PrimaryTumor.tsv.gz
    │   ├── protein_activity-PCPG-PrimaryTumor.tsv.gz
    │   ├── protein_activity-PRAD-PrimaryTumor.tsv.gz
    │   ├── protein_activity-READ-PrimaryTumor.tsv.gz
    │   ├── protein_activity-SARC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-SKCM-PrimaryTumor.tsv.gz
    │   ├── protein_activity-STAD-PrimaryTumor.tsv.gz
    │   ├── protein_activity-TGCT-PrimaryTumor.tsv.gz
    │   ├── protein_activity-THCA-PrimaryTumor.tsv.gz
    │   ├── protein_activity-THYM-PrimaryTumor.tsv.gz
    │   ├── protein_activity-UCEC-PrimaryTumor.tsv.gz
    │   ├── protein_activity-UCS-PrimaryTumor.tsv.gz
    │   ├── protein_activity-UVM-PrimaryTumor.tsv.gz
    │   ├── READ-PrimaryTumor.tsv.gz
    │   ├── SARC-PrimaryTumor.tsv.gz
    │   ├── SKCM-PrimaryTumor.tsv.gz
    │   ├── STAD-PrimaryTumor.tsv.gz
    │   ├── TGCT-PrimaryTumor.tsv.gz
    │   ├── THCA-PrimaryTumor.tsv.gz
    │   ├── THYM-PrimaryTumor.tsv.gz
    │   ├── UCEC-PrimaryTumor.tsv.gz
    │   ├── UCS-PrimaryTumor.tsv.gz
    │   └── UVM-PrimaryTumor.tsv.gz
    ├── signatures
    │   ├── ACC-PrimaryTumor-EX.tsv.gz
    │   ├── Bian2018-EX.tsv.gz
    │   ├── BLCA-PrimaryTumor-EX.tsv.gz
    │   ├── BLCA-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── BLCA-SolidTissueNormal-EX.tsv.gz
    │   ├── BRCA-Metastatic-EX.tsv.gz
    │   ├── BRCA-Metastatic_vs_PrimaryTumor-EX.tsv.gz
    │   ├── BRCA-PrimaryTumor-EX.tsv.gz
    │   ├── BRCA-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── BRCA-SolidTissueNormal-EX.tsv.gz
    │   ├── CCLE-EX.tsv.gz
    │   ├── CESC-PrimaryTumor-EX.tsv.gz
    │   ├── CESC-SolidTissueNormal-EX.tsv.gz
    │   ├── CHOL-PrimaryTumor-EX.tsv.gz
    │   ├── CHOL-SolidTissueNormal-EX.tsv.gz
    │   ├── COAD-PrimaryTumor-EX.tsv.gz
    │   ├── COAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── COAD-SolidTissueNormal-EX.tsv.gz
    │   ├── DLBC-PrimaryTumor-EX.tsv.gz
    │   ├── driver_mutations-EX.tsv.gz
    │   ├── ESCA-PrimaryTumor-EX.tsv.gz
    │   ├── ESCA-SolidTissueNormal-EX.tsv.gz
    │   ├── GBM-PrimaryTumor-EX.tsv.gz
    │   ├── GBM-RecurrentTumor-EX.tsv.gz
    │   ├── GBM-SolidTissueNormal-EX.tsv.gz
    │   ├── HNSC-PrimaryTumor-EX.tsv.gz
    │   ├── HNSC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── HNSC-SolidTissueNormal-EX.tsv.gz
    │   ├── KICH-PrimaryTumor-EX.tsv.gz
    │   ├── KICH-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── KICH-SolidTissueNormal-EX.tsv.gz
    │   ├── KIRC-PrimaryTumor-EX.tsv.gz
    │   ├── KIRC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── KIRC-SolidTissueNormal-EX.tsv.gz
    │   ├── KIRP-PrimaryTumor-EX.tsv.gz
    │   ├── KIRP-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── KIRP-SolidTissueNormal-EX.tsv.gz
    │   ├── LAML-PrimaryBloodDerivedCancerPeripheralBlood-EX.tsv.gz
    │   ├── LGG-PrimaryTumor-EX.tsv.gz
    │   ├── LGG-RecurrentTumor-EX.tsv.gz
    │   ├── LIHC-PrimaryTumor-EX.tsv.gz
    │   ├── LIHC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── LIHC-RecurrentTumor-EX.tsv.gz
    │   ├── LIHC-SolidTissueNormal-EX.tsv.gz
    │   ├── LUAD-PrimaryTumor-EX.tsv.gz
    │   ├── LUAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── LUAD-SolidTissueNormal-EX.tsv.gz
    │   ├── LUSC-PrimaryTumor-EX.tsv.gz
    │   ├── LUSC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── LUSC-SolidTissueNormal-EX.tsv.gz
    │   ├── MESO-PrimaryTumor-EX.tsv.gz
    │   ├── OV-PrimaryTumor-EX.tsv.gz
    │   ├── OV-RecurrentTumor-EX.tsv.gz
    │   ├── PAAD-PrimaryTumor-EX.tsv.gz
    │   ├── PAAD-SolidTissueNormal-EX.tsv.gz
    │   ├── PANCAN-SolidTissueNormal-EX.tsv.gz
    │   ├── PCPG-PrimaryTumor-EX.tsv.gz
    │   ├── PCPG-SolidTissueNormal-EX.tsv.gz
    │   ├── PRAD-PrimaryTumor-EX.tsv.gz
    │   ├── PRAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── PRAD-SolidTissueNormal-EX.tsv.gz
    │   ├── READ-PrimaryTumor-EX.tsv.gz
    │   ├── READ-SolidTissueNormal-EX.tsv.gz
    │   ├── Riaz2017-NIV3-NAIVE-EX.tsv.gz
    │   ├── Riaz2017-NIV3-PROG-EX.tsv.gz
    │   ├── Riaz2017-ON_CombPD1_CTLA4-EX.tsv.gz
    │   ├── Riaz2017-On-EX.tsv.gz
    │   ├── Riaz2017-ON-EX.tsv.gz
    │   ├── Riaz2017-ON_PD1-EX.tsv.gz
    │   ├── Riaz2017-PRE_CombPD1_CTLA4-EX.tsv.gz
    │   ├── Riaz2017-Pre-EX.tsv.gz
    │   ├── Riaz2017-PRE-EX.tsv.gz
    │   ├── Riaz2017-PRE_PD1-EX.tsv.gz
    │   ├── SARC-PrimaryTumor-EX.tsv.gz
    │   ├── SARC-RecurrentTumor-EX.tsv.gz
    │   ├── SKCM-Metastatic-EX.tsv.gz
    │   ├── SKCM-Metastatic_vs_PrimaryTumor-EX.tsv.gz
    │   ├── SKCM-PrimaryTumor-EX.tsv.gz
    │   ├── STAD-PrimaryTumor-EX.tsv.gz
    │   ├── STAD-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── STAD-SolidTissueNormal-EX.tsv.gz
    │   ├── TGCT-PrimaryTumor-EX.tsv.gz
    │   ├── THCA-Metastatic-EX.tsv.gz
    │   ├── THCA-Metastatic_vs_PrimaryTumor-EX.tsv.gz
    │   ├── THCA-PrimaryTumor-EX.tsv.gz
    │   ├── THCA-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── THCA-SolidTissueNormal-EX.tsv.gz
    │   ├── THYM-PrimaryTumor-EX.tsv.gz
    │   ├── tumorigenesis-EX.tsv.gz
    │   ├── UCEC-PrimaryTumor-EX.tsv.gz
    │   ├── UCEC-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz
    │   ├── UCEC-SolidTissueNormal-EX.tsv.gz
    │   ├── UCS-PrimaryTumor-EX.tsv.gz
    │   └── UVM-PrimaryTumor-EX.tsv.gz
    ├── survival_analysis
    │   ├── genexpr_tpm-ACC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-ACC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-BLCA-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-BLCA-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-BRCA-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-BRCA-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-CESC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-CESC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-CHOL-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-CHOL-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-COAD-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-COAD-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-DLBC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-DLBC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-ESCA-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-ESCA-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-GBM-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-GBM-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-HNSC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-HNSC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-KICH-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-KICH-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-KIRC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-KIRC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-KIRP-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-KIRP-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-LAML-PrimaryBloodDerivedCancerPeripheralBlood-cat.tsv.gz
    │   ├── genexpr_tpm-LAML-PrimaryBloodDerivedCancerPeripheralBlood-surv.tsv.gz
    │   ├── genexpr_tpm-LGG-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-LGG-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-LIHC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-LIHC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-LUAD-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-LUAD-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-LUSC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-LUSC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-MESO-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-MESO-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-OV-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-OV-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-PAAD-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-PAAD-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-PCPG-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-PCPG-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-PRAD-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-PRAD-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-READ-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-READ-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-SARC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-SARC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-SKCM-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-SKCM-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-STAD-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-STAD-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-TGCT-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-TGCT-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-THCA-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-THCA-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-THYM-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-THYM-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-UCEC-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-UCEC-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-UCS-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-UCS-PrimaryTumor-surv.tsv.gz
    │   ├── genexpr_tpm-UVM-PrimaryTumor-cat.tsv.gz
    │   ├── genexpr_tpm-UVM-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-ACC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-ACC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-BLCA-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-BLCA-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-BRCA-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-BRCA-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-CESC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-CESC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-CHOL-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-CHOL-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-COAD-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-COAD-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-DLBC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-DLBC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-ESCA-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-ESCA-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-GBM-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-GBM-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-HNSC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-HNSC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-KICH-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-KICH-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-KIRC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-KIRC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-KIRP-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-KIRP-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-LAML-PrimaryBloodDerivedCancerPeripheralBlood-cat.tsv.gz
    │   ├── protein_activity-LAML-PrimaryBloodDerivedCancerPeripheralBlood-surv.tsv.gz
    │   ├── protein_activity-LGG-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-LGG-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-LIHC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-LIHC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-LUAD-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-LUAD-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-LUSC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-LUSC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-MESO-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-MESO-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-OV-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-OV-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-PAAD-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-PAAD-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-PCPG-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-PCPG-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-PRAD-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-PRAD-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-READ-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-READ-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-SARC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-SARC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-SKCM-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-SKCM-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-STAD-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-STAD-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-TGCT-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-TGCT-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-THCA-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-THCA-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-THYM-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-THYM-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-UCEC-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-UCEC-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-UCS-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-UCS-PrimaryTumor-surv.tsv.gz
    │   ├── protein_activity-UVM-PrimaryTumor-cat.tsv.gz
    │   ├── protein_activity-UVM-PrimaryTumor-surv.tsv.gz
    │   ├── splicing-EX-Riaz2017-ON-cat.tsv.gz
    │   ├── splicing-EX-Riaz2017-ON-surv.tsv.gz
    │   ├── splicing-EX-Riaz2017-PRE-cat.tsv.gz
    │   └── splicing-EX-Riaz2017-PRE-surv.tsv.gz
    └── survival_analysis_with_confounders
        ├── genexpr_tpm-ACC-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-ACC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-BLCA-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-BLCA-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-BRCA-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-BRCA-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-CESC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-CHOL-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-CHOL-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-COAD-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-COAD-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-ESCA-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-ESCA-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-GBM-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-HNSC-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-HNSC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-KICH-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-KICH-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-KIRC-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-KIRC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-KIRP-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-KIRP-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-LGG-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-LIHC-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-LIHC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-LUAD-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-LUAD-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-LUSC-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-LUSC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-MESO-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-MESO-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-OV-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-PAAD-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-PAAD-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-PCPG-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-PRAD-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-READ-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-READ-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-SARC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-SKCM-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-SKCM-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-STAD-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-STAD-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-TGCT-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-TGCT-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-THCA-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-THCA-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-UCEC-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-UCS-PrimaryTumor-surv.tsv.gz
        ├── genexpr_tpm-UVM-PrimaryTumor-cat.tsv.gz
        ├── genexpr_tpm-UVM-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-ACC-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-ACC-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-BLCA-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-BLCA-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-BRCA-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-BRCA-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-CHOL-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-CHOL-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-COAD-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-COAD-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-ESCA-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-ESCA-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-HNSC-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-HNSC-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-KICH-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-KICH-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-KIRC-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-KIRC-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-KIRP-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-KIRP-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-LIHC-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-LIHC-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-LUAD-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-LUAD-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-LUSC-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-LUSC-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-MESO-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-MESO-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-PAAD-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-PAAD-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-READ-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-READ-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-SKCM-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-SKCM-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-STAD-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-STAD-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-TGCT-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-TGCT-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-THCA-PrimaryTumor-cat.tsv.gz
        ├── protein_activity-THCA-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-UCEC-PrimaryTumor-surv.tsv.gz
        ├── protein_activity-UVM-PrimaryTumor-cat.tsv.gz
        └── protein_activity-UVM-PrimaryTumor-surv.tsv.gz
```