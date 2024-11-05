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

## Expected outputs
```{shell}
$ tree data/prep/
data/prep/
├── clip_peaks_mapped
│   └── POSTAR3.tsv.gz
├── demeter2
│   └── CCLE.tsv.gz
├── doubling_times
│   └── CCLE.tsv.gz
├── event_psi
│   ├── ACC-ALTA.tsv.gz
│   ├── ACC-ALTD.tsv.gz
│   ├── ACC-EX.tsv.gz
│   ├── ACC-INT.tsv.gz
│   ├── ACC-PrimaryTumor-EX.tsv.gz
│   ├── Bian2018-ALTA.tsv.gz
│   ├── Bian2018-ALTD.tsv.gz
│   ├── Bian2018-EX.tsv.gz
│   ├── Bian2018-INT.tsv.gz
│   ├── BLCA-ALTA.tsv.gz
│   ├── BLCA-ALTD.tsv.gz
│   ├── BLCA-EX.tsv.gz
│   ├── BLCA-INT.tsv.gz
│   ├── BLCA-PrimaryTumor-EX.tsv.gz
│   ├── BLCA-SolidTissueNormal-EX.tsv.gz
│   ├── BRCA-ALTA.tsv.gz
│   ├── BRCA-ALTD.tsv.gz
│   ├── BRCA-EX.tsv.gz
│   ├── BRCA-INT.tsv.gz
│   ├── BRCA-Metastatic-EX.tsv.gz
│   ├── BRCA-PrimaryTumor-EX.tsv.gz
│   ├── BRCA-SolidTissueNormal-EX.tsv.gz
│   ├── CardosoMoreira2020-ALTA.tsv.gz
│   ├── CardosoMoreira2020-ALTD.tsv.gz
│   ├── CardosoMoreira2020-EX.tsv.gz
│   ├── CardosoMoreira2020-INT.tsv.gz
│   ├── CCLE-ALTA.tsv.gz
│   ├── CCLE-ALTD.tsv.gz
│   ├── CCLE-EX.tsv.gz
│   ├── CCLE-INT.tsv.gz
│   ├── CESC-ALTA.tsv.gz
│   ├── CESC-ALTD.tsv.gz
│   ├── CESC-EX.tsv.gz
│   ├── CESC-INT.tsv.gz
│   ├── CESC-Metastatic-EX.tsv.gz
│   ├── CESC-PrimaryTumor-EX.tsv.gz
│   ├── CESC-SolidTissueNormal-EX.tsv.gz
│   ├── CHOL-ALTA.tsv.gz
│   ├── CHOL-ALTD.tsv.gz
│   ├── CHOL-EX.tsv.gz
│   ├── CHOL-INT.tsv.gz
│   ├── CHOL-PrimaryTumor-EX.tsv.gz
│   ├── CHOL-SolidTissueNormal-EX.tsv.gz
│   ├── COAD-ALTA.tsv.gz
│   ├── COAD-ALTD.tsv.gz
│   ├── COAD-EX.tsv.gz
│   ├── COAD-INT.tsv.gz
│   ├── COAD-Metastatic-EX.tsv.gz
│   ├── COAD-PrimaryTumor-EX.tsv.gz
│   ├── COAD-SolidTissueNormal-EX.tsv.gz
│   ├── DLBC-ALTA.tsv.gz
│   ├── DLBC-ALTD.tsv.gz
│   ├── DLBC-EX.tsv.gz
│   ├── DLBC-INT.tsv.gz
│   ├── DLBC-PrimaryTumor-EX.tsv.gz
│   ├── ENASFS-ALTA.tsv.gz
│   ├── ENASFS-ALTD.tsv.gz
│   ├── ENASFS-EX.tsv.gz
│   ├── ENASFS-INT.tsv.gz
│   ├── ENCORE-ALTA.tsv.gz
│   ├── ENCORE-ALTD.tsv.gz
│   ├── ENCORE-EX.tsv.gz
│   ├── ENCORE-INT.tsv.gz
│   ├── ENCOREKD-ALTA.tsv.gz
│   ├── ENCOREKD-ALTD.tsv.gz
│   ├── ENCOREKD-EX.tsv.gz
│   ├── ENCOREKD-INT.tsv.gz
│   ├── ENCOREKO-ALTA.tsv.gz
│   ├── ENCOREKO-ALTD.tsv.gz
│   ├── ENCOREKO-EX.tsv.gz
│   ├── ENCOREKO-INT.tsv.gz
│   ├── ESCA-ALTA.tsv.gz
│   ├── ESCA-ALTD.tsv.gz
│   ├── ESCA-EX.tsv.gz
│   ├── ESCA-INT.tsv.gz
│   ├── ESCA-Metastatic-EX.tsv.gz
│   ├── ESCA-PrimaryTumor-EX.tsv.gz
│   ├── ESCA-SolidTissueNormal-EX.tsv.gz
│   ├── GBM-ALTA.tsv.gz
│   ├── GBM-ALTD.tsv.gz
│   ├── GBM-EX.tsv.gz
│   ├── GBM-INT.tsv.gz
│   ├── GBM-PrimaryTumor-EX.tsv.gz
│   ├── GBM-RecurrentTumor-EX.tsv.gz
│   ├── GBM-SolidTissueNormal-EX.tsv.gz
│   ├── Hafner2019-ALTA.tsv.gz
│   ├── Hafner2019-ALTD.tsv.gz
│   ├── Hafner2019-EX.tsv.gz
│   ├── Hafner2019-INT.tsv.gz
│   ├── HNSC-ALTA.tsv.gz
│   ├── HNSC-ALTD.tsv.gz
│   ├── HNSC-EX.tsv.gz
│   ├── HNSC-INT.tsv.gz
│   ├── HNSC-Metastatic-EX.tsv.gz
│   ├── HNSC-PrimaryTumor-EX.tsv.gz
│   ├── HNSC-SolidTissueNormal-EX.tsv.gz
│   ├── ipsc_differentiation-ALTA.tsv.gz
│   ├── ipsc_differentiation-ALTD.tsv.gz
│   ├── ipsc_differentiation-EX.tsv.gz
│   ├── ipsc_differentiation-INT.tsv.gz
│   ├── KICH-ALTA.tsv.gz
│   ├── KICH-ALTD.tsv.gz
│   ├── KICH-EX.tsv.gz
│   ├── KICH-INT.tsv.gz
│   ├── KICH-PrimaryTumor-EX.tsv.gz
│   ├── KICH-SolidTissueNormal-EX.tsv.gz
│   ├── KIRC-ALTA.tsv.gz
│   ├── KIRC-ALTD.tsv.gz
│   ├── KIRC-EX.tsv.gz
│   ├── KIRC-INT.tsv.gz
│   ├── KIRC-PrimaryTumor-EX.tsv.gz
│   ├── KIRC-SolidTissueNormal-EX.tsv.gz
│   ├── KIRP-ALTA.tsv.gz
│   ├── KIRP-ALTD.tsv.gz
│   ├── KIRP-EX.tsv.gz
│   ├── KIRP-INT.tsv.gz
│   ├── KIRP-PrimaryTumor-EX.tsv.gz
│   ├── KIRP-SolidTissueNormal-EX.tsv.gz
│   ├── LAML-ALTA.tsv.gz
│   ├── LAML-ALTD.tsv.gz
│   ├── LAML-EX.tsv.gz
│   ├── LAML-INT.tsv.gz
│   ├── LAML-PrimaryBloodDerivedCancerPeripheralBlood-EX.tsv.gz
│   ├── LGG-ALTA.tsv.gz
│   ├── LGG-ALTD.tsv.gz
│   ├── LGG-EX.tsv.gz
│   ├── LGG-INT.tsv.gz
│   ├── LGG-PrimaryTumor-EX.tsv.gz
│   ├── LGG-RecurrentTumor-EX.tsv.gz
│   ├── LIHC-ALTA.tsv.gz
│   ├── LIHC-ALTD.tsv.gz
│   ├── LIHC-EX.tsv.gz
│   ├── LIHC-INT.tsv.gz
│   ├── LIHC-PrimaryTumor-EX.tsv.gz
│   ├── LIHC-RecurrentTumor-EX.tsv.gz
│   ├── LIHC-SolidTissueNormal-EX.tsv.gz
│   ├── Lu2021-ALTA.tsv.gz
│   ├── Lu2021-ALTD.tsv.gz
│   ├── Lu2021-EX.tsv.gz
│   ├── Lu2021-INT.tsv.gz
│   ├── LUAD-ALTA.tsv.gz
│   ├── LUAD-ALTD.tsv.gz
│   ├── LUAD-EX.tsv.gz
│   ├── LUAD-INT.tsv.gz
│   ├── LUAD-PrimaryTumor-EX.tsv.gz
│   ├── LUAD-SolidTissueNormal-EX.tsv.gz
│   ├── LUSC-ALTA.tsv.gz
│   ├── LUSC-ALTD.tsv.gz
│   ├── LUSC-EX.tsv.gz
│   ├── LUSC-INT.tsv.gz
│   ├── LUSC-PrimaryTumor-EX.tsv.gz
│   ├── LUSC-SolidTissueNormal-EX.tsv.gz
│   ├── MESO-ALTA.tsv.gz
│   ├── MESO-ALTD.tsv.gz
│   ├── MESO-EX.tsv.gz
│   ├── MESO-INT.tsv.gz
│   ├── MESO-PrimaryTumor-EX.tsv.gz
│   ├── Nijhuis2020-ALTA.tsv.gz
│   ├── Nijhuis2020-ALTD.tsv.gz
│   ├── Nijhuis2020-EX.tsv.gz
│   ├── Nijhuis2020-INT.tsv.gz
│   ├── OV-ALTA.tsv.gz
│   ├── OV-ALTD.tsv.gz
│   ├── OV-EX.tsv.gz
│   ├── OV-INT.tsv.gz
│   ├── OV-PrimaryTumor-EX.tsv.gz
│   ├── OV-RecurrentTumor-EX.tsv.gz
│   ├── PAAD-ALTA.tsv.gz
│   ├── PAAD-ALTD.tsv.gz
│   ├── PAAD-EX.tsv.gz
│   ├── PAAD-INT.tsv.gz
│   ├── PAAD-Metastatic-EX.tsv.gz
│   ├── PAAD-PrimaryTumor-EX.tsv.gz
│   ├── PAAD-SolidTissueNormal-EX.tsv.gz
│   ├── PANCAN-PrimaryTumor-EX.tsv.gz
│   ├── PANCAN-SolidTissueNormal-EX.tsv.gz
│   ├── PCPG-ALTA.tsv.gz
│   ├── PCPG-ALTD.tsv.gz
│   ├── PCPG-EX.tsv.gz
│   ├── PCPG-INT.tsv.gz
│   ├── PCPG-Metastatic-EX.tsv.gz
│   ├── PCPG-PrimaryTumor-EX.tsv.gz
│   ├── PCPG-SolidTissueNormal-EX.tsv.gz
│   ├── PRAD-ALTA.tsv.gz
│   ├── PRAD-ALTD.tsv.gz
│   ├── PRAD-EX.tsv.gz
│   ├── PRAD-INT.tsv.gz
│   ├── PRAD-Metastatic-EX.tsv.gz
│   ├── PRAD-PrimaryTumor-EX.tsv.gz
│   ├── PRAD-SolidTissueNormal-EX.tsv.gz
│   ├── READ-ALTA.tsv.gz
│   ├── READ-ALTD.tsv.gz
│   ├── READ-EX.tsv.gz
│   ├── READ-INT.tsv.gz
│   ├── READ-PrimaryTumor-EX.tsv.gz
│   ├── READ-SolidTissueNormal-EX.tsv.gz
│   ├── Riaz2017-ALTA.tsv.gz
│   ├── Riaz2017-ALTD.tsv.gz
│   ├── Riaz2017-EX.tsv.gz
│   ├── Riaz2017-INT.tsv.gz
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
│   ├── SARC-ALTA.tsv.gz
│   ├── SARC-ALTD.tsv.gz
│   ├── SARC-EX.tsv.gz
│   ├── SARC-INT.tsv.gz
│   ├── SARC-Metastatic-EX.tsv.gz
│   ├── SARC-PrimaryTumor-EX.tsv.gz
│   ├── SARC-RecurrentTumor-EX.tsv.gz
│   ├── SARC-SolidTissueNormal-EX.tsv.gz
│   ├── sf_drugs-ALTA.tsv.gz
│   ├── sf_drugs-ALTD.tsv.gz
│   ├── sf_drugs-EX.tsv.gz
│   ├── sf_drugs-INT.tsv.gz
│   ├── sf_ptms-ALTA.tsv.gz
│   ├── sf_ptms-ALTD.tsv.gz
│   ├── sf_ptms-EX.tsv.gz
│   ├── sf_ptms-INT.tsv.gz
│   ├── SKCM-ALTA.tsv.gz
│   ├── SKCM-ALTD.tsv.gz
│   ├── SKCM-EX.tsv.gz
│   ├── SKCM-INT.tsv.gz
│   ├── SKCM-Metastatic-EX.tsv.gz
│   ├── SKCM-PrimaryTumor-EX.tsv.gz
│   ├── SKCM-SolidTissueNormal-EX.tsv.gz
│   ├── STAD-ALTA.tsv.gz
│   ├── STAD-ALTD.tsv.gz
│   ├── STAD-EX.tsv.gz
│   ├── STAD-INT.tsv.gz
│   ├── STAD-PrimaryTumor-EX.tsv.gz
│   ├── STAD-SolidTissueNormal-EX.tsv.gz
│   ├── TGCT-ALTA.tsv.gz
│   ├── TGCT-ALTD.tsv.gz
│   ├── TGCT-EX.tsv.gz
│   ├── TGCT-INT.tsv.gz
│   ├── TGCT-PrimaryTumor-EX.tsv.gz
│   ├── THCA-ALTA.tsv.gz
│   ├── THCA-ALTD.tsv.gz
│   ├── THCA-EX.tsv.gz
│   ├── THCA-INT.tsv.gz
│   ├── THCA-Metastatic-EX.tsv.gz
│   ├── THCA-PrimaryTumor-EX.tsv.gz
│   ├── THCA-SolidTissueNormal-EX.tsv.gz
│   ├── THYM-ALTA.tsv.gz
│   ├── THYM-ALTD.tsv.gz
│   ├── THYM-EX.tsv.gz
│   ├── THYM-INT.tsv.gz
│   ├── THYM-PrimaryTumor-EX.tsv.gz
│   ├── THYM-SolidTissueNormal-EX.tsv.gz
│   ├── tumorigenesis-ALTA.tsv.gz
│   ├── tumorigenesis-ALTD.tsv.gz
│   ├── tumorigenesis-EX.tsv.gz
│   ├── tumorigenesis-INT.tsv.gz
│   ├── UCEC-ALTA.tsv.gz
│   ├── UCEC-ALTD.tsv.gz
│   ├── UCEC-EX.tsv.gz
│   ├── UCEC-INT.tsv.gz
│   ├── UCEC-PrimaryTumor-EX.tsv.gz
│   ├── UCEC-SolidTissueNormal-EX.tsv.gz
│   ├── UCS-ALTA.tsv.gz
│   ├── UCS-ALTD.tsv.gz
│   ├── UCS-EX.tsv.gz
│   ├── UCS-INT.tsv.gz
│   ├── UCS-PrimaryTumor-EX.tsv.gz
│   ├── UVM-ALTA.tsv.gz
│   ├── UVM-ALTD.tsv.gz
│   ├── UVM-EX.tsv.gz
│   ├── UVM-INT.tsv.gz
│   └── UVM-PrimaryTumor-EX.tsv.gz
├── event_psi_imputation_benchmark
│   ├── CardosoMoreira2020-EX-k100.tsv.gz
│   ├── CardosoMoreira2020-EX-k10.tsv.gz
│   ├── CardosoMoreira2020-EX-k2.tsv.gz
│   ├── CardosoMoreira2020-EX-k50.tsv.gz
│   ├── CardosoMoreira2020-EX-k5.tsv.gz
│   └── CardosoMoreira2020-EX.tsv.gz
├── event_psi_imputed
│   ├── ACC-EX.tsv.gz
│   ├── BLCA-EX.tsv.gz
│   ├── BRCA-EX.tsv.gz
│   ├── CardosoMoreira2020-EX.tsv.gz
│   ├── CCLE-EX.tsv.gz
│   ├── CESC-EX.tsv.gz
│   ├── CHOL-EX.tsv.gz
│   ├── COAD-EX.tsv.gz
│   ├── DLBC-EX.tsv.gz
│   ├── ESCA-EX.tsv.gz
│   ├── GBM-EX.tsv.gz
│   ├── HNSC-EX.tsv.gz
│   ├── KICH-EX.tsv.gz
│   ├── KIRC-EX.tsv.gz
│   ├── KIRP-EX.tsv.gz
│   ├── LAML-EX.tsv.gz
│   ├── LGG-EX.tsv.gz
│   ├── LIHC-EX.tsv.gz
│   ├── LUAD-EX.tsv.gz
│   ├── LUSC-EX.tsv.gz
│   ├── MESO-EX.tsv.gz
│   ├── OV-EX.tsv.gz
│   ├── PAAD-EX.tsv.gz
│   ├── PANCAN-PrimaryTumor-EX.tsv.gz
│   ├── PANCAN-SolidTissueNormal-EX.tsv.gz
│   ├── PCPG-EX.tsv.gz
│   ├── PRAD-EX.tsv.gz
│   ├── READ-EX.tsv.gz
│   ├── SARC-EX.tsv.gz
│   ├── SKCM-EX.tsv.gz
│   ├── STAD-EX.tsv.gz
│   ├── TGCT-EX.tsv.gz
│   ├── THCA-EX.tsv.gz
│   ├── THYM-EX.tsv.gz
│   ├── UCEC-EX.tsv.gz
│   ├── UCS-EX.tsv.gz
│   └── UVM-EX.tsv.gz
├── event_psi_imputed_discretized_gaussian
│   ├── ACC-EX.tsv.gz
│   ├── BLCA-EX.tsv.gz
│   ├── BRCA-EX.tsv.gz
│   ├── CESC-EX.tsv.gz
│   ├── CHOL-EX.tsv.gz
│   ├── COAD-EX.tsv.gz
│   ├── DLBC-EX.tsv.gz
│   ├── ESCA-EX.tsv.gz
│   ├── GBM-EX.tsv.gz
│   ├── HNSC-EX.tsv.gz
│   ├── KICH-EX.tsv.gz
│   ├── KIRC-EX.tsv.gz
│   ├── KIRP-EX.tsv.gz
│   ├── LAML-EX.tsv.gz
│   ├── LGG-EX.tsv.gz
│   ├── LIHC-EX.tsv.gz
│   ├── LUAD-EX.tsv.gz
│   ├── LUSC-EX.tsv.gz
│   ├── MESO-EX.tsv.gz
│   ├── OV-EX.tsv.gz
│   ├── PAAD-EX.tsv.gz
│   ├── PCPG-EX.tsv.gz
│   ├── PRAD-EX.tsv.gz
│   ├── READ-EX.tsv.gz
│   ├── SARC-EX.tsv.gz
│   ├── SKCM-EX.tsv.gz
│   ├── STAD-EX.tsv.gz
│   ├── TGCT-EX.tsv.gz
│   ├── THCA-EX.tsv.gz
│   ├── THYM-EX.tsv.gz
│   ├── UCEC-EX.tsv.gz
│   ├── UCS-EX.tsv.gz
│   └── UVM-EX.tsv.gz
├── event_psi_imputed_discretized_qep
│   ├── CCLE-EX.tsv.gz
│   ├── LAML-EX.tsv.gz
│   └── LIHC-EX.tsv.gz
├── genexpr_counts
│   ├── ACC.tsv.gz
│   ├── BLCA.tsv.gz
│   ├── CESC.tsv.gz
│   ├── CHOL.tsv.gz
│   ├── COAD.tsv.gz
│   ├── DLBC.tsv.gz
│   ├── ESCA.tsv.gz
│   ├── GBM.tsv.gz
│   ├── HNSC.tsv.gz
│   ├── KICH.tsv.gz
│   ├── KIRC.tsv.gz
│   ├── KIRP.tsv.gz
│   ├── LAML.tsv.gz
│   ├── LGG.tsv.gz
│   ├── LIHC.tsv.gz
│   ├── LUAD.tsv.gz
│   ├── LUSC.tsv.gz
│   ├── MESO.tsv.gz
│   ├── PAAD.tsv.gz
│   ├── PCPG.tsv.gz
│   ├── PRAD.tsv.gz
│   ├── READ.tsv.gz
│   ├── SARC.tsv.gz
│   ├── SKCM.tsv.gz
│   ├── STAD.tsv.gz
│   ├── TGCT.tsv.gz
│   ├── THCA.tsv.gz
│   ├── THYM.tsv.gz
│   ├── UCEC.tsv.gz
│   ├── UCS.tsv.gz
│   └── UVM.tsv.gz
├── genexpr_tpm
│   ├── ACC-PrimaryTumor.tsv.gz
│   ├── ACC.tsv.gz
│   ├── Bian2018.tsv.gz
│   ├── BLCA-PrimaryTumor.tsv.gz
│   ├── BLCA-SolidTissueNormal.tsv.gz
│   ├── BLCA.tsv.gz
│   ├── BRCA-Metastatic.tsv.gz
│   ├── BRCA-PrimaryTumor.tsv.gz
│   ├── BRCA-SolidTissueNormal.tsv.gz
│   ├── BRCA.tsv.gz
│   ├── CardosoMoreira2020.tsv.gz
│   ├── CCLE.tsv.gz
│   ├── CESC-PrimaryTumor.tsv.gz
│   ├── CESC-SolidTissueNormal.tsv.gz
│   ├── CESC.tsv.gz
│   ├── CHOL-PrimaryTumor.tsv.gz
│   ├── CHOL-SolidTissueNormal.tsv.gz
│   ├── CHOL.tsv.gz
│   ├── COAD-PrimaryTumor.tsv.gz
│   ├── COAD-SolidTissueNormal.tsv.gz
│   ├── COAD.tsv.gz
│   ├── DLBC-PrimaryTumor.tsv.gz
│   ├── DLBC.tsv.gz
│   ├── ENASFS.tsv.gz
│   ├── ENCOREKD.tsv.gz
│   ├── ENCOREKO.tsv.gz
│   ├── ENCORE.tsv.gz
│   ├── ESCA-PrimaryTumor.tsv.gz
│   ├── ESCA-SolidTissueNormal.tsv.gz
│   ├── ESCA.tsv.gz
│   ├── GBM-PrimaryTumor.tsv.gz
│   ├── GBM-RecurrentTumor.tsv.gz
│   ├── GBM-SolidTissueNormal.tsv.gz
│   ├── GBM.tsv.gz
│   ├── Hafner2019.tsv.gz
│   ├── HNSC-PrimaryTumor.tsv.gz
│   ├── HNSC-SolidTissueNormal.tsv.gz
│   ├── HNSC.tsv.gz
│   ├── ipsc_differentiation.tsv.gz
│   ├── KICH-PrimaryTumor.tsv.gz
│   ├── KICH-SolidTissueNormal.tsv.gz
│   ├── KICH.tsv.gz
│   ├── KIRC-PrimaryTumor.tsv.gz
│   ├── KIRC-SolidTissueNormal.tsv.gz
│   ├── KIRC.tsv.gz
│   ├── KIRP-PrimaryTumor.tsv.gz
│   ├── KIRP-SolidTissueNormal.tsv.gz
│   ├── KIRP.tsv.gz
│   ├── LAML-PrimaryBloodDerivedCancerPeripheralBlood.tsv.gz
│   ├── LAML.tsv.gz
│   ├── LGG-PrimaryTumor.tsv.gz
│   ├── LGG-RecurrentTumor.tsv.gz
│   ├── LGG.tsv.gz
│   ├── LIHC-PrimaryTumor.tsv.gz
│   ├── LIHC-RecurrentTumor.tsv.gz
│   ├── LIHC-SolidTissueNormal.tsv.gz
│   ├── LIHC.tsv.gz
│   ├── Lu2021.tsv.gz
│   ├── LUAD-PrimaryTumor.tsv.gz
│   ├── LUAD-SolidTissueNormal.tsv.gz
│   ├── LUAD.tsv.gz
│   ├── LUSC-PrimaryTumor.tsv.gz
│   ├── LUSC-SolidTissueNormal.tsv.gz
│   ├── LUSC.tsv.gz
│   ├── MESO-PrimaryTumor.tsv.gz
│   ├── MESO.tsv.gz
│   ├── Nijhuis2020.tsv.gz
│   ├── OV-PrimaryTumor.tsv.gz
│   ├── OV-RecurrentTumor.tsv.gz
│   ├── OV.tsv.gz
│   ├── PAAD-PrimaryTumor.tsv.gz
│   ├── PAAD-SolidTissueNormal.tsv.gz
│   ├── PAAD.tsv.gz
│   ├── PANCAN-PrimaryTumor.tsv.gz
│   ├── PANCAN-SolidTissueNormal.tsv.gz
│   ├── PCPG-PrimaryTumor.tsv.gz
│   ├── PCPG-SolidTissueNormal.tsv.gz
│   ├── PCPG.tsv.gz
│   ├── PRAD-PrimaryTumor.tsv.gz
│   ├── PRAD-SolidTissueNormal.tsv.gz
│   ├── PRAD.tsv.gz
│   ├── READ-PrimaryTumor.tsv.gz
│   ├── READ-SolidTissueNormal.tsv.gz
│   ├── READ.tsv.gz
│   ├── Riaz2017-ON.tsv.gz
│   ├── Riaz2017-PRE.tsv.gz
│   ├── Riaz2017.tsv.gz
│   ├── SARC-PrimaryTumor.tsv.gz
│   ├── SARC-RecurrentTumor.tsv.gz
│   ├── SARC.tsv.gz
│   ├── sf_drugs.tsv.gz
│   ├── sf_ptms.tsv.gz
│   ├── SKCM-Metastatic.tsv.gz
│   ├── SKCM-PrimaryTumor.tsv.gz
│   ├── SKCM.tsv.gz
│   ├── STAD-PrimaryTumor.tsv.gz
│   ├── STAD-SolidTissueNormal.tsv.gz
│   ├── STAD.tsv.gz
│   ├── TGCT-PrimaryTumor.tsv.gz
│   ├── TGCT.tsv.gz
│   ├── THCA-Metastatic.tsv.gz
│   ├── THCA-PrimaryTumor.tsv.gz
│   ├── THCA-SolidTissueNormal.tsv.gz
│   ├── THCA.tsv.gz
│   ├── THYM-PrimaryTumor.tsv.gz
│   ├── THYM.tsv.gz
│   ├── tumorigenesis.tsv.gz
│   ├── UCEC-PrimaryTumor.tsv.gz
│   ├── UCEC-SolidTissueNormal.tsv.gz
│   ├── UCEC.tsv.gz
│   ├── UCS-PrimaryTumor.tsv.gz
│   ├── UCS.tsv.gz
│   ├── UVM-PrimaryTumor.tsv.gz
│   └── UVM.tsv.gz
├── genexpr_tpm_discretized_qep
│   ├── CCLE.tsv.gz
│   ├── LAML.tsv.gz
│   └── LIHC.tsv.gz
├── ground_truth_kd
│   └── ENCORE
│       ├── HepG2
│       │   ├── delta_psi-EX-masked.tsv.gz
│       │   ├── delta_psi-EX.tsv.gz
│       │   ├── delta_psi_rel-EX-masked.tsv.gz
│       │   ├── delta_psi_rel-EX.tsv.gz
│       │   └── log2_fold_change_tpm.tsv.gz
│       └── K562
│           ├── delta_psi-EX-masked.tsv.gz
│           ├── delta_psi-EX.tsv.gz
│           ├── delta_psi_rel-EX-masked.tsv.gz
│           ├── delta_psi_rel-EX.tsv.gz
│           └── log2_fold_change_tpm.tsv.gz
├── ground_truth_pert
│   ├── ENASFS
│   │   ├── delta_psi-EX.tsv.gz
│   │   └── log2_fold_change_tpm.tsv.gz
│   ├── ENCOREKD
│   │   ├── HepG2
│   │   │   ├── delta_psi-EX.tsv.gz
│   │   │   └── log2_fold_change_tpm.tsv.gz
│   │   └── K562
│   │       ├── delta_psi-EX.tsv.gz
│   │       └── log2_fold_change_tpm.tsv.gz
│   ├── ENCOREKO
│   │   ├── HepG2
│   │   │   ├── delta_psi-EX.tsv.gz
│   │   │   └── log2_fold_change_tpm.tsv.gz
│   │   └── K562
│   │       ├── delta_psi-EX.tsv.gz
│   │       └── log2_fold_change_tpm.tsv.gz
│   └── SplicingLore
│       └── delta_psi-EX.tsv.gz
├── kd_transcriptomes
│   └── ENCORE
│       ├── HepG2
│       │   ├── delta_psi-EX.tsv.gz
│       │   ├── delta_psi_rel-EX.tsv.gz
│       │   └── log2_fold_change_tpm.tsv.gz
│       └── K562
│           ├── delta_psi-EX.tsv.gz
│           ├── delta_psi_rel-EX.tsv.gz
│           └── log2_fold_change_tpm.tsv.gz
├── metadata
│   ├── ACC.tsv.gz
│   ├── Bian2018.tsv.gz
│   ├── BLCA.tsv.gz
│   ├── BRCA.tsv.gz
│   ├── CardosoMoreira2020.tsv.gz
│   ├── CCLE.tsv.gz
│   ├── CESC.tsv.gz
│   ├── CHOL.tsv.gz
│   ├── COAD.tsv.gz
│   ├── DLBC.tsv.gz
│   ├── ENASFS.tsv.gz
│   ├── ENCOREKD.tsv.gz
│   ├── ENCOREKO.tsv.gz
│   ├── ESCA.tsv.gz
│   ├── GBM.tsv.gz
│   ├── Hafner2019.tsv.gz
│   ├── HNSC.tsv.gz
│   ├── ipsc_differentiation.tsv.gz
│   ├── KICH.tsv.gz
│   ├── KIRC.tsv.gz
│   ├── KIRP.tsv.gz
│   ├── LAML.tsv.gz
│   ├── LGG.tsv.gz
│   ├── LIHC.tsv.gz
│   ├── Lu2021.tsv.gz
│   ├── LUAD.tsv.gz
│   ├── LUSC.tsv.gz
│   ├── MESO.tsv.gz
│   ├── Nijhuis2020.tsv.gz
│   ├── OV.tsv.gz
│   ├── PAAD.tsv.gz
│   ├── PANCAN.tsv.gz
│   ├── PCPG.tsv.gz
│   ├── PRAD.tsv.gz
│   ├── READ.tsv.gz
│   ├── Riaz2017.tsv.gz
│   ├── SARC.tsv.gz
│   ├── sf_drugs.tsv.gz
│   ├── sf_ptms.tsv.gz
│   ├── SKCM.tsv.gz
│   ├── SplicingLore.tsv.gz
│   ├── STAD.tsv.gz
│   ├── TGCT.tsv.gz
│   ├── THCA.tsv.gz
│   ├── THYM.tsv.gz
│   ├── tumorigenesis.tsv.gz
│   ├── UCEC.tsv.gz
│   ├── UCS.tsv.gz
│   └── UVM.tsv.gz
├── metmap
│   └── CCLE.tsv.gz
├── mutations
│   ├── ACC.tsv.gz
│   ├── BLCA.tsv.gz
│   ├── BRCA.tsv.gz
│   ├── CCLE.tsv.gz
│   ├── CESC.tsv.gz
│   ├── CHOL.tsv.gz
│   ├── COAD.tsv.gz
│   ├── DLBC.tsv.gz
│   ├── ESCA.tsv.gz
│   ├── GBM.tsv.gz
│   ├── HNSC.tsv.gz
│   ├── KICH.tsv.gz
│   ├── KIRC.tsv.gz
│   ├── KIRP.tsv.gz
│   ├── LAML.tsv.gz
│   ├── LGG.tsv.gz
│   ├── LIHC.tsv.gz
│   ├── LUAD.tsv.gz
│   ├── LUSC.tsv.gz
│   ├── MESO.tsv.gz
│   ├── OV.tsv.gz
│   ├── PAAD.tsv.gz
│   ├── PCPG.tsv.gz
│   ├── PRAD.tsv.gz
│   ├── READ.tsv.gz
│   ├── SARC.tsv.gz
│   ├── SKCM.tsv.gz
│   ├── STAD.tsv.gz
│   ├── TGCT.tsv.gz
│   ├── THCA.tsv.gz
│   ├── THYM.tsv.gz
│   ├── UCEC.tsv.gz
│   ├── UCS.tsv.gz
│   └── UVM.tsv.gz
├── pert_transcriptomes
│   ├── ENASFS
│   │   ├── delta_psi-EX.tsv.gz
│   │   └── log2_fold_change_tpm.tsv.gz
│   ├── ENCOREKD
│   │   ├── HepG2
│   │   │   ├── delta_psi-EX.tsv.gz
│   │   │   ├── delta_psi_rel-EX.tsv.gz
│   │   │   └── log2_fold_change_tpm.tsv.gz
│   │   └── K562
│   │       ├── delta_psi-EX.tsv.gz
│   │       ├── delta_psi_rel-EX.tsv.gz
│   │       └── log2_fold_change_tpm.tsv.gz
│   └── ENCOREKO
│       ├── HepG2
│       │   ├── delta_psi-EX.tsv.gz
│       │   ├── delta_psi_rel-EX.tsv.gz
│       │   └── log2_fold_change_tpm.tsv.gz
│       └── K562
│           ├── delta_psi-EX.tsv.gz
│           ├── delta_psi_rel-EX.tsv.gz
│           └── log2_fold_change_tpm.tsv.gz
├── phosphoproteomics
│   └── Hafner2019-log2_fold_changes.tsv.gz
├── ppi
│   └── STRINGDB.tsv.gz
├── references
│   └── splicing_factors
│       ├── splicing_factors-ensembl.txt
│       ├── splicing_factors-symbol.txt
│       └── splicing_factors.tsv
└── summary_stats
    ├── event_psi
    │   ├── ACC-EX.tsv.gz
    │   ├── BLCA-EX.tsv.gz
    │   ├── BRCA-EX.tsv.gz
    │   ├── CCLE-EX.tsv.gz
    │   ├── CESC-EX.tsv.gz
    │   ├── CHOL-EX.tsv.gz
    │   ├── COAD-EX.tsv.gz
    │   ├── DLBC-EX.tsv.gz
    │   ├── GBM-EX.tsv.gz
    │   ├── HNSC-EX.tsv.gz
    │   ├── KICH-EX.tsv.gz
    │   ├── KIRC-EX.tsv.gz
    │   ├── KIRP-EX.tsv.gz
    │   ├── LGG-EX.tsv.gz
    │   ├── LIHC-EX.tsv.gz
    │   ├── LUAD-EX.tsv.gz
    │   ├── LUSC-EX.tsv.gz
    │   ├── MESO-EX.tsv.gz
    │   ├── PAAD-EX.tsv.gz
    │   ├── PCPG-EX.tsv.gz
    │   ├── PRAD-EX.tsv.gz
    │   ├── READ-EX.tsv.gz
    │   ├── SARC-EX.tsv.gz
    │   ├── SKCM-EX.tsv.gz
    │   ├── TGCT-EX.tsv.gz
    │   ├── THCA-EX.tsv.gz
    │   ├── THYM-EX.tsv.gz
    │   ├── UCEC-EX.tsv.gz
    │   ├── UCS-EX.tsv.gz
    │   └── UVM-EX.tsv.gz
    ├── event_psi_imputed
    │   ├── ACC-EX.tsv.gz
    │   ├── BLCA-EX.tsv.gz
    │   ├── BRCA-EX.tsv.gz
    │   ├── CCLE-EX.tsv.gz
    │   ├── CESC-EX.tsv.gz
    │   ├── CHOL-EX.tsv.gz
    │   ├── COAD-EX.tsv.gz
    │   ├── DLBC-EX.tsv.gz
    │   ├── ESCA-EX.tsv.gz
    │   ├── GBM-EX.tsv.gz
    │   ├── HNSC-EX.tsv.gz
    │   ├── KICH-EX.tsv.gz
    │   ├── KIRC-EX.tsv.gz
    │   ├── KIRP-EX.tsv.gz
    │   ├── LAML-EX.tsv.gz
    │   ├── LGG-EX.tsv.gz
    │   ├── LIHC-EX.tsv.gz
    │   ├── LUAD-EX.tsv.gz
    │   ├── LUSC-EX.tsv.gz
    │   ├── MESO-EX.tsv.gz
    │   ├── OV-EX.tsv.gz
    │   ├── PAAD-EX.tsv.gz
    │   ├── PCPG-EX.tsv.gz
    │   ├── PRAD-EX.tsv.gz
    │   ├── READ-EX.tsv.gz
    │   ├── SARC-EX.tsv.gz
    │   ├── SKCM-EX.tsv.gz
    │   ├── STAD-EX.tsv.gz
    │   ├── TGCT-EX.tsv.gz
    │   ├── THCA-EX.tsv.gz
    │   ├── THYM-EX.tsv.gz
    │   ├── UCEC-EX.tsv.gz
    │   ├── UCS-EX.tsv.gz
    │   └── UVM-EX.tsv.gz
    └── genexpr_tpm
        ├── ACC.tsv.gz
        ├── BLCA.tsv.gz
        ├── BRCA.tsv.gz
        ├── CCLE.tsv.gz
        ├── CESC.tsv.gz
        ├── CHOL.tsv.gz
        ├── COAD.tsv.gz
        ├── DLBC.tsv.gz
        ├── ESCA.tsv.gz
        ├── GBM.tsv.gz
        ├── HNSC.tsv.gz
        ├── KICH.tsv.gz
        ├── KIRC.tsv.gz
        ├── KIRP.tsv.gz
        ├── LAML.tsv.gz
        ├── LGG.tsv.gz
        ├── LIHC.tsv.gz
        ├── LUAD.tsv.gz
        ├── LUSC.tsv.gz
        ├── MESO.tsv.gz
        ├── OV.tsv.gz
        ├── PAAD.tsv.gz
        ├── PCPG.tsv.gz
        ├── PRAD.tsv.gz
        ├── READ.tsv.gz
        ├── SARC.tsv.gz
        ├── SKCM.tsv.gz
        ├── STAD.tsv.gz
        ├── TGCT.tsv.gz
        ├── THCA.tsv.gz
        ├── THYM.tsv.gz
        ├── UCEC.tsv.gz
        ├── UCS.tsv.gz
        └── UVM.tsv.gz
```