# Prepare supplementary tables for submission

Once all workflows have been run, this workflow generates the supplementary tables for the publication.

## Expected outputs
```{shell}
$ tree results/prepare_submission/
results/prepare_submission/
└── files
    ├── ENASFS.tsv.gz
    └── supplementary_tables
        ├── suptab01_consensus_list_splicing_factors.tsv.gz
        ├── suptab02_pert_datasets.tsv.gz
        ├── suptab03_cancer_splicing_programs.tsv.gz
        └── suptab04_program_targets_enrichments.tsv.gz
```