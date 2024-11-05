## workflows
### Execute them in the following order
1. `obtain_data`: Download data
2. `preprocess_data`: Preprocess downloaded data for the project
3. `regulon_inference`: Inference of SF-exon networks
4. `sf_activity_validation`: Validation of SF activity estimation with VIPER and SF-exon networks
5. `cancer_splicing_program`: Identifying a recurrent cancer splicing program
6. `prepare_submission`: Prepare supplementary tables for submission

### Stucture of each workflow
Inside each workflow folder you'll find a `README.md` file explaining how to run the corresponding workflow as well as Snakefile(s) to run the workflows. In most cases there is also a `<workflow_dir>/scripts` folder containing helper scripts.

```{shell}
workflows/<workflow_dir>/
├── <workflow_snakefile>.smk
├── README.md
└── scripts
    └── <helper_script>
```