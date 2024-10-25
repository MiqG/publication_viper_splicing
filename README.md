# Disentangling the splicing factor programs underlying complex molecular phenotypes

## Overview
Repurposing VIPER to estimate differential splicing factor activity.

## Structure
- `data`: where data will be downloaded.
- `results`: where results will be outputed.
- `src`: project-wide scripts.
- `support`: supporting files which creation is hard to automate.
- `workflows`
    1. `obtain_data`: data download
    2. `preprocess_data`:
        - quantification of gene expression and exon inclusion
    3. `regulon_inference`: benchmark of SF-exon network inference approaches
        - experimentally derived regulons
        - aracne regulons
        - regulon evaluation
        - regulon eda
    4. `sf_activity_validation`: validation of splicing factor activity estimation
        - protein perturbations
        - combinatorial perturbations
        - SF3b complex
    5. `cancer_splicing_program`: recurrent cancer-driver-like splicing program
        - definition using TCGA
        - patient prognosis
        - hallmarks
            - proliferation
            - immune evasion
        - tumorigenesis
            - from fibroblasts to cancer
    6. `prepare_submission`: prepare supplementary tables

## Usage
### 1. Install environment
We recommend installing mamba miniforge rather than conda in your corresponding OS: https://github.com/conda-forge/miniforge.
```shell
mamba env create -f environment.yml
```
This may take a while

### 2. Run snakemake with the `publication_viper_splicing` environment
All workflows were written as `snakemake` pipelines. Refer to each specific workflow to run them.

```shell
# run locally
snakemake -s {snakefile path} --cores {n_cpus} --use-conda

# cluster job submission example (SGE scheduler)
snakemake --cluster "qsub -cwd -pe smp {threads} -l virtual_free={resources.memory}G,h_rt={resources.runtime} -q long,short" --use-conda --jobs 30 -s {snakefile path}
```

## Authors
- [Miquel Anglada Girotto](https://orcid.org/0000-0003-1885-8649)
- [Daniel F. Moackley](https://orcid.org/0000-0002-3279-4189)
- [Chaolin Zhang](https://orcid.org/0000-0002-8310-7537)
- [Samuel Miravet-Verde](https://orcid.org/0000-0002-1542-5912)
- [Andrea Califano](https://orcid.org/0000-0003-4742-3679)
- [Luis Serrano](https://orcid.org/0000-0002-5276-1392)

## License

