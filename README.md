# Exon inclusion signatures enable accurate estimation of splicing factor activity

Repurposing VIPER to estimate differential splicing factor activity. This repository contains all the code to reproduce our analyses.

## Structure
- `config.yaml`: important paths used througout the workflows to be set by the user.
- `environment.yaml`: conda/mamba environment specification file.
- `data`: where data will be downloaded.
- `results`: where results will be outputed.
- `src`: project-wide scripts.
- `support`: supporting files which creation is hard to automate.
- `workflows`
    1. `obtain_data` ([Details](https://github.com/MiqG/publication_viper_splicing/tree/main/workflows/obtain_data)): raw data download
    
    2. `preprocess_data` ([Details](https://github.com/MiqG/publication_viper_splicing/tree/main/workflows/preprocess_data)):
        - quantification of gene expression and exon inclusion
    
    3. `regulon_inference` ([Details](https://github.com/MiqG/publication_viper_splicing/tree/main/workflows/regulon_inference)): benchmark of SF-exon network inference approaches
        - experimentally derived regulons
        - aracne regulons
        - regulon evaluation
        - regulon eda
    
    4. `sf_activity_validation` ([Details](https://github.com/MiqG/publication_viper_splicing/tree/main/workflows/sf_activity_validation)): validation of splicing factor activity estimation
        - protein perturbations
        - combinatorial perturbations
        - SF3b complex
    
    5. `cancer_splicing_program` ([Details](https://github.com/MiqG/publication_viper_splicing/tree/main/workflows/cancer_splicing_program)): recurrent cancer-driver-like splicing program
        - definition using TCGA
        - patient prognosis
        - hallmarks
            - proliferation
            - immune evasion
        - tumorigenesis
            - from fibroblasts to cancer
    
    6. `prepare_submission` ([Details](https://github.com/MiqG/publication_viper_splicing/tree/main/workflows/prepare_submission)): prepare supplementary tables


## Installation and requirements
### Environment
We recommend installing mamba miniforge rather than conda in your corresponding OS: https://github.com/conda-forge/miniforge.
```shell
mamba env create -f environment.yaml
```

### `vast-tools` (requires manual installation)
Based on https://github.com/vastgroup/vast-tools?tab=readme-ov-file#installation.

Clone vast-tools repository in your desired path:
```{shell}
cd ~/repositories # or your desired path
git clone https://github.com/vastgroup/vast-tools.git
```

Download VastDB Homo sapiens (Hs2) genome assembly in your desired path:
```{shell}
cd ~/projects/publication_viper/data/raw/VastDB/assemblies # or your desired path
wget https://vastdb.crg.eu/libs/vastdb.hs2.20.12.19.tar.gz
tar -xvzf vastdb.hs2.20.12.19.tar.gz
```

Update `config.yaml` file accordingly with the path to the vast-tools directory and to the VastDB genome assembly.

#### TCGA restricted data access
Obtain a GDC token file to download data (see: https://docs.gdc.cancer.gov/Data/Data_Security/Data_Security/) and place its path in `config.yaml`.


## Usage
All workflows were written as `snakemake` pipelines. To execute each workflow:

1. activate the project's environment
```{shell}
mamba activate publication_viper_splicing
```

2. use the following command format:
```{shell}
snakemake -s <workflow_name>.smk --cores <number_of_cores>
```

In case you want to run the workflows on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) to adapt the command according to your job submission scheduler.

## Issues
Please, report any issues here: https://github.com/MiqG/publication_viper_splicing/issues

## Authors
- [Miquel Anglada Girotto](https://orcid.org/0000-0003-1885-8649)
- [Daniel F. Moackley](https://orcid.org/0000-0002-3279-4189)
- [Chaolin Zhang](https://orcid.org/0000-0002-8310-7537)
- [Samuel Miravet-Verde](https://orcid.org/0000-0002-1542-5912)
- [Andrea Califano](https://orcid.org/0000-0003-4742-3679)
- [Luis Serrano](https://orcid.org/0000-0002-5276-1392)

## Citation
Anglada-Girotto, M., Moakley, D. F., Zhang, C., Miravet-Verde, S., Califano, A., & Serrano, L. (2024). Disentangling the splicing factor programs underlying complex molecular phenotypes. bioRxiv, 2024-06.