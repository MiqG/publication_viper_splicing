import os
import pandas as pd
import numpy as np

# unpack config
configfile: "../../config.yaml"
PATHS = config["PATHS"]
PATHS = {k: os.path.expanduser(v) for k, v in PATHS.items()} # make sure to have full paths (without ~)
VASTDB_DIR = PATHS["VAST_TOOLS"]["VASTDB"]
VAST_TOOLS_DIR = PATHS["VAST_TOOLS"]["BIN"]

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
ARTICLES_DIR = os.path.join(RAW_DIR,'articles')
ARTICLE_DIR = os.path.join(ARTICLES_DIR,'Nijhuis2020')

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJNA673205-Nijhuis2020.tsv'))
metadata = metadata.loc[metadata["library_source"]=="TRANSCRIPTOMIC"]
metadata = metadata.loc[metadata["fastq_ftp"].str.split(';').apply(len)==2]

## URLS to download
URLS = metadata['fastq_ftp'].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS = {os.path.basename(url): url for url in URLS}
SAMPLES = list(URLS.keys())
N_SAMPLES = len(SAMPLES)

## fastq sizes
SIZES = metadata.set_index("run_accession")["fastq_bytes"].astype(str).str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()
SIZE_THRESH = 5e9

ENDS = ["1","2"]

##### RULES #####
rule all:
    input:
        # PROTEOMICS
        # download supplementary tables
        os.path.join(ARTICLE_DIR,"supplementary_data","rnaseq_geo"),
        os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S1-12_search_results_txt_folder"),
        os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S13-16_search_results_txt_folder"),
        
        # preprocess supplementary tables
        os.path.join(ARTICLE_DIR,"supplementary_data","proteomics_lfq_intensity.tsv.gz"),
        
        # TRANSCRIPTOMICS
        # download metadata
        os.path.join(ARTICLE_DIR,"metadata.tsv"),
        
        # download .fastq
        expand(os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}'), end=ENDS, sample=SAMPLES),
        
        # quantify event PSI and gene expression TPM
        expand(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES),
        
        # combine
        os.path.join(ARTICLE_DIR,'vast_out','.done/vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES),
        
        # tidy PSI
        os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        '.done/Nijhuis2020.done'
        
        
        
rule download_supdata:
    params:
        rnaseq_geo = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE160446&format=file",
        proteomics1 = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/01/PXD022164/MQ_S1-12_search_results_txt_folder.zip",
        proteomics2 = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/01/PXD022164/MQ_S13-16_search_results_txt_folder.zip"
    output:
        rnaseq_geo = directory(os.path.join(ARTICLE_DIR,"supplementary_data","rnaseq_geo")),
        proteomics1 = directory(os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S1-12_search_results_txt_folder")),
        proteomics2 = directory(os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S13-16_search_results_txt_folder"))
    shell:
        """
        set -eo pipefail

        # rna seq counts
        wget '{params.rnaseq_geo}' -O {output.rnaseq_geo}
        # proteomics
        wget '{params.proteomics1}' -O {output.proteomics1}.zip
        wget '{params.proteomics2}' -O {output.proteomics2}.zip

        # prep RNA seq
        mkdir rnaseq_geo
        tar xvf {output.rnaseq_geo} -C rnaseq_geo/

        # prep proteomics
        unzip {output.proteomics1}.zip -d {output.proteomics1}
        unzip {output.proteomics2}.zip -d {output.proteomics2}

        echo "Done!"
        """

        
rule prep_supdata:
    input:
        proteomics1 = os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S1-12_search_results_txt_folder"),
        proteomics2 = os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S13-16_search_results_txt_folder")
    output:
        proteomics = os.path.join(ARTICLE_DIR,"supplementary_data","proteomics_lfq_intensity.tsv.gz")
    run:
        import pandas as pd
        import numpy as np
        
        # load
        ## proteomics
        proteomics1 = pd.read_table(os.path.join(input.proteomics1,"txt","proteinGroups.txt"), low_memory=False)
        proteomics2 = pd.read_table(os.path.join(input.proteomics2,"txt","proteinGroups.txt"), low_memory=False)

        # prep proteomics
        common_cols = set(proteomics1.columns).intersection(proteomics2.columns)
        samples1 = set(proteomics1.columns) - set(common_cols)
        samples2 = set(proteomics2.columns) - set(common_cols)
        proteomics1 = proteomics1[["Gene names"]+list(samples1)]
        proteomics2 = proteomics2[["Gene names"]+list(samples2)]
        proteomics1["Gene names"] = proteomics1["Gene names"].str.split(";")
        proteomics2["Gene names"] = proteomics2["Gene names"].str.split(";")
        proteomics = pd.merge(
            proteomics1.explode("Gene names"),
            proteomics2.explode("Gene names"),
            on="Gene names", how="outer"
        )
        cols_oi = [c for c in proteomics.columns if "LFQ" in c]
        proteomics = proteomics[["Gene names"] + sorted(cols_oi)].copy()
        proteomics = proteomics.loc[~proteomics["Gene names"].isnull()].set_index("Gene names").copy()
        proteomics.index.name = "GENE"
        
        # sample identifiers
        proteomics.columns = [c.replace("LFQ intensity ","sample") for c in proteomics.columns]
        
        # save
        proteomics.reset_index().to_csv(output.proteomics, **SAVE_PARAMS)
        
        print("Done!")

        
rule download_metadata:
    params:
        url = "https://www.ebi.ac.uk/ena/portal/api/search?query=tax_eq(9606)%20AND%20study_accession=%22PRJNA673205%22&result=read_run&fields=accession,altitude,assembly_quality,assembly_software,base_count,binning_software,bio_material,broker_name,cell_line,cell_type,center_name,checklist,collected_by,collection_date,collection_date_submitted,completeness_score,contamination_score,country,cram_index_aspera,cram_index_ftp,cram_index_galaxy,cultivar,culture_collection,depth,description,dev_stage,ecotype,elevation,environment_biome,environment_feature,environment_material,environmental_package,environmental_sample,experiment_accession,experiment_alias,experiment_title,experimental_factor,fastq_aspera,fastq_bytes,fastq_ftp,fastq_galaxy,fastq_md5,first_created,first_public,germline,host,host_body_site,host_genotype,host_gravidity,host_growth_conditions,host_phenotype,host_sex,host_status,host_tax_id,identified_by,instrument_model,instrument_platform,investigation_type,isolate,isolation_source,last_updated,lat,library_construction_protocol,library_layout,library_name,library_selection,library_source,library_strategy,location,lon,mating_type,nominal_length,nominal_sdev,parent_study,ph,project_name,protocol_label,read_count,run_accession,run_alias,salinity,sample_accession,sample_alias,sample_capture_status,sample_collection,sample_description,sample_material,sample_title,sampling_campaign,sampling_platform,sampling_site,scientific_name,secondary_sample_accession,secondary_study_accession,sequencing_method,serotype,serovar,sex,specimen_voucher,sra_aspera,sra_bytes,sra_ftp,sra_galaxy,sra_md5,strain,study_accession,study_alias,study_title,sub_species,sub_strain,submission_accession,submission_tool,submitted_aspera,submitted_bytes,submitted_format,submitted_ftp,submitted_galaxy,submitted_host_sex,submitted_md5,submitted_sex,target_gene,tax_id,taxonomic_classification,taxonomic_identity_marker,temperature,tissue_lib,tissue_type,variety&limit=0&download=true&format=tsv"
    output:
        os.path.join(RAW_DIR,"articles","Nijhuis2020","metadata.tsv")
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 2
    shell:
        """
        set -eo pipefail
        
        wget --user-agent="Chrome" \
             --no-check-certificate \
             "{params.url}" \
             -O {output}
             
        echo "Done!"
        """
        
        
rule download:
    params:
        sample = '{sample}',
        end = "{end}",
        url = lambda wildcards: URLS[wildcards.sample],
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs')
    output:
        download_done = os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}')
    threads: 1
    resources:
        runtime = 7200, # 2h
        memory = 2
    shell:
        """
        set -eo pipefail
        
        # download
        echo "Downloading {params.sample}..."
        
        wget --user-agent="Chrome" \
             --no-check-certificate \
             {params.url}/{params.sample}_{params.end}.fastq.gz \
             -O {params.fastqs_dir}/{params.sample}_{params.end}.fastq.gz
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule align:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir=VAST_TOOLS_DIR,
        vast_out = directory(os.path.join(ARTICLE_DIR,'vast_out','{sample}'))
    input:
        dbDir = VASTDB_DIR,
        download_done = [os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS]
    output:
        align_done = touch(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'))
    threads: 12
    resources:
        runtime = lambda wildcards: 86400 if SIZES[wildcards.sample]>SIZE_THRESH else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    shell:
        """
        set -eo pipefail
        
        # align paired reads
        echo "Aligning {params.sample}..."
        {params.bin_dir}/vast-tools align \
                    {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                    {params.fastqs_dir}/{params.sample}_2.fastq.gz \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --expr \
                    --EEJ_counts \
                    --cores {threads} \
                    --output {params.vast_out}
        echo "Finished aligning {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule vasttools_combine:
    input:
        done = [os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}').format(sample=sample) for sample in SAMPLES],
        dbDir = VASTDB_DIR
    output:
        touch(os.path.join(ARTICLE_DIR,'vast_out','.done','vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES)),
        tpm = os.path.join(ARTICLE_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES),
        psi = os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES)
    params:
        bin_dir=VAST_TOOLS_DIR,
        folder = os.path.join(ARTICLE_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 86400, # 24h
        memory = 20
    shell:
        """
        set -eo pipefail
        
        # group results
        echo "Grouping results..."
        mkdir -p {params.folder}/to_combine
        ln -s {params.folder}/*/to_combine/* {params.folder}/to_combine/
        mkdir -p {params.folder}/expr_out
        ln -s {params.folder}/*/expr_out/* {params.folder}/expr_out/
        
        # combine runs
        echo "Combining runs..."
        {params.bin_dir}/vast-tools combine \
                    --cores {threads} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --keep_raw_reads \
                    --keep_raw_incl \
                    --output {params.folder} \
                    --TPM
        
        # compress outputs
        echo "Compressing outputs..."
        gzip -f {params.folder}/raw_incl/*
        gzip -f {params.folder}/raw_reads/*
        gzip -f {params.folder}/*.tab
        
        # remove grouped results
        echo "Removing grouped results..."
        rm -r {params.folder}/to_combine
        rm -r {params.folder}/expr_out
        
        echo "Done!"
        """
    
    
rule vasttools_tidy:
    input:
        os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES) ## combined table
    output:
        touch('.done/Nijhuis2020.done'),
        tidy = os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir=VAST_TOOLS_DIR
    threads: 1
    resources:
        runtime = 3600*12, # 12h
        memory = 20
    shell:
        """
        set -eo pipefail
        
        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """