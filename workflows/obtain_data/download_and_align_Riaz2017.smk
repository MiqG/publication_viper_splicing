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
DATASET_DIR = os.path.join(RAW_DIR,"articles","Riaz2017")

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJNA356761-Riaz2017.tsv'))
metadata = metadata.loc[metadata["library_source"]=="TRANSCRIPTOMIC"]
metadata_paired = metadata.loc[metadata["library_layout"]=="PAIRED"]
metadata_single = metadata.loc[metadata["library_layout"]=="SINGLE"]

## URLS to download
URLS_PAIRED = metadata_paired['fastq_ftp'].astype("str").str.split(';').str[0].apply(os.path.dirname).to_list()
URLS_PAIRED = {os.path.basename(url): url for url in URLS_PAIRED}
SAMPLES_PAIRED = list(URLS_PAIRED.keys())
URLS_SINGLE = metadata_single['fastq_ftp'].astype("str").str.split(';').str[0].apply(os.path.dirname).to_list()
URLS_SINGLE = {os.path.basename(url): url for url in URLS_SINGLE}
SAMPLES_SINGLE = list(URLS_SINGLE.keys())

N_SAMPLES = len(SAMPLES_PAIRED) + len(SAMPLES_SINGLE)

## fastq sizes
SIZES_PAIRED = metadata_paired.set_index("run_accession")["fastq_bytes"].astype("str").str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()
SIZES_SINGLE = metadata_single.set_index("run_accession")["fastq_bytes"].astype("str").str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()

SIZE_THRESH = 5e9

##### RULES #####
rule all:
    input:
        # download metadata
        os.path.join(DATASET_DIR,"metadata.tsv"),
        
        # download .fastq
        ## paired single end
        expand(os.path.join(DATASET_DIR,'fastqs','.done','{sample}_{end}_paired'), end=["1","2"], sample=SAMPLES_PAIRED),
        ## single end
        expand(os.path.join(DATASET_DIR,'fastqs','.done','{sample}_{end}_single'), end=["1"], sample=SAMPLES_SINGLE),
        
        # quantify event PSI and gene expression TPM
        ## paired end
        expand(os.path.join(DATASET_DIR,'vast_out','.done','{sample}_paired'), sample=SAMPLES_PAIRED),
        ## single end
        expand(os.path.join(DATASET_DIR,'vast_out','.done','{sample}_single'), sample=SAMPLES_SINGLE),

        # combine
        os.path.join(DATASET_DIR,'vast_out','.done/vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES),

        # tidy PSI
        os.path.join(DATASET_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        '.done/Riaz2017.done'

        
rule download_metadata:
    params:
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA356761&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created,bam_ftp,bam_bytes,bam_md5&format=tsv&download=true&limit=0"
    output:
        os.path.join(DATASET_DIR,"metadata.tsv")
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
        
        
rule download_paired:
    params:
        sample = '{sample}',
        end = "{end}",
        url = lambda wildcards: URLS_PAIRED[wildcards.sample],
        fastqs_dir = os.path.join(DATASET_DIR,'fastqs')
    output:
        download_done = os.path.join(DATASET_DIR,'fastqs','.done','{sample}_{end}_paired')
    threads: 1
    resources:
        runtime = 7200, # 2h
        memory = 2
    shell:
        """
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
        
        
rule download_single:
    params:
        sample = '{sample}',
        end = "{end}",
        url = lambda wildcards: URLS_SINGLE[wildcards.sample],
        fastqs_dir = os.path.join(DATASET_DIR,'fastqs'),
        bin_dir=VAST_TOOLS_DIR
    output:
        download_done = os.path.join(DATASET_DIR,'fastqs','.done','{sample}_{end}_single')
    threads: 1
    resources:
        runtime = 7200, # 2h
        memory = 2
    shell:
        """
        # download
        echo "Downloading {params.sample}..."
        
        wget --user-agent="Chrome" \
             --no-check-certificate \
             {params.url}/{params.sample}.fastq.gz \
             -O {params.fastqs_dir}/{params.sample}_{params.end}.fastq.gz
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule align_paired:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(DATASET_DIR,'fastqs'),
        bin_dir=VAST_TOOLS_DIR,
        vast_out = directory(os.path.join(DATASET_DIR,'vast_out','{sample}'))
    input:
        dbDir = VASTDB_DIR,
        download_done = [os.path.join(DATASET_DIR,'fastqs','.done','{sample}_{end}_paired').format(end=end, sample='{sample}') for end in ["1","2"]]
    output:
        align_done = touch(os.path.join(DATASET_DIR,'vast_out','.done','{sample}_paired'))
    threads: 16
    resources:
        runtime = lambda wildcards: 86400 if SIZES_PAIRED[wildcards.sample]>SIZE_THRESH else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    shell:
        """
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
        
        
rule align_single:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(DATASET_DIR,'fastqs'),
        bin_dir=VAST_TOOLS_DIR,
        vast_out = directory(os.path.join(DATASET_DIR,'vast_out','{sample}'))
    input:
        dbDir = VASTDB_DIR,
        download_done = [os.path.join(DATASET_DIR,'fastqs','.done','{sample}_{end}_single').format(end=end, sample='{sample}') for end in ["1"]]
    output:
        align_done = touch(os.path.join(DATASET_DIR,'vast_out','.done','{sample}_single'))
    threads: 16
    resources:
        runtime = lambda wildcards: 86400 if SIZES_SINGLE[wildcards.sample]>SIZE_THRESH else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    shell:
        """
        # align paired reads
        echo "Aligning {params.sample}..."
        {params.bin_dir}/vast-tools align \
                    {params.fastqs_dir}/{params.sample}_1.fastq.gz \
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
        done = [os.path.join(DATASET_DIR,'vast_out','.done','{sample}_paired').format(sample=sample) for sample in SAMPLES_PAIRED] + [os.path.join(DATASET_DIR,'vast_out','.done','{sample}_single').format(sample=sample) for sample in SAMPLES_SINGLE],
        dbDir = VASTDB_DIR
    output:
        touch(os.path.join(DATASET_DIR,'vast_out','.done','vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES)),
        tpm = os.path.join(DATASET_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES),
        psi = os.path.join(DATASET_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES)
    params:
        bin_dir=VAST_TOOLS_DIR,
        folder = os.path.join(DATASET_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 86400, # 24h
        memory = 80
    shell:
        """
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
        os.path.join(DATASET_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES) ## combined table
    output:
        touch('.done/Riaz2017.done'),
        tidy = os.path.join(DATASET_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir=VAST_TOOLS_DIR
    threads: 1
    resources:
        runtime = 3600*12, # 12h
        memory = 20
    shell:
        """
        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """
