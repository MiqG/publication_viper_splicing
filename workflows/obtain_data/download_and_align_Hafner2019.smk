import os
import pandas as pd
import numpy as np

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
DATASET_DIR = os.path.join(RAW_DIR,"articles","Hafner2019")
VASTDB_DIR = os.path.join(RAW_DIR,'VastDB')

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJNA387311_PRJNA514039-Hafner2019.tsv'))
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
        '.done/Hafner2019.done'

        
rule download_metadata:
    params:
        url = "https://www.ebi.ac.uk/ena/portal/api/search?query=study_accession=%22PRJNA514039%22%20OR%20study_accession=%22PRJNA387311%22&result=read_run&fields=run_accession,experiment_title,study_accession,aligned,altitude,assembly_quality,assembly_software,bam_aspera,bam_bytes,bam_ftp,bam_galaxy,bam_md5,base_count,binning_software,bio_material,bisulfite_protocol,broad_scale_environmental_context,broker_name,cage_protocol,cell_line,cell_type,center_name,checklist,chip_ab_provider,chip_protocol,chip_target,collected_by,collection_date,collection_date_end,collection_date_start,completeness_score,contamination_score,control_experiment,country,cultivar,culture_collection,depth,description,dev_stage,dnase_protocol,ecotype,elevation,environment_biome,environment_feature,environment_material,environmental_medium,environmental_sample,experiment_accession,experiment_alias,experiment_target,experimental_factor,experimental_protocol,extraction_protocol,faang_library_selection,fastq_aspera,fastq_bytes,fastq_ftp,fastq_galaxy,fastq_md5,file_location,first_created,first_public,germline,hi_c_protocol,host,host_body_site,host_genotype,host_gravidity,host_growth_conditions,host_phenotype,host_scientific_name,host_sex,host_status,host_tax_id,identified_by,instrument_model,instrument_platform,investigation_type,isolate,isolation_source,last_updated,lat,library_construction_protocol,library_gen_protocol,library_layout,library_max_fragment_size,library_min_fragment_size,library_name,library_pcr_isolation_protocol,library_prep_date,library_prep_date_format,library_prep_latitude,library_prep_location,library_prep_longitude,library_selection,library_source,library_strategy,local_environmental_context,location,location_end,location_start,lon,marine_region,mating_type,ncbi_reporting_standard,nominal_length,nominal_sdev,pcr_isolation_protocol,ph,project_name,protocol_label,read_count,read_strand,restriction_enzyme,restriction_enzyme_target_sequence,restriction_site,rna_integrity_num,rna_prep_3_protocol,rna_prep_5_protocol,rna_purity_230_ratio,rna_purity_280_ratio,rt_prep_protocol,run_alias,run_date,salinity,sample_accession,sample_alias,sample_capture_status,sample_collection,sample_description,sample_material,sample_prep_interval,sample_prep_interval_units,sample_storage,sample_storage_processing,sample_title,sampling_campaign,sampling_platform,sampling_site,scientific_name,secondary_project,secondary_sample_accession,secondary_study_accession,sequencing_date,sequencing_date_format,sequencing_location,sequencing_longitude,sequencing_method,sequencing_primer_catalog,sequencing_primer_lot,sequencing_primer_provider,serotype,serovar,sex,specimen_voucher,sra_aspera,sra_bytes,sra_ftp,sra_galaxy,sra_md5,status,strain,study_alias,study_title,sub_species,sub_strain,submission_accession,submission_tool,submitted_aspera,submitted_bytes,submitted_format,submitted_ftp,submitted_galaxy,submitted_host_sex,submitted_md5,submitted_read_type,tag,target_gene,tax_id,taxonomic_classification,taxonomic_identity_marker,temperature,tissue_lib,tissue_type,transposase_protocol,variety&limit=0&download=true&format=tsv"
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
        fastqs_dir = os.path.join(DATASET_DIR,'fastqs'),
        bin_dir="~/repositories/vast-tools/"
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
        bin_dir="~/repositories/vast-tools/"
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
        bin_dir="~/repositories/vast-tools/",
        vast_out = directory(os.path.join(DATASET_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
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
        bin_dir="~/repositories/vast-tools/",
        vast_out = directory(os.path.join(DATASET_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
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
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(DATASET_DIR,'vast_out','.done','vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES)),
        tpm = os.path.join(DATASET_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES),
        psi = os.path.join(DATASET_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES)
    params:
        bin_dir="~/repositories/vast-tools/",
        folder = os.path.join(DATASET_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 86400, # 24h
        memory = 20
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
        touch('.done/Hafner2019.done'),
        tidy = os.path.join(DATASET_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir="~/repositories/vast-tools/"
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
