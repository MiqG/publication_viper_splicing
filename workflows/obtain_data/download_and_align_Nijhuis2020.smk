import os
import pandas as pd
import numpy as np

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
ARTICLES_DIR = os.path.join(RAW_DIR,'articles')
ARTICLE_DIR = os.path.join(ARTICLES_DIR,'Nijhuis2020')
VASTDB_DIR = os.path.join(RAW_DIR,'VastDB')

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJNA673205-Nijhuis2020.tsv'))
metadata = metadata.loc[metadata["library_source"]=="TRANSCRIPTOMIC"]

## URLS to download
URLS = metadata['fastq_ftp'].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS = {os.path.basename(url): url for url in URLS}
SAMPLES = list(URLS.keys())
N_SAMPLES = len(SAMPLES)

## fastq sizes
SIZES = metadata.set_index("run_accession")["fastq_bytes"].astype(str).str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()
SIZE_THRESH = 5e9

##### RULES #####
rule all:
    input:
        # download supplementary tables
        os.path.join(ARTICLE_DIR,"supplementary_data","rnaseq_geo"),
        os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S1-12_search_results_txt_folder"),
        os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S13-16_search_results_txt_folder"),
        
        # download metadata
        os.path.join(ARTICLE_DIR,"metadata.tsv"),
        
        # download .fastq
        expand(os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}'), end=["1"], sample=SAMPLES),
        
        # quantify event PSI and gene expression TPM
        expand(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES),
        
        # combine
        #os.path.join(ARTICLE_DIR,'vast_out','.done/vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES),
        
        # tidy PSI
        #os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        #'.done/Nijhuis2020.done'


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
        wget "{params.rnaseq_geo}" -O {output.rnaseq_geo}
        # proteomics
        wget "{params.proteomics1}" -O {output.proteomics1}.zip
        wget "{params.proteomics2}" -O {output.proteomics2}.zip

        # prep RNA seq
        mkdir rnaseq_geo
        tar xvf {} -C rnaseq_geo/

        # prep proteomics
        unzip {output.proteomics1}.zip -d {output.proteomics1}
        unzip {output.proteomics2}.zip -d {output.proteomics2}

        echo "Done!"
        """

rule prep_supdata:
    input:
        rnaseq_geo = os.path.join(ARTICLE_DIR,"supplementary_data","rnaseq_geo"),
        proteomics1 = os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S1-12_search_results_txt_folder"),
        proteomics2 = os.path.join(ARTICLE_DIR,"supplementary_data","MQ_S13-16_search_results_txt_folder"),
        gene_annotation = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        transcriptomics = ,
        proteomics = 
    run:
        import pandas as pd
        import numpy as np
        
        # load
        ## transcriptomics
        transcriptomics = pd.concat([
            pd.read_table(os.path.join(input.rnaseq_geo,f), index_col=0) for f in os.listdir(input.rnaseq_geo)
        ], axis=1)
        ## proteomics
        proteomics1 = pd.read_table(os.path.join(input.proteomics1,"txt","proteinGroups.txt"), low_memory=False)
        proteomics2 = pd.read_table(os.path.join(input.proteomics2,"txt","proteinGroups.txt"), low_memory=False)
        ## gene annotaion
        gene_annotation = pd.read_table(input.gene_annotation)
        
        # prep transcriptomics (translate Entrez to Ensembl)
        gene_annotation = gene_annotation.loc[~gene_annotation["NCBI Gene ID"].isnull()].copy()
        gene_annotation["NCBI Gene ID"] = gene_annotation["NCBI Gene ID"].astype(int)
        transcriptomics = pd.merge(
            transcriptomics.reset_index(), 
            gene_annotation[["NCBI Gene ID","Ensembl gene ID"]], 
            how="left", left_on="index", right_on="NCBI Gene ID"
        ).dropna().set_index("Ensembl gene ID").drop(columns=["index","NCBI Gene ID"])
        transcriptomics.index.name = "ID"
        
        ## log-transform raw counts to counts per million
        transcriptomics = np.log2(transcriptomics / transcriptomics.sum() * 10**6 + 1)
        #df = df - df[["S1.06hVC.sorted.bam","S4.16hVC.sorted.bam"]].mean(axis=1).values.reshape(-1,1)
        #df.loc["ENSG00000131051"].sort_values().hist()
        
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
        ## change sample names
        #proteomics.loc["RBM39"].sort_values().hist()
        
        # save
        transcriptomics.reset_index().to_csv(output.transcriptomics, **SAVE_PARAMS)
        proteomics.to_csv(output.proteomics, **SAVE_PARAMS)
        
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
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir="~/repositories/vast-tools/"
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
             {params.url}/{params.sample}.fastq.gz \
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
        bin_dir="~/repositories/vast-tools/",
        vast_out = directory(os.path.join(ARTICLE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = [os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ["1"]]
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
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(ARTICLE_DIR,'vast_out','.done','vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES)),
        tpm = os.path.join(ARTICLE_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES),
        psi = os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES)
    params:
        bin_dir="~/repositories/vast-tools/",
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
        bin_dir="~/repositories/vast-tools/"
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