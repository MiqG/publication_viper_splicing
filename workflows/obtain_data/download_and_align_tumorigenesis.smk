import os
import pandas as pd
import numpy as np

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
DATASET_DIR = os.path.join(RAW_DIR,"ENA","tumorigenesis")
VASTDB_DIR = os.path.join(RAW_DIR,'VastDB')

# parameters
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-tumorigenesis.tsv'))
metadata = metadata.loc[(metadata["library_source"]=="TRANSCRIPTOMIC")].copy()
metadata = metadata.loc[~metadata["run_accession"].isin(["SRR13291149","SRR12101841"])].copy() # no fastq file link found and very few reads
metadata_paired = metadata.loc[metadata["library_layout"]=="PAIRED"]
metadata_single = metadata.loc[metadata["library_layout"]=="SINGLE"]

## URLS to download
URLS_PAIRED = metadata_paired['fastq_ftp'].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS_PAIRED = {os.path.basename(url): url for url in URLS_PAIRED}
SAMPLES_PAIRED = list(URLS_PAIRED.keys())
URLS_SINGLE = metadata_single['fastq_ftp'].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS_SINGLE = {os.path.basename(url): url for url in URLS_SINGLE}
SAMPLES_SINGLE = list(URLS_SINGLE.keys())

## fastq sizes
SIZES_PAIRED = metadata_paired.set_index("run_accession")["fastq_bytes"].str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()
SIZES_SINGLE = metadata_single.set_index("run_accession")["fastq_bytes"].str.split(";").apply(lambda x: max(np.array(x, dtype=int))).to_dict()

SIZE_THRESH = 5e9

# metadata only for samples to merge
## get samples
idx_tomerge = metadata["study_accession"].isin(["PRJNA769402"])
SAMPLES_PAIRED_TOMERGE = metadata.loc[idx_tomerge & (metadata["library_layout"]=="PAIRED"),"run_accession"].tolist()
SAMPLES_SINGLE_TOMERGE = metadata.loc[idx_tomerge & (metadata["library_layout"]=="SINGLE"),"run_accession"].tolist()

SAMPLES_PAIRED_NOTMERGE = metadata.loc[~idx_tomerge & (metadata["library_layout"]=="PAIRED"),"run_accession"].tolist()
SAMPLES_SINGLE_NOTMERGE = metadata.loc[~idx_tomerge & (metadata["library_layout"]=="SINGLE"),"run_accession"].tolist()

## subset
metadata_to_merge = metadata.loc[idx_tomerge].copy()
metadata_to_merge["group_label"] = np.nan
## study PRJNA769402
idx = metadata_to_merge["study_accession"]=="PRJNA769402"
cols_oi = [
    "study_accession","cell_line_name","condition","pert_time",
    "pert_time_units","pert_concentration","pert_concentration_units","replicate"
]
metadata_to_merge.loc[idx,"group_label"] = metadata_to_merge[cols_oi].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)

metadata_to_merge["run_id"] = metadata_to_merge["run_accession"] + "_1"
metadata_to_merge = metadata_to_merge[["run_accession","run_id","group_label"]]
#metadata_to_merge.to_csv(os.path.join(SUPPORT_DIR,"ENA_filereport-tumorigenesis-low_read_count_groups.tsv"), index=False, header=True, sep="\t")

MERGE_GROUPS = {
    g: metadata_to_merge.loc[metadata_to_merge["group_label"]==g,"run_accession"].tolist() 
    for g in metadata_to_merge["group_label"].unique()
}

SAMPLES_TO_COMBINE = SAMPLES_PAIRED_NOTMERGE + SAMPLES_SINGLE_NOTMERGE + list(MERGE_GROUPS.keys())
N_SAMPLES = len(SAMPLES_TO_COMBINE)

##### RULES #####
rule all:
    input:
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
        
        # merge low coverage
        ## make metadata
        [os.path.join(DATASET_DIR,"vast_out","to_merge","metadata","{group}.tsv").format(group=g) for g in MERGE_GROUPS.keys()],
        ## merge
        expand(os.path.join(DATASET_DIR,"vast_out","to_merge",".done","{group}"), group=MERGE_GROUPS.keys()),
        
        # combine
        os.path.join(DATASET_DIR,'vast_out','.done/vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES),
        
        # tidy PSI
        os.path.join(DATASET_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        '.done/tumorigenesis.done'
        
        
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
        set -eo pipefail
        
        # download
        echo "Downloading {params.sample}..."
        
        nice wget --user-agent="Chrome" \
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
        set -eo pipefail
        
        # download
        echo "Downloading {params.sample}..."
        
        nice wget --user-agent="Chrome" \
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

        
rule groups_merge_low_coverage:
    input:
        done = [os.path.join(DATASET_DIR,'vast_out','.done','{sample}_paired').format(sample=sample) 
                for sample in SAMPLES_PAIRED_TOMERGE
               ] + [
                os.path.join(DATASET_DIR,'vast_out','.done','{sample}_single').format(sample=sample) 
                for sample in SAMPLES_SINGLE_TOMERGE
               ],
        metadata_to_merge = os.path.join(SUPPORT_DIR,"ENA_filereport-tumorigenesis-low_read_count_groups.tsv")
    output:
        [os.path.join(DATASET_DIR,"vast_out","to_merge","metadata","{group}.tsv").format(group=g) for g in MERGE_GROUPS.keys()]
    params:
        output_dir = os.path.join(DATASET_DIR,"vast_out","to_merge","metadata")
    resources:
        runtime = 3600 * 1,
        memory = 15
    run:
        import pandas as pd
        
        metadata_to_merge = pd.read_table(input.metadata_to_merge)
        output_dir = params.output_dir
        
        for g in metadata_to_merge["group_label"].unique():
            # subset
            idx = metadata_to_merge["group_label"]==g
            group_metadata = metadata_to_merge.loc[idx,["run_id","group_label"]]
            # save
            filename = os.path.join(output_dir,g+".tsv")
            group_metadata.to_csv(filename, sep="\t", index=False, header=False)
        
        print("Done!")
        
    
rule merge_low_coverage:
    input:
        group = os.path.join(DATASET_DIR,"vast_out","to_merge","metadata","{group}.tsv"),
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        merge_done = touch(os.path.join(DATASET_DIR,"vast_out","to_merge",".done","{group}"))
    params:
        bin_dir="~/repositories/vast-tools/",
        folder = os.path.join(DATASET_DIR,'vast_out'),
        group = "{group}",
        samples = lambda wildcards: MERGE_GROUPS[wildcards.group]
    threads: 1
    resources:
        runtime = 3600*1, # h
        memory = 10 # GB
    shell:
        """
        set -eo pipefail
        
        # create merge directory
        mkdir -p {params.folder}/to_merge/to_combine/
        mkdir -p {params.folder}/to_merge/expr_out/
        
        # link every sample to combine in the corresponding folders
        for SAMPLE in {params.samples}
        do
            ln -nsf {params.folder}/$SAMPLE/to_combine/* {params.folder}/to_merge/to_combine/
            ln -nsf {params.folder}/$SAMPLE/expr_out/* {params.folder}/to_merge/expr_out/
        done
        
        # merge samples
        {params.bin_dir}/vast-tools merge \
                    --groups {input.group} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --expr \
                    -o {params.folder}/to_merge/
                    
        # move merged files
        mkdir -p {params.folder}/merged/to_combine/
        mkdir -p {params.folder}/merged/expr_out/
        
        mv {params.folder}/to_merge/to_combine/{params.group}* {params.folder}/merged/to_combine/
        mv {params.folder}/to_merge/expr_out/{params.group}* {params.folder}/merged/expr_out/
        
        echo "Done!"
        """
    
    
rule vasttools_combine:
    input:
        done = [os.path.join(DATASET_DIR,'vast_out','.done','{sample}_paired').format(sample=sample) for sample in SAMPLES_PAIRED_NOTMERGE] + [os.path.join(DATASET_DIR,'vast_out','.done','{sample}_single').format(sample=sample) for sample in SAMPLES_SINGLE_NOTMERGE] + [os.path.join(DATASET_DIR,"vast_out","to_merge",".done","{sample}").format(sample=sample) for sample in list(MERGE_GROUPS.keys())],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(DATASET_DIR,'vast_out','.done','vasttools_combine-{n_samples}').format(n_samples=N_SAMPLES)),
        tpm = os.path.join(DATASET_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES),
        psi = os.path.join(DATASET_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=N_SAMPLES)
    params:
        bin_dir="~/repositories/vast-tools/",
        folder = os.path.join(DATASET_DIR,'vast_out'),
        samples_merged = list(MERGE_GROUPS.keys()),
        samples_notmerged = SAMPLES_SINGLE_NOTMERGE + SAMPLES_PAIRED_NOTMERGE
    threads: 16
    resources:
        runtime = 3600*72, # 72h
        memory = 150 # GB
    shell:
        """
        set -eo pipefail
        
        # group results
        echo "Grouping results..."
        mkdir -p {params.folder}/to_combine
        mkdir -p {params.folder}/expr_out
        
        # link every sample to combine in the corresponding folders
        for SAMPLE in {params.samples_merged}
        do
            ln -nsf {params.folder}/merged/to_combine/$SAMPLE* {params.folder}/to_combine/ 
            ln -nsf {params.folder}/merged/expr_out/$SAMPLE* {params.folder}/expr_out/ 
        done
        
        for SAMPLE in {params.samples_notmerged}
        do
            ln -nsf {params.folder}/$SAMPLE/to_combine/* {params.folder}/to_combine/
            ln -nsf {params.folder}/$SAMPLE/expr_out/* {params.folder}/expr_out/
        done
        
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
        touch('.done/tumorigenesis.done'),
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