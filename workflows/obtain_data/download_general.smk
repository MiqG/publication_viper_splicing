import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        # download BIOMART human-mouse homologies
        os.path.join(RAW_DIR,"BIOMART","human2mouse.tsv"),
        
        # DepMap - DEMETER2
        os.path.join(RAW_DIR,"DepMap",'demeter2.zip'),
        
        # HGNC - gene ids
        os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz"),
        
        # VastDB
        os.path.join(RAW_DIR,"VastDB",'EXPRESSION_TABLE-hg38.tab.gz'),
        os.path.join(RAW_DIR,"VastDB",'PSI_TABLE-hg38.tab.gz'),
        os.path.join(RAW_DIR,"VastDB",'EVENT_INFO-hg38.tab.gz'),
        os.path.join(RAW_DIR,"VastDB",'EVENT_INFO-hg38_noseqs.tsv'),
        os.path.join(RAW_DIR,"VastDB",'PROT_IMPACT-hg38-v3.tab.gz'),
        os.path.join(RAW_DIR,"VastDB",'PROT_ISOFORMS-hg38-v3.tab.gz'),        
        os.path.join(RAW_DIR,"VastDB",'assemblies','Hs2'),
        
        # UCSCXena
        os.path.join(RAW_DIR,"UCSCXena","GDC","snv","GDC-PANCAN.mutect2_snv.tsv.gz"),
        os.path.join(RAW_DIR,"UCSCXena",'TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz'),
        os.path.join(RAW_DIR,"UCSCXena",'TCGA','phenotype','TCGASubtype.20170308.tsv.gz'),
        os.path.join(RAW_DIR,"UCSCXena",'TCGA','phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz'),
        
        # STRINGDB
        os.path.join(RAW_DIR,"STRINGDB",'9606.protein.links.full.v11.5.txt.gz'),
        os.path.join(RAW_DIR,"STRINGDB",'9606.protein.aliases.v11.5.txt.gz'),
        
        # harmonizome
        os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz")
        
        
rule download_biomart_orthologues:
    params:
        human2mouse = """http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /><Attribute name = "mmusculus_homolog_ensembl_gene" /><Attribute name = "mmusculus_homolog_associated_gene_name" /></Dataset></Query>"""
    output:
        human2mouse = os.path.join(RAW_DIR,"BIOMART","human2mouse.tsv")
    shell:
        """
        wget --user-agent="Chrome" --no-clobber --no-check-certificate '{params.human2mouse}' -O {output.human2mouse}        
        
        echo Done!
        """
        
rule download_depmap:
    params:
        demeter2 = 'https://ndownloader.figshare.com/articles/6025238/versions/6',
    output:
        demeter2 = os.path.join(RAW_DIR,"DepMap",'demeter2.zip'),
        readme = os.path.join(RAW_DIR,"DepMap",'README.md')
    shell:
        """
        # download
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.demeter2} -O {output.demeter2}
        
        # README
        echo "Downloaded on $(date)." > {output.readme}
        echo Done!
        """
        
rule download_hgnc:
    params:
        annot = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_ensembl_id&col=gd_aliases&col=gd_enz_ids&col=family.id&col=gd_app_name&col=gd_status&col=gd_locus_type&col=gd_locus_group&col=gd_prev_sym&col=gd_prev_name&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_date2app_or_res&col=gd_date_mod&col=gd_date_sym_change&col=gd_date_name_change&col=gd_pub_acc_ids&col=gd_pub_eg_id&col=gd_mgd_id&col=gd_other_ids&col=gd_other_ids_list&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=family.name&col=gd_ccds_ids&col=gd_vega_ids&col=gd_lsdb_links&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    output:
        annot = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    shell:
        """
        wget --user-agent="Chrome" --no-check-certificate "{params.annot}" -O {output.annot}
        
        gzip --force {output.annot}
        
        mv {output.annot}.gz {output.annot}
        """
        
rule download_vastdb_hg38:
    params:
        event_genexpr = "https://vastdb.crg.eu/downloads/hg38/EXPRESSION_TABLE-hg38.tab.gz",
        event_psi = 'https://vastdb.crg.eu/downloads/hg38/PSI_TABLE-hg38.tab.gz',
        event_info = 'https://vastdb.crg.eu/downloads/hg38/EVENT_INFO-hg38.tab.gz',
        event_impact = 'https://vastdb.crg.eu/downloads/hg38/PROT_IMPACT-hg38-v3.tab.gz',
        protein_isoforms = 'https://vastdb.crg.eu/downloads/hg38/PROT_ISOFORMS-hg38.tab.gz',
        cols_todrop = 'Seq_C1,Seq_A,Seq_C2'
    output:
        event_genexpr = os.path.join(RAW_DIR,"VastDB",'EXPRESSION_TABLE-hg38.tab.gz'),
        event_psi = os.path.join(RAW_DIR,"VastDB",'PSI_TABLE-hg38.tab.gz'),
        event_info = os.path.join(RAW_DIR,"VastDB",'EVENT_INFO-hg38.tab.gz'),
        event_info_clean = os.path.join(RAW_DIR,"VastDB",'EVENT_INFO-hg38_noseqs.tsv'),
        event_impact = os.path.join(RAW_DIR,"VastDB",'PROT_IMPACT-hg38-v3.tab.gz'),
        protein_isoforms = os.path.join(RAW_DIR,"VastDB",'PROT_ISOFORMS-hg38-v3.tab.gz'),
        readme = os.path.join(RAW_DIR,"VastDB",'README.md')
    shell:
        """
        # download
        wget --user-agent="Chrome" --no-check-certificate {params.event_genexpr} -O {output.event_genexpr}
        wget --user-agent="Chrome" --no-check-certificate {params.event_psi} -O {output.event_psi}
        wget --user-agent="Chrome" --no-check-certificate {params.event_info} -O {output.event_info}
        wget --user-agent="Chrome" --no-check-certificate {params.event_impact} -O {output.event_impact}
        wget --user-agent="Chrome" --no-check-certificate {params.protein_isoforms} -O {output.protein_isoforms}
        
        # remove sequence columns
        csvcut --tabs --maxfieldsize=10000000 --not-columns {params.cols_todrop} {output.event_info} | csvformat --out-tabs > {output.event_info_clean}
        
        echo "Downloaded on $(date)." > {output.readme}
        echo Done!
        """
        
rule download_vastdb_assemblies:
    params:
        Hs2 = 'https://vastdb.crg.eu/libs/vastdb.hs2.20.12.19.tar.gz'
    output:
        Hs2 = directory(os.path.join(RAW_DIR,"VastDB",'assemblies','Hs2'))
    shell:
        """
        # Hs2
        wget --user-agent="Chrome" --no-check-certificate {params.Hs2} -O {output.Hs2}.tar.gz
        tar xzvf {output.Hs2}.tar.gz
        rm {output.Hs2}.tar.gz
        
        echo "Done!"
        """
        
rule download_ucscxena_gdc_pancan:
    params:
        snv = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/GDC-PANCAN.mutect2_snv.tsv.gz"
    output:
        snv = os.path.join(RAW_DIR,"UCSCXena","GDC","snv","GDC-PANCAN.mutect2_snv.tsv.gz"),
        readme = os.path.join(RAW_DIR,"UCSCXena",'GDC','README.md')
    shell:
        """
        # SNV
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.snv} -O {output.snv}
        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """        
        
rule download_ucscxena_tcga:
    output:
        clinical = os.path.join(RAW_DIR,"UCSCXena",'TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz'),
        molecular_subtype = os.path.join(RAW_DIR,"UCSCXena",'TCGA','phenotype','TCGASubtype.20170308.tsv.gz'),
        sample_type = os.path.join(RAW_DIR,"UCSCXena",'TCGA','phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz'),
        readme = os.path.join(RAW_DIR,"UCSCXena",'TCGA','README.md')
    params:
        clinical = 'https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/latest/Survival_SupplementalTable_S1_20171025_xena_sp.gz',
        molecular_subtype = 'https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/latest/TCGASubtype.20170308.tsv.gz',
        sample_type = 'https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/latest/TCGA_phenotype_denseDataOnlyDownload.tsv.gz'
    shell:
        """
        # Phenotype
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.clinical} -O {output.clinical}
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.molecular_subtype} -O {output.molecular_subtype}
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.sample_type} -O {output.sample_type}        
        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """
        
rule download_stringdb:
    params:
        links = 'https://stringdb-static.org/download/protein.links.full.v11.5/9606.protein.links.full.v11.5.txt.gz',
        aliases = 'https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz'
    output:
        links = os.path.join(RAW_DIR,"STRINGDB",'9606.protein.links.full.v11.5.txt.gz'),
        aliases = os.path.join(RAW_DIR,"STRINGDB",'9606.protein.aliases.v11.5.txt.gz'),
        readme = os.path.join(RAW_DIR,"STRINGDB",'README.md')
    shell:
        """
        # download
        wget --no-check-certificate {params.links} -O {output.links}
        wget --no-check-certificate {params.aliases} -O {output.aliases}
        
        # add readme
        echo "Downloaded on $(date)." > {output.readme}
        echo Done!
        """
        
rule download_harmonizome_gene_sets:
    params:
        chea = "https://maayanlab.cloud/static/hdfs/harmonizome/data/cheappi/gene_set_library_crisp.gmt.gz"
    output:
        chea = os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz"),
        readme = os.path.join(RAW_DIR,"Harmonizome",'README.md')
    shell:
        """
        # aneuploidy
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.chea} -O {output.chea}        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """ 