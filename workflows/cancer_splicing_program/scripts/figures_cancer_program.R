#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - Validate whether protein activity inferrence indicates that indisulam targets RBM39
# - Investigate the effect of MS023, inhibitor of Type I PRMT enzymes
#     - no cell growth effect in vitro, but strong suppression of tumor growth in vivo

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)
require(clusterProfiler)
require(org.Hs.eg.db)
require(pROC)
require(readxl)
require(ComplexHeatmap)
require(ggplotify)
require(ggbeeswarm)

# variables
THRESH_FDR = 0.05
ORGDB = org.Hs.eg.db
THRESH_N_SUM = 5.5
DRIVER_TYPES = c("Oncogenic","Tumor suppressor")
SPLICEOSOME_STAGES = c('B complex','Bact complex','C complex','P complex','A complex','U5 snRNP','tri-snRNP','17S U2 snRNP','U1 snRNP')

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = c(
    "Random Genes"="darkgreen",
    "Not Driver-like"="grey",
    "Oncogenic"="#F6AE2D",
    "Tumor suppressor"="#6C98B3"
)

PAL_CANCER_TYPES = setNames(
    c("#A6CEE3","#3385BB","#84BF96","#6DBD57","#7F9D55","#F57C7C","#E42622",
      "#FBB268","#FE8D19","#DE9E83","#9D7BBA","#977899","#F3E587","#B15928"),
    c("BLCA","BRCA","COAD","HNSC","KICH","KIRC","KIRP",
      "LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
)

PAL_IMMUNE_SCREEN = setNames(c("#EB9486","#7E7F9A"), c(TRUE, FALSE))

PAL_SAMPLE_TYPE = setNames(
    get_palette(palette = "npg", k=3),
    c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood","Solid Tissue Normal")
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","cancer_splicing_program")
# REGINF_DIR = file.path(ROOT,"results","regulon_inference")
# diff_activity_file = file.path(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# diff_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN','genexpr_tpm-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# diff_genexpr_deseq_file = file.path(RESULTS_DIR,'files','PANCAN','genexpr_counts_deseq2-deseq2-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","cancer_program")
# demeter2_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# ccle_metadata_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# survival_activity_file = file.path(RESULTS_DIR,'files','PANCAN',"protein_activity-survival_analysis-surv.tsv.gz")
# survival_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-survival_analysis-surv.tsv.gz")
# survival_activity_conf_file = file.path(RESULTS_DIR,'files','PANCAN',"protein_activity-survival_analysis_with_confounders-surv.tsv.gz")
# survival_genexpr_conf_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-survival_analysis_with_confounders-surv.tsv.gz")
# sf_crossreg_activity_file = file.path(RESULTS_DIR,'files','PANCAN',"protein_activity-sf_cross_regulation.tsv.gz")
# sf_crossreg_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-sf_cross_regulation.tsv.gz")
# ontology_chea_file = file.path(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz")
# sf_activity_vs_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm_vs_activity.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","PANCAN.tsv.gz")
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# regulons_path = file.path(REGINF_DIR,"files","experimentally_derived_regulons_pruned-EX")
# regulons_jaccard_file = file.path(REGINF_DIR,"files","regulons_eda_jaccard","experimentally_derived_regulons_pruned-EX.tsv.gz")
# msigdb_dir = file.path(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# immune_screen_file = file.path(SUPPORT_DIR,"supplementary_tables_literature","Dubrot2022-suptabs-41590_2022_1315_MOESM2_ESM.xlsx") # Sup. Tab. 13
# human2mouse_file = file.path(RAW_DIR,"BIOMART","human2mouse.tsv")
# event_prior_knowledge_file = file.path(SUPPORT_DIR,"supplementary_tables_literature","AngladaGirotto2024-supdata01_event_prior_knowledge.txt")
# rbpdb_file = file.path(SUPPORT_DIR,"supplementary_tables_literature","Cook2011-RBPDB_v1.3.1_proteins_human_2012-11-21.tsv")
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")

##### FUNCTIONS #####
load_regulons = function(regulons_path, patt=NULL){
    if (file.exists(regulons_path) && !dir.exists(regulons_path)){
        # regulons_path is a file, we load only that regulon (we'll tun regular VIPER)
        regulon_files = list(regulons_path)
    }else if (dir.exists(regulons_path)){
        # regulons_path is a directory, we load all regulons contained (we'll run metaVIPER)
        regulon_files = list.files(regulons_path, pattern=patt, full.names=TRUE)
    }else {
        stop("Invalid regulons_path.")
    }
    
    regulons = sapply(regulon_files, function(regulon_file){
        regulon = read_tsv(regulon_file)
        return(regulon)
    }, simplify=FALSE)
    
    return(regulons)
}


plot_comparison = function(diff_activity, diff_genexpr, survival_activity, survival_genexpr, sf_activity_vs_genexpr){
    plts = list()    
    
    X = diff_activity %>%
        filter(is_significant) %>%
        dplyr::rename(activity = "condition_a-median") %>%
        distinct(GENE, activity, cancer_type) %>%
        left_join(
            diff_genexpr %>%
            filter(is_significant) %>%
            dplyr::rename(genexpr_log2FC = "median_log2FC") %>%
            distinct(GENE, genexpr_log2FC, cancer_type),
            by=c("GENE","cancer_type")
        )
    
    # do activities and changes in gene expression correlate?
    plts[["comparison-diff_analysis-scatter"]] = X %>%
        ggscatter(x="activity", y="genexpr_log2FC", alpha=0.5, size=1, color="brown") +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Protein Activity", y="Gene Expression log2FC")
    
    # do coxph coefficients correlate?
    X = survival_activity %>%
        distinct(GENE, cancer_type, coxph_coef) %>%
        left_join(
            survival_genexpr %>%
            distinct(GENE, cancer_type, coxph_coef),
            by=c("GENE","cancer_type"),
            suffix=c("_activity","_genexpr")
        )
    
    plts[["comparison-survival-scatter"]] = X %>%
        ggscatter(x="coxph_coef_activity", y="coxph_coef_genexpr", alpha=0.5, size=1) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="CoxPH Protein Activity", y="CoxPH Gene Expression")
    
    # are there splicing factors whose activity is highly correlated with their expression?
    X = sf_activity_vs_genexpr
    plts[["comparison-correlation_by_cancer-violin"]] = X %>%
        ggplot(aes(x=cancer_type, y=correlation)) +
        geom_violin(aes(fill=cancer_type), color=NA, trim=TRUE) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1) +
        fill_palette(get_palette("Paired", 33)) + 
        theme_pubr(x.text.angle=75) +
        guides(fill="none") +
        labs(x="Cancer Type", y="Spearman Correlation")
    
    plts[["comparison-correlation_median_pancan-violin"]] = X %>%
        group_by(GENE) %>%
        summarize(
            correlation = median(correlation),
            comparison = "Protein Activity vs Gene Expression"
        ) %>%
        ungroup() %>%
        ggplot(aes(x=comparison, y=correlation)) +
        geom_violin(aes(fill=comparison), color=NA, trim=TRUE) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1) +
        fill_palette("lightgreen") + 
        geom_text_repel(
            aes(label=GENE),
            . %>% slice_max(correlation, n=5),
            size=FONT_SIZE, family=FONT_FAMILY,
            segment.size=0.1, max.overlaps=50
        ) +
        theme_pubr(x.text.angle=0) +
        guides(fill="none") +
        labs(x="Comparison", y="Median Spearman Correlation")
    
    return(plts)
}


plot_driver_selection = function(driver_activity, driver_genexpr, diff_activity, diff_genexpr){
    plts = list()
    
    # SF activity
    X = driver_activity
    x = X %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    
    sf_order = x %>%
        pivot_wider(id_cols="GENE", names_from="driver_type", values_from="n_sign", values_fill=0) %>%
        rowwise() %>%
        mutate(diff=sum(`Tumor suppressor`,`Oncogenic`)) %>%
        ungroup() %>%
        arrange(diff,`Tumor suppressor`,`Oncogenic`) %>%
        pull(GENE)
    
    sfs_oi = x %>% 
        slice_max(n_sum, n=5) %>% 
        filter(driver_type=="Oncogenic") %>%
        bind_rows(
            x %>% 
            slice_min(n_sum, n=5) %>% 
            filter(driver_type=="Tumor suppressor")
        )

    plts[["driver_selection-n_signif_vs_driver_type-activity-bar"]] = x %>%
        mutate(GENE = factor(GENE, levels=sf_order)) %>%
        ggbarplot(x="GENE", y="n", fill="driver_type", color=NA) +
        geom_text(
            aes(label=GENE),
            sfs_oi %>% filter(driver_type=="Oncogenic"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=1, vjust=1, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE),
            sfs_oi %>% filter(driver_type=="Tumor suppressor"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=45, hjust=0, vjust=0.5, nudge_y=0.25
        ) +
        geom_hline(yintercept=THRESH_N_SUM, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_text(
            aes(x=x, label=label),
            . %>% 
            filter(abs(n_sum)>THRESH_N_SUM) %>%
            group_by(GENE) %>%
            slice_max(n, n=1) %>%
            ungroup() %>%
            count(driver_type) %>% 
                mutate(
                    label=paste0("n=",n),
                    n=c(15,15),
                    x=c(120,40)
                ),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~driver_type, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="Count", fill="Driver Type")
    
    mat = X %>%
        distinct(driver_type, cancer_type, GENE) %>%
        mutate(detected = ifelse(driver_type=="Oncogenic",1,2)) %>%
        pivot_wider(id_cols="cancer_type", names_from="GENE", values_from="detected", 
                    values_fill=NA, values_fn = ~ sum(.x, na.rm=TRUE)) %>%
        column_to_rownames("cancer_type") %>%
        as.matrix()
    
    colors = structure(c('#6C98B3','#F6AE2D','green','gray'), names = c(2,1,3,NA))
    plts[["driver_selection-sf_vs_cancer_type-heatmap"]] = mat[,sf_order] %>%
        Heatmap(
            col=colors, 
            name="Is activated",
            row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            heatmap_legend_param = list(legend_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY)),
            cluster_rows=TRUE,
            cluster_columns=FALSE
        ) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    # are activities different between groups?
    Y = x %>% 
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>%
        left_join(diff_activity, by="GENE") %>%
        drop_na(driver_type, cancer_type)
    plts[["driver_selection-drivers_vs_cancer_type-activity_drivers_vs_activity-violin"]] = Y %>%
        ggplot(aes(x=cancer_type, y=`condition_a-median`, group=interaction(cancer_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Cancer Type", y="median(Protein Activity)", fill="Driver Type")

    Y = x %>% 
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>%
        left_join(diff_genexpr, by="GENE") %>%
        drop_na(driver_type, cancer_type)
    plts[["driver_selection-drivers_vs_cancer_type-activity_drivers_vs_genexpr-violin"]] = Y %>%
        ggplot(aes(x=cancer_type, y=`condition_a-median`, group=interaction(cancer_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Cancer Type", y="median(Protein Activity)", fill="Driver Type")
    
    # SF genexpr
    X = driver_genexpr
    x = X %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    
    plts[["driver_selection-n_signif_vs_driver_type-genexpr-bar"]] = x %>%
        mutate(GENE = factor(GENE, levels=sf_order)) %>%
        ggbarplot(x="GENE", y="n", fill="driver_type", color=NA) +
        geom_text(
            aes(label=GENE),
            . %>% filter(driver_type=="Oncogenic" & GENE%in%(sfs_oi %>% filter(driver_type=="Oncogenic") %>% pull(GENE))),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=1, vjust=1, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE),
            . %>% filter(driver_type=="Tumor suppressor" & GENE%in%(sfs_oi %>% filter(driver_type=="Tumor suppressor") %>% pull(GENE))),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=0, vjust=1, nudge_y=-0.25
        ) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~driver_type, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="Count", fill="Driver Type")

    sf_order = x %>%
        pivot_wider(id_cols="GENE", names_from="driver_type", values_from="n_sign", values_fill=0) %>%
        rowwise() %>%
        mutate(diff=sum(`Tumor suppressor`,`Oncogenic`)) %>%
        ungroup() %>%
        arrange(diff,`Tumor suppressor`,`Oncogenic`) %>%
        pull(GENE)
    
    plts[["driver_selection-n_signif_vs_driver_type_by_genexpr-genexpr-bar"]] = x %>%
        mutate(GENE = factor(GENE, levels=sf_order)) %>%
        ggbarplot(x="GENE", y="n", fill="driver_type", color=NA) +
        geom_text(
            aes(label=GENE),
            . %>% filter(driver_type=="Oncogenic" & GENE%in%(sfs_oi %>% filter(driver_type=="Oncogenic") %>% pull(GENE))),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=1, vjust=1, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE),
            . %>% filter(driver_type=="Tumor suppressor" & GENE%in%(sfs_oi %>% filter(driver_type=="Tumor suppressor") %>% pull(GENE))),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=45, hjust=0, vjust=0.5, nudge_y=0.25
        ) +
        geom_hline(yintercept=THRESH_N_SUM, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_text(
            aes(x=x, label=label),
            . %>% 
            filter(abs(n_sum)>THRESH_N_SUM) %>%
            group_by(GENE) %>%
            slice_max(n, n=1) %>%
            ungroup() %>%
            count(driver_type) %>% 
                mutate(
                    label=paste0("n=",n),
                    n=c(15,15),
                    x=c(120,40)
                ),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~driver_type, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="Count", fill="Driver Type")
    
    # are expressions different between oncogenic and tumor suppressor selected via gene expression?
    Y = x %>% 
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>%
        left_join(diff_genexpr, by="GENE") %>%
        drop_na(driver_type, cancer_type)
    plts[["driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_genexpr-violin"]] = Y %>%
        ggplot(aes(x=cancer_type, y=`condition_a-median`, group=interaction(cancer_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Cancer Type", y="median(Gene Expression)", fill="Driver Type")
    
    # are activities different between oncogenic and tumor suppressor selected via gene expression?
    Y = x %>% 
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>%
        left_join(diff_activity, by="GENE") %>%
        drop_na(driver_type, cancer_type)
    plts[["driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_activity-violin"]] = Y %>%
        ggplot(aes(x=cancer_type, y=`condition_a-median`, group=interaction(cancer_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Cancer Type", y="median(Protein Activity)", fill="Driver Type")
    
    return(plts)
}


plot_eda_programs = function(driver_activity, splicing_factors, rbpdb){
    plts = list()
    
    X = driver_activity %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>%
        left_join(
                splicing_factors %>% 
                    distinct(GENE, spliceosome_db_complex) %>%
                    drop_na(), 
            by="GENE"
        ) %>%
        left_join(
            rbpdb %>% 
                distinct(X5,X9) %>%
                dplyr::rename(GENE=X5, rbp_family=X9) %>%
                drop_na(), 
            by="GENE"
        )
    
    # stages of the spliceosome, only 11 splicing factors could be considered
    x = X %>%
        separate_rows(spliceosome_db_complex, sep=", ") %>%
        drop_na() %>% 
        distinct(GENE, spliceosome_db_complex, driver_type) %>%
        group_by(spliceosome_db_complex, driver_type) %>%
        mutate(sfs = paste(GENE, collapse=", ")) %>%
        ungroup() %>%
        count(spliceosome_db_complex, driver_type, sfs)
        
    plts[["eda_programs-spliceosome_db_complex-bar"]] = x %>%
        mutate(spliceosome_db_complex = factor(spliceosome_db_complex, levels=SPLICEOSOME_STAGES)) %>%
        ggbarplot(x="spliceosome_db_complex", y="n", fill="driver_type", color=NA, position=position_dodge(0.9)) +
        geom_text(aes(label=sfs, group=driver_type), size=FONT_SIZE, family=FONT_FAMILY, angle=90, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 45) +
        labs(x="Spliceosome Complex", y="Count", subtitle="11 SFs out of 509")
    
    # RBP families
    x = X %>%
        separate_rows(rbp_family, sep="; ") %>%
        drop_na() %>% 
        distinct(GENE, rbp_family, driver_type) %>%
        group_by(rbp_family, driver_type) %>%
        mutate(sfs = paste(GENE, collapse=", ")) %>%
        ungroup() %>%
        count(rbp_family, driver_type, sfs)
        
    plts[["eda_programs-rbpdb_family-bar"]] = x %>%
        ggbarplot(x="rbp_family", y="n", fill="driver_type", color=NA, position=position_dodge(0.9)) +
        geom_text(aes(label=sfs, group=driver_type), size=FONT_SIZE, family=FONT_FAMILY, angle=90, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 45) +
        labs(x="RBPDB Family (7 out of 49)", y="Count", subtitle="Total SFs = 11")
    
    return(plts)
}

plot_prolif_driver = function(diff_activity, demeter2){
    plts = list()
    
    X = diff_activity %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    
    # Do inferred tumor suppressor and oncogenic activities match gene dependencies?
    gene_dependencies = demeter2 %>%
        mutate(
            cell_type = case_when(
                primary_disease=="Engineered" ~ "Engineered",
                primary_disease=="Fibroblast" ~ "Fibroblast",
                str_detect(primary_disease, "Cancer") | primary_disease%in%c("Leukemia","Myeloma","Lymphoma","Rhabdoid","Neuroblastoma") ~ "Cancerous"
            )
        )
    
    set.seed(1234)
    random_genes = gene_dependencies %>% 
        distinct(index, is_pan_essential) %>%
        group_by(is_pan_essential) %>%
        slice_sample(n=50) %>%
        ungroup()
    dependencies_random_genes = random_genes %>%
        left_join(gene_dependencies, by=c("index","is_pan_essential")) %>%
        mutate(
            driver_type = "Random Genes",
            GENE = index
        )
    
    x = X %>%
        mutate(like_driver = abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1, with_ties=FALSE) %>% 
        ungroup() %>%
        mutate(
            driver_type = ifelse(like_driver, driver_type, "Not Driver-like")
        ) %>%
        left_join(gene_dependencies, by=c("GENE"="index")) %>%
        bind_rows(
            dependencies_random_genes
        )
    
    n_cell_types = x %>% 
        distinct(DepMap_ID, cell_type) %>%
        count(cell_type) %>%
        mutate(cell_type_lab = sprintf("%s (n=%s)", cell_type, n))
    
    x = x %>%
        drop_na(cell_type, driver_type, gene_dependency) %>%
        group_by(GENE, driver_type, cell_type, is_pan_essential) %>%
        summarize(
            gene_dependency = median(gene_dependency)
        ) %>%
        ungroup() %>%
        left_join(n_cell_types, by="cell_type") %>%
        mutate(is_pan_essential = ifelse(is_pan_essential, "Essential", "Not essential"))
    
    comparisons = list(c("Random Genes","Not Driver-like"),c("Random Genes","Oncogenic"),c("Random Genes","Tumor suppressor"))
    plts[["prolif_driver-driver_type_vs_demeter2-essential-box"]] = x %>%
        drop_na(driver_type, gene_dependency) %>%
        filter(is_pan_essential=="Essential" & driver_type%in%DRIVER_TYPES) %>%
        mutate(driver_type=factor(driver_type, levels=DRIVER_TYPES)) %>%
        ggplot(aes(x=driver_type, y=gene_dependency)) +
        geom_quasirandom(aes(color=driver_type), size=0.1, varwidth=0.5) +
        geom_boxplot(width=0.2, outlier.shape=NA, color="black", fill=NA) +
        color_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text(
            aes(y=-3, label=label),
            . %>% count(cell_type_lab, driver_type, is_pan_essential) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~is_pan_essential+cell_type_lab) + 
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Driver Type", y="median(Demeter2 Gene Dependency)")
    
    plts[["prolif_driver-driver_type_vs_demeter2-not_essential-box"]] = x %>%
        drop_na(driver_type, gene_dependency) %>%
        filter(is_pan_essential=="Not essential" & driver_type%in%DRIVER_TYPES) %>%
        mutate(driver_type=factor(driver_type, levels=DRIVER_TYPES)) %>%
        ggplot(aes(x=driver_type, y=gene_dependency)) +
        geom_quasirandom(aes(color=driver_type), size=0.1, varwidth=0.5) +
        geom_boxplot(width=0.2, outlier.shape=NA, color="black", fill=NA) +
        color_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text(
            aes(y=-1.1, label=label),
            . %>% count(cell_type_lab, driver_type, is_pan_essential) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~is_pan_essential+cell_type_lab) + 
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Driver Type", y="median(Demeter2 Gene Dependency)")
    
    return(plts)
}


make_roc_analysis = function(driver_classif, survival_analysis_surv){
    # roc curve based on thresholds of n_sum for each type of cancer
    # We expect to classify Tumor Suppressor as Low Risk and Oncogenic as High Risk
    driver_types = c("Tumor suppressor", "Oncogenic")
    driver_levels = list(
        "Tumor suppressor" = c("High Risk", "Low Risk"),
        "Oncogenic" = c("Low Risk", "High Risk")
    )
    
    driver_classif = driver_classif %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    roc_analysis = lapply(driver_types, function(driver_type_oi){
        result = driver_classif %>%
            filter(driver_type==driver_type_oi) %>%
            left_join(
                survival_analysis_surv,
                by="GENE"
            ) %>%
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")) %>%
            count(GENE, driver_type, surv_type) %>%
            roc(response=surv_type, predictor=n, 
                levels=driver_levels[[driver_type_oi]], direction="<")
        
        result = result %>%
            coords(transpose=FALSE, re="all") %>% 
            mutate(
                fpr = 1 - specificity,
                auc = result[["auc"]],
                driver_type = driver_type_oi
            ) %>%
            ungroup()
    }) %>% do.call(rbind, .)
    
    return(roc_analysis)
}


plot_survival_analysis = function(survival_roc, survival_omic, survival_omic_conf, driver_omic, patt=""){
    plts = list()
    
    # ROC curves
    X = survival_roc
    plts[["survival_analysis-cancers_all-roc_curves"]] = X %>% 
        filter(cancers_subset=="full") %>%
        mutate(auc_lab = sprintf("AUC(%s)=%s",driver_type,round(auc,2))) %>%
        ggscatter(x="fpr", y="sensitivity", size=1,
                  color="driver_type", palette=PAL_DRIVER_TYPE) +
        geom_path(aes(color=driver_type), linetype="dashed", size=LINE_SIZE) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab, color=driver_type),
            . %>% distinct(auc_lab, driver_type) %>%
                mutate(x=0.12, y=c(0,0.075)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(
            x="False Positive Rate (1 - Specificity)", 
            y="True Positive Rate (Sensitivity)",
            color="Driver Type"
        )
    
    plts[["survival_analysis-cancers_differential-roc_curves"]] = X %>% 
        filter(cancers_subset=="differential") %>%
        mutate(auc_lab = sprintf("AUC(%s)=%s",driver_type,round(auc,2))) %>%
        ggscatter(x="fpr", y="sensitivity", size=1,
                  color="driver_type", palette=PAL_DRIVER_TYPE) +
        geom_path(aes(color=driver_type), linetype="dashed", size=LINE_SIZE) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab, color=driver_type),
            . %>% distinct(auc_lab, driver_type) %>%
                mutate(x=0.12, y=c(0,0.075)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text_repel(aes(label=threshold), size=FONT_SIZE, family=FONT_FAMILY,
                        segment.size=0.1, max.overlaps=50) +
        theme(aspect.ratio=1) +
        labs(
            x="False Positive Rate (1 - Specificity)", 
            y="True Positive Rate (Sensitivity)",
            color="Driver Type"
        )
    
    
    # Frequency of risk association
    X = driver_omic %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    plts[["survival_analysis-cancer_all-examples-bar"]] = X %>%
        left_join(
            survival_omic %>%            
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")),
            by="GENE"
        ) %>%
        count(GENE, driver_type, surv_type) %>%
        drop_na() %>%
        filter(GENE %in% c("SNRNP200","HNRNPD")) %>%
        #mutate(GENE = factor(GENE, levels=c("SNRNP200","HNRNPD"))) %>%
        ggbarplot(
            x="surv_type", y="n", fill="driver_type", color=NA, 
            palette=PAL_DRIVER_TYPE, position=position_dodge(0.9), 
            label=TRUE, lab.family=FONT_FAMILY, lab.size=FONT_SIZE
        ) +
        facet_wrap(~GENE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="N. Cancer Types", fill="Driver Type")
        
    
    plts[["survival_analysis-cancers_all-violin"]] = X %>%
        left_join(
            survival_omic %>%            
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")),
            by="GENE"
        ) %>%
        count(GENE, driver_type, surv_type) %>%
        drop_na() %>%
        ggplot(aes(x=surv_type, y=n, group=interaction(surv_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        geom_text(
            aes(label=label, y=3), 
            . %>% count(surv_type,driver_type) %>% mutate(label=paste0("n=",n)),
            position=position_dodge(0.9),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        labs(x="Risk with High Protein Activity", y="N. Cancers", fill="Driver Type")
    
    
    plts[["survival_analysis_conf-cancer_all-examples-bar"]] = X %>%
        left_join(
            survival_omic_conf %>%            
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")),
            by="GENE"
        ) %>%
        count(GENE, driver_type, surv_type) %>%
        drop_na() %>%
        filter(GENE %in% c("SNRNP200","HNRNPD")) %>%
        #mutate(GENE = factor(GENE, levels=c("SNRNP200","HNRNPD"))) %>%
        ggbarplot(
            x="surv_type", y="n", fill="driver_type", color=NA, 
            palette=PAL_DRIVER_TYPE, position=position_dodge(0.9), 
            label=TRUE, lab.family=FONT_FAMILY, lab.size=FONT_SIZE
        ) +
        facet_wrap(~GENE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Splicing Factor", y="N. Cancer Types", fill="Driver Type")
        
    
    plts[["survival_analysis_conf-cancers_all-violin"]] = X %>%
        left_join(
            survival_omic_conf %>%            
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")),
            by="GENE"
        ) %>%
        count(GENE, driver_type, surv_type) %>%
        drop_na() %>%
        ggplot(aes(x=surv_type, y=n, group=interaction(surv_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        geom_text(
            aes(label=label, y=3), 
            . %>% count(surv_type,driver_type) %>% mutate(label=paste0("n=",n)),
            position=position_dodge(0.9),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        labs(x="Risk with High Protein Activity", y="N. Cancers", fill="Driver Type")    
    
    # different n_sum thresholds
    thresholds = c(8,10,12)
    plts[["survival_analysis-cancers_all_vs_thresh-violin"]] = lapply(thresholds, function(thresh){
            X %>%
                filter(abs(n_sum)>=thresh) %>%
                left_join(
                    survival_omic %>%            
                    mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")),
                    by="GENE"
                ) %>%
                count(GENE, driver_type, surv_type) %>%
                drop_na() %>%
                mutate(n_sum_threshold = sprintf("thresh=%s",thresh))
        }) %>% 
        bind_rows() %>%
        mutate(n_sum_threshold = factor(n_sum_threshold, levels=sprintf("thresh=%s",thresholds))) %>%
        ggplot(aes(x=surv_type, y=n, group=interaction(surv_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        geom_text(
            aes(label=label, y=3), 
            . %>% count(n_sum_threshold,surv_type,driver_type) %>% mutate(label=paste0("n=",n)),
            position=position_dodge(0.9),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~n_sum_threshold) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Risk with High Protein Activity", y="N. Cancers", fill="Driver Type")
    
    
    plts[["survival_analysis-cancers_differential-violin"]] = X %>%
        left_join(
            survival_omic %>%
            filter(cancer_type %in% driver_omic[["cancer_type"]]) %>%
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")),
            by="GENE"
        ) %>%
        count(GENE, driver_type, surv_type) %>%
        drop_na() %>%
        ggplot(aes(x=surv_type, y=n, group=interaction(surv_type,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        geom_text(
            aes(label=label, y=0), 
            . %>% count(surv_type,driver_type) %>% mutate(label=paste0("n=",n)),
            position=position_dodge(0.9),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        labs(x="Risk with High Protein Activity", y="N. Cancers", fill="Driver Type")
    
    
    names(plts) = paste0(names(plts),patt)
    
    return(plts)
}


plot_sf_crossreg = function(driver_omic, sf_crossreg, regulons_jaccard, patt=""){
    plts = list()
    
    lab_pairs = c(
        "Tumor suppressor\n&\nTumor suppressor",
        "Oncogenic\n&\nTumor suppressor",
        "Oncogenic\n&\nOncogenic"
    )
    
    driver_classif = driver_omic %>%
        count(GENE, ENSEMBL, driver_type) %>%
        group_by(GENE, ENSEMBL) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE, ENSEMBL) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    oncogenic = driver_classif %>% filter(driver_type=="Oncogenic") %>% pull(ENSEMBL)
    suppressor = driver_classif %>% filter(driver_type=="Tumor suppressor") %>% pull(ENSEMBL)
    
    X = sf_crossreg %>%
        mutate(
            driver_type_a = case_when(
                regulator_a %in% oncogenic ~ "Oncogenic",
                regulator_a %in% suppressor ~ "Tumor suppressor"
            ),
            driver_type_b = case_when(
                regulator_b %in% oncogenic ~ "Oncogenic",
                regulator_b %in% suppressor ~ "Tumor suppressor"
            ),
            same_driver_type = driver_type_a==driver_type_b
        ) %>%
        drop_na(driver_type_a, driver_type_b) %>%
        rowwise() %>%
        mutate(lab_pair = paste(sort(c(driver_type_a,driver_type_b)), collapse="\n&\n")) %>%
        ungroup() %>%
        group_by(regulator_a,regulator_b,lab_pair) %>%
        summarize(correlation = median(correlation)) %>%
        ungroup()
    
    # are oncogenic vs suppressor corregulated?
    plts[["sf_cross_regulation-correlations-violin"]] = X %>%
        mutate(lab_pair=factor(lab_pair, levels=lab_pairs)) %>%
        ggviolin(x="lab_pair", y="correlation", fill="lab_pair", color=NA, trim=TRUE,
                 palette=c(PAL_DRIVER_TYPE[[1]],"lightgray",PAL_DRIVER_TYPE[[2]])) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, position=position_dodge(0.9)) +
        stat_compare_means(
            method="wilcox.test", label="p.format", ref.group="Oncogenic\n&\nTumor suppressor",
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-0.5, label=label),
            . %>% count(lab_pair) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="", y="SF-SF Protein Activity Correlation")
    
    # are regulator targets shared among cancer-driver genes?
    X = regulons_jaccard %>%
        mutate(
            driver_type_a = case_when(
                regulator_a %in% oncogenic ~ "Oncogenic",
                regulator_a %in% suppressor ~ "Tumor suppressor"
            ),
            driver_type_b = case_when(
                regulator_b %in% oncogenic ~ "Oncogenic",
                regulator_b %in% suppressor ~ "Tumor suppressor"
            ),
            same_driver_type = driver_type_a==driver_type_b,
            similarity = 1 - distance
        ) %>%
        drop_na(driver_type_a, driver_type_b) %>%
        rowwise() %>%
        mutate(lab_pair = paste(sort(c(driver_type_a,driver_type_b)), collapse="\n&\n")) %>%
        ungroup() %>%
        group_by(regulator_a,regulator_b,lab_pair) %>%
        summarize(similarity = median(similarity)) %>%
        ungroup()
    
    
    plts[["sf_cross_regulation-regulon_similarity-violin"]] = X %>%
        mutate(lab_pair=factor(lab_pair, levels=lab_pairs)) %>%
        ggviolin(x="lab_pair", y="similarity", fill="lab_pair", color=NA, trim=TRUE,
                 palette=c(PAL_DRIVER_TYPE[[1]],"lightgray",PAL_DRIVER_TYPE[[2]])) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, position=position_dodge(0.9)) +
        stat_compare_means(
            method="wilcox.test", label="p.format", ref.group="Oncogenic\n&\nTumor suppressor",
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-0.005, label=label),
            . %>% count(lab_pair) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="", y="SF-SF Regulons Jaccard Similarity")
    
    
    names(plts) = paste0(names(plts),patt)
    
    return(plts)
}


make_enrichments = function(driver_classif, regulons, ontologies){
    
    X = driver_classif %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, ENSEMBL, driver_type) %>%
        group_by(GENE, ENSEMBL) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE, ENSEMBL) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    oncogenics = X %>% filter(driver_type=="Oncogenic") %>% distinct(GENE,ENSEMBL)
    suppressors = X %>% filter(driver_type=="Tumor suppressor") %>% distinct(GENE,ENSEMBL)
    enrichments = list(
        "CHEA" = list(
            "oncogenics" = enricher(oncogenics[["GENE"]], TERM2GENE=ontologies[["CHEA"]]),
            "suppressors" = enricher(suppressors[["GENE"]], TERM2GENE=ontologies[["CHEA"]])
        ),
        "reactome" = list(
            "oncogenics" = enricher(
                regulons %>% filter(regulator %in% oncogenics[["ENSEMBL"]]) %>% pull(GENE),
                TERM2GENE=ontologies[["reactome"]]
            ),
            "suppressors" = enricher(
                regulons %>% filter(regulator %in% suppressors[["ENSEMBL"]]) %>% pull(GENE),
                TERM2GENE=ontologies[["reactome"]]
            )
        )
    )
    
    return(enrichments)
}


plot_enrichments = function(enrichments_reactome, immune_screen){
    plts = list()
    
    X = enrichments_reactome
    
    terms_oi = X %>%
        group_by(gene_set) %>%
        slice_max(GeneRatio, n=10) %>%
        ungroup() %>%
        pull(Description) %>%
        unique()
    
    plts[["reactome_enrichments-bar"]] = X %>%
        filter(Description %in% terms_oi) %>%
        group_by(Description) %>%
        mutate(ratio_sums = sum(GeneRatio)) %>%
        ungroup() %>%
        arrange(GeneRatio) %>%
        ggbarplot(x="Description", y="GeneRatio", fill="driver_type", color=NA,
                  palette=PAL_DRIVER_TYPE, position=position_dodge(0.9)) +
        geom_text(aes(label=Count, group=driver_type), 
                  size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9), hjust=-0.1) +
        labs(x="Description", y="GeneRatio", fill="Driver Type") +
        coord_flip()
    
    # antigen processing related genes found in immunogenesis screen
    x = immune_screen %>%
        left_join(
            X %>%
                filter(str_detect(ID, "ANTIGEN")) %>%
                separate_rows(geneID),
            by=c("human_symbol"="geneID"),
            relationship="many-to-many"
        ) %>%
        mutate(
            in_reactome = !is.na(ID),
            score = Sign*`Average Score`
        ) %>% 
        distinct(Gene, score, in_reactome, Dataset, Comparison) %>% 
        group_by(Dataset, Comparison) %>% 
        arrange(score) %>% 
        mutate(ranking = row_number()) %>% 
        ungroup()
    
    genes_oi = x %>%
        group_by(Dataset, Comparison, in_reactome) %>%
        slice_max(abs(score), n=6) %>%
        ungroup() %>%
        left_join(immune_screen %>% distinct(Gene, human_symbol), by="Gene") %>%
        mutate(gene_lab = sprintf("%s (%s)", Gene, human_symbol))
    
    plts[["reactome_enrichments-immune_screen-scatter"]] = x %>% 
        ggplot(aes(x=ranking, y=score, color=in_reactome)) + 
        #geom_point(data = . %>% filter(!in_reactome)) + 
        geom_point(data = . %>% filter(in_reactome)) + 
        color_palette(PAL_IMMUNE_SCREEN) +
        geom_text_repel(
            aes(label=gene_lab),
            genes_oi %>% filter(in_reactome),
            color="black",
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1
        ) +
        facet_wrap(~Dataset+Comparison) + 
        theme_pubr() +
        labs(x="Ranking", y="Score")
    
    return(plts)
}


plot_n_samples = function(n_samples){
    plts = list()
    
    X = n_samples
    
    plts[["n_samples-balloon"]] = X %>%
        arrange(n) %>%
        ggballoonplot(x="sample_type", y="cancer_type", size="n",
                      fill="sample_type", color="sample_type", palette=PAL_SAMPLE_TYPE) +
        coord_radial(theta="y", inner.radius=0.3, rotate_angle=TRUE) +
        guides(color="none", theta=guide_axis_theta(angle = 90)) +
        labs(size="N. Samples", fill="Sample Type")
    
    return(plts)
}


make_plots = function(
    diff_activity, diff_genexpr,
    demeter2, 
    survival_roc_activity, survival_roc_genexpr, 
    survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
    survival_activity, survival_genexpr, 
    survival_activity_conf, survival_genexpr_conf, 
    driver_activity, driver_genexpr, 
    sf_crossreg_activity, sf_crossreg_genexpr, 
    enrichments_reactome, immune_screen, sf_activity_vs_genexpr, regulons_jaccard,
    n_samples,
    splicing_factors, rbpdb
){
    plts = list(
        plot_driver_selection(driver_activity, driver_genexpr, diff_activity, diff_genexpr),
        plot_prolif_driver(diff_activity, demeter2),
        plot_survival_analysis(survival_roc_activity, survival_activity, survival_activity_conf, driver_activity, "-activity"),
        plot_survival_analysis(survival_roc_genexpr, survival_genexpr, survival_genexpr_conf, driver_genexpr, "-genexpr"),
        plot_survival_analysis(survival_roc_genexpr_w_activity_labs, survival_genexpr, survival_genexpr_conf, driver_activity, "-genexpr_w_activity_labs"),
        plot_survival_analysis(survival_roc_activity_w_genexpr_labs, survival_activity, survival_activity_conf, driver_genexpr, "-activity_w_genexpr_labs"),
        plot_sf_crossreg(driver_activity, sf_crossreg_activity, regulons_jaccard, "-activity"),
        plot_sf_crossreg(driver_genexpr, sf_crossreg_genexpr, regulons_jaccard, "-genexpr"),
        plot_enrichments(enrichments_reactome, immune_screen),
        plot_comparison(diff_activity, diff_genexpr, survival_activity, survival_genexpr, sf_activity_vs_genexpr),
        plot_n_samples(n_samples),
        plot_eda_programs(driver_activity, splicing_factors, rbpdb)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    diff_activity, diff_genexpr,
    demeter2, 
    survival_roc_activity, survival_roc_genexpr, 
    survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
    survival_activity, survival_genexpr, 
    survival_activity_conf, survival_genexpr_conf, 
    driver_activity, driver_genexpr, 
    sf_crossreg_activity, sf_crossreg_genexpr, 
    enrichments_reactome, immune_screen, sf_activity_vs_genexpr, regulons_jaccard,
    n_samples,
    splicing_factors, rbpdb
){
    figdata = list(
        "cancer_program" = list(
            "enrichments_reactome" = enrichments_reactome
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)   
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    # cancer splicing program definition
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-activity-bar", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-genexpr-bar", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "driver_selection-n_signif_vs_driver_type_by_genexpr-genexpr-bar", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-activity_drivers_vs_activity-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-activity_drivers_vs_genexpr-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_genexpr-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_activity-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-sf_vs_cancer_type-heatmap", '.pdf', figs_dir, width=8, height=5, format=FALSE)
    
    # gene dependencies
    save_plt(plts, "prolif_driver-driver_type_vs_demeter2-essential-box", '.pdf', figs_dir, width=9, height=7)
    save_plt(plts, "prolif_driver-driver_type_vs_demeter2-not_essential-box", '.pdf', figs_dir, width=9, height=7)
    
    # survival - activity
    save_plt(plts, 'survival_analysis-cancers_all-roc_curves-activity', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_differential-roc_curves-activity', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancer_all-examples-bar-activity', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_all-violin-activity', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancer_all-examples-bar-activity', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancers_all-violin-activity', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all_vs_thresh-violin-activity", '.pdf', figs_dir, width=9, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-violin-activity", '.pdf', figs_dir, width=5, height=6)
    
    # survival - genexpr
    save_plt(plts, 'survival_analysis-cancers_all-roc_curves-genexpr', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_differential-roc_curves-genexpr', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancer_all-examples-bar-genexpr', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_all-violin-genexpr', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancer_all-examples-bar-genexpr', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancers_all-violin-genexpr', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all_vs_thresh-violin-genexpr", '.pdf', figs_dir, width=9, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-violin-genexpr", '.pdf', figs_dir, width=5, height=6)

    # survival - activity with genexpr cancer splicing program labels
    save_plt(plts, 'survival_analysis-cancers_all-roc_curves-activity_w_genexpr_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_differential-roc_curves-activity_w_genexpr_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancer_all-examples-bar-activity_w_genexpr_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_all-violin-activity_w_genexpr_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancer_all-examples-bar-activity_w_genexpr_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancers_all-violin-activity_w_genexpr_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all_vs_thresh-violin-activity_w_genexpr_labs", '.pdf', figs_dir, width=9, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-violin-activity_w_genexpr_labs", '.pdf', figs_dir, width=5, height=6)
    
    # survival - genexpr with activity cancer splicing program labels
    save_plt(plts, 'survival_analysis-cancers_all-roc_curves-genexpr_w_activity_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_differential-roc_curves-genexpr_w_activity_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancer_all-examples-bar-genexpr_w_activity_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis-cancers_all-violin-genexpr_w_activity_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancer_all-examples-bar-genexpr_w_activity_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'survival_analysis_conf-cancers_all-violin-genexpr_w_activity_labs', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all_vs_thresh-violin-genexpr_w_activity_labs", '.pdf', figs_dir, width=9, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-violin-genexpr_w_activity_labs", '.pdf', figs_dir, width=5, height=6)
    
    # EDA programs
    save_plt(plts, "eda_programs-spliceosome_db_complex-bar", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "eda_programs-rbpdb_family-bar", '.pdf', figs_dir, width=5, height=7)
    
    # cross regulation splicing factors
    save_plt(plts, "sf_cross_regulation-correlations-violin-activity", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-regulon_similarity-violin-activity", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-correlations-violin-genexpr", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-regulon_similarity-violin-genexpr", '.pdf', figs_dir, width=5, height=5)
    
    save_plt(plts, "reactome_enrichments-bar", '.pdf', figs_dir, width=16, height=8.5)
    save_plt(plts, "reactome_enrichments-immune_screen-scatter", '.pdf', figs_dir, width=7, height=8)
    
    save_plt(plts, "comparison-diff_analysis-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "comparison-survival-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "comparison-correlation_by_cancer-violin", '.pdf', figs_dir, width=12, height=5)
    save_plt(plts, "comparison-correlation_median_pancan-violin", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "n_samples-balloon", '.pdf', figs_dir, width=23, height=10, format=FALSE)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        d = file.path(dir,'figdata',x)
        dir.create(d, recursive=TRUE)
        lapply(names(figdata[[x]]), function(nm){
            df = figdata[[x]][[nm]]
            filename = file.path(d, paste0(nm,'.tsv.gz'))
            write_tsv(df, filename)
            
            print(filename)
        })
    })
}


parseargs = function(){
    
    option_list = list( 
        make_option("--diff_activity_file", type="character"),
        make_option("--diff_genexpr_file", type="character"),
        make_option("--diff_genexpr_deseq_file", type="character"),
        make_option("--survival_activity_file", type="character"),
        make_option("--survival_genexpr_file", type="character"),
        make_option("--survival_activity_conf_file", type="character"),
        make_option("--survival_genexpr_conf_file", type="character"),
        make_option("--sf_crossreg_activity_file", type="character"),
        make_option("--sf_crossreg_genexpr_file", type="character"),
        make_option("--demeter2_file", type="character"),
        make_option("--ccle_metadata_file", type="character"),
        make_option("--ontology_chea_file", type="character"),
        make_option("--sf_activity_vs_genexpr_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--regulons_path", type="character"),
        make_option("--regulons_jaccard_file", type="character"),
        make_option("--annotation_file", type="character"),
        make_option("--gene_annotation_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--immune_screen_file", type="character"),
        make_option("--human2mouse_file", type="character"),
        make_option("--event_prior_knowledge_file", type="character"),
        make_option("--rbpdb_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--figs_dir", type="character")
        make_option("--random_seed", type="integer", default=1234)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    diff_activity_file = args[["diff_activity_file"]]
    diff_genexpr_file = args[["diff_genexpr_file"]]
    diff_genexpr_deseq_file = args[["diff_genexpr_deseq_file"]]
    survival_activity_file = args[["survival_activity_file"]]
    survival_genexpr_file = args[["survival_genexpr_file"]]
    survival_activity_conf_file = args[["survival_activity_conf_file"]]
    survival_genexpr_conf_file = args[["survival_genexpr_conf_file"]]
    sf_crossreg_activity_file = args[["sf_crossreg_activity_file"]]
    sf_crossreg_genexpr_file = args[["sf_crossreg_genexpr_file"]]
    demeter2_file = args[["demeter2_file"]]
    ccle_metadata_file = args[["ccle_metadata_file"]]
    sf_activity_vs_genexpr_file = args[["sf_activity_vs_genexpr_file"]]
    metadata_file = args[["metadata_file"]]
    regulons_path = args[["regulons_path"]]
    annotation_file = args[["annotation_file"]]
    regulons_jaccard_file = args[["regulons_jaccard_file"]]
    gene_annotation_file = args[["gene_annotation_file"]]
    ontology_chea_file = args[["ontology_chea_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    immune_screen_file = args[["immune_screen_file"]]
    human2mouse_file = args[["human2mouse_file"]]
    event_prior_knowledge_file = args[["event_prior_knowledge_file"]]
    rbpdb_file = args[["rbpdb_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    figs_dir = args[["figs_dir"]]
    random_seed = args[["random_seed"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_activity = read_tsv(diff_activity_file)
    diff_genexpr = read_tsv(diff_genexpr_file)
    diff_genexpr_deseq = read_tsv(diff_genexpr_deseq_file)
    survival_activity = read_tsv(survival_activity_file)
    survival_genexpr = read_tsv(survival_genexpr_file)
    survival_activity_conf = read_tsv(survival_activity_conf_file)
    survival_genexpr_conf = read_tsv(survival_genexpr_conf_file)    
    sf_crossreg_activity = read_tsv(sf_crossreg_activity_file)
    sf_crossreg_genexpr = read_tsv(sf_crossreg_genexpr_file)
    demeter2 = read_tsv(demeter2_file)
    ccle_metadata = read_tsv(ccle_metadata_file)
    sf_activity_vs_genexpr = read_tsv(sf_activity_vs_genexpr_file)
    metadata = read_tsv(metadata_file)
    regulons = load_regulons(regulons_path)
    annot = read_tsv(annotation_file)
    regulons_jaccard = read_tsv(regulons_jaccard_file)
    gene_annotation = read_tsv(gene_annotation_file) %>%
        dplyr::rename(
            GENE = `Approved symbol`,
            ENSEMBL = `Ensembl gene ID`
        )
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt")),
        "CHEA" = read.gmt(ontology_chea_file)
    )
    immune_screen = read_excel(immune_screen_file, sheet="Supplementary Table 13")
    human2mouse = read_tsv(
        human2mouse_file, col_names=c("human_ensembl","human_symbol","mouse_ensembl","mouse_symbol")
    )
    event_prior_knowledge = read_tsv(event_prior_knowledge_file)
    rbpdb = read_tsv(rbpdb_file, col_names=FALSE)
    splicing_factors = read_tsv(splicing_factors_file)
    
    # prep
    n_samples = metadata %>%
        filter(sample_type %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood","Solid Tissue Normal")) %>%
        count(cancer_type, sample_type) %>%
        mutate(sample_type = factor(sample_type, levels=c("Solid Tissue Normal","Primary Blood Derived Cancer - Peripheral Blood","Primary Tumor")))
    
    diff_activity = diff_activity %>%
        group_by(cancer_type) %>% # correct p-values for each type of cancer
        mutate(
            padj = p.adjust(pvalue, method="fdr"),
            is_significant = padj < THRESH_FDR
        ) %>%
        ungroup() %>%
        dplyr::rename("ENSEMBL"="regulator") %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = "ENSEMBL"
        )
    
    diff_genexpr = diff_genexpr %>%
        dplyr::rename("ENSEMBL"="ID") %>%
        filter(ENSEMBL %in% diff_activity[["ENSEMBL"]]) %>%
        group_by(cancer_type) %>% # correct p-values for each type of cancer
        mutate(
            padj = p.adjust(pvalue, method="fdr"),
            is_significant = padj < THRESH_FDR
        ) %>%
        ungroup() %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = "ENSEMBL"
        ) 
    
    diff_genexpr_deseq = diff_genexpr_deseq %>%
        dplyr::rename("ENSEMBL"="gene", "median_log2FC"="log2FoldChange") %>%
        mutate(ENSEMBL = gsub("*.\\.","",ENSEMBL)) %>% # remove gene version
        filter(ENSEMBL %in% diff_activity[["ENSEMBL"]]) %>%
        group_by(cancer_type) %>% # correct p-values for each type of cancer
        mutate(
            padj = p.adjust(pvalue, method="fdr"),
            is_significant = padj < THRESH_FDR
        ) %>%
        ungroup() %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = "ENSEMBL"
        )
    
    driver_activity = diff_activity %>%
        mutate(driver_type = ifelse(
            `condition_a-median`>0, "Oncogenic", "Tumor suppressor")
        ) %>%
        filter(is_significant)    
    
    driver_genexpr = diff_genexpr %>%
        filter(GENE%in%diff_activity[["GENE"]] & is_significant) %>%
        mutate(driver_type=ifelse(median_log2FC>0, "Oncogenic", "Tumor suppressor")) %>%
        drop_na(driver_type)
    
    driver_genexpr_deseq = diff_genexpr_deseq %>%
        filter(GENE%in%diff_activity[["GENE"]] & is_significant) %>%
        mutate(driver_type=ifelse(median_log2FC>0, "Oncogenic", "Tumor suppressor")) %>%
        drop_na(driver_type)    
    
    survival_activity = survival_activity %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("feature"="ENSEMBL")
        )
    
    survival_activity_conf = survival_activity_conf %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("feature"="ENSEMBL")
        )    
    
    survival_genexpr = survival_genexpr %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("feature"="ENSEMBL")
        )
    
    survival_genexpr_conf = survival_genexpr_conf %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("feature"="ENSEMBL")
        )
    
    sf_activity_vs_genexpr = sf_activity_vs_genexpr %>%
        filter(sf_genexpr==sf_activity) %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("sf_activity"="ENSEMBL")
        )
    
    immune_screen = immune_screen %>%
        left_join(
            human2mouse %>% drop_na(human_symbol, mouse_symbol), 
            by=c("Gene"="mouse_symbol"),
            relationship = "many-to-many" # WARNING! There are duplicates!
        )
    
    pan_essential_genes = event_prior_knowledge %>% filter(is_pan_essential) %>% pull(GENE)
    demeter2 = demeter2 %>%
        pivot_longer(-index, names_to="DepMap_ID", values_to="gene_dependency") %>%
        drop_na() %>%
        left_join(ccle_metadata, by="DepMap_ID") %>%
        mutate(is_pan_essential = index %in% pan_essential_genes)
    
    # enrichment
    regulons = regulons %>% 
        bind_rows() %>%
        distinct(regulator, target) %>%
        left_join(annot %>% distinct(EVENT, GENE), by=c("target"="EVENT"))
    
    enrichments = make_enrichments(driver_activity, regulons, ontologies)
    
    enrichments_reactome = lapply(names(enrichments[["reactome"]]), function(gene_set_oi){
            x = enrichments[["reactome"]][[gene_set_oi]] %>%
                as.data.frame() %>%
                mutate(gene_set = gene_set_oi)
            return(x)
        }) %>%
        bind_rows() %>%
        filter(p.adjust < THRESH_FDR) %>%
        rowwise() %>%
        mutate(GeneRatio = eval(parse(text=GeneRatio))) %>%
        ungroup() %>%
        mutate(
            driver_type = case_when(
                gene_set=="suppressors" ~ "Tumor suppressor",
                gene_set=="oncogenics" ~ "Oncogenic"
            )
        )
    
    # roc analysis
    survival_roc_activity = make_roc_analysis(driver_activity, survival_activity) %>%
        mutate(cancers_subset="full") %>%
        rbind(
            make_roc_analysis(
                driver_activity %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]]), 
                survival_activity %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]])
            ) %>%
            mutate(cancers_subset="differential")
        )
    
    survival_roc_genexpr = make_roc_analysis(driver_genexpr, survival_genexpr) %>%
        mutate(cancers_subset="full") %>%
        rbind(
            make_roc_analysis(
                driver_genexpr %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]]), 
                survival_genexpr %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]])
            ) %>%
            mutate(cancers_subset="differential")
        )
    
    survival_roc_genexpr_w_activity_labs = make_roc_analysis(driver_activity, survival_genexpr) %>%
        mutate(cancers_subset="full") %>%
        rbind(
            make_roc_analysis(
                driver_activity %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]]), 
                survival_genexpr %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]])
            ) %>%
            mutate(cancers_subset="differential")
        )
    
    survival_roc_activity_w_genexpr_labs = make_roc_analysis(driver_genexpr, survival_activity) %>%
        mutate(cancers_subset="full") %>%
        rbind(
            make_roc_analysis(
                driver_genexpr %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]]), 
                survival_activity %>% 
                filter(cancer_type %in% diff_activity[["cancer_type"]])
            ) %>%
            mutate(cancers_subset="differential")
        )
    
    # plot
    plts = make_plots(
        diff_activity, diff_genexpr_deseq,
        demeter2, 
        survival_roc_activity, survival_roc_genexpr, 
        survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
        survival_activity, survival_genexpr, 
        survival_activity_conf, survival_genexpr_conf, 
        driver_activity, driver_genexpr_deseq, 
        sf_crossreg_activity, sf_crossreg_genexpr, 
        enrichments_reactome, immune_screen, sf_activity_vs_genexpr, regulons_jaccard,
        n_samples,
        splicing_factors, rbpdb
    )
    
    # make figdata
    figdata = make_figdata(
        diff_activity, diff_genexpr_deseq,
        demeter2, 
        survival_roc_activity, survival_roc_genexpr, 
        survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
        survival_activity, survival_genexpr, 
        survival_activity_conf, survival_genexpr_conf, 
        driver_activity, driver_genexpr_deseq,
        sf_crossreg_activity, sf_crossreg_genexpr, 
        enrichments_reactome, immune_screen, sf_activity_vs_genexpr, regulons_jaccard,
        n_samples,
        splicing_factors, rbpdb
    )

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}