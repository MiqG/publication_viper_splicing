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
library(pROC)

# variables
THRESH_FDR = 0.05
ORGDB = org.Hs.eg.db
THRESH_N_SUM = 5.5

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = setNames(
    c("#6C98B3","#F6AE2D"),
    c("Tumor suppressor","Oncogenic")
)

PAL_CANCER_TYPES = setNames(
    c("#A6CEE3","#3385BB","#84BF96","#6DBD57","#7F9D55","#F57C7C","#E42622",
      "#FBB268","#FE8D19","#DE9E83","#9D7BBA","#977899","#F3E587","#B15928"),
    c("BLCA","BRCA","COAD","HNSC","KICH","KIRC","KIRP",
      "LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_tcga")
# REGINF_DIR = file.path(ROOT,"results","regulon_inference")
# diff_activity_file = file.path(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# diff_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN','genexpr_tpm-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","cancer_program")
# assocs_gene_dependency_file = file.path(RESULTS_DIR,"..","sf_activity_ccle","files","protein_activity_vs_demeter2","CCLE.tsv.gz")
# survival_activity_file = file.path(RESULTS_DIR,'files','PANCAN',"protein_activity-survival_analysis-surv.tsv.gz")
# survival_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-survival_analysis-surv.tsv.gz")
# sf_crossreg_activity_file = file.path(RESULTS_DIR,'files','PANCAN',"protein_activity-sf_cross_regulation.tsv.gz")
# sf_crossreg_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-sf_cross_regulation.tsv.gz")
# ontology_chea_file = file.path(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz")
# sf_activity_vs_genexpr_file = file.path(RESULTS_DIR,'files','PANCAN',"genexpr_tpm_vs_activity.tsv.gz")
# protein_activity_stn_file = file.path(RESULTS_DIR,"files","protein_activity","PANCAN-SolidTissueNormal-EX.tsv.gz")
# genexpr_tpm_stn_file = file.path(PREP_DIR,"genexpr_tpm","PANCAN-SolidTissueNormal.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","PANCAN.tsv.gz")
# regulons_jaccard_file = file.path(REGINF_DIR,"files","regulons_eda_jaccard","experimentally_derived_regulons_pruned-EX.tsv.gz")

##### FUNCTIONS #####
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
        ggscatter(x="activity", y="genexpr_log2FC", alpha=0.5, size=1) +
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
        ggbarplot(x="GENE", y="n_sign", fill="driver_type", color=NA) +
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
            angle=-45, hjust=0, vjust=1, nudge_y=-0.25
        ) +
        geom_hline(yintercept=c(-THRESH_N_SUM,THRESH_N_SUM), linetype="dashed", color="black", size=LINE_SIZE) +
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
                    n_sign=c(10,-10),
                    x=c(120,40)
                ),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=70) +
        labs(x="Splicing Factor", y="Count", fill="Driver Type")
    
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
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
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
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
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
        ggbarplot(x="GENE", y="n_sign", fill="driver_type", color=NA) +
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
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
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
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Cancer Type", y="median(Protein Activity)", fill="Driver Type")
    
    return(plts)
}


plot_prolif_driver = function(diff_activity, assocs_gene_dependency){
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
    x = X %>%
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>% 
        left_join(
            assocs_gene_dependency,
            by="GENE"
        )
    
    plts[["prolif_driver-driver_type_vs_demeter2-violin"]] = x %>%
        ggviolin(x="driver_type", y="pearson_coef", color=NA, fill="driver_type", palette=PAL_DRIVER_TYPE, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text_repel(
            aes(label=GENE),
            x %>% slice_max(pearson_coef, n=5) %>% bind_rows(x %>% slice_min(pearson_coef, n=5)),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1
        ) +
        geom_text(
            aes(y=-0.3, label=label),
            . %>% count(driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="Driver Type", y="Demeter2 Pearson Correlation")
    
    
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


plot_survival_analysis = function(survival_roc, survival_omic, driver_omic, patt=""){
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
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
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
        stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
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
            method="wilcox.test", label="p.signif", ref.group="Oncogenic\n&\nTumor suppressor",
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
            method="wilcox.test", label="p.signif", ref.group="Oncogenic\n&\nTumor suppressor",
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


make_enrichments = function(driver_classif, ontology_chea){
    
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
    
    oncogenics = X %>% filter(driver_type=="Oncogenic") %>% pull(GENE)
    suppressors = X %>% filter(driver_type=="Tumor suppressor") %>% pull(GENE)
    tf_enrichments = list(
        "oncogenics" = enricher(oncogenics, TERM2GENE=ontology_chea),
        "suppressors" = enricher(suppressors, TERM2GENE=ontology_chea)
    )
    
    return(tf_enrichments)
}


plot_tf_enrichments = function(tf_enrichments){
    plts = list()
    
    X = tf_enrichments
    
    plts[["tf_enrichments-oncogenic-cnet"]] = X[["oncogenics"]] %>% 
        cnetplot(cex_label_category=0.5, cex_label_gene=0.5, 
                 cex_family_category=FONT_FAMILY, cex_family_category=FONT_FAMILY) +
        guides(size="none") +
        theme(aspect.ratio=1)
    
    plts[["tf_enrichments-suppressor-cnet"]] = X[["suppressors"]] %>% 
        cnetplot(cex_label_category=0.5, cex_label_gene=0.5, 
                 cex_family_category=FONT_FAMILY, cex_family_category=FONT_FAMILY) +
        guides(size="none") +
        theme(aspect.ratio=1)
    
    return(plts)
}


plot_examples = function(){
    # SRSF1 is a driver of breast cancer and pancreatic cancer
    diff_activity %>% filter(GENE == "SRSF1")
    diff_genexpr %>% filter(GENE == "SRSF1")
    
    # 
}


plot_drivers_activity_stn = function(protein_activity_stn){
    plts = list()
    
    X = protein_activity_stn
    
    # median activity of oncogenic and tumor suppressor SFs per sample across cancers
    cancer_order = X %>%
        group_by(cancer_type) %>%
        summarize(MKI67 = median(MKI67, na.rm=TRUE)) %>%
        ungroup() %>%
        arrange(MKI67) %>%
        pull(cancer_type)
    
    plts[["drivers_activity_stn-cancer_vs_MKI67-line"]] = X %>%
        mutate(cancer_type = factor(cancer_type, levels=cancer_order)) %>%
        ggplot(aes(x=cancer_type, y=MKI67, group=1)) +
        geom_smooth(color="black", fill="lightgrey", size=LINE_SIZE, linetype="dashed") + 
        theme_pubr(x.text.angle = 70) +
        labs(x="STN Tissue Type", y="log2(TPM+1) of MKI67")

    plts[["drivers_activity_stn-cancer_vs_oncogenic-line"]] = X %>%
        group_by(cancer_type, driver_type, GENE) %>%
        summarize(
            activity = median(activity, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        mutate(cancer_type = factor(cancer_type, levels=cancer_order)) %>%
        ggplot(aes(x=cancer_type, y=activity, group=driver_type)) +
        geom_smooth(aes(color=driver_type), size=LINE_SIZE, linetype="dashed", fill="lightgrey") +
        color_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        labs(x="STN Tissue Type", y="median(Protein Activity per SF)", color="Driver Type")
    
    return(plts)
}


make_plots = function(
    diff_activity, diff_genexpr,
    assocs_gene_dependency, 
    survival_roc_activity, survival_roc_genexpr, 
    survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
    survival_activity, survival_genexpr, 
    driver_activity, driver_genexpr, 
    sf_crossreg_activity, sf_crossreg_genexpr, 
    tf_enrichments, sf_activity_vs_genexpr,
    protein_activity_stn, regulons_jaccard
){
    plts = list(
        plot_driver_selection(driver_activity, driver_genexpr, diff_activity, diff_genexpr),
        plot_prolif_driver(diff_activity, assocs_gene_dependency),
        plot_survival_analysis(survival_roc_activity, survival_activity, driver_activity, "-activity"),
        plot_survival_analysis(survival_roc_genexpr, survival_genexpr, driver_genexpr, "-genexpr"),
        plot_survival_analysis(survival_roc_genexpr_w_activity_labs, survival_genexpr, driver_activity, "-genexpr_w_activity_labs"),
        plot_survival_analysis(survival_roc_activity_w_genexpr_labs, survival_activity, driver_genexpr, "-activity_w_genexpr_labs"),
        plot_sf_crossreg(driver_activity, sf_crossreg_activity, regulons_jaccard, "-activity"),
        plot_sf_crossreg(driver_genexpr, sf_crossreg_genexpr, regulons_jaccard, "-genexpr"),
        plot_tf_enrichments(tf_enrichments),
        plot_comparison(diff_activity, diff_genexpr, survival_activity, survival_genexpr, sf_activity_vs_genexpr),
        plot_drivers_activity_stn(protein_activity_stn)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    diff_activity, diff_genexpr,
    assocs_gene_dependency, 
    survival_roc_activity, survival_roc_genexpr, 
    survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
    survival_activity, survival_genexpr, 
    driver_activity, driver_genexpr, 
    sf_crossreg_activity, sf_crossreg_genexpr, 
    tf_enrichments, sf_activity_vs_genexpr,
    protein_activity_stn, regulons_jaccard
){
    figdata = list(
        "tcga_tumorigenesis" = list(
            "diff_protein_activity" = diff_activity,
            "assocs_gene_dependency" = assocs_gene_dependency,
            "survival_roc_activity" = survival_roc_activity,
            "sf_cross_regulation" = sf_crossreg_activity
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
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-genexpr-bar", '.pdf', figs_dir, width=12, height=9)
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-genexpr-bar", '.pdf', figs_dir, width=12, height=9)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-activity_drivers_vs_activity-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-activity_drivers_vs_genexpr-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_genexpr-violin", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "driver_selection-drivers_vs_cancer_type-genexpr_drivers_vs_activity-violin", '.pdf', figs_dir, width=8, height=6)
    
    save_plt(plts, "prolif_driver-driver_type_vs_demeter2-violin", '.pdf', figs_dir, width=5, height=5)
    
    save_plt(plts, "survival_analysis-cancers_all-roc_curves-activity", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-roc_curves-activity", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all-violin-activity", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-violin-activity", '.pdf', figs_dir, width=5, height=6)
    
    save_plt(plts, "survival_analysis-cancers_all-roc_curves-genexpr", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-roc_curves-genexpr", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all-violin-genexpr", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-violin-genexpr", '.pdf', figs_dir, width=5, height=6)

    save_plt(plts, "survival_analysis-cancers_all-roc_curves-genexpr_w_activity_labs", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-roc_curves-genexpr_w_activity_labs", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_all-roc_curves-activity_w_genexpr_labs", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-roc_curves-activity_w_genexpr_labs", '.pdf', figs_dir, width=5, height=6)
    
    save_plt(plts, "sf_cross_regulation-correlations-violin-activity", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-regulon_similarity-violin-activity", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-correlations-violin-genexpr", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-regulon_similarity-violin-genexpr", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "tf_enrichments-oncogenic-cnet", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "tf_enrichments-suppressor-cnet", '.pdf', figs_dir, width=4, height=4)
    
    save_plt(plts, "comparison-diff_analysis-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "comparison-survival-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "comparison-correlation_by_cancer-violin", '.pdf', figs_dir, width=12, height=5)
    save_plt(plts, "comparison-correlation_median_pancan-violin", '.pdf', figs_dir, width=4, height=4)
    
    save_plt(plts, "drivers_activity_stn-cancer_vs_MKI67-line", '.pdf', figs_dir, width=6, height=3)
    save_plt(plts, "drivers_activity_stn-cancer_vs_oncogenic-line", '.pdf', figs_dir, width=6, height=6)
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
        make_option("--survival_activity_file", type="character"),
        make_option("--survival_genexpr_file", type="character"),
        make_option("--sf_crossreg_activity_file", type="character"),
        make_option("--sf_crossreg_genexpr_file", type="character"),
        make_option("--assocs_gene_dependency_file", type="character"),
        make_option("--ontology_chea_file", type="character"),
        make_option("--sf_activity_vs_genexpr_file", type="character"),
        make_option("--genexpr_tpm_stn_file", type="character"),
        make_option("--protein_activity_stn_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--regulons_jaccard_file", type="character"),
        make_option("--gene_annotation_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    diff_activity_file = args[["diff_activity_file"]]
    diff_genexpr_file = args[["diff_genexpr_file"]]
    survival_activity_file = args[["survival_activity_file"]]
    survival_genexpr_file = args[["survival_genexpr_file"]]
    sf_crossreg_activity_file = args[["sf_crossreg_activity_file"]]
    sf_crossreg_genexpr_file = args[["sf_crossreg_genexpr_file"]]
    assocs_gene_dependency_file = args[["assocs_gene_dependency_file"]]
    ontology_chea_file = args[["ontology_chea_file"]]
    sf_activity_vs_genexpr_file = args[["sf_activity_vs_genexpr_file"]]
    genexpr_tpm_stn_file = args[["genexpr_tpm_stn_file"]]
    protein_activity_stn_file = args[["protein_activity_stn_file"]]
    metadata_file = args[["metadata_file"]]
    regulons_jaccard_file = args[["regulons_jaccard_file"]]
    gene_annotation_file = args[["gene_annotation_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_activity = read_tsv(diff_activity_file)
    diff_genexpr = read_tsv(diff_genexpr_file)
    survival_activity = read_tsv(survival_activity_file)
    survival_genexpr = read_tsv(survival_genexpr_file)
    sf_crossreg_activity = read_tsv(sf_crossreg_activity_file)
    sf_crossreg_genexpr = read_tsv(sf_crossreg_genexpr_file)
    assocs_gene_dependency = read_tsv(assocs_gene_dependency_file)
    ontology_chea = read.gmt(ontology_chea_file)
    sf_activity_vs_genexpr = read_tsv(sf_activity_vs_genexpr_file)
    genexpr_tpm_stn = read_tsv(genexpr_tpm_stn_file)
    protein_activity_stn = read_tsv(protein_activity_stn_file)
    metadata = read_tsv(metadata_file)
    regulons_jaccard = read_tsv(regulons_jaccard_file)
    gene_annotation = read_tsv(gene_annotation_file) %>%
        dplyr::rename(
            GENE = `Approved symbol`,
            ENSEMBL = `Ensembl gene ID`
        )
    
    # prep
    diff_activity = diff_activity %>%
        mutate(
            is_significant = padj < THRESH_FDR
        ) %>%
        dplyr::rename("ENSEMBL"="regulator") %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = "ENSEMBL"
        )
    
    diff_genexpr = diff_genexpr %>%
        dplyr::rename("ENSEMBL"="ID") %>%
        filter(ENSEMBL %in% diff_activity[["ENSEMBL"]]) %>%
        mutate(
            padj = p.adjust(pvalue, method="fdr"),
            is_significant = padj < THRESH_FDR
        ) %>%
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
    
    survival_activity = survival_activity %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("feature"="ENSEMBL")
        )    
    
    survival_genexpr = survival_genexpr %>%
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
    
    genexpr_mki67 = genexpr_tpm_stn %>%
        filter(ID == "ENSG00000148773") %>%
        pivot_longer(-ID, names_to="sampleID", values_to="MKI67")
    
    protein_activity_stn = driver_activity %>%
        count(GENE, ENSEMBL, driver_type) %>%
        group_by(GENE, ENSEMBL) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup() %>% 
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE, ENSEMBL) %>%
        slice_max(n, n=1) %>%
        ungroup() %>%
        left_join(
            protein_activity_stn %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity"),
            by=c("ENSEMBL"="regulator")
        ) %>%
        left_join(
            genexpr_mki67,
            by="sampleID"
        ) %>%
        left_join(
            metadata,
            by="sampleID"
        )
    
    # enrichment
    tf_enrichments = make_enrichments(driver_activity, ontology_chea)
    
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
        diff_activity, diff_genexpr,
        assocs_gene_dependency, 
        survival_roc_activity, survival_roc_genexpr, 
        survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
        survival_activity, survival_genexpr, 
        driver_activity, driver_genexpr, 
        sf_crossreg_activity, sf_crossreg_genexpr, 
        tf_enrichments, sf_activity_vs_genexpr,        
        protein_activity_stn, regulons_jaccard
    )
    
    # make figdata
    figdata = make_figdata(
        diff_activity, diff_genexpr,
        assocs_gene_dependency, 
        survival_roc_activity, survival_roc_genexpr,
        survival_roc_genexpr_w_activity_labs, survival_roc_activity_w_genexpr_labs,
        survival_activity, survival_genexpr, 
        driver_activity, driver_genexpr, 
        sf_crossreg_activity, sf_crossreg_genexpr, 
        tf_enrichments, sf_activity_vs_genexpr,
        protein_activity_stn, regulons_jaccard
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