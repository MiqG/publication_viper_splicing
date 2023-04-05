#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables
RANDOM_SEED = 1234

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "darkgrey"
PAL_ACCENT = "orange"
PAL_DUAL = c(PAL_DARK, PAL_ACCENT)
PAL_CONTRAST = c("darkgrey","darkred")

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","regulon_evaluation")

# assocs_mi_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","aracne.tsv.gz")
# assocs_spear_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","correlation_spearman.tsv.gz")
# assocs_lm_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","linear_model.tsv.gz")

# REGINF_DIR = file.path(ROOT,"results","regulon_inference")
# regulons_clip_file = file.path(REGINF_DIR,"files","regulons","CLIP","POSTAR3.tsv.gz")
# regulons_pert_dpsi_file = file.path(REGINF_DIR,"files","regulons","pert_rnaseq","delta_psi-EX-merged.tsv.gz")
# regulons_pert_dpsi_rel_file = file.path(REGINF_DIR,"files","regulons","pert_rnaseq","delta_psi_rel-EX-merged.tsv.gz")

# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","LIHC.tsv.gz")
# regulators_file = file.path(REGINF_DIR,"files","inputs","LIHC","genexpr_tpm-splicing_factors.tsv")
# targets_file = file.path(REGINF_DIR,"files","inputs","LIHC","event_psi_imputed-exons.tsv")

# summary_stats_genexpr_file = file.path(PREP_DIR,"summary_stats","genexpr_tpm","LIHC.tsv.gz")
# summary_stats_splicing_file = file.path(PREP_DIR,"summary_stats","event_psi_imputed","LIHC-EX.tsv.gz")
# event_annot_file = file.path(RAW_DIR,"VastDB","event_annotation-Hs2.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,'figures','assocs_evaluation','LIHC')

##### FUNCTIONS #####
plot_assocs = function(assocs, summary_stats){
    plts = list()
    
    X = assocs %>%
        left_join(
            summary_stats,
            by=c("target"="EVENT")
        )
    
    # distributions of p-values
    plts[["assocs-pvalues-lm_pvalue-distr"]] = X %>%
        gghistogram(x="lm_pvalue", bins=50, fill=PAL_DARK, color=NA) +
        theme(aspect.ratio=1) +
        labs(x="LM p-value", y="Count")

    plts[["assocs-pvalues-spear_pvalue-distr"]] = X %>%
        gghistogram(x="spear_pvalue", bins=50, fill=PAL_DARK, color=NA) +
        theme(aspect.ratio=1) +
        labs(x="Spearman p-value", y="Count")
        
    # distributions of associations
    plts[["assocs-coefs-mi-distr"]] = X %>%
        gghistogram(x="mutual_information", bins=50, fill=PAL_DARK, color=NA) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="Count")
    
    plts[["assocs-coefs-lm_coef-distr"]] = X %>%
        drop_na(lm_coef) %>%
        gghistogram(x="lm_coef", bins=50, fill=PAL_DARK, color=NA) +
        theme(aspect.ratio=1) +
        labs(x="LM Coefficient", y="Count")

    plts[["assocs-coefs-spear_coef-distr"]] = X %>%
        gghistogram(x="spear_coef", bins=50, fill=PAL_DARK, color=NA) +
        theme(aspect.ratio=1) +
        labs(x="Spearman Coef.", y="Count")
    
    # how do associations compare to each other?
    plts[["assocs-mi_vs_lm_coef"]] = X %>%
        ggplot(aes(x=mutual_information, y=abs(lm_coef))) +
        geom_scattermore(color=PAL_DARK, pixels=c(1000,1000), pointsize=1, alpha=0.5) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="|LM Coefficient|")

    plts[["assocs-mi_vs_spear_coef"]] = X %>%
        ggplot(aes(x=mutual_information, y=abs(spear_coef))) +
        geom_scattermore(color=PAL_DARK, pixels=c(1000,1000), pointsize=1, alpha=0.5) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="|Spearman Coef.|")
    
    plts[["assocs-lm_vs_spear_coef"]] = X %>%
        ggplot(aes(x=lm_coef, y=spear_coef)) +
        geom_scattermore(color=PAL_DARK, pixels=c(1000,1000), pointsize=1, alpha=0.5) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="LM Coef.", y="Spearman Coef.")
    
    # how do p-values compare to each other?
    plts[["assocs-mi_vs_lm_pvalue-scatter"]] = X %>%
        drop_na(log_lm_pvalue) %>%
        ggplot(aes(x=mutual_information, y=log_lm_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="-log10(LM p-value)")

    plts[["assocs-mi_vs_spear_pvalue-scatter"]] = X %>%
        drop_na(log_lm_pvalue) %>%
        ggplot(aes(x=mutual_information, y=log_spear_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="-log10(Spearman p-value)")
    
    # is the association coefficient related to the target gene expression?
    ## scatters
    plts[["assocs-genexpr_vs_mi-scatter"]] = X %>%
        group_by(target) %>%
        slice_max(abs(mutual_information), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_mean, y=mutual_information)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="Mutual Information")
    
    plts[["assocs-genexpr_vs_spear_coef-scatter"]] = X %>%
        group_by(target) %>%
        slice_max(abs(spear_coef), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_mean, y=abs(spear_coef))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="|Spearman Coef.|")
    
    plts[["assocs-genexpr_vs_spear_pvalue-scatter"]] = X %>%
        group_by(target) %>%
        slice_max(abs(log_spear_pvalue), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_mean, y=log_spear_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="-log10(Spearman p-value)")
    
    plts[["assocs-genexpr_vs_lm_coef-scatter"]] = X %>%
        group_by(target) %>%
        slice_max(abs(lm_coef), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_mean, y=abs(lm_coef))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="|LM Coef.|")
    
    plts[["assocs-genexpr_vs_lm_pvalue-scatter"]] = X %>%
        group_by(target) %>%
        slice_max(abs(log_lm_pvalue), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_mean, y=log_lm_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="-log10(LM p-value)")
    
    ### boxplots
    plts[["assocs-genexpr_vs_mi-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(mutual_information, n=1) %>%
        ungroup() %>%
        mutate(assoc_bins = cut(
            mutual_information, include.lowest=TRUE, 
            breaks=seq(0,max(mutual_information, na.rm=TRUE)+0.1,length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=ENSEMBL_mean)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="Mean log2(TPM+1)")
    
    plts[["assocs-genexpr_vs_spear_coef-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(abs(spear_coef), n=1) %>%
        ungroup() %>%
        mutate(assoc_bins = cut(
            abs(spear_coef), include.lowest=TRUE, 
            breaks=seq(0,max(abs(spear_coef), na.rm=TRUE)+0.1,length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=ENSEMBL_mean)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="|Spearman Coef.|", y="Mean log2(TPM+1)")
    
    
    plts[["assocs-genexpr_vs_lm_coef-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(abs(lm_coef), n=1) %>%
        ungroup() %>%
        mutate(assoc_bins = cut(
            abs(lm_coef), include.lowest=TRUE, 
            breaks=seq(0,max(abs(lm_coef), na.rm=TRUE)+0.1,length.out=8)
        )) %>%
        ggplot(aes(x=assoc_bins, y=ENSEMBL_mean)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="|LM Coef.|", y="Mean log2(TPM+1)")

    
    # splicing variation vs associations
    ## scatters
    plts[["assocs-splicing_vs_mi-scatter"]] = X %>%
        ggplot(aes(x=log10(iqr), y=mutual_information)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="log10(IQR PSI)", y="Mutual Information")
    
    plts[["assocs-splicing_vs_spear_coef-scatter"]] = X %>%
        ggplot(aes(x=log10(iqr), y=abs(spear_coef))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="log10(IQR PSI)", y="|Spearman Coef.|")
    
    plts[["assocs-splicing_vs_spear_pvalue-scatter"]] = X %>%
        ggplot(aes(x=log10(iqr), y=log_spear_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="log10(IQR PSI)", y="-log10(Spearman p-value)")
    
    plts[["assocs-splicing_vs_lm_coef-scatter"]] = X %>%
        ggplot(aes(x=log10(iqr), y=abs(lm_coef))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="log10(IQR PSI)", y="|LM Coef.|")
    
    plts[["assocs-splicing_vs_lm_pvalue-scatter"]] = X %>%
        ggplot(aes(x=log10(iqr), y=log_lm_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="log10(IQR PSI)", y="-log10(LM p-value)")
    
    ## boxplots
    plts[["assocs-splicing_vs_mi-box"]] = X %>% 
        drop_na(mutual_information) %>%
        mutate(assoc_bins = cut(
            mutual_information, include.lowest=TRUE, 
            breaks=seq(0,max(mutual_information, na.rm=TRUE),length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=log10(iqr))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="log10(IQR PSI)")
    
    plts[["assocs-splicing_vs_spear_coef-box"]] = X %>% 
        drop_na(spear_coef) %>%
        mutate(assoc_bins = cut(
            abs(spear_coef), include.lowest=TRUE, 
            breaks=seq(0,max(abs(spear_coef), na.rm=TRUE),length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=log10(iqr))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="|Spearman Coef.|", y="log10(IQR PSI)")
    
    plts[["assocs-splicing_vs_spear_pvalue-box"]] = X %>% 
        drop_na(log_spear_pvalue) %>%
        mutate(assoc_bins = cut(
            abs(log_spear_pvalue), include.lowest=TRUE, 
            breaks=seq(0,max(abs(log_spear_pvalue), na.rm=TRUE),length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=log10(iqr))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="-log10(Spearman p-value)", y="log10(IQR PSI)")
    
    plts[["assocs-splicing_vs_lm_coef-box"]] = X %>% 
        drop_na(lm_coef) %>%
        mutate(assoc_bins = cut(
            abs(lm_coef), include.lowest=TRUE, 
            breaks=seq(0,max(abs(lm_coef), na.rm=TRUE),length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=log10(iqr))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(y="log10(IQR PSI)", x="|LM Coef.|")
    
    plts[["assocs-splicing_vs_lm_pvalue-box"]] = X %>% 
        drop_na(log_lm_pvalue) %>%
        mutate(assoc_bins = cut(
            abs(log_lm_pvalue), include.lowest=TRUE, 
            breaks=seq(0,max(abs(log_lm_pvalue), na.rm=TRUE),length.out=11)
        )) %>%
        ggplot(aes(x=assoc_bins, y=log10(iqr))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(assoc_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 70) +
        theme(aspect.ratio=1) +
        labs(x="-log10(LM p-value)", y="log10(IQR PSI)")
    
    # splicing variation vs gene expression
    plts[["assocs-splicing_iqr_vs_genexpr_mean-scatter"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(iqr, n=1) %>%
        ungroup() %>%
        drop_na(iqr) %>%
        ggplot(aes(x=log10(iqr), y=ENSEMBL_mean)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="log10(IQR PSI)", y="Mean log2(TPM+1)")  
    
    # splicing average vs gene average (sanity check)
    plts[["assocs-splicing_mean_vs_genexpr_mean-scatter"]] = X %>%
        drop_na(ENSEMBL_mean, EVENT_mean) %>%
        ggplot(aes(x=EVENT_mean, y=ENSEMBL_mean)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean PSI", y="Mean log2(TPM+1)")
    
    plts[["assocs-splicing_mean_vs_genexpr_mean-box"]] = X %>%
        drop_na(ENSEMBL_mean, EVENT_mean) %>%
        mutate(psi_bins = cut(
            EVENT_mean, include.lowest=TRUE, 
            breaks=seq(0, 100, length.out=8)
        )) %>%
        ggplot(aes(x=psi_bins, y=ENSEMBL_mean)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1, width=0.5) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(psi_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean PSI", y="Mean log2(TPM+1)")
    
    # pull out some examples SF genexpr vs exon PSI
    ## low MI and high spearman
    #     edge_oi = X %>%
    #         mutate(diff = abs(spear_coef) - mutual_information) %>%
    #         slice_max(diff, n=1) 
    #     sf_oi = edge_oi %>% pull(regulator)
    #     exon_oi = edge_oi %>% pull(target)

    #     x = regulators %>% filter(regulator %in% sf_oi) %>% pivot_longer(-regulator, values_to="reg_genexpr")
    #     y = targets %>% filter(target %in% exon_oi) %>% pivot_longer(-target, values_to="target_psi")

    #     df = x %>%
    #         left_join(y, by="name") %>%
    #         left_join(edge_oi, by=c("regulator","target")) %>%
    #         mutate(
    #             label = sprintf(
    #                 "MI=%s | Spear. Coef.=%s | Spear. p-value=%s | LM Coef.=%s | LM p-value=%s",
    #                 signif(mutual_information,2), signif(spear_coef,2), signif(spear_pvalue,2), 
    #                 signif(lm_coef,2), signif(lm_pvalue,2)),
    #             x = min(reg_genexpr), 
    #             y = max(target_psi)
    #         )

    #     plts[["assocs-example-low_mi_vs_high_spear-scatter"]] = df %>% 
    #         ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
    #         geom_text(
    #             aes(x=x, y=y, label=label),
    #             . %>% distinct(x,y,label),
    #             hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
    #         ) +
    #         theme(aspect.ratio=1) +
    #         labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))  
    
    return(plts)
}


compute_precision = function(labels, values, len=11){
    threshs = seq(10, length(values), length.out=len) # they must be ordered

    precisions = sapply(threshs, function(thresh){
        preds = values > values[thresh]
        # how many of predicted TRUE, are TRUE
        TP = sum( labels[preds] )
        # how many of predicted TRUE, are FALSE
        FP = sum( !labels[preds] )
        precision = TP / (TP + FP)
        return(precision)
    })
    
    return(precisions)
}


compute_recall = function(labels, values, len=11){
    threshs = seq(10, length(values), length.out=len) # they must be ordered
    
    recalls = sapply(threshs, function(thresh){
        preds = values > values[thresh]
        # how many of predicted TRUE, are TRUE
        TP = sum( labels[preds] )
        # how many of predicted FALSE, are TRUE
        FN = sum( labels[!preds] )
        recall = TP / (TP + FN)
        return(recall)
    })
    
    return(recalls)
    
}


plot_eval_clip = function(assocs, regulons_clip){
    plts = list()
    
    X = assocs
    
    print(sprintf("CLIP total interactions: %s", sum(X[["in_clip"]])))
    
    eval_clip = X
    eval_clip = eval_clip %>%
        arrange(-abs(mutual_information)) %>%
        reframe(
            precision = compute_precision(in_clip, -abs(mutual_information)),
            recall = compute_recall(in_clip, -abs(mutual_information)),
            ranking_var = "mutual_information"
        ) %>% 
        ungroup() %>%
        bind_rows(
            eval_clip %>%
            arrange(-abs(lm_coef)) %>%
            reframe(
                precision = compute_precision(in_clip, -abs(lm_coef)),
                recall = compute_recall(in_clip, -abs(lm_coef)),
                ranking_var = "lm_coef"
            )
        ) %>%
        bind_rows(
            eval_clip %>%
            arrange(lm_pvalue) %>%
            reframe(
                precision = compute_precision(in_clip, lm_pvalue),
                recall = compute_recall(in_clip, lm_pvalue),
                ranking_var = "lm_pvalue"
            )
        ) %>%
        bind_rows(
            eval_clip %>%
            arrange(-abs(spear_coef)) %>%
            reframe(
                precision = compute_precision(in_clip, -abs(spear_coef)),
                recall = compute_recall(in_clip, -abs(spear_coef)),
                ranking_var = "spear_coef"
            )
        ) %>%
        bind_rows(
            eval_clip %>%
            arrange(spear_pvalue) %>%
            reframe(
                precision = compute_precision(in_clip, spear_pvalue),
                recall = compute_recall(in_clip, spear_pvalue),
                ranking_var = "spear_pvalue"
            )
        ) %>%
        drop_na()
    gc()
    
    # do clip interactions tend to have large association values?   
    plts[["eval_clip-recall_vs_precision-scatter"]] = eval_clip %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("simpsons") +
        theme_pubr() +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Recall", y="Precision", color="Association")
    
    return(plts)
}


plot_eval_pert = function(assocs, regulons_pert){
    plts = list()
    
    X = assocs %>%
        left_join(
            regulons_pert,
            by = c("regulator", "target")
        ) %>%
        drop_na(delta_psi)
    
    print(sprintf("Total interactions: %s (Thresh = %s)", 
          X %>% filter(is_inter) %>% count(cell_line, is_inter) %>% pull(n), 10))

    eval_pert = X
    eval_pert = eval_pert %>%
        group_by(cell_line) %>%
        arrange(-abs(mutual_information)) %>%
        reframe(
            precision = compute_precision(is_inter, -abs(mutual_information)),
            recall = compute_recall(is_inter, -abs(mutual_information)),
            ranking_var = "mutual_information"
        ) %>% 
        ungroup %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(-abs(lm_coef)) %>%
            reframe(
                precision = compute_precision(is_inter, -abs(lm_coef)),
                recall = compute_recall(is_inter, -abs(lm_coef)),
                ranking_var = "lm_coef"
            ) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(lm_pvalue) %>%
            reframe(
                precision = compute_precision(is_inter, lm_pvalue),
                recall = compute_recall(is_inter, lm_pvalue),
                ranking_var = "lm_pvalue"
            ) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(-abs(spear_coef)) %>%
            reframe(
                precision = compute_precision(is_inter, -abs(spear_coef)),
                recall = compute_recall(is_inter, -abs(spear_coef)),
                ranking_var = "spear_coef"
            ) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(spear_pvalue) %>%
            reframe(
                precision = compute_precision(is_inter, spear_pvalue),
                recall = compute_recall(is_inter, spear_pvalue),
                ranking_var = "spear_pvalue"
            ) %>%
            ungroup()
        ) %>%
        drop_na()
    gc()
    
    # do clip interactions tend to have large association values?
    plts[["eval_pert-recall_vs_precision-scatter"]] = eval_pert %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("simpsons") +
        theme_pubr() +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Recall", y="Precision", color="Association")
    
    # are associations informative of perturbations magnitude?
    ## delta PSI
    #     plts[["eval_pert-delta_psi_vs_mi-scatter"]] = X %>%
    #         ggplot(aes(x=abs(delta_psi), y=mutual_information)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         guides(color="none") +
    #         labs(x="|Delta PSI|", y="Mutual Information")

    #     plts[["eval_pert-delta_psi_vs_spear_coef-scatter"]] = X %>%
    #         ggplot(aes(x=delta_psi, y=spear_coef)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         guides(color="none") +
    #         labs(x="Delta PSI", y="Spearman Coef.")

    #     plts[["eval_pert-delta_psi_vs_spear_pvalue-scatter"]] = X %>%
    #         ggplot(aes(x=abs(delta_psi), y=log_spear_pvalue)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         guides(color="none") +
    #         labs(x="|Delta PSI|", y="-log10(Spearman p-value)")

    #     plts[["eval_pert-delta_psi_vs_lm_coef-scatter"]] = X %>%
    #         ggplot(aes(x=delta_psi, y=log_lm_coef)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         guides(color="none") +
    #         labs(x="Delta PSI", y="log10(LM Coefficient+1)")

    #     plts[["eval_pert-delta_psi_vs_lm_pvalue-scatter"]] = X %>%
    #         ggplot(aes(x=abs(delta_psi), y=log_lm_pvalue)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         guides(color="none") +
    #         labs(x="|Delta PSI|", y="-log10(LM p-value)")


    #     ## relative delta PSI 
    #     plts[["eval_pert-delta_psi_rel_vs_mi-scatter"]] = X %>%
    #         ggplot(aes(x=abs(delta_psi_rel), y=mutual_information)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         guides(color="none") +
    #         labs(x="|Delta PSI Rel.|", y="Mutual Information")

    #     plts[["eval_pert-delta_psi_rel_vs_spear_coef-scatter"]] = X %>%
    #         ggplot(aes(x=delta_psi_rel, y=spear_coef)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         guides(color="none") +
    #         labs(x="Delta PSI Rel.", y="Spearman Coef.")

    #     plts[["eval_pert-delta_psi_rel_vs_spear_pvalue-scatter"]] = X %>%
    #         ggplot(aes(x=abs(delta_psi_rel), y=log_spear_pvalue)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         guides(color="none") +
    #         labs(x="|Delta PSI Rel.|", y="-log10(Spearman p-value)")

    #     plts[["eval_pert-delta_psi_rel_vs_lm_coef-scatter"]] = X %>%
    #         ggplot(aes(x=delta_psi_rel, y=log_lm_coef)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         guides(color="none") +
    #         labs(x="Delta PSI Rel.", y="log10(LM Coefficient+1)")

    #     plts[["eval_pert-delta_psi_rel_vs_lm_coef_custom-scatter"]] = X %>%
    #         filter(abs(delta_psi)>10) %>%
    #         mutate(log_lm_coef = ifelse(delta_psi_rel<0, -log_lm_coef, log_lm_coef)) %>%
    #         ggplot(aes(x=abs(delta_psi_rel), y=log_lm_coef)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
    #         guides(color="none") +
    #         labs(x="Delta PSI Rel.", y="log10(LM Coefficient+1)")

    #     plts[["eval_pert-delta_psi_rel_vs_lm_pvalue-scatter"]] = X %>%
    #         ggplot(aes(x=abs(delta_psi_rel), y=log_lm_pvalue)) +
    #         geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
    #         color_palette("Dark2") +
    #         theme_pubr() +
    #         stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
    #         facet_wrap(~cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         guides(color="none") +
    #         labs(x="|Delta PSI Rel.|", y="-log10(LM p-value)")
    
    #     # do signs coincide?
    #     eval_signs_pert = X %>% 
    #         filter(abs(delta_psi)>10) %>%
    #         mutate(
    #             delta_psi = -sign(delta_psi),
    #             lm_coef = sign(lm_coef),
    #             spear_coef = sign(spear_coef)
    #         ) %>% 
    #         pivot_longer(c(lm_coef, spear_coef))

    #     tests = eval_signs_pert %>%
    #         group_by(cell_line, name) %>%
    #         summarize(
    #             contingency = list(table(delta_psi, value))
    #         ) %>%
    #         mutate(
    #             fisher_test = map(contingency, fisher.test),
    #             p_value = map_dbl(fisher_test, pluck, "p.value"),
    #             odds_ratio = map_dbl(fisher_test, pluck, "estimate")
    #         ) %>% 
    #         ungroup()

    #     plts[["eval_pert-delta_psi_sign_vs_coefs_sign-bar"]] = eval_signs_pert %>%
    #         count(cell_line, name, delta_psi, value) %>%
    #         mutate(label = sprintf("Delta PSI = %s & Coef. = %s", delta_psi, value)) %>%
    #         ggbarplot(x="label", y="n", fill="name", position=position_dodge(0.9)) +
    #         facet_wrap(~cell_line) +
    #         theme_pubr(x.text.angle = 70) +
    #         labs(x="Sign Combination", y="Count")

    #     # do changes coincide?
    #     eval_changes = X %>% 
    #         pivot_longer(c(lm_pvalue, spear_pvalue)) %>%
    #         mutate(
    #             delta_psi = abs(delta_psi) > 10,
    #             value = value < 0.05
    #         )

    #     tests = eval_changes %>%
    #         group_by(cell_line, name) %>%
    #         summarize(
    #             contingency = list(table(delta_psi, value))
    #         ) %>%
    #         mutate(
    #             fisher_test = map(contingency, fisher.test, alternative="less"), # TRUE & TRUE
    #             p_value = map_dbl(fisher_test, pluck, "p.value"),
    #             odds_ratio = map_dbl(fisher_test, pluck, "estimate")
    #         ) %>% 
    #         ungroup()

    #     plts[["eval_pert-delta_psi_vs_coefs_pvalue-bar"]] = eval_changes %>%
    #         count(cell_line, name, delta_psi, value) %>%
    #         mutate(label = sprintf("Delta PSI High = %s & p-value = %s", delta_psi, value)) %>%
    #         ggbarplot(x="label", y="n", fill="name", position=position_dodge(0.9)) +
    #         facet_wrap(~cell_line) +
    #         theme_pubr(x.text.angle = 70) +
    #         labs(x="Significant Combination", y="Count")    
        
    
    return(plts)
}


plot_clip_vs_pert = function(regulons_pert, regulons_clip){
    plts = list()
    
    X = regulons_pert %>% 
        filter(cell_line != "merged") %>%
        left_join(
            regulons_clip %>%
            mutate(in_clip = TRUE),
            by=c("regulator"="ENSEMBL","target")
        ) %>%
        mutate(in_clip=replace_na(in_clip, FALSE)) %>%
        drop_na(cell_line) %>%
        distinct(cell_line, regulator, target, delta_psi, delta_psi_rel, in_clip)
    
    # do exons with CLIP signal change more?
    #     plts[["clip_vs_pert-delta_psi-boxplot"]] = X %>%
    #         ggplot(aes(x=in_clip, y=abs(delta_psi))) +
    #         geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
    #         stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
    #         theme_pubr() +
    #         fill_palette(PAL_DUAL) +
    #         guides(fill="none") +
    #         theme(aspect.ratio=1) +
    #         labs(x="CLIP Interaction", y="|Delta PSI|") 

    #     plts[["clip_vs_pert-delta_psi_rel-boxplot"]] = X %>%
    #         ggplot(aes(x=in_clip, y=abs(delta_psi_rel))) +
    #         geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
    #         stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
    #         theme_pubr() +
    #         fill_palette(PAL_DUAL) +
    #         guides(fill="none") +
    #         theme(aspect.ratio=1) +
    #         labs(x="CLIP Interaction", y="|Delta PSI Rel.|")
    
    # predictive power of delta PSI and delta PSI rel.
    evaluation = X %>%
        group_by(cell_line) %>%
        arrange(-abs(delta_psi)) %>%
        reframe(
            precision = compute_precision(in_clip, -abs(delta_psi), 25),
            recall = compute_recall(in_clip, -abs(delta_psi), 25),
            ranking_var = "delta_psi"
        ) %>% 
        ungroup() %>%
        bind_rows(
            X %>%
            group_by(cell_line) %>%
            arrange(-abs(delta_psi_rel)) %>%
            reframe(
                precision = compute_precision(in_clip, -abs(delta_psi_rel), 25),
                recall = compute_recall(in_clip, -abs(delta_psi_rel), 25),
                ranking_var = "delta_psi_rel"
            ) %>%
            ungroup()
        ) %>% drop_na()
    
    # do clip interactions tend to have large association values?
    plts[["eval_clip_vs_pert-recall_vs_precision-scatter"]] = evaluation %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("simpsons") +
        theme_pubr() +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Recall", y="Precision", color="Association")

    return(plts)
}


make_plots = function(assocs, summary_stats, regulons_clip, regulons_pert){
    plts = list(
        plot_assocs(assocs, summary_stats),
        plot_eval_clip(assocs, regulons_clip),
        plot_eval_pert(assocs, regulons_pert),
        plot_clip_vs_pert(regulons_pert, regulons_clip)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(assocs, summary_stats, regulons_clip, regulons_pert){
    figdata = list(
        "assocs_evaluation" = list(
            "associations" = assocs,
            "summary_statisitcs" = summary_stats,
            "ground_truth_clip" = regulons_clip,
            "ground_truth_pert" = regulons_pert
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
    save_plt(plts, "assocs-pvalues-lm_pvalue-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-pvalues-spear_pvalue-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-coefs-mi-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-coefs-lm_coef-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-coefs-spear_coef-distr", '.pdf', figs_dir, width=4, height=4)
    
    save_plt(plts, "assocs-mi_vs_lm_coef", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-mi_vs_spear_coef", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-lm_vs_spear_coef", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-mi_vs_lm_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-mi_vs_spear_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    
    save_plt(plts, "assocs-genexpr_vs_mi-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_spear_coef-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_spear_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_lm_coef-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_lm_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_mi-box", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_spear_coef-box", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-genexpr_vs_lm_coef-box", '.pdf', figs_dir, width=4, height=4)
    
    save_plt(plts, "assocs-splicing_vs_mi-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_spear_coef-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_spear_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_lm_coef-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_lm_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_mi-box", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_spear_coef-box", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_lm_coef-box", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_spear_pvalue-box", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_vs_lm_pvalue-box", '.pdf', figs_dir, width=4, height=4)

    save_plt(plts, "assocs-splicing_iqr_vs_genexpr_mean-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_mean_vs_genexpr_mean-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "assocs-splicing_mean_vs_genexpr_mean-box", '.pdf', figs_dir, width=4, height=4)
    
    save_plt(plts, "eval_clip-recall_vs_precision-scatter", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "eval_pert-recall_vs_precision-scatter", '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, "eval_clip_vs_pert-recall_vs_precision-scatter", '.pdf', figs_dir, width=10, height=6)
    
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
        make_option("--assocs_mi_file", type="character"),
        make_option("--assocs_spear_file", type="character"),
        make_option("--assocs_lm_file", type="character"),
        make_option("--regulons_clip_file", type="character"),
        make_option("--regulons_pert_dpsi_file", type="character"),
        make_option("--regulons_pert_dpsi_rel_file", type="character"),
        make_option("--regulators_file", type="character"),
        make_option("--targets_file", type="character"),
        make_option("--summary_stats_genexpr_file", type="character"),
        make_option("--summary_stats_splicing_file", type="character"),
        make_option("--event_annot_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    assocs_mi_file = args[["assocs_mi_file"]]
    assocs_spear_file = args[["assocs_spear_file"]]
    assocs_lm_file = args[["assocs_lm_file"]]
    regulons_clip_file = args[["regulons_clip_file"]]
    regulons_pert_dpsi_file = args[["regulons_pert_dpsi_file"]]
    regulons_pert_dpsi_rel_file = args[["regulons_pert_dpsi_rel_file"]]
    #regulators_file = args[["regulators_file"]]
    #targets_file = args[["targets_file"]]
    summary_stats_genexpr_file = args[["summary_stats_genexpr_file"]]
    summary_stats_splicing_file = args[["summary_stats_splicing_file"]]
    event_annot_file = args[["event_annot_file"]]
    figs_dir = args[["figs_dir"]]
    
    set.seed(RANDOM_SEED)
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    assocs_mi = read_tsv(assocs_mi_file)
    assocs_spear = read_tsv(assocs_spear_file)
    assocs_lm = read_tsv(assocs_lm_file)
    regulons_clip = read_tsv(regulons_clip_file)
    regulons_pert_dpsi = read_tsv(regulons_pert_dpsi_file)
    regulons_pert_dpsi_rel = read_tsv(regulons_pert_dpsi_rel_file)
    #regulators = read_tsv(regulators_file)
    #targets = read_tsv(targets_file)
    
    summary_stats_genexpr = read_tsv(summary_stats_genexpr_file)
    summary_stats_splicing = read_tsv(summary_stats_splicing_file)
    event_annot = read_tsv(event_annot_file)
    
    gc()
    
    # prep
    ## summary stats
    summary_stats_splicing = summary_stats_splicing %>% 
        mutate(cv = EVENT_std / EVENT_mean, 
               iqr = EVENT_q75 - EVENT_q25) %>%
        drop_na(iqr)
    
    summary_stats = summary_stats_splicing %>%
        left_join(
            event_annot,
            by="EVENT"
        ) %>%
        left_join(
            summary_stats_genexpr,
            by=c("ENSEMBL"="ID")
        )
    
    ## high-variant exons (TO REMOVE)
    events_oi = summary_stats_splicing %>% 
        #filter(iqr >= 1) %>% 
        pull(EVENT)
    
    ## merge sf-exon associations
    assocs = assocs_mi %>%
        rename(mutual_information = association) %>%
        # linear model
        left_join(
            assocs_lm %>%
            rename(
                lm_coef = target_coefficient_mean,
                lm_pvalue = lr_pvalue
            ) %>%
            dplyr::select(regulator, target, lm_coef, lm_pvalue),
            by = c("regulator","target")
        ) %>%
        # spearman coefficient
        left_join(
            assocs_spear %>%
            rename(
                spear_coef = statistic,
                spear_pvalue = pvalue
            ) %>%
            dplyr::select(regulator, target, spear_coef, spear_pvalue),
            by = c("regulator","target")
        ) %>%
        filter(target %in% events_oi)
    gc()
    
    ## merge perturbation regulons
    regulons_pert = regulons_pert_dpsi %>%
        mutate(delta_psi = likelihood*tfmode) %>%
        dplyr::select(cell_line, regulator, target, delta_psi) %>%
        left_join(
            regulons_pert_dpsi_rel %>%
            mutate(delta_psi_rel = likelihood*tfmode) %>%
            dplyr::select(cell_line, regulator, target, delta_psi_rel),
            by = c("regulator","target","cell_line")
        ) %>%
        mutate(is_inter = abs(delta_psi) > 10) %>% # define interactions based on dpsi threshold
        filter(target %in% events_oi)

    cell_lines = regulons_pert %>% pull(cell_line) %>% unique()
    
    merged_regulons_pert = regulons_pert %>% 
            pivot_wider(names_from="cell_line", values_from="delta_psi")
    merged_regulons_pert[["delta_psi"]] = apply(
        merged_regulons_pert[, cell_lines], 1, function(x){
            x[which.max(abs(x))]
        })
    merged_regulons_pert = merged_regulons_pert %>%
            mutate(
                cell_line = "merged",
                is_inter = abs(delta_psi) > 10
            ) %>%
            distinct(regulator, target, delta_psi, is_inter, cell_line)
    
    regulons_pert = regulons_pert %>% 
        bind_rows(merged_regulons_pert)
    
    gc()
    
    ## add regulons to assocs
    assocs = assocs %>%
        left_join(
            regulons_clip %>% 
            mutate(
                in_clip = TRUE,
                regulator = ENSEMBL
            ) %>%
            dplyr::select(regulator, target, in_clip),
            by=c("regulator","target")
        ) %>%
        mutate(in_clip = replace_na(in_clip, FALSE)) %>%
        mutate(
            log_lm_coef = sign(lm_coef)*log10(abs(lm_coef)+1),
            log_lm_pvalue = -log10(lm_pvalue),
            log_spear_pvalue = -log10(spear_pvalue)
        ) %>%
        filter(target %in% events_oi)
    gc()
    
    # plot
    plts = make_plots(assocs, summary_stats, regulons_clip, regulons_pert)
    gc()
    
    # make figdata
    figdata = make_figdata(assocs, summary_stats, regulons_clip, regulons_pert)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}