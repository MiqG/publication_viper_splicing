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
# REGINF_DIR = file.path(ROOT,"results","regulon_inference")

# assocs_mi_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","aracne.tsv.gz")
# assocs_spear_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","correlation_spearman.tsv.gz")
# assocs_lm_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","linear_model.tsv.gz")

# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","LIHC.tsv.gz")
# regulators_file = file.path(REGINF_DIR,"files","inputs","LIHC","genexpr_tpm-splicing_factors.tsv")
# targets_file = file.path(REGINF_DIR,"files","inputs","LIHC","event_psi_imputed-exons.tsv")

# summary_stats_genexpr_file = file.path(PREP_DIR,"summary_stats","genexpr_tpm","LIHC.tsv.gz")
# summary_stats_splicing_file = file.path(PREP_DIR,"summary_stats","event_psi_imputed","LIHC-EX.tsv.gz")
# event_annot_file = file.path(RAW_DIR,"VastDB","event_annotation-Hs2.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,'figures','assocs_eda','LIHC')

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
    
    return(plts)
}


plot_examples = function(assocs, regulators, targets, genexpr){
    plts = list()
    
    # pull out some examples SF genexpr vs exon PSI
    ## average gene expression of regulators
    mean_regs = regulators %>%
        filter(regulator %in% assocs[["regulator"]]) %>%
        pivot_longer(-regulator) %>%
        group_by(name) %>%
        summarize(avg_regs = mean(value, na.rm=TRUE)) %>%
        ungroup()

    ## MI vs spearman
    ### low MI and high spearman
    edge_oi = assocs %>%
        mutate(diff = abs(spear_coef) - mutual_information) %>%
        slice_max(diff, n=1) 
    sf_oi = edge_oi %>% pull(regulator)
    exon_oi = edge_oi %>% pull(target)

    x = regulators %>% filter(regulator %in% sf_oi) %>% pivot_longer(-regulator, values_to="reg_genexpr")
    y = targets %>% filter(target %in% exon_oi) %>% pivot_longer(-target, values_to="target_psi")

    df = x %>%
        left_join(y, by="name") %>%
        left_join(edge_oi, by=c("regulator","target")) %>%
        mutate(
            label = sprintf(
                "MI=%s | Spear. Coef.=%s | Spear. p-value=%s | LM Coef.=%s | LM p-value=%s",
                signif(mutual_information,2), signif(spear_coef,2), signif(spear_pvalue,2), 
                signif(lm_coef,2), signif(lm_pvalue,2)),
            x = min(reg_genexpr), 
            y = max(target_psi)
        )

    plts[["examples-low_mi_vs_high_spear-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
        geom_text(
            aes(x=x, y=y, label=label),
            . %>% distinct(x,y,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))
    
    ### high MI and low spearman
    edge_oi = assocs %>%
        mutate(diff = mutual_information - abs(spear_coef)) %>%
        slice_max(diff, n=1) 
    sf_oi = edge_oi %>% pull(regulator)
    exon_oi = edge_oi %>% pull(target)

    x = regulators %>% filter(regulator %in% sf_oi) %>% pivot_longer(-regulator, values_to="reg_genexpr")
    y = targets %>% filter(target %in% exon_oi) %>% pivot_longer(-target, values_to="target_psi")

    df = x %>%
        left_join(y, by="name") %>%
        left_join(edge_oi, by=c("regulator","target")) %>%
        mutate(
            label = sprintf(
                "MI=%s | Spear. Coef.=%s | Spear. p-value=%s | LM Coef.=%s | LM p-value=%s",
                signif(mutual_information,2), signif(spear_coef,2), signif(spear_pvalue,2), 
                signif(lm_coef,2), signif(lm_pvalue,2)),
            x = min(reg_genexpr), 
            y = max(target_psi)+1
        )

    plts[["examples-high_mi_vs_low_spear-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
        geom_text(
            aes(x=x, y=y, label=label),
            . %>% distinct(x,y,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))
    
    
    # ENSG00000011304 or PTBP1 controls inclusion of HsaEX0054252 (exon 2 in PTBP3 or ENSG00000119314)
    ## genexpr SF (regulator) vs exon inclusion (target)
    sf_oi = "ENSG00000011304"
    exon_oi = "HsaEX0054252"
    gene_oi = "ENSG00000119314"
    edge_oi = assocs %>% filter(regulator==sf_oi & target==exon_oi)
    
    x = regulators %>% filter(regulator %in% sf_oi) %>% pivot_longer(-regulator, values_to="reg_genexpr")
    y = targets %>% filter(target %in% exon_oi) %>% pivot_longer(-target, values_to="target_psi")
    z = genexpr %>% filter(ID %in% gene_oi) %>% pivot_longer(-ID, values_to="target_genexpr")
    
    df = x %>%
        left_join(y, by="name") %>%
        left_join(z, by="name") %>%
        left_join(edge_oi, by=c("regulator","target")) %>%
        left_join(mean_regs, by="name") %>%
        mutate(
            target_genexpr_reg_ratio = target_genexpr - avg_regs,
            label = sprintf(
                "MI=%s | Spear. Coef.=%s | Spear. p-value=%s | LM Coef.=%s | LM p-value=%s",
                signif(mutual_information,2), signif(spear_coef,2), signif(spear_pvalue,2), 
                signif(lm_coef,2), signif(lm_pvalue,2)),
            x = min(reg_genexpr), 
            y_psi = max(target_psi)+1,
            y_genexpr = max(target_genexpr)+1
        )
    
    plts[["examples-ptbp1-regulator_genexpr_vs_target_psi-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
        geom_text(
            aes(x=x, y=y_psi, label=label),
            . %>% distinct(x,y_psi,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))

    plts[["examples-ptbp1-regulator_genexpr_vs_target_genexpr-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_genexpr", size=1) + 
        geom_text(
            aes(x=x, y=y_genexpr, label=label),
            . %>% distinct(x,y_genexpr,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s log2(TPM+1)", gene_oi))
    
    plts[["examples-ptbp1-regulator_genexpr_vs_avgreg_genexpr-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="avg_regs", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y="Regulators Avg. log2(TPM+1)")
    
    
    plts[["examples-ptbp1-target_genexpr_vs_avgreg_genexpr-scatter"]] = df %>% 
        ggscatter(x="target_genexpr", y="avg_regs", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", gene_oi), y="Regulators Avg. log2(TPM+1)")
    
    
    # ENSG00000182872 or RBM10 controls inclusion of HsaEX0044216 (exon ? in NUMB or ENSG00000133961)
    ## genexpr SF (regulator) vs exon inclusion (target)
    sf_oi = "ENSG00000182872"
    exon_oi = "HsaEX0044216"
    gene_oi = "ENSG00000133961"
    edge_oi = assocs %>% filter(regulator==sf_oi & target==exon_oi)
    
    x = regulators %>% filter(regulator %in% sf_oi) %>% pivot_longer(-regulator, values_to="reg_genexpr")
    y = targets %>% filter(target %in% exon_oi) %>% pivot_longer(-target, values_to="target_psi")
    z = genexpr %>% filter(ID %in% gene_oi) %>% pivot_longer(-ID, values_to="target_genexpr")
    
    df = x %>%
        left_join(y, by="name") %>%
        left_join(z, by="name") %>%
        left_join(edge_oi, by=c("regulator","target")) %>%
        left_join(mean_regs, by="name") %>%
        mutate(
            target_genexpr_reg_ratio = target_genexpr - avg_regs,
            label = sprintf(
                "MI=%s | Spear. Coef.=%s | Spear. p-value=%s | LM Coef.=%s | LM p-value=%s",
                signif(mutual_information,2), signif(spear_coef,2), signif(spear_pvalue,2), 
                signif(lm_coef,2), signif(lm_pvalue,2)),
            x = min(reg_genexpr), 
            y_psi = max(target_psi)+1,
            y_genexpr = max(target_genexpr)+1
        )
    
    plts[["examples-rbm10-regulator_genexpr_vs_target_psi-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
        geom_text(
            aes(x=x, y=y_psi, label=label),
            . %>% distinct(x,y_psi,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))

    plts[["examples-rbm10-regulator_genexpr_vs_target_genexpr-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_genexpr", size=1) + 
        geom_text(
            aes(x=x, y=y_genexpr, label=label),
            . %>% distinct(x,y_genexpr,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s log2(TPM+1)", gene_oi))
    
    plts[["examples-rbm10-regulator_genexpr_vs_avgreg_genexpr-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="avg_regs", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y="Regulators Avg. log2(TPM+1)")
    
    
    plts[["examples-rbm10-target_genexpr_vs_avgreg_genexpr-scatter"]] = df %>% 
        ggscatter(x="target_genexpr", y="avg_regs", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", gene_oi), y="Regulators Avg. log2(TPM+1)")
    
    # ENSG00000131051 or RBM39 controls inclusion of HsaEX0034998 (exon ? in KRAS or ENSG00000133703)
    ## genexpr SF (regulator) vs exon inclusion (target)
    sf_oi = "ENSG00000131051"
    exon_oi = "HsaEX6071800"
    gene_oi = "ENSG00000258289"
    edge_oi = assocs %>% filter(regulator==sf_oi & target==exon_oi)
    
    x = regulators %>% filter(regulator %in% sf_oi) %>% pivot_longer(-regulator, values_to="reg_genexpr")
    y = targets %>% filter(target %in% exon_oi) %>% pivot_longer(-target, values_to="target_psi")
    z = genexpr %>% filter(ID %in% gene_oi) %>% pivot_longer(-ID, values_to="target_genexpr")
    
    df = x %>%
        left_join(y, by="name") %>%
        left_join(z, by="name") %>%
        left_join(edge_oi, by=c("regulator","target")) %>%
        left_join(mean_regs, by="name") %>%
        mutate(
            target_genexpr_reg_ratio = target_genexpr - avg_regs,
            label = sprintf(
                "MI=%s | Spear. Coef.=%s | Spear. p-value=%s | LM Coef.=%s | LM p-value=%s",
                signif(mutual_information,2), signif(spear_coef,2), signif(spear_pvalue,2), 
                signif(lm_coef,2), signif(lm_pvalue,2)),
            x = min(reg_genexpr), 
            y_psi = max(target_psi)+1,
            y_genexpr = max(target_genexpr)+1
        )
    
    plts[["examples-rbm39-regulator_genexpr_vs_target_psi-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
        geom_text(
            aes(x=x, y=y_psi, label=label),
            . %>% distinct(x,y_psi,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))

    plts[["examples-rbm39-regulator_genexpr_vs_target_genexpr-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_genexpr", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s log2(TPM+1)", gene_oi))
    
    plts[["examples-rbm39-regulator_genexpr_vs_avgreg_genexpr-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="avg_regs", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y="Regulators Avg. log2(TPM+1)")
    
    plts[["examples-rbm39-target_genexpr_vs_avgreg_genexpr-scatter"]] = df %>% 
        ggscatter(x="target_genexpr", y="avg_regs", size=1) + 
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", gene_oi), y="Regulators Avg. log2(TPM+1)")
    
    return(plts)
}


make_plots = function(assocs, summary_stats){
    plts = list(
        plot_assocs(assocs, summary_stats),
        plot_examples()
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(assocs, summary_stats){
    figdata = list(
        "assocs_evaluation" = list(
            "associations" = assocs,
            "summary_statisitcs" = summary_stats,
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
    
    save_plt(plts, "examples-low_mi_vs_high_spear-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "examples-high_mi_vs_low_spear-scatter", '.pdf', figs_dir, width=4, height=4)
    
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
    regulators_file = args[["regulators_file"]]
    targets_file = args[["targets_file"]]
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
    regulators = read_tsv(regulators_file)
    targets = read_tsv(targets_file)
    genexpr = read_tsv(genexpr_file)
    
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
                lm_pvalue = lr_pvalue,
                lm_pearson = pearson_correlation_mean
            ) %>%
            dplyr::select(regulator, target, lm_coef, lm_pvalue, lm_pearson),
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
    
    # plot
    plts = make_plots(assocs, summary_stats)
    gc()
    
    # make figdata
    figdata = make_figdata(assocs, summary_stats)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}