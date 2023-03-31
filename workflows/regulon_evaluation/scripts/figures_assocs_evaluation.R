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

# figs_dir = file.path(RESULTS_DIR,'figures','assocs_evaluation')

##### FUNCTIONS #####
plot_assocs = function(assocs){
    plts = list()
    
    X = assocs %>%
        left_join(
            summary_stats,
            by=c("target"="EVENT")
        )
    
    # distributions of p-values
    plts[["assocs-coefs-lm_pvalue-distr"]] = X %>%
        gghistogram(x="lm_pvalue", bins=50, fill=PAL_DARK, color=NA) +
        theme(aspect.ratio=1) +
        labs(x="LM p-value", y="Count")

    plts[["assocs-coefs-spear_pvalue-distr"]] = X %>%
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
    plts[["assocs-coefs-mi_vs_lm"]] = X %>%
        ggplot(aes(x=mutual_information, y=abs(lm_coef))) +
        geom_scattermore(color=PAL_DARK, pixels=c(1000,1000), pointsize=1, alpha=0.5) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="|LM Coefficient|")

    plts[["assocs-coefs-mi_vs_spear"]] = X %>%
        ggplot(aes(x=mutual_information, y=abs(spear_coef))) +
        geom_scattermore(color=PAL_DARK, pixels=c(1000,1000), pointsize=1, alpha=0.5) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="|Spearman Coef.|")
    
    plts[["assocs-coefs-lm_vs_spear"]] = X %>%
        ggplot(aes(x=lm_coef, y=spear_coef)) +
        geom_scattermore(color=PAL_DARK, pixels=c(1000,1000), pointsize=1, alpha=0.5) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="LM Coef.", y="Spearman Coef.")
    
    # how do p-values compare to each other?
    plts[["assocs-pvalues-mi_vs_lm"]] = X %>%
        drop_na(log_lm_pvalue) %>%
        ggplot(aes(x=mutual_information, y=log_lm_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="-log10(LM p-value)")

    plts[["assocs-pvalues-mi_vs_spear"]] = X %>%
        drop_na(log_lm_pvalue) %>%
        ggplot(aes(x=mutual_information, y=log_spear_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(aspect.ratio=1) +
        labs(x="Mutual Information", y="-log10(Spearman p-value)")
    
    # is the association coefficient related to the target gene expression?
    plts[["assocs-genexpr_vs_mi"]] = X %>%
        group_by(target) %>%
        slice_max(abs(mutual_information), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_median, y=mutual_information)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="Mutual Information")
    
    plts[["assocs-genexpr_vs_spear_coef"]] = X %>%
        group_by(target) %>%
        slice_max(abs(spear_coef), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_median, y=abs(spear_coef))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="|Spearman Coef.|")
    
    plts[["assocs-genexpr_vs_spear_pvalue"]] = X %>%
        group_by(target) %>%
        slice_max(abs(log_spear_pvalue), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_median, y=log_spear_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="-log10(Spearman p-value)")
    
    plts[["assocs-genexpr_vs_lm_coef"]] = X %>%
        group_by(target) %>%
        slice_max(abs(lm_coef), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_median, y=abs(lm_coef))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.5, color=PAL_DARK) +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="|LM Coef.|")
    
    plts[["assocs-genexpr_vs_lm_pvalue"]] = X %>%
        group_by(target) %>%
        slice_max(abs(log_lm_pvalue), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=ENSEMBL_median, y=log_lm_pvalue)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Median log2(TPM+1)", y="-log10(LM p-value)")
    
    plts[["assocs-genexpr_vs_mi-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(abs(mutual_information), n=1) %>%
        ungroup() %>%
        mutate(ENSEMBL_mean_bins = cut(
            ENSEMBL_mean, include.lowest=TRUE, 
            breaks=seq(0,max(ENSEMBL_mean, na.rm=TRUE),length.out=8)
        )) %>%
        ggplot(aes(x=ENSEMBL_mean_bins, y=mutual_information)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(ENSEMBL_mean_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean log2(TPM+1)", y="Mutual Information")
    
    plts[["assocs-genexpr_vs_spear_coef-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(abs(spear_coef), n=1) %>%
        ungroup() %>%
        mutate(ENSEMBL_mean_bins = cut(
            ENSEMBL_mean, include.lowest=TRUE, 
            breaks=seq(0,max(ENSEMBL_mean, na.rm=TRUE)+0.1,length.out=8)
        )) %>%
        ggplot(aes(x=ENSEMBL_mean_bins, y=abs(spear_coef))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(ENSEMBL_mean_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean log2(TPM+1)", y="|Spearman Coef.|")
    
    
    plts[["assocs-genexpr_vs_lm_coef-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(abs(lm_coef), n=1) %>%
        ungroup() %>%
        mutate(ENSEMBL_mean_bins = cut(
            ENSEMBL_mean, include.lowest=TRUE, 
            breaks=seq(0,max(ENSEMBL_mean, na.rm=TRUE)+0.1,length.out=8)
        )) %>%
        ggplot(aes(x=ENSEMBL_mean_bins, y=abs(lm_coef))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(ENSEMBL_mean_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean log2(TPM+1)", y="|LM Coef.|")

    
    # splicing variation vs associations
    plts[["assocs-splicing_vs_mi-box"]] = X %>%
        drop_na(iqr) %>%
        mutate(iqr_bins = cut(
            iqr, include.lowest=TRUE, 
            breaks=seq(0,max(iqr, na.rm=TRUE),length.out=8)
        )) %>%
        ggplot(aes(x=iqr_bins, y=mutual_information)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(iqr_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="IQR PSI", y="Mutual Information")
    
    plts[["assocs-splicing_vs_spear_coef-box"]] = X %>%
        drop_na(iqr) %>%
        mutate(iqr_bins = cut(
            iqr, include.lowest=TRUE, 
            breaks=seq(0,max(iqr, na.rm=TRUE),length.out=8)
        )) %>%
        ggplot(aes(x=iqr_bins, y=abs(spear_coef))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(iqr_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="IQR PSI", y="|Spearman Coef.|")
    
    plts[["assocs-splicing_vs_lm_coef-box"]] = X %>%
        drop_na(iqr) %>%
        mutate(iqr_bins = cut(
            iqr, include.lowest=TRUE, 
            breaks=seq(0,max(iqr, na.rm=TRUE),length.out=8)
        )) %>%
        ggplot(aes(x=iqr_bins, y=abs(lm_coef))) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(iqr_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="IQR PSI", y="|LM Coef.|")
    
    
    # splicing variation vs gene expression
    plts[["assocs-splicing_iqr_vs_genexpr_mean-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(iqr, n=1) %>%
        ungroup() %>%
        drop_na(iqr) %>%
        mutate(iqr_bins = cut(
            iqr, include.lowest=TRUE, 
            breaks=seq(0,max(iqr, na.rm=TRUE),length.out=6)
        )) %>%
        ggplot(aes(x=iqr_bins, y=ENSEMBL_mean)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(iqr_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="IQR PSI", y="Mean log2(TPM+1)")
    
    plts[["assocs-genexpr_mean_vs_splicing_iqr-box"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(iqr, n=1) %>%
        ungroup() %>%
        drop_na(iqr) %>%
        mutate(ENSEMBL_mean_bins = cut(
            ENSEMBL_mean, include.lowest=TRUE, 
            breaks=seq(0,max(ENSEMBL_mean, na.rm=TRUE)+0.1,length.out=6)
        )) %>%
        ggplot(aes(x=ENSEMBL_mean_bins, y=iqr)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(ENSEMBL_mean_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean log2(TPM+1)", y="IQR PSI")
    
    plts[["assocs-genexpr_mean_vs_splicing_iqr-scatter"]] = X %>%
        drop_na(ENSEMBL_mean) %>%
        group_by(target) %>%
        slice_max(iqr, n=1) %>%
        ungroup() %>%
        drop_na(iqr) %>%
        ggplot(aes(x=ENSEMBL_mean, y=iqr)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=1, alpha=0.1, color=PAL_DARK) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_density_2d(size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean log2(TPM+1)", y="IQR PSI")
    
    
    # splicing average vs gene average (sanity check)
    plts[["assocs-splicing_mean_vs_genexpr_mean-box"]] = X %>%
        drop_na(ENSEMBL_mean, EVENT_mean) %>%
        mutate(psi_bins = cut(
            EVENT_mean, include.lowest=TRUE, 
            breaks=seq(0, 100, length.out=8)
        )) %>%
        ggplot(aes(x=psi_bins, y=ENSEMBL_mean)) +
        geom_boxplot(fill=PAL_ACCENT, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(psi_bins) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Mean PSI", y="Mean log2(TPM+1)")
    
    # pull out some examples SF genexpr vs exon PSI (TODO)
    ## low MI and high spearman
    edge_oi = X %>%
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
    
    plts[["assocs-example-low_mi_vs_high_spear-scatter"]] = df %>% 
        ggscatter(x="reg_genexpr", y="target_psi", size=1) + 
        geom_text(
            aes(x=x, y=y, label=label),
            . %>% distinct(x,y,label),
            hjust=0, vjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x=sprintf("%s log2(TPM+1)", sf_oi), y=sprintf("%s PSI", exon_oi))  
    
    return(plts)
}


plot_eval_clip = function(assocs, regulons_clip){
    plts = list()
    
    set.seed(RANDOM_SEED)
    
    X = assocs
    
    print(sprintf("CLIP total interactions: %s", sum(X[["in_clip"]])))
    
    eval_clip = X %>%
        arrange(-abs(mutual_information)) %>%
        mutate(
            ranking = row_number(),
            ranking_ratio = ranking / n(),
            cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
            ranking_var = "mutual_information"
        ) %>%
        filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
        bind_rows(
            X %>%
            arrange(-abs(lm_coef)) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "lm_coef"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1))))
        ) %>%
        bind_rows(
            X %>%
            arrange(lm_pvalue) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "lm_pvalue"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1))))
        ) %>%
        bind_rows(
            X %>%
            arrange(-abs(spear_coef)) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "spear_coef"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1))))
        ) %>%
        bind_rows(
            X %>%
            arrange(spear_pvalue) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "spear_pvalue"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1))))
        ) %>%
        bind_rows(
            X %>%
            arrange(-in_clip) %>% 
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "ground_truth"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1))))
        ) %>%
        bind_rows(
            X %>%
            sample_frac(1L) %>% 
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "random"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1))))
        )
    
    # calculate auc
    eval_clip = eval_clip %>%
        group_by(ranking_var) %>%
        mutate(auc = sum(diff(ranking_ratio) * (head(cumsum_in_clip,-1)+tail(cumsum_in_clip,-1)))/2) %>%
        ungroup()
    
    # do clip interactions tend to have large association values?
    plts[["eval_clip-rankings_by_var-scatter"]] = eval_clip %>%
        mutate(auc_lab = sprintf("AUC(%s)=%s",ranking_var,round(auc,2))) %>%
        ggplot(aes(x=ranking_ratio, y=cumsum_in_clip)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("simpsons") +
        theme_pubr() +
        geom_text(
            aes(x=x, y=y, label=auc_lab, color=ranking_var),
            . %>% distinct(auc, auc_lab, ranking_var) %>% 
                arrange(auc) %>%
                mutate(x=0.65, y=seq(0,0.25,length.out=7)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(x="Ranking Ratio Association", y="TPR CLIP Interactions", color="Association")
    
    plts[["eval_clip-in_clip_vs_mi-box"]] = X %>%
        ggplot(aes(x=in_clip, y=mutual_information)) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.05),
            . %>% count(in_clip) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="Mutual Information")

    plts[["eval_clip-in_clip_vs_spear_coef-box"]] = X %>%
        ggplot(aes(x=in_clip, y=abs(spear_coef))) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(in_clip) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="|Spearman Coef.|")

    plts[["eval_clip-in_clip_vs_spear_pvalue-box"]] = X %>%
        ggplot(aes(x=in_clip, y=-log10(spear_pvalue))) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-1),
            . %>% count(in_clip) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="-log10(Spearman p-value)")

    plts[["eval_clip-in_clip_vs_lm_coef-box"]] = X %>%
        ggplot(aes(x=in_clip, y=abs(log_lm_coef))) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.1),
            . %>% count(in_clip) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="|log10(LM Coefficient+1)|")
        # here it is higher in the FALSE class!

    plts[["eval_clip-in_clip_vs_lm_pvalue-box"]] = X %>%
        ggplot(aes(x=in_clip, y=-log10(lm_pvalue))) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-1),
            . %>% count(in_clip) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="-log10(LM p-value)")
    
    return(plts)
}


plot_eval_pert = function(){
    plts = list()
    
    X = assocs %>%
        left_join(
            regulons_pert,
            by = c("regulator", "target")
        ) %>%
        drop_na()
    
    # are associations informative of perturbations magnitude?
    ## delta PSI
    plts[["eval_pert-delta_psi_vs_mi-scatter"]] = X %>%
        ggplot(aes(x=abs(delta_psi), y=mutual_information)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="|Delta PSI|", y="Mutual Information")
        
    plts[["eval_pert-delta_psi_vs_spear_coef-scatter"]] = X %>%
        ggplot(aes(x=delta_psi, y=spear_coef)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        guides(color="none") +
        labs(x="Delta PSI", y="Spearman Coef.")
        
    plts[["eval_pert-delta_psi_vs_spear_pvalue-scatter"]] = X %>%
        ggplot(aes(x=abs(delta_psi), y=log_spear_pvalue)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="|Delta PSI|", y="-log10(Spearman p-value)")
         
    plts[["eval_pert-delta_psi_vs_lm_coef-scatter"]] = X %>%
        ggplot(aes(x=delta_psi, y=log_lm_coef)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        guides(color="none") +
        labs(x="Delta PSI", y="log10(LM Coefficient+1)")
        
    plts[["eval_pert-delta_psi_vs_lm_pvalue-scatter"]] = X %>%
        ggplot(aes(x=abs(delta_psi), y=log_lm_pvalue)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="|Delta PSI|", y="-log10(LM p-value)")
       
    
    ## relative delta PSI 
    plts[["eval_pert-delta_psi_rel_vs_mi-scatter"]] = X %>%
        ggplot(aes(x=abs(delta_psi_rel), y=mutual_information)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="|Delta PSI Rel.|", y="Mutual Information")
        
    plts[["eval_pert-delta_psi_rel_vs_spear_coef-scatter"]] = X %>%
        ggplot(aes(x=delta_psi_rel, y=spear_coef)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        guides(color="none") +
        labs(x="Delta PSI Rel.", y="Spearman Coef.")
        
    plts[["eval_pert-delta_psi_rel_vs_spear_pvalue-scatter"]] = X %>%
        ggplot(aes(x=abs(delta_psi_rel), y=log_spear_pvalue)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="|Delta PSI Rel.|", y="-log10(Spearman p-value)")
         
    plts[["eval_pert-delta_psi_rel_vs_lm_coef-scatter"]] = X %>%
        ggplot(aes(x=delta_psi_rel, y=log_lm_coef)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        guides(color="none") +
        labs(x="Delta PSI Rel.", y="log10(LM Coefficient+1)")
    
    plts[["eval_pert-delta_psi_rel_vs_lm_coef_custom-scatter"]] = X %>%
        filter(abs(delta_psi)>10) %>%
        mutate(log_lm_coef = ifelse(delta_psi_rel<0, -log_lm_coef, log_lm_coef)) %>%
        ggplot(aes(x=abs(delta_psi_rel), y=log_lm_coef)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_hline(yintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype="dashed", color="black", size=LINE_SIZE) +
        guides(color="none") +
        labs(x="Delta PSI Rel.", y="log10(LM Coefficient+1)")
        
    plts[["eval_pert-delta_psi_rel_vs_lm_pvalue-scatter"]] = X %>%
        ggplot(aes(x=abs(delta_psi_rel), y=log_lm_pvalue)) +
        geom_scattermore(aes(color=cell_line), pixels=c(1000,1000), pointsize=1, alpha=0.1) +
        color_palette("Dark2") +
        theme_pubr() +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="|Delta PSI Rel.|", y="-log10(LM p-value)")
    
    dpsi_thresh = 10
    set.seed(RANDOM_SEED)
    
    X = X %>%
        mutate(is_inter = abs(delta_psi) > dpsi_thresh) # define interactions based on dpsi threshold

    print(sprintf("Total interactions: %s (Thresh = %s)", 
                  X %>% filter(is_inter) %>% count(cell_line, is_inter) %>% pull(n), dpsi_thresh))

    eval_pert = X
    eval_pert = eval_pert %>%
        group_by(cell_line) %>%
        arrange(-abs(mutual_information)) %>%
        mutate(
            ranking = row_number(),
            ranking_ratio = ranking / n(),
            cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
            ranking_var = "mutual_information"
        ) %>%
        filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
        ungroup() %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(abs(lm_coef)) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
                ranking_var = "lm_coef"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(lm_pvalue) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
                ranking_var = "lm_pvalue"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(-abs(spear_coef)) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
                ranking_var = "spear_coef"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(spear_pvalue) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
                ranking_var = "spear_pvalue"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            arrange(-is_inter) %>% 
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
                ranking_var = "ground_truth"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        ) %>%
        bind_rows(
            eval_pert %>%
            group_by(cell_line) %>%
            sample_frac(1L) %>% 
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_is_inter = cumsum(is_inter) / sum(is_inter),
                ranking_var = "random"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        )

    # calculate auc
    eval_pert = eval_pert %>%
        group_by(cell_line, ranking_var) %>%
        mutate(auc = sum(diff(ranking_ratio) * (head(cumsum_is_inter,-1)+tail(cumsum_is_inter,-1)))/2) %>%
        ungroup()

    # do clip interactions tend to have large association values?
    plts[[sprintf("eval_pert-rankings_by_var-scatter-thresh%s",dpsi_thresh)]] = eval_pert %>%
        ggplot(aes(x=ranking_ratio, y=cumsum_is_inter)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("simpsons") +
        theme_pubr() +
        geom_text(
            aes(x=x, y=y, label=auc_lab, color=ranking_var),
            . %>% 
                group_by(cell_line) %>%
                distinct(cell_line, auc, ranking_var) %>% 
                arrange(auc) %>%
                mutate(
                    auc_lab = sprintf("AUC(%s)=%s",ranking_var,round(auc,2)),
                    x=0.35, 
                    y=seq(0,0.25,length.out=7)
                ) %>%
                ungroup(),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Ranking Ratio Association", 
             y=sprintf("Cumulative TPR (|dPSI|>%s)", dpsi_thresh), color="Association")

    plts[[sprintf("eval_pert-is_inter_vs_mi-box-thresh%s",dpsi_thresh)]] = X %>%
        ggplot(aes(x=is_inter, y=mutual_information)) +
        geom_boxplot(aes(fill=is_inter), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.05),
            . %>% count(is_inter) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="Mutual Information")

    plts[[sprintf("eval_pert-is_inter_vs_spear_coef-box-thresh%s",dpsi_thresh)]] = X %>%
        ggplot(aes(x=is_inter, y=abs(spear_coef))) +
        geom_boxplot(aes(fill=is_inter), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.01),
            . %>% count(is_inter) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="|Spearman Coef.|")

    plts[[sprintf("eval_pert-is_inter_vs_spear_pvalue-box-thresh%s",dpsi_thresh)]] = X %>%
        ggplot(aes(x=is_inter, y=-log10(spear_pvalue))) +
        geom_boxplot(aes(fill=is_inter), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-1),
            . %>% count(is_inter) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="-log10(Spearman p-value)")

    plts[[sprintf("eval_pert-is_inter_vs_lm_coef-box-thresh%s",dpsi_thresh)]] = X %>%
        ggplot(aes(x=is_inter, y=abs(log_lm_coef))) +
        geom_boxplot(aes(fill=is_inter), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-0.1),
            . %>% count(is_inter) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="|log10(LM Coefficient+1)|")
        # here it is higher in the FALSE class!

    plts[[sprintf("eval_pert-is_inter_vs_lm_pvalue-box-thresh%s",dpsi_thresh)]] = X %>%
        ggplot(aes(x=is_inter, y=-log10(lm_pvalue))) +
        geom_boxplot(aes(fill=is_inter), outlier.size=0.1) +
        geom_text(
            aes(label=label, y=-1),
            . %>% count(is_inter) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="-log10(LM p-value)")        
    
    # do signs coincide?
    eval_signs = X %>% 
        filter(abs(delta_psi)>10) %>%
        mutate(
            delta_psi = -sign(delta_psi),
            lm_coef = sign(lm_coef),
            spear_coef = sign(spear_coef)
        ) %>% 
        pivot_longer(c(lm_coef, spear_coef))
    
    tests = eval_signs %>%
        group_by(cell_line, name) %>%
        summarize(
            contingency = list(table(delta_psi, value))
        ) %>%
        mutate(
            fisher_test = map(contingency, fisher.test),
            p_value = map_dbl(fisher_test, pluck, "p.value"),
            odds_ratio = map_dbl(fisher_test, pluck, "estimate")
        ) %>% 
        ungroup()
        
    plts[["eval_pert-delta_psi_sign_vs_coefs_sign-bar"]] = eval_signs %>%
        count(cell_line, name, delta_psi, value) %>%
        mutate(label = sprintf("Delta PSI = %s & Coef. = %s", delta_psi, value)) %>%
        ggbarplot(x="label", y="n", fill="name", position=position_dodge(0.9)) +
        facet_wrap(~cell_line) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Sign Combination", y="Count")
    
    # do changes coincide?
    eval_changes = X %>% 
        pivot_longer(c(lm_pvalue, spear_pvalue)) %>%
        mutate(
            delta_psi = abs(delta_psi) > 10,
            value = value < 0.05
        )
    
    tests = eval_changes %>%
        group_by(cell_line, name) %>%
        summarize(
            contingency = list(table(delta_psi, value))
        ) %>%
        mutate(
            fisher_test = map(contingency, fisher.test, alternative="less"), # TRUE & TRUE
            p_value = map_dbl(fisher_test, pluck, "p.value"),
            odds_ratio = map_dbl(fisher_test, pluck, "estimate")
        ) %>% 
        ungroup()
    
    plts[["eval_pert-delta_psi_vs_coefs_pvalue-bar"]] = eval_changes %>%
        count(cell_line, name, delta_psi, value) %>%
        mutate(label = sprintf("Delta PSI High = %s & p-value = %s", delta_psi, value)) %>%
        ggbarplot(x="label", y="n", fill="name", position=position_dodge(0.9)) +
        facet_wrap(~cell_line) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Significant Combination", y="Count")    
        
    
    return(plts)
}


plot_clip_vs_pert = function(){
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
    plts[["clip_vs_pert-delta_psi-boxplot"]] = X %>%
        ggplot(aes(x=in_clip, y=abs(delta_psi))) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="|Delta PSI|") 
    
    plts[["clip_vs_pert-delta_psi_rel-boxplot"]] = X %>%
        ggplot(aes(x=in_clip, y=abs(delta_psi_rel))) +
        geom_boxplot(aes(fill=in_clip), outlier.size=0.1) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        fill_palette(PAL_DUAL) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="CLIP Interaction", y="|Delta PSI Rel.|")
    
    # predictive power of delta PSI and delta PSI rel.
    evaluation = X %>%
        group_by(cell_line) %>%
        arrange(-abs(delta_psi)) %>%
        mutate(
            ranking = row_number(),
            ranking_ratio = ranking / n(),
            cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
            ranking_var = "delta_psi"
        ) %>%
        filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
        ungroup() %>%
        bind_rows(
            X %>%
            group_by(cell_line) %>%
            arrange(-abs(delta_psi_rel)) %>%
            mutate(
                ranking = row_number(),
                ranking_ratio = ranking / n(),
                cumsum_in_clip = cumsum(in_clip) / sum(in_clip),
                ranking_var = "delta_psi_rel"
            ) %>%
            filter(ranking %in% round(quantile(ranking, seq(0,1,0.1)))) %>%
            ungroup()
        )
    
    # calculate auc
    evaluation = evaluation %>%
        group_by(cell_line, ranking_var) %>%
        mutate(auc = sum(diff(ranking_ratio) * (head(cumsum_in_clip,-1)+tail(cumsum_in_clip,-1)))/2) %>%
        ungroup()
    
    # do clip interactions tend to have large association values?
    plts[["clip_vs_pert-tpr-scatter"]] = evaluation %>%
        ggplot(aes(x=ranking_ratio, y=cumsum_in_clip)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("simpsons") +
        theme_pubr() +
        geom_text(
            aes(x=x, y=y, label=auc_lab, color=ranking_var),
            . %>% 
                group_by(cell_line) %>%
                distinct(cell_line, auc, ranking_var) %>% 
                arrange(auc) %>%
                mutate(
                    auc_lab = sprintf("AUC(%s)=%s",ranking_var,signif(auc,2)),
                    x=0.58, 
                    y=seq(0,0.05,length.out=2)
                ) %>%
                ungroup(),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        facet_wrap(~cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Ranking Ratio Association", 
             y=sprintf("Cumulative TPR (|dPSI|>%s)", dpsi_thresh), color="Association")
    
    return(plts)
}


# plot_evaluation = function(evaluation){
#     plts = list()
    
#     X = evaluation
    
#     n_kds = X %>%
#         distinct(sf_target_inference_method, kd_cell_line, KD_ENSEMBL) %>%
#         count(sf_target_inference_method, kd_cell_line) %>%
#         mutate(label=sprintf("%s (n=%s)",kd_cell_line,n))
    
#     plts[["evaluation-thresh_vs_prop_correct-line"]] = X %>%
#         left_join(n_kds, by=c("sf_target_inference_method","kd_cell_line")) %>%
#         ggplot(aes(x=threshold_classification, y=prop_correct, group=KD_ENSEMBL)) +
#         geom_line(size=0.1, color="grey", alpha=0.5) +
#         geom_smooth(aes(color=kd_cell_line, fill=kd_cell_line, group=kd_cell_line), 
#                     se=FALSE, span=0.2, size=LINE_SIZE, linetype="dashed", alpha=0.5, method="loess") +
#         color_palette("Dark2") +
#         theme_pubr(legend="none") +
#         facet_wrap(~label+sf_target_inference_method) +
#         theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="Classification Threshold", y="Proportion Correct") +
#         lims(x=c(0,1), y=c(0,1))
    
#     x = X %>%
#         left_join(n_kds, by=c("kd_cell_line","sf_target_inference_method")) %>% 
#         group_by(sf_target_inference_method, label, kd_cell_line, threshold_classification) %>%
#         summarize(
#             med_tpr=median(tpr, na.rm=TRUE), 
#             med_fpr=median(fpr, na.rm=TRUE),
#             med_recall=median(recall, na.rm=TRUE),
#             med_precision=median(precision, na.rm=TRUE)
#         )
#     plts[["evaluation-fpr_vs_tpr-line"]] = x %>%
#         arrange(threshold_classification) %>%
#         ggplot(aes(x=med_fpr, y=med_tpr)) +
#         geom_line(aes(color=sf_target_inference_method), size=LINE_SIZE, linetype="dashed") +
#         geom_point(aes(color=sf_target_inference_method), size=1) +
#         color_palette("Dark2") +
#         facet_wrap(~label) +
#         theme_pubr() + 
#         theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="FPR", y="TPR", color="Inference Method")
    
#     plts[["evaluation-recall_vs_precision-line"]] = x %>%
#         ggplot(aes(x=med_recall, y=med_precision)) +
#         geom_line(aes(color=sf_target_inference_method), size=LINE_SIZE, linetype="dashed") +
#         geom_point(aes(color=sf_target_inference_method), size=1) +
#         color_palette("Dark2") +
#         facet_wrap(~label) +
#         theme_pubr() + 
#         theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="Recall", y="Precision", color="Inference Method")
    
#     return(plts)
# }


make_plots = function(viper_result, encore_logfc, regulons, sf_info){
    plts = list(
        plot_viper_activities(viper_result, encore_logfc, sf_info),
        plot_regulons(regulons, sf_info)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(viper_result, encore_logfc, regulons, sf_info){
    figdata = list(
        "target_inference" = list(
            "viper_result" = viper_result
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
#     save_plt(plts, "evaluation-thresh_vs_prop_correct-line", '.pdf', figs_dir, width=8, height=8)
#     save_plt(plts, "evaluation-fpr_vs_tpr-line", '.pdf', figs_dir, width=8, height=6)
#     save_plt(plts, "evaluation-recall_vs_precision-line", '.pdf', figs_dir, width=8, height=6)
    
    save_plt(plts, "viper_activities-activity_vs_is_validated-violin", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "viper_activities-activity_ranking_vs_is_validated-violin", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "viper_activities-fc_kds_vs_activity-scatter", '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, "viper_activities-fc_kds_vs_activity_ranking-scatter", '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, "regulons-n_targets_per_reg-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "regulons-pleiotropy-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "regulons-microexons-violin", '.pdf', figs_dir, width=4, height=4)
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
        make_option("--viper_result_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    viper_result_file = args[["viper_result_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    assocs_mi = read_tsv(assocs_mi_file)
    assocs_spear = read_tsv(assocs_spear_file)
    assocs_lm = read_tsv(assocs_lm_file)
    regulons_clip = read_tsv(regulons_clip_file)
    regulons_pert_dpsi = read_tsv(regulons_pert_dpsi_file)
    regulons_pert_dpsi_rel = read_tsv(regulons_pert_dpsi_rel_file)
    regulators = read_tsv(regulators_file)
    targets = read_tsv(targets_file)
    
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
        filter(target %in% events_oi)
    
    regulons_pert = regulons_pert %>%
        bind_rows(
            regulons_pert %>% 
            group_by(cell_line) %>%
            ungroup() %>%
            mutate(cell_line = "merged") %>%
            group_by(regulator, target) %>%
            filter(abs(delta_psi) == max(abs(delta_psi))) %>%
            filter(abs(delta_psi_rel) == max(abs(delta_psi_rel))) %>%
            ungroup() %>%
            distinct()
        )
    
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
    
    # plot
    plts = make_plots(viper_result, encore_logfc, regulons, sf_info)
    
    # make figdata
    figdata = make_figdata(viper_result, encore_logfc, regulons)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}