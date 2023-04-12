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
# assocs_lm2_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","linear_model2.tsv.gz")

# REGINF_DIR = file.path(ROOT,"results","regulon_inference")
# regulons_clip_file = file.path(REGINF_DIR,"files","regulons","CLIP","POSTAR3.tsv.gz")
# regulons_pert_dpsi_file = file.path(REGINF_DIR,"files","regulons","pert_rnaseq","delta_psi-EX-merged.tsv.gz")
# regulons_pert_dpsi_rel_file = file.path(REGINF_DIR,"files","regulons","pert_rnaseq","delta_psi_rel-EX-merged.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,'figures','assocs_evaluation','LIHC')

##### FUNCTIONS #####
compute_precision = function(labels, values, len=11){
    threshs = round(seq(10, length(values), length.out=len)) # they must be ordered

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
    threshs = round(seq(10, length(values), length.out=len)) # they must be ordered
    
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


compute_prec_rec = function(.data, true_var, pred_var, label, len=11){
    df = .data %>%
        arrange({{pred_var}}) %>%
        drop_na({{pred_var}}) %>%
        reframe(
            precision = compute_precision({{true_var}}, {{pred_var}}, len),
            recall = compute_recall({{true_var}}, {{pred_var}}, len),
            ranking_var = label,
            n_obs = n()
        )
    return(df)
}


compute_prec_rec_mult = function(.data, true_var, vars_oi, len=11){
    
    result = lapply(
        vars_oi, function(var_oi){
            result = .data %>% compute_prec_rec(.data[[true_var]], .data[[var_oi]], var_oi)
            
            return(result)
        }) %>%
        do.call(rbind, .) %>%
        drop_na()
    
    return(result)
}


plot_eval_clip = function(assocs, regulons_clip){
    plts = list()
    
    X = assocs
    
    print(sprintf("CLIP total interactions: %s", sum(X[["in_clip"]])))
    
    vars_oi = c("mutual_information",
                "abs_lm_coef","abs_lm_pearson","lm_pearson","lm_pvalue",
                "abs_lm2_coef","abs_lm2_pearson","lm2_pearson","lm2_pvalue",
                "abs_spear_coef","spear_pvalue")

    eval_clip = X %>% compute_prec_rec_mult("in_clip", vars_oi)
    gc()
    
    # do clip interactions tend to have large association values?   
    # FDR thresholds
    thresholds = c(1e-30, 1e-25, 1e-20, 1e-15, 1e-10, 1)
    eval_clip_threshs_fdr = lapply(
        thresholds, function(thresh){
            # spearman
            eval_spear = X %>% 
                filter(spear_padj <= thresh) %>%
                compute_prec_rec_mult("in_clip", "abs_spear_coef") %>%
                mutate(thresh_fdr = thresh)
            
            # linear model
            eval_lm = X %>% 
                filter(lm_padj <= thresh) %>%
                compute_prec_rec_mult("in_clip", c("abs_lm_coef","lm_pearson")) %>%
                mutate(thresh_fdr = thresh)
            
            # linear model 2
            eval_lm2 = X %>% 
                filter(lm2_padj <= thresh) %>%
                compute_prec_rec_mult("in_clip", c("abs_lm2_coef","lm2_pearson")) %>%
                mutate(thresh_fdr = thresh)

            
            evals = rbind(eval_spear, eval_lm, eval_lm2) %>%
                mutate(n_obs_lab = sprintf("n(%s)=%s",ranking_var,n_obs))
            return(evals)
    }) %>%
    do.call(rbind, .)
    
    plts[["eval_clip-recall_vs_precision-fdr_threshs-scatter"]] = eval_clip_threshs_fdr %>%
        #filter(ranking_var %in% c("lm_pearson","lm2_pearson")) %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("Paired") +
        theme_pubr() +
        facet_wrap(~thresh_fdr) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        ylim(0,NA) +
        labs(x="Recall", y="Precision", color="Association")
    
    
    plts[["eval_clip-recall_vs_precision-fdr_threshs-bar"]] = eval_clip_threshs_fdr %>%
        distinct(thresh_fdr, ranking_var, n_obs) %>%
        ggbarplot(x="thresh_fdr", y="n_obs", fill="ranking_var", 
                  color=NA, palette="Paired", position=position_dodge(0.9)) +
        yscale("log10", .format=TRUE) +
        labs(x="FDR Threshold", y="Total SF-exon interactions", fill="Association")
        
    ## correlation thresholds
    thresholds = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
    eval_clip_threshs_corr = lapply(
        thresholds, function(thresh){
            # spearman
            eval_spear = X %>%
                filter(abs(abs_spear_coef) >= thresh) %>%
                compute_prec_rec_mult("in_clip", "spear_padj") %>%
                mutate(thresh_corr = thresh)
            
            # linear model
            eval_lm = X %>% 
                filter(lm_pearson >= thresh) %>%
                compute_prec_rec_mult("in_clip", c("abs_lm_coef","lm_padj")) %>%
                mutate(thresh_corr = thresh)
            
            # linear model 2
            eval_lm2 = X %>% 
                filter(lm2_pearson >= thresh) %>%
                compute_prec_rec_mult("in_clip", c("abs_lm2_coef","lm2_padj")) %>%
                mutate(thresh_corr = thresh)
            
            evals = rbind(eval_spear, eval_lm, eval_lm2)
            return(evals)
    }) %>%
    do.call(rbind, .)
    
    plts[["eval_clip-recall_vs_precision-corr_threshs-scatter"]] = eval_clip_threshs_corr %>%
        #filter(ranking_var %in% c("lm_pearson","lm2_pearson")) %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("Paired") +
        theme_pubr() +
        facet_wrap(~thresh_corr, ncol=3) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        ylim(0, NA) +
        labs(x="Recall", y="Precision", color="Association")
   
   plts[["eval_clip-recall_vs_precision-corr_threshs-bar"]] = eval_clip_threshs_corr %>%
        distinct(thresh_corr, ranking_var, n_obs) %>%
        ggbarplot(x="thresh_corr", y="n_obs", fill="ranking_var", 
                  color=NA, palette="Paired", position=position_dodge(0.9)) +
        yscale("log10", .format=TRUE) +
        labs(x="Correlation Threshold", y="Total SF-exon interactions", fill="Association")
    
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

    vars_oi = c("mutual_information",
                "abs_lm_coef","abs_lm_pearson","lm_pearson","lm_pvalue",
                "abs_lm2_coef","abs_lm2_pearson","lm2_pearson","lm2_pvalue",
                "abs_spear_coef","spear_pvalue")
    
    eval_pert = X %>%
        group_by(cell_line) %>%
        compute_prec_rec_mult("is_inter", vars_oi) %>%
        ungroup()
    gc()
    
    # do clip interactions tend to have large association values?
    # evaluate precision and recall using additional thresholds
    ## FDR thresholds
    thresholds = c(1e-30, 1e-25, 1e-20, 1e-15, 1e-10, 1)
    eval_pert_threshs_fdr = lapply(
        thresholds, function(thresh){
            # spearman
            eval_spear = X %>% 
                filter(spear_padj <= thresh) %>%
                group_by(cell_line) %>%
                compute_prec_rec_mult("is_inter", "abs_spear_coef") %>%
                ungroup() %>%
                mutate(thresh_fdr = thresh)
            
            # linear model
            eval_lm = X %>% 
                filter(lm_padj <= thresh) %>%
                group_by(cell_line) %>%
                compute_prec_rec_mult("is_inter", c("abs_lm_coef","lm_pearson")) %>%
                ungroup() %>%
                mutate(thresh_fdr = thresh)
            
            # linear model 2
            eval_lm2 = X %>% 
                filter(lm2_padj <= thresh) %>%
                group_by(cell_line) %>%
                compute_prec_rec_mult("is_inter", c("abs_lm2_coef","lm2_pearson")) %>%
                ungroup() %>%
                mutate(thresh_fdr = thresh)
            
            evals = rbind(eval_spear, eval_lm, eval_lm2)
            return(evals)
    }) %>%
    do.call(rbind, .)
    
    plts[["eval_pert-recall_vs_precision-fdr_threshs-scatter"]] = eval_pert_threshs_fdr %>%
        #filter(ranking_var %in% c("lm_pearson","lm2_pearson")) %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("Paired") +
        theme_pubr() +
        facet_wrap(~thresh_fdr+cell_line, ncol=3) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        ylim(0, NA) +
        labs(x="Recall", y="Precision", color="Association")
   
   plts[["eval_pert-recall_vs_precision-fdr_threshs-bar"]] = eval_pert_threshs_fdr %>%
        distinct(cell_line, thresh_fdr, ranking_var, n_obs) %>%
        ggbarplot(x="thresh_fdr", y="n_obs", fill="ranking_var", 
                  color=NA, palette="Paired", position=position_dodge(0.9)) +
        yscale("log10", .format=TRUE) +
        facet_wrap(~cell_line) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="FDR Threshold", y="Total SF-exon interactions", fill="Association")
    
    ## correlation thresholds
    thresholds = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
    eval_pert_threshs_corr = lapply(
        thresholds, function(thresh){
            # spearman
            eval_spear = X %>%
                filter(abs(abs_spear_coef) >= thresh) %>%
                group_by(cell_line) %>%
                compute_prec_rec_mult("is_inter", "spear_padj") %>%
                ungroup() %>%
                mutate(thresh_corr = thresh)
            
            # linear model
            eval_lm = X %>% 
                filter(lm_pearson >= thresh) %>%
                group_by(cell_line) %>%
                compute_prec_rec_mult("is_inter", c("abs_lm_coef","lm_padj")) %>%
                ungroup() %>%
                mutate(thresh_corr = thresh)
            
            # linear model 2
            eval_lm2 = X %>% 
                filter(lm2_pearson >= thresh) %>%
                group_by(cell_line) %>%
                compute_prec_rec_mult("is_inter", c("abs_lm2_coef","lm2_padj")) %>%
                ungroup() %>%
                mutate(thresh_corr = thresh)
            
            evals = rbind(eval_spear, eval_lm, eval_lm2)
            return(evals)
    }) %>%
    do.call(rbind, .)
    
    plts[["eval_pert-recall_vs_precision-corr_threshs-scatter"]] = eval_pert_threshs_corr %>%
        #filter(ranking_var %in% c("lm_pearson","lm2_pearson")) %>%
        ggplot(aes(x=recall, y=precision)) +
        geom_line(aes(color=ranking_var), size=LINE_SIZE) +
        geom_point(aes(color=ranking_var), size=1) +
        color_palette("Paired") +
        theme_pubr() +
        facet_wrap(~thresh_corr+cell_line, ncol=3) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        ylim(0, NA) +
        labs(x="Recall", y="Precision", color="Association")
   
   plts[["eval_pert-recall_vs_precision-corr_threshs-bar"]] = eval_pert_threshs_corr %>%
        distinct(cell_line, thresh_corr, ranking_var, n_obs) %>%
        ggbarplot(x="thresh_corr", y="n_obs", fill="ranking_var", 
                  color=NA, palette="Paired", position=position_dodge(0.9)) +
        yscale("log10", .format=TRUE) +
        facet_wrap(~cell_line) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Correlation Threshold", y="Total SF-exon interactions", fill="Association")
    
    return(plts)
}


# plot_clip_vs_pert = function(regulons_pert, regulons_clip){
#     plts = list()
    
#     X = regulons_pert %>% 
#         filter(cell_line != "merged") %>%
#         left_join(
#             regulons_clip %>%
#             mutate(in_clip = TRUE),
#             by=c("regulator"="ENSEMBL","target")
#         ) %>%
#         mutate(in_clip=replace_na(in_clip, FALSE)) %>%
#         drop_na(cell_line) %>%
#         distinct(cell_line, regulator, target, delta_psi, delta_psi_rel, in_clip)
    
#     # predictive power of delta PSI and delta PSI rel.
#     evaluation = X %>%
#         group_by(cell_line) %>%
#         arrange(-abs(delta_psi)) %>%
#         reframe(
#             precision = compute_precision(in_clip, -abs(delta_psi), 25),
#             recall = compute_recall(in_clip, -abs(delta_psi), 25),
#             ranking_var = "delta_psi"
#         ) %>% 
#         ungroup() %>%
#         bind_rows(
#             X %>%
#             group_by(cell_line) %>%
#             arrange(-abs(delta_psi_rel)) %>%
#             reframe(
#                 precision = compute_precision(in_clip, -abs(delta_psi_rel), 25),
#                 recall = compute_recall(in_clip, -abs(delta_psi_rel), 25),
#                 ranking_var = "delta_psi_rel"
#             ) %>%
#             ungroup()
#         ) %>% drop_na()
    
#     # do clip interactions tend to have large association values?
#     plts[["eval_clip_vs_pert-recall_vs_precision-scatter"]] = evaluation %>%
#         ggplot(aes(x=recall, y=precision)) +
#         geom_line(aes(color=ranking_var), size=LINE_SIZE) +
#         geom_point(aes(color=ranking_var), size=1) +
#         color_palette("simpsons") +
#         theme_pubr() +
#         facet_wrap(~cell_line) +
#         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="Recall", y="Precision", color="Association")

#     return(plts)
# }


make_plots = function(assocs, regulons_clip, regulons_pert){
    plts = list(
        plot_eval_clip(assocs, regulons_clip),
        plot_eval_pert(assocs, regulons_pert)#,
        #plot_clip_vs_pert(regulons_pert, regulons_clip)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(assocs, regulons_clip, regulons_pert){
    figdata = list(
        "assocs_evaluation" = list(
            "associations" = assocs,
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
    save_plt(plts, "eval_clip-recall_vs_precision-fdr_threshs-scatter", '.pdf', figs_dir, width=15, height=12)
    save_plt(plts, "eval_clip-recall_vs_precision-fdr_threshs-bar", '.pdf', figs_dir, width=5, height=7)
    save_plt(plts, "eval_clip-recall_vs_precision-corr_threshs-scatter", '.pdf', figs_dir, width=15, height=12)
    save_plt(plts, "eval_clip-recall_vs_precision-corr_threshs-bar", '.pdf', figs_dir, width=5, height=7)
    
    save_plt(plts, "eval_pert-recall_vs_precision-fdr_threshs-scatter", '.pdf', figs_dir, width=15, height=35)
    save_plt(plts, "eval_pert-recall_vs_precision-fdr_threshs-bar", '.pdf', figs_dir, width=15, height=8)
    save_plt(plts, "eval_pert-recall_vs_precision-corr_threshs-scatter", '.pdf', figs_dir, width=15, height=35)
    save_plt(plts, "eval_pert-recall_vs_precision-corr_threshs-bar", '.pdf', figs_dir, width=15, height=8)
    
    #save_plt(plts, "eval_clip_vs_pert-recall_vs_precision-scatter", '.pdf', figs_dir, width=10, height=6)
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
        make_option("--assocs_lm2_file", type="character"),
        make_option("--regulons_clip_file", type="character"),
        make_option("--regulons_pert_dpsi_file", type="character"),
        make_option("--regulons_pert_dpsi_rel_file", type="character"),
        make_option("--regulators_file", type="character"),
        make_option("--targets_file", type="character"),
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
    assocs_lm2_file = args[["assocs_lm2_file"]]
    regulons_clip_file = args[["regulons_clip_file"]]
    regulons_pert_dpsi_file = args[["regulons_pert_dpsi_file"]]
    regulons_pert_dpsi_rel_file = args[["regulons_pert_dpsi_rel_file"]]
    figs_dir = args[["figs_dir"]]
    
    set.seed(RANDOM_SEED)
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    assocs_mi = read_tsv(assocs_mi_file)
    assocs_spear = read_tsv(assocs_spear_file)
    assocs_lm = read_tsv(assocs_lm_file)
    assocs_lm2 = read_tsv(assocs_lm2_file)
    regulons_clip = read_tsv(regulons_clip_file)
    regulons_pert_dpsi = read_tsv(regulons_pert_dpsi_file)
    regulons_pert_dpsi_rel = read_tsv(regulons_pert_dpsi_rel_file)
    
    gc()
    
    # prep
    ## merge sf-exon associations
    assocs = assocs_mi %>%
        rename(mutual_information = association) %>%
        # linear model
        left_join(
            assocs_lm %>%
            rename(
                lm_coef = target_coefficient_mean,
                lm_pvalue = lr_pvalue,
                lm_padj = lr_padj,
                lm_pearson = pearson_correlation_mean
            ) %>%
            dplyr::select(regulator, target, lm_coef, lm_pvalue, lm_pearson, lm_padj),
            by = c("regulator","target")
        ) %>%
        # spearman coefficient
        left_join(
            assocs_spear %>% 
            rename(
                spear_coef = statistic,
                spear_pvalue = pvalue,
                spear_padj = padj
            ) %>%
            dplyr::select(regulator, target, spear_coef, spear_pvalue, spear_padj),
            by = c("regulator","target")
        ) %>%
        # linear model2
        left_join(
            assocs_lm2 %>%
            rename(
                lm2_coef = target_splicing_coefficient_mean,
                lm2_pvalue = lr_pvalue,
                lm2_padj = lr_padj,
                lm2_pearson = pearson_correlation_mean
            ) %>%
            dplyr::select(regulator, target, lm2_coef, lm2_pvalue, lm2_pearson, lm2_padj),
            by = c("regulator","target")
        ) 
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
        mutate(is_inter = abs(delta_psi) > 10) # define interactions based on dpsi threshold

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
            log_lm_pvalue = -log10(lm_pvalue),
            log_lm2_pvalue = -log10(lm2_pvalue),
            log_spear_pvalue = -log10(spear_pvalue)
        ) %>%
        mutate(
            abs_mutual_information = -abs(mutual_information),
            abs_lm_coef = -abs(lm_coef),
            abs_lm2_coef = -abs(lm2_coef),
            abs_lm_pearson = -abs(lm_pearson),
            abs_lm2_pearson = -abs(lm2_pearson),
            abs_spear_coef = -abs(spear_coef)
        )
    gc()
    
    # plot
    plts = make_plots(assocs, regulons_clip, regulons_pert)
    gc()
    
    # make figdata
    figdata = make_figdata(assocs, regulons_clip, regulons_pert)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}