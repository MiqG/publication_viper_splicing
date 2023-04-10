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

# figs_dir = file.path(RESULTS_DIR,'figures','assocs_evaluation','LIHC')

##### FUNCTIONS #####
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


make_plots = function(assocs, regulons_clip, regulons_pert){
    plts = list(
        plot_eval_clip(assocs, regulons_clip),
        plot_eval_pert(assocs, regulons_pert),
        plot_clip_vs_pert(regulons_pert, regulons_clip)
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
            log_lm_coef = sign(lm_coef)*log10(abs(lm_coef)+1),
            log_lm_pvalue = -log10(lm_pvalue),
            log_spear_pvalue = -log10(spear_pvalue)
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