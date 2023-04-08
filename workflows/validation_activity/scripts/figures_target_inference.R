#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# make figures of evaluate target inference algorithms.
# median accuracy of 0.973 and 0.975 in HepG2 and K562 with threshold of 0.2

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "darkgreen"
PAL_DUAL = c("darkgreen","#7570B3") # '#1B9E77''#7570B3'

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","sf_targets_inference")
# SUPPORT_DIR = file.path(ROOT,"support")

#evaluation_file = file.path(RESULTS_DIR,"files","inference_evaluation","merged.tsv.gz")

# viper_result_file = file.path(RESULTS_DIR,"files","validations","aracne","event_psi","ENCORE","HepG2-LIHC_regulons.tsv.gz")
# regulons_file = file.path(RESULTS_DIR,"files","target_inference","aracne","event_psi","LIHC","regulons.tsv.gz")
# encore_logfc_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"HepG2",'log2_fold_change_tpm.tsv.gz')
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# sf_info_file = file.path(SUPPORT_DIR,"Rodolska2022-suptab1-splicing_factors.tsv")

# figs_dir = file.path(RESULTS_DIR,'figures','validations')

##### FUNCTIONS #####
prep_encore_logfc = function(encore_logfc, thresh_sign=0.1, thresh_log2fc=1){
    
    encore_pvalues = encore_logfc %>%
        column_to_rownames("ID") %>%
        mutate_all(function(x){
            z_score = (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
            p_value = 2*pnorm(-abs(z_score))
            fdr = p.adjust(p_value, method="fdr")
            return(fdr)
        })
        
    is_significant = sapply(colnames(encore_pvalues), function(KD){
        logfc = encore_logfc %>%
            filter(ID==KD) %>%
            dplyr::select(all_of(KD)) %>%
            pull() %>%
            abs()
        pval = encore_pvalues[KD,KD]
        is_significant = (pval < thresh_sign) & (logfc > thresh_log2fc)
        if (is.na(is_significant)){
            is_significant = FALSE
        }
        return(is_significant)
    })
    
    encore_logfc = encore_logfc[,c(TRUE,is_significant)]
    
    return(encore_logfc)
}


plot_regulons = function(regulons){
    plts = list()
    
    X = regulons
    
    total_upreg = X[["regulator"]] %>% unique() %>% length()
    total_targets = X[["target"]] %>% unique() %>% length()
    
    # number of targets per upstream regulator?
    plts[["regulons-n_targets_per_reg-distr"]] = X %>%
        count(regulator) %>%
        gghistogram(x="n", fill=PAL_DARK, color=NA) +
        labs(x="# Targets per Upstream Regulator", y="Count")
    
    # pleiotropy? how many upstream regulators per target?
    plts[["regulons-pleiotropy-distr"]] = X %>%
        count(target) %>%
        gghistogram(x="n", fill=PAL_DARK, color=NA) +
        labs(x="# Upstream Regulators per Target", y="Count")
    
    # target length of SRRM3 and SRRM4?
    plts[["regulons-microexons-violin"]] = X %>%
        mutate(
            sf_type = ifelse(regulator %in% c("ENSG00000177679","ENSG00000139767"), 
                            "Short", "Long")
        ) %>% # SRRM3, SRRM4
        distinct(sf_type, target, target_event_length) %>%
        ggviolin(x="sf_type", y="target_event_length", fill="sf_type", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, fill=NA, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text(
            aes(label=label, y=1),
            . %>% count(sf_type) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY,
        ) +
        guides(fill="none") +
        yscale("log10", .format=TRUE) +
        labs(x="SF Specificity", y="Target Event Length")
    
    return(plts)
}


plot_viper_activities = function(viper_result, encore_logfc, sf_info){
    plts = list()
    
    X = viper_result %>%
        pivot_longer(
            -regulator, 
            names_to="KD", 
            values_to="viper_activity"
        ) %>%
        left_join(
            encore_logfc %>%
            pivot_longer(-ID, names_to="KD", values_to="log2_fc"),
            by = c("regulator"="ID","KD")
        ) %>%
        left_join(
            encore_logfc %>%
            mutate_if(is.numeric, rank) %>%
            pivot_longer(-ID, names_to="KD", values_to="log2_fc_ranking"),
            by = c("regulator"="ID","KD")
        ) %>%
        drop_na() %>%
        group_by(KD) %>%
        mutate(
            viper_activity_ranking = rank(viper_activity)
        ) %>%
        ungroup() %>%
        left_join(
            sf_info %>% 
            distinct(ENSEMBL, CLASS, SYMBOL),
            by = c("KD"="ENSEMBL")
        ) %>%
        drop_na(CLASS)
        
    
    # randomize selection of the same N
    x_validated = X %>% filter(is_validated)
    x_random_non_validated = X %>%
        filter(!is_validated) %>%
        sample_n(size=nrow(x_validated), seed=1234)
    x = rbind(x_validated, x_random_non_validated)
    
    # are the rankings of validated SFs higher than non-validated?
    plts[["viper_activities-activity_vs_is_validated-violin"]] = x %>%
        ggviolin(x="is_validated", y="viper_activity", 
                 fill="is_validated", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black") +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        fill_palette(PAL_DUAL) +
        geom_text(
            aes(label=label, y=-2),
            . %>% count(is_validated) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY,
        ) +
        guides(fill="none") +
        labs(x="Is Validated Experimentally", y="VIPER Activity")
    
    plts[["viper_activities-activity_ranking_vs_is_validated-violin"]] = x %>%
        ggviolin(x="is_validated", y="viper_activity_ranking", 
                 fill="is_validated", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black") +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        fill_palette(PAL_DUAL) +
        geom_text(
            aes(label=label, y=-5),
            . %>% count(is_validated) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="Is Validated Experimentally", y="VIPER Activity Ranking")
    
    # does predicted activity correlate with experimental logFC?
    plts[["viper_activities-fc_kds_vs_activity-scatter"]] = X %>%
        filter(is_validated) %>%
        ggplot(aes(x=log2_fc, y=viper_activity)) +
        geom_scattermore(
            aes(color=CLASS),
            pixels = c(1000,1000), pointsize=10, alpha=0.8) +
        color_palette(get_palette("Paired",20)) +
        theme_pubr() +
        facet_wrap(~is_validated, scales="free") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        labs(x="Log2FC KD vs WT", y="Viper Activity")
    
    plts[["viper_activities-fc_kds_vs_activity_ranking-scatter"]] = X %>%
        ggplot(aes(x=log2_fc, y=viper_activity_ranking)) +
        geom_scattermore(
            aes(color=is_validated),
            pixels = c(1000,1000), pointsize=10, alpha=0.8) +
        color_palette(PAL_DUAL) +
        theme_pubr() +
        facet_wrap(~is_validated, scales="free_x") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        guides(color="none") +
        labs(x="Log2FC KD vs WT", y="Viper Activity Ranking")

    return(plts)
}


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
        make_option("--encore_logfc_file", type="character"),
        make_option("--regulons_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    viper_result_file = args[["viper_result_file"]]
    encore_logfc_file = args[["encore_logfc_file"]]
    regulons_file = args[["regulons_file"]]
    event_info_file = args[["event_info_file"]]
    sf_info_file = args[["sf_info_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    viper_result = read_tsv(viper_result_file)
    encore_logfc = read_tsv(encore_logfc_file)
    regulons = read_tsv(regulons_file)
    event_info = read_tsv(event_info_file)
    sf_info = read_tsv(sf_info_file)
    
    # prep
    ## select significant KDs
    encore_logfc = prep_encore_logfc(encore_logfc)
    ## add event length info to targets
    regulons = regulons %>%
        left_join(
            event_info %>% 
                dplyr::select(EVENT, LE_o) %>%
                rename(target_event_length = LE_o),
            by = c("target"="EVENT")
        )
    ## filter regulons with few targets
    regulators_oi = regulons %>%
        count(regulator) %>%
        filter(n > 20) %>%
        pull(regulator)
    
    regulons = regulons %>%
        filter(regulator %in% regulators_oi)
    
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