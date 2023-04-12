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
# RESULTS_DIR = file.path(ROOT,"results","validation_activity")
# SUPPORT_DIR = file.path(ROOT,"support")

# protein_activities_file = file.path(RESULTS_DIR,'files','ground_truth','viper',"merged-genexpr_vs_psi_imputed-delta_psi.tsv.gz")
# encore_logfc_hepg2_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"HepG2",'log2_fold_change_tpm.tsv.gz')
# encore_logfc_k562_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"K562",'log2_fold_change_tpm.tsv.gz')

# figs_dir = file.path(RESULTS_DIR,'figures','validations')

##### FUNCTIONS #####
prep_encore_logfc = function(encore_logfc, thresh_sign=0.05, thresh_log2fc=1){
    
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


plot_viper_activities = function(protein_activities, encore_logfc){
    plts = list()
    
    X = protein_activities %>%
        left_join(encore_logfc, by = c("regulator"="ID","KD","cell_line")
        ) %>%
        drop_na() %>%
        mutate(
            is_validated = regulator == KD
        )
    
    # does predicted activity correlate with experimental logFC of knocked down genes?
    plts[["viper_activities-fc_kds_vs_protein_activity-scatter"]] = X %>%
        filter(is_validated) %>%
        ggplot(aes(x=log2_fc, y=protein_activity)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=10, alpha=0.8) +
        color_palette(get_palette("Paired",20)) +
        theme_pubr() +
        facet_wrap(~assoc_method+regulon+cell_line, scales="free", ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        labs(x="Log2FC KD vs WT", y="Viper Activity")
    
    # for those genes with largest predicted changes in activity?
    plts[["viper_activities-fc_kds_vs_max_protein_activity-scatter"]] = X %>%
        group_by(assoc_method, regulon, cell_line, KD) %>%
        slice_max(abs(protein_activity), n=1) %>%
        ungroup() %>%
        ggplot(aes(x=log2_fc, y=protein_activity)) +
        geom_scattermore(aes(color=is_validated), pixels = c(1000,1000), pointsize=10, alpha=1) +
        color_palette("lancet") + 
        theme_pubr() +
        facet_wrap(~assoc_method+regulon+cell_line, scales="free", ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        labs(x="Log2FC KD vs WT", y="Viper Activity", color="Is KD")
    
    # what is the ranking of each splicing factor when validated across all experiments?
    set.seed(1234)
    x = X %>%
        group_by(assoc_method, regulon, cell_line, regulator) %>%
        arrange(abs(protein_activity)) %>% 
        mutate(ranking = row_number()) %>%
        ungroup() %>%
        filter(is_validated) %>%
        mutate(ordering_type = "real") %>%
        bind_rows(
            X %>%
            group_by(assoc_method, regulon, cell_line) %>%
            mutate(is_validated = sample(is_validated)) %>%
            ungroup() %>%
            group_by(assoc_method, regulon, cell_line, regulator) %>%
            arrange(abs(protein_activity)) %>% 
            mutate(ranking = row_number()) %>%
            ungroup() %>%
            filter(is_validated) %>%
            mutate(ordering_type = "random")
        )
    
    plts[["viper_activities-specificity-violin"]] = x %>%
        ggviolin(x="ordering_type", y="ranking", fill="ordering_type", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, fill=NA, outlier.size=0.1) + 
        facet_wrap(~assoc_method+regulon+cell_line, scales="free", ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(fill="none") + 
        labs(x="Sample Ranking", y="Viper Activity Ranking Across Experiments")

    
    return(plts)
}


make_plots = function(protein_activities, encore_logfc){
    plts = list(
        plot_viper_activities(protein_activities, encore_logfc)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activities, encore_logfc){
    
    X = protein_activities %>%
        left_join(encore_logfc, by = c("regulator"="ID","KD","cell_line")
        ) %>%
        drop_na() %>%
        mutate(
            is_validated = regulator == KD
        )
    
    figdata = list(
        "target_inference" = list(
            "protein_activities_vs_kd_fc" = X
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
    save_plt(plts, "viper_activities-fc_kds_vs_protein_activity-scatter", '.pdf', figs_dir, width=12, height=35)
    save_plt(plts, "viper_activities-fc_kds_vs_max_protein_activity-scatter", '.pdf', figs_dir, width=12, height=35)
    save_plt(plts, "viper_activities-specificity-violin", '.pdf', figs_dir, width=12, height=35)
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
        make_option("--protein_activities_file", type="character"),
        make_option("--encore_logfc_k562_file", type="character"),
        make_option("--encore_logfc_hepg2_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    protein_activities_file = args[["protein_activities_file"]]
    encore_logfc_hepg2_file = args[["encore_logfc_hepg2_file"]]
    encore_logfc_k562_file = args[["encore_logfc_k562_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activities = read_tsv(protein_activities_file)
    encore_logfc_hepg2 = read_tsv(encore_logfc_hepg2_file)
    encore_logfc_k562 = read_tsv(encore_logfc_k562_file)
    
    # prep
    ## select significant KDs
    encore_logfc_hepg2 = prep_encore_logfc(encore_logfc_hepg2)
    encore_logfc_k562 = prep_encore_logfc(encore_logfc_k562)
    ## combine
    encore_logfc = encore_logfc_hepg2 %>%
        pivot_longer(-ID, names_to="KD", values_to="log2_fc") %>%
        mutate(cell_line = "HepG2") %>%
        bind_rows(
            encore_logfc_k562 %>%
            pivot_longer(-ID, names_to="KD", values_to="log2_fc") %>%
            mutate(cell_line = "K562")
        )
    
    # plot
    plts = make_plots(protein_activities, encore_logfc)
    
    # make figdata
    figdata = make_figdata(protein_activities, encore_logfc)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}