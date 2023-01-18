#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA datasets used as ground truth.
# 
# Outline
# -------
# - using relative delta PSI rather than raw delta PSI

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(extrafont)

# variables

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#6AC2BF"
PAL_SINGLE_DARK = "#716454"
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","eda")

# dpsi_k562_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE','K562','delta_psi-EX.tsv.gz')
# dpsi_hepg2_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE','HepG2','delta_psi-EX.tsv.gz')
# rel_k562_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE','K562','delta_psi_rel-EX.tsv.gz')
# rel_hepg2_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE','HepG2','delta_psi_rel-EX.tsv.gz')

# figs_dir = file.path(RESULTS_DIR,'figures','ground_truth')

##### FUNCTIONS #####
compute_correlations = function(mat_a, mat_b){
    mat_a = mat_a %>% column_to_rownames("EVENT") %>% as.data.frame()
    mat_b = mat_b %>% column_to_rownames("EVENT") %>% as.data.frame()
        
    columns = intersect(colnames(mat_a), colnames(mat_b))
    rows = intersect(rownames(mat_a), rownames(mat_b))
    
    correls = lapply(columns, function(col_oi){
        a = mat_a[rows,col_oi]
        b = mat_b[rows,col_oi]
        
        idx = !is.na(a) & !is.na(b)
        
        corr = data.frame(
            KD_ENSEMBL = col_oi,
            correlation_pearson = cor(a[idx], b[idx], method="pearson"),
            correlation_spearman = cor(a[idx], b[idx], method="spearman"),
            n_obs = sum(idx)
        )
        
        return(corr)
    })
    correls = do.call(rbind, correls)
    
    return(correls)
}


plot_reproducibility = function(correls, dpsi_hepg2, dpsi_k562, rel_hepg2, rel_k562){
    plts = list()
    
    X = correls %>%
        pivot_longer(
            c(correlation_pearson, correlation_spearman), 
            names_to="correlation_method", values_to="correlation"
        )
    
    # does reproducibility improve with relative delta PSI?
    plts[["reproducibility-correls-distr"]] = X %>%
        ggviolin(x="data_type", y="correlation", fill=PAL_SINGLE_ACCENT, color=NA) +
        geom_boxplot(width=0.1, fill=NA, outlier.size=0.1) +
        facet_wrap(~correlation_method, scales="free_y") +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Type of Delta PSI", y="Reproducibility K562 and HepG2")
    
    # plot best/worst examples
    x = X %>%
        filter(correlation_method=="correlation_pearson") %>%
        pivot_wider(names_from="data_type", values_from="correlation") %>%
        mutate(corr_diff = dpsi - rel_dpsi)
    
    kds_oi = rbind(
        x %>% group_by(correlation_method) %>% slice_min(corr_diff, n=1),
        x %>% group_by(correlation_method) %>% slice_max(corr_diff, n=1)
    )
    
    ## best
    kd_oi = kds_oi[["KD_ENSEMBL"]][1]
    df = data.frame(
        dpsi_hepg2 = dpsi_hepg2 %>% pull(kd_oi),
        dpsi_k562 = dpsi_k562 %>% pull(kd_oi),
        rel_hepg2 = rel_hepg2 %>% pull(kd_oi),
        rel_k562 = rel_k562 %>% pull(kd_oi)
    )
    plts[["reproducibility-example-best-dpsi-scatter"]] = df %>%
        ggplot(aes(x=dpsi_hepg2, y=dpsi_k562)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y.npc=0.9) +
        labs(x="deltaPSI HepG2", y="deltaPSI K562") +
        theme_pubr()
    
    plts[["reproducibility-example-best-rel_dpsi-scatter"]] = df %>%
        ggplot(aes(x=rel_hepg2, y=rel_k562)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y.npc=0.9) +
        labs(x="Rel. deltaPSI HepG2", y="Rel. deltaPSI K562") +
        theme_pubr()
    
    ## worst
    kd_oi = kds_oi[["KD_ENSEMBL"]][2]
    df = data.frame(
        dpsi_hepg2 = dpsi_hepg2 %>% pull(kd_oi),
        dpsi_k562 = dpsi_k562 %>% pull(kd_oi),
        rel_hepg2 = rel_hepg2 %>% pull(kd_oi),
        rel_k562 = rel_k562 %>% pull(kd_oi)
    )
    plts[["reproducibility-example-worst-dpsi-scatter"]] = df %>%
        ggplot(aes(x=dpsi_hepg2, y=dpsi_k562)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y.npc=0.9) +
        labs(x="deltaPSI HepG2", y="deltaPSI K562") +
        theme_pubr()
    
    plts[["reproducibility-example-worst-rel_dpsi-scatter"]] = df %>%
        ggplot(aes(x=rel_hepg2, y=rel_k562)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY, label.y.npc=0.9) +
        labs(x="Rel. deltaPSI HepG2", y="Rel. deltaPSI K562") +
        theme_pubr()
    
    return(plts)
}


make_plots = function(correls, dpsi_hepg2, dpsi_k562, rel_hepg2, rel_k562){
    plts = list(
        plot_reproducibility(correls, dpsi_hepg2, dpsi_k562, rel_hepg2, rel_k562)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(correls, dpsi_hepg2, dpsi_k562, rel_hepg2, rel_k562){
    figdata = list(
        "reproducibility" = list(
            "regular_dpsi-K562" = dpsi_k562,
            "regular_dpsi-HepG2" = dpsi_hepg2,
            "relative_dpsi-K562" = rel_k562,
            "relative_dpsi-HepG2" = rel_hepg2,
            "correlations" = correls
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
    save_plt(plts, "reproducibility-correls-distr", '.pdf', figs_dir, width=8, height=5)
    
    save_plt(plts, "reproducibility-example-best-dpsi-scatter", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "reproducibility-example-best-rel_dpsi-scatter", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "reproducibility-example-worst-dpsi-scatter", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "reproducibility-example-worst-rel_dpsi-scatter", '.pdf', figs_dir, width=5, height=5)
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
        make_option("--dpsi_k562_file", type="character"),
        make_option("--dpsi_hepg2_file", type="character"),
        make_option("--rel_k562_file", type="character"),
        make_option("--rel_hepg2_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    dpsi_k562_file = args[["dpsi_k562_file"]]
    dpsi_hepg2_file = args[["dpsi_hepg2_file"]]
    rel_k562_file = args[["rel_k562_file"]]
    rel_hepg2_file = args[["rel_hepg2_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    dpsi_k562 = read_tsv(dpsi_k562_file)
    dpsi_hepg2 = read_tsv(dpsi_hepg2_file)
    rel_k562 = read_tsv(rel_k562_file)
    rel_hepg2 = read_tsv(rel_hepg2_file)
    
    # correlate genes together
    correls = rbind(
        compute_correlations(dpsi_k562[,1:25], dpsi_hepg2[,1:25]) %>% mutate(data_type="dpsi"),
        compute_correlations(rel_k562[,1:25], rel_hepg2[,1:25]) %>% mutate(data_type="rel_dpsi")
    )
    
    # plot
    plts = make_plots(correls, dpsi_hepg2, dpsi_k562, rel_hepg2, rel_k562)
    
    # make figdata
    figdata = make_figdata(correls, dpsi_hepg2, dpsi_k562, rel_hepg2, rel_k562)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}