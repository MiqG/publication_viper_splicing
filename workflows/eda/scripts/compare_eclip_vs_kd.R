#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - consider on exon or cryptic 

require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(optparse)

ROOT = here::here()

# variables

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# eclip_file = file.path(PREP_DIR,"eclip_peaks_mapped","merged.tsv.gz")
# delta_psi_file = file.path(PREP_DIR,"kd_transcriptomes","ENCORE","delta_psi-EX.tsv.gz")
# delta_psi_file = file.path(PREP_DIR,"kd_transcriptomes","ENCORE","delta_psi_rel-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,'metadata','ENCORE.tsv.gz')
# figs_dir = file.path(ROOT,'results','eda','figures','eda')

##### FUNCTIONS #####
plot_comparison = function(eclip, dpsi){
    plts = list()
    
    # do strong splicing changes occur where the RBP binds?
    plts[["comparison-fc_vs_pvalue"]] = eclip %>% 
        group_by(cell_line, EVENT) %>%
        slice_max(log10_pvalue, n=1) %>%
        ungroup() %>%
        ggplot(aes(x=log10_pvalue, y=log2_fc)) +
        geom_scattermore(pointsize=8, alpha=0.5) +
        facet_wrap(~cell_line) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="log2(FC)", y="log10(p-value)")
    
    # distributions of delta PSI vs binds RBP or not
    dpsi %>%
        ggboxplot(x="binds_rbp", y="abs_deltaPSI") +
        stat_compare_means(method="wilcox.test") +
        facet_wrap(~KD_GENE+cell_line)
    
    dpsi %>% 
        mutate(deltaPSI_binned = cut(abs_deltaPSI, breaks=seq(0,100,5))) %>%
        drop_na() %>%
        count(KD_GENE, cell_line, deltaPSI_binned, binds_rbp) %>%
        group_by(KD_GENE, cell_line, binds_rbp) %>%
        mutate(perc = n / sum(n)) %>%
        ungroup() %>%
        group_by(KD_GENE, cell_line, deltaPSI_binned) %>%
        mutate(rel_perc = perc / sum(perc)) %>%
        ungroup() %>%        
        ggbarplot(x="deltaPSI_binned", y="rel_perc", fill="binds_rbp") +
        facet_wrap(~KD_GENE+cell_line) +
        theme_pubr(x.text.angle = 70)
    
    # KD_GENE PSI changes vs eCLIP peaks
    X = eclip %>%
        left_join(dpsi, by=c("EVENT","cell_line","rbp"="KD_GENE")) %>%
        drop_na(deltaPSI)
    
    # does change in PSI correlate with eCLIP features?
    plts[["comparison-dpsi_vs_fc"]] = X %>% 
        ggplot(aes(x=deltaPSI, y=log2_fc)) +
        geom_scattermore(pointsize=8, alpha=0.5) +
        facet_wrap(~cell_line+bin_abs_distance_event_to_peak) +
        stat_cor(method="spearman") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Delta PSI", y="log2(FC)")

    plts[["comparison-dpsi_vs_pvalue"]] = X %>%
        ggplot(aes(x=deltaPSI, y=log10_pvalue)) +
        geom_scattermore(pointsize=8, alpha=0.5) +
        facet_wrap(~cell_line) +
        stat_cor(method="spearman") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Delta PSI", y="log10(p-value)")
    
    plts[["comparison-abs_dpsi_vs_fc"]] = X %>%
        ggplot(aes(x=abs(deltaPSI), y=log2_fc)) +
        geom_scattermore(pointsize=8, alpha=0.5) +
        facet_wrap(~cell_line) +
        stat_cor(method="spearman") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="|Delta PSI|", y="log2(FC)")

    plts[["comparison-abs_dpsi_vs_pvalue"]] = X %>%
        ggplot(aes(x=abs(deltaPSI), y=log10_pvalue)) +
        geom_scattermore(pointsize=8, alpha=0.5) +
        facet_wrap(~cell_line) +
        stat_cor(method="spearman") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="|Delta PSI|", y="log10(p-value)")

    
    # should we consider the distances to event splice sites?
    X %>%
        count(peak_position) %>%
        ggbarplot(x="peak_position", y="n")
    
    plts[["comparison-distances"]] = X %>%
        gghistogram(x="distance_event_to_peak") +
        facet_wrap(~peak_position+cell_line, scales="free")
    
    # correlation between delta PSI and FC or p-value
    x = X %>%
        filter(abs(deltaPSI)>5) %>%
        group_by(cell_line, rbp) %>%
        summarize(correlation = cor(deltaPSI, log2_fc, method="spearman")) %>%
        ungroup()
    
    plts[[""]] = x %>%
        gghistogram(x="correlation") +
        facet_wrap(~cell_line)
    

    
    return(plts)
}


make_plots = function(){
    plts = list(
        plot_model_selection(),
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(){

    figdata = list(
        "model_selection" = list(
            "model_summaries"= models,
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension=".pdf", 
                    directory="", dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # model selection
    save_plt(plts, "model_sel-deps_sorted_vs_std_ctl_neg", ".pdf", figs_dir, width=5, height=5)
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
        make_option("--eclip_file", type="character"),
        make_option("--delta_psi_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    eclip_file = args[["eclip_file"]]
    delta_psi_file = args[["delta_psi_file"]]
    metadata_file = args[["metadata_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    eclip = read_tsv(eclip_file)
    delta_psi = read_tsv(delta_psi_file)
    metadata = read_tsv(metadata_file)
    
    gc()
    
    # preprocess
    eclip = eclip %>%
        filter(str_detect(sample_id,"HepG2") | str_detect(sample_id,"K562")) %>%
        separate(sample_id, into = c("rbp","cell_line","idr"), remove=FALSE)
    
    # measure distances
    ## event coordinates vs peak
    eclip[["distance_event_start_to_peak_start"]] = eclip[["EVENT_start"]] - eclip[["peak_start"]]
    eclip[["distance_event_start_to_peak_end"]] = eclip[["EVENT_start"]] - eclip[["peak_end"]]
    eclip[["distance_event_end_to_peak_start"]] = eclip[["EVENT_end"]] - eclip[["peak_start"]]
    eclip[["distance_event_end_to_peak_end"]] = eclip[["EVENT_end"]] - eclip[["peak_end"]]
    ## get peak coordinate closest to event coordinate
    cols_oi = grep("distance_event", colnames(eclip), value=TRUE)
    eclip[["distance_event_to_peak"]] = apply(eclip[,cols_oi], 1, function(x){ x[which.min(abs(x))] })
    eclip[["abs_distance_event_to_peak"]] = abs(eclip[["distance_event_to_peak"]])
    eclip[["peak_on_intron_upstream"]] = (
        (eclip[["peak_start"]] < eclip[["EVENT_start"]]) & (eclip[["peak_end"]] < eclip[["EVENT_start"]])
    )
    eclip[["peak_on_intron_downstream"]] = (
        (eclip[["EVENT_end"]] < eclip[["peak_start"]]) & (eclip[["EVENT_end"]] < eclip[["peak_end"]])
    )
    eclip[["peak_on_intron"]] = eclip[["peak_on_intron_upstream"]] | eclip[["peak_on_intron_downstream"]]
    eclip[["peak_on_exon"]] = (
        (eclip[["EVENT_start"]] <= eclip[["peak_start"]]) & (eclip[["EVENT_end"]] >= eclip[["peak_end"]])
    )
    eclip[["peak_overlaps_3ss"]] = (
        (eclip[["EVENT_start"]] >= eclip[["peak_start"]]) & (eclip[["EVENT_end"]] >= eclip[["peak_end"]])
    )
    eclip[["peak_overlaps_5ss"]] = (
        (eclip[["EVENT_start"]] <= eclip[["peak_start"]]) & (eclip[["EVENT_end"]] <= eclip[["peak_end"]])
    )
    eclip[["peak_on_ss"]] = eclip[["peak_overlaps_3ss"]] | eclip[["peak_overlaps_5ss"]]
    eclip[["peak_covers_exon"]] = (
        (eclip[["peak_start"]] < eclip[["EVENT_start"]]) & (eclip[["EVENT_end"]] < eclip[["peak_end"]])
    )
    eclip = eclip %>%
        mutate(
            peak_position = case_when(
                peak_on_intron ~ "peak_on_intron",
                peak_on_exon ~ "peak_on_exon",
                peak_on_ss ~ "peak_on_ss",
                peak_covers_exon ~ "peak_covers_exon"
            ),
            bin_abs_distance_event_to_peak = cut(abs_distance_event_to_peak, breaks=seq(0,250,50))
        )
    
    events_oi = eclip %>% distinct(EVENT) %>% pull()
    samples_oi = metadata %>%
        filter(KD_GENE %in% (eclip %>% pull(rbp))) %>%
        filter(KD_GENE %in% c("AQR","SF3B1")) %>%
        distinct(sampleID) %>%
        pull(sampleID)
    dpsi = delta_psi %>%
        # only measured eCLIPs
        dplyr::select(all_of(c("EVENT",samples_oi))) %>%
        # drop all NA
        filter(rowSums(!is.na(across(all_of(samples_oi)))) > 0) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="deltaPSI") %>%
        drop_na() %>%
        left_join(metadata %>% distinct(sampleID,cell_line,KD_GENE), by="sampleID") %>%
        group_by(EVENT, cell_line, KD_GENE) %>%
        summarize(deltaPSI = median(deltaPSI, na.rm=TRUE)) %>%
        ungroup() %>%
        drop_na(deltaPSI) %>%
        left_join(
            eclip %>% 
            distinct(EVENT, cell_line, rbp) %>%
            mutate(binds_rbp = TRUE), 
            by=c("EVENT","cell_line","KD_GENE"="rbp")
        ) %>%
        mutate(
            binds_rbp = replace_na(binds_rbp, FALSE),
            abs_deltaPSI = abs(deltaPSI)
        )
    
    # plot
    plts = make_plots()

    # make figdata
    figdata = make_figdata()
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}