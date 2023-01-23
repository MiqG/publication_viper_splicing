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
require(ggrepel)
require(extrafont)


# variables

# formatting
FONT_SIZE = 2
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,'data','prep')
# eclip_file = file.path(PREP_DIR,"eclip_peaks_mapped","merged.tsv.gz")
# delta_psi_hepg2_file = file.path(PREP_DIR,"ground_truth_kd","ENCORE","HepG2","delta_psi-EX.tsv.gz")
# delta_psi_k562_file = file.path(PREP_DIR,"ground_truth_kd","ENCORE","K562","delta_psi-EX.tsv.gz")
# delta_psi_hepg2_file = file.path(PREP_DIR,"ground_truth_kd","ENCORE","HepG2","delta_psi_rel-EX.tsv.gz")
# delta_psi_k562_file = file.path(PREP_DIR,"ground_truth_kd","ENCORE","K562","delta_psi_rel-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,'metadata','ENCORE.tsv.gz')
# figs_dir = file.path(ROOT,'results','eda','figures','eda')

##### FUNCTIONS #####
plot_comparison = function(eclip, binding_dpsi){
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
    plts[["comparison-dpsi_vs_rbp_binding-scatter"]] = binding_dpsi %>%
        ggplot(aes(x=delta_abs_psi, y=-log10(p.adj))) +
        geom_point(size=1) +
        theme_pubr() +
        facet_wrap(~cell_line) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text_repel(
            aes(label=KD_GENE),
            . %>% group_by(cell_line) %>% slice_max(delta_abs_psi, n=1) %>% ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="|deltaPSI| Binding vs Not Binding", y="-log10(FDR)")
    
    return(plts)
}


make_plots = function(eclip, binding_dpsi){
    plts = list(
        plot_comparison(eclip, binding_dpsi)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(eclip, binding_dpsi, dpsi){

    figdata = list(
        "eclip_vs_kd" = list(
            "eclip" = eclip,
            "binding_dpsi" = binding_dpsi,
            "dpsi" = dpsi
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
    save_plt(plts, "comparison-fc_vs_pvalue", ".pdf", figs_dir, width=8, height=5)
    save_plt(plts, "comparison-dpsi_vs_rbp_binding-scatter", ".pdf", figs_dir, width=8, height=5)    
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
        make_option("--delta_psi_hepg2_file", type="character"),
        make_option("--delta_psi_k562_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    eclip_file = args[["eclip_file"]]
    delta_psi_hepg2_file = args[["delta_psi_hepg2_file"]]
    delta_psi_k562_file = args[["delta_psi_k562_file"]]
    metadata_file = args[["metadata_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    eclip = read_tsv(eclip_file)
    delta_psi_hepg2 = read_tsv(delta_psi_hepg2_file)
    delta_psi_k562 = read_tsv(delta_psi_k562_file)
    metadata = read_tsv(metadata_file)
    
    gc()
    
    # preprocess
    eclip = eclip %>%
        filter(str_detect(sample_id,"HepG2") | str_detect(sample_id,"K562")) %>%
        separate(sample_id, into = c("rbp_GENE","cell_line","idr"), remove=FALSE) %>%
        left_join(metadata %>% distinct(KD_GENE, KD_ENSEMBL), by=c("rbp_GENE"="KD_GENE"))
    
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
    
    # combine
    dpsi = rbind(
            delta_psi_hepg2 %>%
            pivot_longer(-EVENT, names_to="KD_ENSEMBL", values_to="deltaPSI") %>%
            drop_na() %>% mutate(cell_line="HepG2"),
            delta_psi_k562 %>%
            pivot_longer(-EVENT, names_to="KD_ENSEMBL", values_to="deltaPSI") %>%
            drop_na() %>% mutate(cell_line="K562")
        ) %>%  
        left_join(
            eclip %>% 
            distinct(EVENT, cell_line, KD_ENSEMBL) %>%
            mutate(binds_rbp = TRUE), 
            by=c("EVENT","cell_line","KD_ENSEMBL")
        ) %>%
        left_join(
            metadata %>% distinct(KD_GENE, KD_ENSEMBL), by="KD_ENSEMBL"
        ) %>%
        mutate(
            binds_rbp = replace_na(binds_rbp, FALSE),
            abs_deltaPSI = abs(deltaPSI)
        )
    
    # compare |dPSI| of events binding vs not binding
    kds_oi = dpsi %>% 
        count(KD_ENSEMBL, cell_line) %>% 
        count(KD_ENSEMBL) %>% 
        filter(n>1) %>%
        pull(KD_ENSEMBL)
    tests = dpsi %>%
        filter(KD_ENSEMBL %in% kds_oi) %>%
        compare_means(
            formula = abs_deltaPSI ~ binds_rbp,
            method = "wilcox.test",
            group.by = c("KD_ENSEMBL","cell_line"),
            p.adjust.method = "fdr"
        )
    binding_dpsi = dpsi %>%
        group_by(KD_GENE, KD_ENSEMBL, cell_line, binds_rbp) %>%
        summarize(med_abs_dpsi = median(abs_deltaPSI, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(binds_rbp = ifelse(binds_rbp, "abs_psi_binding", "abs_psi_nobinding")) %>%
        pivot_wider(names_from="binds_rbp", values_from="med_abs_dpsi") %>%
        mutate(delta_abs_psi = abs_psi_binding - abs_psi_nobinding) %>%
        left_join(
            tests, by=c("KD_ENSEMBL","cell_line")
        )
    
    # plot
    plts = make_plots(eclip, binding_dpsi)

    # make figdata
    figdata = make_figdata(eclip, binding_dpsi, dpsi)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}