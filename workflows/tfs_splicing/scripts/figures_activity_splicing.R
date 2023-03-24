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
require(ggrepel)

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
# RESULTS_DIR = file.path(ROOT,"results","tfs_splicing")

# correlations_file = file.path(RESULTS_DIR,"files","tf_activity_vs_splicing","spearman_correlation","TCGA","PANCAN.tsv.gz")
# protein_impact_file = file.path(RAW_DIR,"VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# event_annotation_file = file.path(RAW_DIR,"VastDB","event_annotation-Hs2.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,'figures','tf_activity_vs_splicing')

##### FUNCTIONS #####
plot_correlations = function(correlations){
    plts = list()
    
    X = correlations %>%
        distinct(event_gene, median_correlation, impact) %>%
        arrange(median_correlation) %>%
        mutate(ranking = row_number())
    
    # Are there any extreme associations?
    x = X %>%
        slice_max(median_correlation, n=3) %>%
        bind_rows(
            X %>%
            slice_max(-median_correlation, n=3)
        )
    plts[["correlations-ranking-scatter"]] = X %>%
        ggscatter(x="ranking", y="median_correlation", size=0.5, color="black") +
        geom_text_repel(
            aes(label=event_gene, color=impact),
            x,
            max.overlaps=50, segment.size=0.1,
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        color_palette("jco") +
        theme(aspect.ratio=1) +
        labs(x="Ranking", y="median(Spearman Correlation)", color="Splicing Impact")
    
    # what are the correlation distributions of the extreme associations?
    X = correlations %>%
        filter(event_gene %in% x[["event_gene"]]) %>%
        arrange(median_correlation)
    plts[["correlations-top-box"]] = X %>%
        ggstripchart(x="event_gene", y="spearman_correlation", color="cancer_type", size=1) +
        geom_boxplot(width=0.1, outlier.shape=NA, color="black", fill=NA) +
        color_palette(get_palette("Paired", 19)) +
        geom_text(
            aes(label=label, y=1),
            . %>% count(event_gene) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY,
        ) +
        labs(x="Event & Gene", y="Spearman Correlation", color="Cancer Type") +
        theme_pubr(x.text.angle = 70)
    
    return(plts)
}

make_plots = function(correlations){
    plts = list(
        plot_correlations(correlations)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(correlations){
    figdata = list(
        "associations" = list(
            "correlations" = correlations
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
    save_plt(plts, "correlations-ranking-scatter", '.pdf', figs_dir, width=25, height=6.5)
    save_plt(plts, "correlations-top-box", '.pdf', figs_dir, width=5, height=10)
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
        make_option("--correlations_file", type="character"),
        make_option("--protein_impact_file", type="character"),
        make_option("--event_annotation_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    correlations_file = args[["correlations_file"]]
    protein_impact_file = args[["protein_impact_file"]]
    event_annotation_file = args[["event_annotation_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    correlations = read_tsv(correlations_file)
    protein_impact = read_tsv(protein_impact_file) %>%
        rename(EVENT=EventID, impact=ONTO)
    event_annot = read_tsv(event_annotation_file)
    
    # prep
    correlations = correlations %>%
        left_join(event_annot, by=c("ENSEMBL","EVENT")) %>%
        mutate(event_gene = paste0(EVENT,"_",GENE)) %>%
        left_join(protein_impact, by="EVENT") %>%
        filter(spearman_n_obs > 20) %>%
        group_by(ENSEMBL, EVENT) %>%
        mutate(
            median_correlation = median(spearman_correlation),
            n_cancer_types = n()
        ) %>%
        ungroup() %>%
        filter(n_cancer_types > 10)
    
    # plot
    plts = make_plots(correlations)
    
    # make figdata
    figdata = make_figdata(correlations)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}