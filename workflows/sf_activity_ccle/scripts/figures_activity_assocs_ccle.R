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
require(ggrepel)

# variables
RANDOM_SEED = 1234

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "darkgrey"
PAL_ACCENT = "darkred"
PAL_DUAL = c(PAL_DARK, PAL_ACCENT)
PAL_CONTRAST = c("darkgrey","darkred")
PAL_CELL_LINES = "Dark2"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_ccle")
# assocs_gene_dependency_file = file.path(RESULTS_DIR,"files","protein_activity_vs_demeter2","CCLE.tsv.gz")
# assocs_metastasis_file = file.path(RESULTS_DIR,"files","protein_activity_vs_metmap","CCLE.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","activity_ccle")

##### FUNCTIONS #####
plot_assocs_metastasis = function(assocs_metastasis){
    plts = list()
    
    X = assocs_metastasis %>%
        group_by(metastatic_tissue) %>%
        arrange(spearman_coef) %>%
        mutate(assoc_ranking = row_number()) %>%
        ungroup()
    
    plts[["assocs_metastasis-rankings-scatter"]] = X %>%
        ggplot(aes(x=assoc_ranking, y=spearman_coef)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, color=PAL_DARK) +
        theme_pubr() +
        facet_wrap(~metastatic_tissue) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text_repel(
            aes(label=GENE),
            X %>% group_by(metastatic_tissue) %>% slice_max(spearman_coef, n=3) %>% ungroup() %>%
                bind_rows(
                    X %>% group_by(metastatic_tissue) %>% slice_min(spearman_coef, n=3) %>% ungroup()
                ),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1
        ) +
        labs(x="Ranking", y="Spearman Coef.")
    
    return(plts)
}


plot_assocs_gene_dependency = function(assocs_gene_dependency){
    plts = list()
    
    X = assocs_gene_dependency %>%
        arrange(spearman_coef) %>%
        mutate(assoc_ranking = row_number())
    
    plts[["assocs_gene_dependency-rankings-scatter"]] = X %>%
        ggplot(aes(x=assoc_ranking, y=spearman_coef)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, color=PAL_DARK) +
        geom_text_repel(
            aes(label=GENE),
            X %>% slice_max(spearman_coef, n=3) %>% ungroup() %>%
                bind_rows(
                    X %>% slice_min(spearman_coef, n=3) %>% ungroup()
                ),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Ranking", y="Spearman Coef.")
    
    return(plts)
}

make_plots = function(assocs_metastasis, assocs_gene_dependency){
    plts = list(
        plot_assocs_metastasis(assocs_metastasis),
        plot_assocs_gene_dependency(assocs_gene_dependency)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(assocs_metastasis, assocs_gene_dependency){
    figdata = list(
        "activity_ccle" = list(
            "assocs_metastasis" = assocs_metastasis,
            "assocs_gene_dependency" = assocs_gene_dependency
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
    save_plt(plts, "assocs_metastasis-rankings-scatter", '.pdf', figs_dir, width=13, height=9)
    save_plt(plts, "assocs_gene_dependency-rankings-scatter", '.pdf', figs_dir, width=4, height=4)
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
        make_option("--assocs_gene_dependency_file", type="character"),
        make_option("--assocs_metastasis_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    assocs_gene_dependency_file = args[["assocs_gene_dependency_file"]]
    assocs_metastasis_file = args[["assocs_metastasis_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    assocs_gene_dependency = read_tsv(assocs_gene_dependency_file)
    assocs_metastasis = read_tsv(assocs_metastasis_file)
    gc()
    
    # plot
    plts = make_plots(assocs_metastasis, assocs_gene_dependency)
    gc()
    
    # make figdata
    figdata = make_figdata(assocs_gene_dependency, assocs_metastasis)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}