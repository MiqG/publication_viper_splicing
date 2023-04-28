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
PAL_CELL_LINES = "Dark2"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_Lu2021")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","Lu2021.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","Lu2021.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Lu2021.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_drug_targets")
# gene_oi = "ENSG00000131051"

##### FUNCTIONS #####
plot_genexpr_gene_oi = function(genexpr, metadata, gene_oi){
    plts = list()
    
    X = genexpr %>%
        filter(ID == gene_oi) %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_tpm") %>%
        left_join(metadata, by="sampleID")
    
    plts[["genexpr-box"]] = X %>%
        ggplot(aes(x=condition, y=genexpr_tpm)) +
        geom_point(aes(color=cell_line_name), position=position_jitter(0.1)) +
        geom_boxplot(fill=NA, width=0.25, outlier.size=0.1) +
        color_palette(PAL_CELL_LINES) +
        theme_pubr() +
        stat_compare_means(method="t.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY, ref.group="DMSO") + 
        labs(x="Condition", y="log2(TPM+1)", color="Cell Line", subtitle=gene_oi)
    
    names(plts) = sprintf("%s-%s",names(plts),gene_oi)
    
    return(plts)
}


plot_activity_gene_oi = function(protein_activity, metadata, gene_oi){
    plts = list()
    
    X = protein_activity %>%
        filter(regulator == gene_oi) %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID")
    
    plts[["activity-box"]] = X %>%
        ggplot(aes(x=condition, y=activity)) +
        geom_point(aes(color=cell_line_name), position=position_jitter(0.1)) +
        geom_boxplot(fill=NA, width=0.25, outlier.size=0.1) +
        color_palette(PAL_CELL_LINES) +
        theme_pubr() +
        stat_compare_means(method="t.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY, ref.group="DMSO") + 
        labs(x="Condition", y="Protein Activity", color="Cell Line", subtitle=gene_oi)
    
    names(plts) = sprintf("%s-%s",names(plts),gene_oi)
    
    return(plts)
}

make_plots = function(genexpr, protein_activity, metadata){
    plts = list(
        plot_genexpr_gene_oi(genexpr, metadata, "ENSG00000131051"),
        plot_activity_gene_oi(protein_activity, metadata, "ENSG00000131051")
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(genexpr, protein_activity, metadata){
    figdata = list(
        "assocs_evaluation" = list(
            "associations" = assocs,
            "summary_statisitcs" = summary_stats,
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
    save_plt(plts, "examples-high_mi_vs_low_spear-scatter", '.pdf', figs_dir, width=4, height=4)
    
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
        make_option("--genexpr_file", type="character"),
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    genexpr_file = args[["genexpr_file"]]
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    genexpr = read_tsv(genexpr_file)
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    gc()
    
    # prep
    metadata = metadata %>%
        mutate(condition = factor(condition, levels=c("DMSO","INDISULAM","MS023")))
    
    # plot
    plts = make_plots(genexpr, protein_activity, metadata)
    gc()
    
    # make figdata
    figdata = make_figdata(genexpr, protein_activity, metadata)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}