#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# See how indisulam affects RBM39 protein and transcriptomic expression and how it compares with inferred protein activities.
# 

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
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_Nijhuis2020")
# proteomics_file = file.path(RAW_DIR,"articles","Nijhuis2020","supplementary_data","proteomics_lfq_intensity.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","Nijhuis2020.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","Nijhuis2020-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Nijhuis2020.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_drug_targets")
# gene_oi = "ENSG00000131051" # RBM39

##### FUNCTIONS #####
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
        "validation_drug_target_activity" = list(
            "genexpr" = genexpr,
            "protein_activity" = protein_activity,
            "metadata" = metadata
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
    save_plt(plts, "genexpr-indisulam_vs_dmso-box-ENSG00000131051", '.pdf', figs_dir, width=5, height=6)
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
        make_option("--proteomics_file", type="character"),
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
    
    proteomics_file = args[["proteomics_file"]]
    genexpr_file = args[["genexpr_file"]]
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    proteomics = read_tsv(proteomics_file)
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
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}