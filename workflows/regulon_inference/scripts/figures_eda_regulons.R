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
PAL_DUAL = c("grey","orange") # '#1B9E77''#7570B3'

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","validation_activity")
# SUPPORT_DIR = file.path(ROOT,"support")
# regulons_dir = file.path(RESULTS_DIR,"files","subsetted_regulons","regulons_selected")
# figs_dir = file.path(RESULTS_DIR,'figures','eda_networks')

##### FUNCTIONS #####
plot_regulons = function(regulons){
    plts = list()
    
    X = regulons
    
    plts[["regulons-n_targets_per_regulator-box"]] = regulons %>%
        count(dataset, regulator) %>%
        ggboxplot(x="dataset", y="n", width=0.5, outlier.size=0.1, fill="orange") +
        guides(fill="none") +
        theme_pubr(x.text.angle=45) +
        theme(aspect.ratio=1) +
        yscale("log10", .format=TRUE) +
        labs(x="Dataset", y="N. Targets per Regulator")
    
    plts[["regulons-n_regulators_per_target-box"]] = regulons %>%
        count(dataset, target) %>%
        ggboxplot(x="dataset", y="n", width=0.5, outlier.size=0.1, fill="orange") +
        guides(fill="none") +
        theme_pubr(x.text.angle=45) +
        theme(aspect.ratio=1) +
        labs(x="Dataset", y="N. Regulators per Target")
    
    return(plts)
}


make_plots = function(regulons){
    plts = list(
        plot_regulons(regulons)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(regulons){
    
    figdata = list(
        "evaluation" = list(
            "regulons" = regulons
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
        save_plt(plts, "regulons-n_targets_per_regulator-box", '.pdf', figs_dir, width=5, height=5)    
        save_plt(plts, "regulons-n_regulators_per_target-box", '.pdf', figs_dir, width=5, height=5)    
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
        make_option("--regulons_dir", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    regulons_dir = args[["regulons_dir"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    regulons = lapply(list.files(regulons_dir, full.names=TRUE), function(regulons_file){
        dataset = basename(regulons_file) %>% gsub("-dpsi_morethan_15.tsv.gz","",.)
        regulons = read_tsv(regulons_file) %>%
            mutate(dataset = dataset)
        return(regulons)
    }) %>% do.call(rbind, .)
    
    # plot
    plts = make_plots(regulons)
    
    # make figdata
    figdata = make_figdata(regulons)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
