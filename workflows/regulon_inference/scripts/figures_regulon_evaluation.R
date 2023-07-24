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


# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# evaluation_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz")

##### FUNCTIONS #####
plot_evaluation = function(evaluation){
    plts = list()
    
    X = evaluation %>%
        filter(signature_id!=regulon_id)
    
    plts[["evaluation-ranking_perc_vs_regulon_set-violin"]] = X %>%
        ggplot(aes(x=regulon_set_id, y=ranking_perc, group=interaction(regulon_set_id, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE, position=position_dodge(0.9)) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_EVAL_TYPE) + 
        facet_wrap(~eval_direction, ncol=1) +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% count(regulon_set_id, eval_direction, eval_type) %>% mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    return(plts)
}


make_plots = function(evaluation){
    plts = list(
        plot_evaluation(evaluation)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(evaluation){
    figdata = list(
        "regulon_evaluation" = list(
            "evaluation" = evaluation
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
    save_plt(plts, "evaluation-ranking_perc_vs_regulon_set-violin", '.pdf', figs_dir, width=5, height=9)
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
        make_option("--evaluation_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    evaluation_file = args[["evaluation_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    evaluation = read_tsv(evaluation_file)
    
    # prep
    evaluation = evaluation %>%
        mutate(signature_id = gsub("-","_",signature_id))
    
    # plot
    plts = make_plots(evaluation)
    
    # make figdata
    figdata = make_figdata(evaluation)
    
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}