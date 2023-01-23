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

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","sf_targets_inference")

# evaluation_file = file.path(RESULTS_DIR,"files","inference_evaluation","merged.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,'figures','target_inference')

##### FUNCTIONS #####
plot_evaluation = function(evaluation){
    plts = list()
    
    X = evaluation
    
    n_kds = X %>%
        distinct(sf_target_inference_method, kd_cell_line, KD_ENSEMBL) %>%
        count(sf_target_inference_method, kd_cell_line) %>%
        mutate(label=sprintf("%s (n=%s)",kd_cell_line,n))
    
    plts[["evaluation-thresh_vs_prop_correct-line"]] = X %>%
        left_join(n_kds, by=c("sf_target_inference_method","kd_cell_line")) %>%
        ggplot(aes(x=threshold_classification, y=prop_correct, group=KD_ENSEMBL)) +
        geom_line(size=0.1, color="grey", alpha=0.5) +
        geom_smooth(aes(color=kd_cell_line, fill=kd_cell_line, group=kd_cell_line), 
                    se=FALSE, span=0.2, size=LINE_SIZE, linetype="dashed", alpha=0.5, method="loess") +
        color_palette("Dark2") +
        theme_pubr(legend="none") +
        facet_wrap(~label+sf_target_inference_method) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Classification Threshold", y="Proportion Correct") +
        lims(x=c(0,1), y=c(0,1))
    
    x = X %>%
        left_join(n_kds, by=c("kd_cell_line","sf_target_inference_method")) %>% 
        group_by(sf_target_inference_method, label, kd_cell_line, threshold_classification) %>%
        summarize(
            med_tpr=median(tpr, na.rm=TRUE), 
            med_fpr=median(fpr, na.rm=TRUE),
            med_recall=median(recall, na.rm=TRUE),
            med_precision=median(precision, na.rm=TRUE)
        )
    plts[["evaluation-fpr_vs_tpr-line"]] = x %>%
        arrange(threshold_classification) %>%
        ggplot(aes(x=med_fpr, y=med_tpr)) +
        geom_line(aes(color=sf_target_inference_method), size=LINE_SIZE, linetype="dashed") +
        geom_point(aes(color=sf_target_inference_method), size=1) +
        color_palette("Dark2") +
        facet_wrap(~label) +
        theme_pubr() + 
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="FPR", y="TPR", color="Inference Method")
    
    plts[["evaluation-recall_vs_precision-line"]] = x %>%
        ggplot(aes(x=med_recall, y=med_precision)) +
        geom_line(aes(color=sf_target_inference_method), size=LINE_SIZE, linetype="dashed") +
        geom_point(aes(color=sf_target_inference_method), size=1) +
        color_palette("Dark2") +
        facet_wrap(~label) +
        theme_pubr() + 
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Recall", y="Precision", color="Inference Method")
    
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
        "target_inference" = list(
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
    save_plt(plts, "evaluation-thresh_vs_prop_correct-line", '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, "evaluation-fpr_vs_tpr-line", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "evaluation-recall_vs_precision-line", '.pdf', figs_dir, width=8, height=6)
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
    
    # plot
    plts = make_plots(evaluation)
    
    # make figdata
    figdata = make_figdata(evaluation)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}