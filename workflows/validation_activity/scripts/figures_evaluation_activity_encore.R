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

# evaluation_rankings_file = file.path(RESULTS_DIR,"files","subsetted_regulons","evaluation_rankings.tsv.gz")
# evaluation_corrs_file = file.path(RESULTS_DIR,"files","subsetted_regulons","evaluation_corrs.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,'figures','evaluation')

##### FUNCTIONS #####
plot_evaluation = function(evaluation_rankings, evaluation_corrs){
    plts = list()
    
    # number of regulators per threshold (TODO) - few overlaps regulator-KD between ENCORE screens!
    #     plts[["evaluation-thresh_vs_n_regulators-lessthan"]] = evaluation_rankings %>% 
    #         filter(regulator==PERT_GENE & thresh_type=="lessthan" & eval_type=="real") %>%
    #         mutate(thresh = sprintf("<=%s",thresh)) %>%
    #         count(thresh, thresh_type, eval_type, dataset_regulon, dataset_signature, cell_line) %>% 
    #         ggboxplot(x="thresh", y="n") +
    #         theme_pubr(x.text.angle=70) +
    #         facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
    #         theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    #         labs(x="|DeltaPSI| Threshold", y="Ranking Between", fill="Ranking Type")
    
    # rankings between samples
    plts[["evaluation-ranking_between-lessthan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="lessthan") %>%
        mutate(thresh = sprintf("<=%s",thresh)) %>%
        ggboxplot(x="thresh", y="ranking_between", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Between", fill="Ranking Type")
    
    plts[["evaluation-rankperc_between-lessthan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="lessthan") %>%
        mutate(thresh = sprintf("<=%s",thresh)) %>%
        ggboxplot(x="thresh", y="rankperc_between", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Between %", fill="Ranking Type")
    
    plts[["evaluation-ranking_between-morethan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="morethan") %>%
        mutate(thresh = sprintf(">=%s",thresh)) %>%
        ggboxplot(x="thresh", y="ranking_between", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Between", fill="Ranking Type")
    
    plts[["evaluation-rankperc_between-morethan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="morethan") %>%
        mutate(thresh = sprintf(">=%s",thresh)) %>%
        ggboxplot(x="thresh", y="rankperc_between", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Between %", fill="Ranking Type")
    
    # rankings within samples
    plts[["evaluation-ranking_within-lessthan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="lessthan") %>%
        mutate(thresh = sprintf("<=%s",thresh)) %>%
        ggboxplot(x="thresh", y="ranking_within", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Within", fill="Ranking Type")
    
    plts[["evaluation-rankperc_within-lessthan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="lessthan") %>%
        mutate(thresh = sprintf("<=%s",thresh)) %>%
        ggboxplot(x="thresh", y="rankperc_within", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Within %", fill="Ranking Type")
    
    plts[["evaluation-ranking_within-morethan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="morethan") %>%
        mutate(thresh = sprintf(">=%s",thresh)) %>%
        ggboxplot(x="thresh", y="ranking_within", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Within", fill="Ranking Type")
    
    plts[["evaluation-rankperc_within-morethan-box"]] = evaluation_rankings %>%
        filter(regulator==PERT_GENE & thresh_type=="morethan") %>%
        mutate(thresh = sprintf(">=%s",thresh)) %>%
        ggboxplot(x="thresh", y="rankperc_within", fill="eval_type", width=0.5, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Ranking Within %", fill="Ranking Type")
    
    # how good are correlations between dataset signatures with the different thresholded networks?
    plts[["evaluation-corrs-lessthan-bar"]] = evaluation_corrs %>%
        filter(thresh_type=="lessthan") %>%
        mutate(thresh = sprintf("<=%s",thresh)) %>%
        ggbarplot(x="thresh", y="pearson_coef", fill="eval_type", color=NA, position=position_dodge(0.9)) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Pearson Coef.", fill="Ranking Type")
    
    plts[["evaluation-corrs-morethan-bar"]] = evaluation_corrs %>%
        filter(thresh_type=="morethan") %>%
        mutate(thresh = sprintf(">=%s",thresh)) %>%
        ggbarplot(x="thresh", y="pearson_coef", fill="eval_type", color=NA, position=position_dodge(0.9)) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset_regulon+dataset_signature+cell_line, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|DeltaPSI| Threshold", y="Pearson Coef.", fill="Ranking Type")
    
    return(plts)
}


make_plots = function(evaluation_rankings, evaluation_corrs){
    plts = list(
        plot_evaluation(evaluation_rankings, evaluation_corrs)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(evaluation_rankings, evaluation_corrs){
    
    figdata = list(
        "evaluation" = list(
            "evaluation_rankings" = evaluation_rankings,
            "evaluation_corrs" = evaluation_corrs
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
    save_plt(plts, "evaluation-ranking_between-lessthan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-rankperc_between-lessthan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-ranking_between-morethan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-rankperc_between-morethan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-ranking_within-lessthan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-rankperc_within-lessthan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-ranking_within-morethan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-rankperc_within-morethan-box", '.pdf', figs_dir, width=25, height=35)
    save_plt(plts, "evaluation-corrs-lessthan-bar", '.pdf', figs_dir, width=15, height=30)
    save_plt(plts, "evaluation-corrs-morethan-bar", '.pdf', figs_dir, width=15, height=30)
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
        make_option("--evaluation_rankings_file", type="character"),
        make_option("--evaluation_corrs_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    evaluation_rankings_file = args[["evaluation_rankings_file"]]
    evaluation_corrs_file = args[["evaluation_corrs_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    evaluation_rankings = read_tsv(evaluation_rankings_file)
    evaluation_corrs = read_tsv(evaluation_corrs_file)
    
    # plot
    plts = make_plots(evaluation_rankings, evaluation_corrs)
    
    # make figdata
    figdata = make_figdata(evaluation_rankings, evaluation_corrs)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}