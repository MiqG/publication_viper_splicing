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
require(ComplexHeatmap)
require(ggplotify)
require(ggvenn)

# variables
IDS = c(
    "CardosoMoreira2020-EX-k2",
    "CardosoMoreira2020-EX-k5",
    "CardosoMoreira2020-EX-k10",
    "CardosoMoreira2020-EX-k50",
    "CardosoMoreira2020-EX-k100"
)

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_SINGLE_DARK = "darkgreen"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","preprocess_data")
# SUPPORT_DIR = file.path(ROOT,"support")

# original_file = file.path(PREP_DIR,'event_psi','CardosoMoreira2020-EX.tsv.gz')
# synthetic_file = file.path(PREP_DIR,'event_psi_imputation_benchmark','CardosoMoreira2020-EX.tsv.gz')
# imputed_files = "~/projects/publication_viper_splicing/data/prep/event_psi_imputation_benchmark/CardosoMoreira2020-EX-k5.tsv.gz,~/projects/publication_viper_splicing/data/prep/event_psi_imputation_benchmark/CardosoMoreira2020-EX-k10.tsv.gz,~/projects/publication_viper_splicing/data/prep/event_psi_imputation_benchmark/CardosoMoreira2020-EX-k50.tsv.gz,~/projects/publication_viper_splicing/data/prep/event_psi_imputation_benchmark/CardosoMoreira2020-EX-k2.tsv.gz,~/projects/publication_viper_splicing/data/prep/event_psi_imputation_benchmark/CardosoMoreira2020-EX-k100.tsv.gz"
# figs_dir = file.path(RESULTS_DIR,'figures','benchmark_imputation')

##### FUNCTIONS #####
plot_imputation_benchmark = function(imputed){
    plts = list()
    
    X = imputed %>%
        mutate(imputation_id = factor(imputation_id, levels=IDS))
    
    # correlation between synthetically missing and imputed PSI
    plts[["imputation_benchmark-real_vs_imputed-scatter"]] = X %>%
        drop_na(psi_real, psi_imputed) %>%
        ggplot(aes(x=psi_real, y=psi_imputed)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4, color="black", alpha=0.5) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~missing_lab+imputation_id, nrow=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="PSI Real", y="PSI Imputed", subtitle="5 randomly selected samples from CardosoMoreira2020")
    
    return(plts)
}


make_plots = function(imputed){
    plts = list(
        plot_imputation_benchmark(imputed)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(splicing_factors){
    
    figdata = list(
        "eda" = list(
            "splicing_factors" = splicing_factors
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
    save_plt(plts, "imputation_benchmark-real_vs_imputed-scatter", '.pdf', figs_dir, width=12, height=9)
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
        make_option("--original_file", type="character"),
        make_option("--synthetic_file", type="character"),
        make_option("--imputed_files", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    original_file = args[["original_file"]]
    synthetic_file = args[["synthetic_file"]]
    imputed_files = args[["imputed_files"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    set.seed(1234)
    
    # load
    original = read_tsv(original_file) %>% 
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi_real")
    samples_oi = original %>% distinct(sampleID) %>% slice_sample(n=5) %>% pull(sampleID) # randomly select 5 samples
    original = original %>% filter(sampleID%in%samples_oi)
    
    synthetic = read_tsv(synthetic_file) %>% 
        select(all_of(c("EVENT", samples_oi))) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi_synthetic")
    
    gt_labels = synthetic %>%
        left_join(original, by=c("sampleID","EVENT")) %>%
        mutate(missing_lab = case_when(
            is.na(psi_real) & is.na(psi_synthetic) ~ "Totally Missing",
            !is.na(psi_real) & is.na(psi_synthetic) ~ "Synthetically Missing",
            !is.na(psi_real) & !is.na(psi_synthetic) ~ "Real"
        ))
    
    imputed_files = unlist(strsplit(imputed_files, ","))
    imputed = sapply(imputed_files, function(file){
        imputed = read_tsv(file) %>% 
            select(all_of(c("EVENT", samples_oi))) %>%
            pivot_longer(-EVENT, names_to="sampleID", values_to="psi_imputed") %>%
            left_join(gt_labels) %>%
            mutate(imputation_id = gsub(".tsv.gz","",basename(file)) )
        return(imputed)
    }, simplify=FALSE) %>% bind_rows()
    
    # plot
    plts = make_plots(imputed)
    
    # make figdata
    #figdata = make_figdata(imputed)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
