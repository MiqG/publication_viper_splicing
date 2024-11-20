require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)

# variables
THRESH_FDR = 0.05

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = c(
    #"Non-driver"="lightgrey",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

PAL_DARK = "brown"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","cancer_splicing_program")
# driver_types_file = file.path(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
# protein_activity_ccle_file = file.path(RESULTS_DIR,"files","protein_activity","CCLE-EX.tsv.gz")
# genexpr_ccle_file = file.path(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
# metadata_ccle_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","proliferation-EX")

##### FUNCTIONS #####
plot_proliferation_ccle = function(protein_activity_ccle){
    plts = list()
    
    X = protein_activity_ccle
    
    # How do cancer-driver activities relate to proliferation?
    plts[["hallmarks_ccle-mki67_vs_activity_median-scatter"]] = X %>% 
        distinct(DepMap_ID, MKI67, activity, driver_type, GENE) %>%
        group_by(DepMap_ID, driver_type, MKI67) %>%
        summarize(activity_median = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        mutate(driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE))) %>%
        ggscatter(x="MKI67", y="activity_median", color="driver_type", 
                  size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) + 
        geom_density_2d(color="black", size=LINE_SIZE) +
        geom_smooth(method="lm", size=LINE_SIZE, color="black", linetype="dashed") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~driver_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log2(TPM+1) MKI67", y="median(Protein Activity SFs)")
    
    plts[["hallmarks_ccle-mki67_vs_activity_median_correlations-violin"]] = X %>% 
        distinct(DepMap_ID, MKI67, activity, driver_type, GENE, primary_disease) %>%
        group_by(DepMap_ID, driver_type, MKI67, primary_disease) %>%
        summarize(activity_median = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        group_by(driver_type, primary_disease) %>%
        summarize(
            correlation = cor(MKI67, activity_median),
            n_obs = n()
        ) %>%
        ungroup() %>%
        drop_na() %>%
        filter(n_obs > 10) %>%
        filter(!(primary_disease %in% c("Fibroblast"))) %>%
        mutate(driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE))) %>%
        ggviolin(x="driver_type", y="correlation", fill="driver_type", palette=PAL_DRIVER_TYPE, color=NA, trim=TRUE) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(fill="none") + 
        theme(aspect.ratio=1) +
        labs(x="Cancer Splicing Program", y="Pearson correlation\nlog2(TPM+1) MKI67 vs median(Protein Activity SFs)\nby Cell Line Cancer Type")
    
    return(plts)
}


make_plots = function(protein_activity_ccle){
    plts = list(
        plot_proliferation_ccle(protein_activity_ccle)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity_ccle){
    figdata = list(
        "proliferation" = list(
            "protein_activity_ccle" = protein_activity_ccle
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
    save_plt(plts, "hallmarks_ccle-mki67_vs_activity_median-scatter", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "hallmarks_ccle-mki67_vs_activity_median_correlations-violin", '.pdf', figs_dir, width=6, height=5)
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
        make_option("--protein_activity_ccle_file", type="character"),
        make_option("--genexpr_ccle_file", type="character"),
        make_option("--metadata_ccle_file", type="character"),
        make_option("--driver_types_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    protein_activity_ccle_file = args[["protein_activity_ccle_file"]]
    genexpr_ccle_file = args[["genexpr_ccle_file"]]
    metadata_ccle_file = args[["metadata_ccle_file"]]
    driver_types_file = args[["driver_types_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_ccle = read_tsv(protein_activity_ccle_file)
    genexpr_ccle = read_tsv(genexpr_ccle_file)
    metadata_ccle = read_tsv(metadata_ccle_file)
    driver_types = read_tsv(driver_types_file)
    
    # prep
    ## CCLE
    genexpr_mki67_ccle = genexpr_ccle %>%
        filter(ID == "ENSG00000148773") %>%
        pivot_longer(-ID, names_to="DepMap_ID", values_to="MKI67")
    
    protein_activity_ccle = protein_activity_ccle %>%
        pivot_longer(-regulator, names_to="DepMap_ID", values_to="activity") %>%
        left_join(
            genexpr_mki67_ccle,
            by="DepMap_ID"
        ) %>%
        left_join(
            metadata_ccle,
            by="DepMap_ID"
        ) %>%
        left_join(
            driver_types,
            by=c("regulator"="ENSEMBL")
        )
    
    # plot
    plts = make_plots(
        protein_activity_ccle
    )
    
    # make figdata
    figdata = make_figdata(
        protein_activity_ccle
    )

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}