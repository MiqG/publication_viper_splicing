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

# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# metadata_encore_kd_file = file.path(PREP_DIR,"metadata","ENCOREKD.tsv.gz")
# metadata_encore_ko_file = file.path(PREP_DIR,"metadata","ENCOREKO.tsv.gz")
# kd_screen_file = file.path(SUPPORT_DIR,"kd_screen-symbol.txt")
# ena_sfs_file = file.path(SUPPORT_DIR,"ENA_filereport-selected_sf_experiments_handcurated.tsv")
# figs_dir = file.path(RESULTS_DIR,'figures','eda')

##### FUNCTIONS #####
plot_splicing_factors = function(splicing_factors){
    plts = list()
    
    X = splicing_factors
    
    # upset plot - splicing factors considered in this analysis
    ## only SF lists
    sfs_oi = list(
        "Rogalska2022" = X %>% filter(in_rogalska) %>% pull(ENSEMBL),
        "Hegele2012" = X %>% filter(in_hegele) %>% pull(ENSEMBL),
        "Head2021" = X %>% filter(in_head) %>% pull(ENSEMBL),
        "Seiler2018" = X %>% filter(in_seiler) %>% pull(ENSEMBL),
        "Papasaikas2015" = X %>% filter(in_papasaikas) %>% pull(ENSEMBL)
    )
    plts[["splicing_factors-only_sf_lists-venn"]] = sfs_oi %>%
        ggvenn(
            fill_color = get_palette("Dark2",length(sfs_oi)),
            stroke_color = NA,
            set_name_size = FONT_SIZE+0.5,
            text_size = FONT_SIZE
        )
    
    m = sfs_oi %>% 
        list_to_matrix() %>% 
        make_comb_mat()
    plts[["splicing_factors-only_sf_lists-upset"]] = m %>%
        UpSet(comb_order = order(comb_size(m)), 
              comb_col = PAL_SINGLE_DARK,
              top_annotation = upset_top_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA)),
              right_annotation = upset_right_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA))) %>%
            draw() %>%
            grid.grabExpr() %>%
            as.ggplot()
    
    ## covered
    sfs_oi = list(
        "AllSFs" = X %>% pull(ENSEMBL),
        #"KDscreen" = X %>% filter(in_kd_screen) %>% pull(ENSEMBL),
        #"KDsENA" = X %>% filter(in_ena_sfs) %>% pull(ENSEMBL),
        "ENCORE KD" = X %>% filter(in_encore_ko) %>% pull(ENSEMBL),
        "ENCORE KO" = X %>% filter(in_encore_kd) %>% pull(ENSEMBL)
    )
    
    plts[["splicing_factors-all_datasets-venn"]] = sfs_oi %>%
        ggvenn(
            fill_color = get_palette("rickandmorty",length(sfs_oi)),
            stroke_color = NA,
            set_name_size = FONT_SIZE+0.5,
            text_size = FONT_SIZE
        )
    
    m = sfs_oi %>% 
        list_to_matrix() %>% 
        make_comb_mat()
    plts[["splicing_factors-all_datasets-upset"]] = m %>%
        UpSet(comb_order = order(comb_size(m)), 
              comb_col = PAL_SINGLE_DARK,
              top_annotation = upset_top_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA)),
              right_annotation = upset_right_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA))) %>%
            draw() %>%
            grid.grabExpr() %>%
            as.ggplot()
    
    return(plts)
}


make_plots = function(splicing_factors){
    plts = list(
        plot_splicing_factors(splicing_factors)
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
    save_plt(plts, "splicing_factors-only_sf_lists-venn", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "splicing_factors-only_sf_lists-upset", '.pdf', figs_dir, width=11.5, height=7)
    save_plt(plts, "splicing_factors-all_datasets-venn", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "splicing_factors-all_datasets-upset", '.pdf', figs_dir, width=10, height=7)
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
        make_option("--splicing_factors_file", type="character"),
        make_option("--metadata_encore_kd_file", type="character"),
        make_option("--metadata_encore_ko_file", type="character"),
        make_option("--kd_screen_file", type="character"),
        make_option("--ena_sfs_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    splicing_factors_file = args[["splicing_factors_file"]]
    metadata_encore_kd_file = args[["metadata_encore_kd_file"]]
    metadata_encore_ko_file = args[["metadata_encore_ko_file"]]
    kd_screen_file = args[["kd_screen_file"]]
    ena_sfs_file = args[["ena_sfs_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    splicing_factors = read_tsv(splicing_factors_file)
    metadata_encore_kd = read_tsv(metadata_encore_kd_file)
    metadata_encore_ko = read_tsv(metadata_encore_ko_file)
    kd_screen = readLines(kd_screen_file)
    ena_sfs = read_tsv(ena_sfs_file)
    
    # prep
    ena_sfs = ena_sfs %>% filter(IS_USEFUL & PERT_TYPE%in%c("KNOCKDOWN","KNOCKOUT","OVEREXPRESSION"))
    
    splicing_factors = splicing_factors %>%
        mutate(
            in_encore_kd = GENE %in% metadata_encore_kd[["PERT_GENE"]],
            in_encore_ko = GENE %in% metadata_encore_ko[["PERT_GENE"]],
            in_kd_screen = GENE %in% kd_screen,
            in_ena_sfs = GENE %in% ena_sfs[["PERT_GENE"]]
        )
    
    # stats
    ## overview
    n_total_sfs = splicing_factors %>% nrow()
    n_sfs_screened = splicing_factors %>% filter(in_encore_kd | in_encore_ko | in_ena_sfs) %>% nrow()
    print(sprintf("Screened %s SFs out of %s", n_sfs_screened, n_total_sfs))
    ## features
    splicing_factors %>% filter(in_encore_kd | in_encore_ko) %>% nrow()
    splicing_factors %>% filter(in_ena_sfs) %>% nrow()
    ## ENA
    ena_sfs %>% distinct(study_accession) %>% nrow()
    
    # plot
    plts = make_plots(splicing_factors)
    
    # make figdata
    figdata = make_figdata(splicing_factors)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}