#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - Validate whether protein activity inferrence indicates that indisulam targets RBM39
# - Investigate the effect of MS023, inhibitor of Type I PRMT enzymes
#     - no cell growth effect in vitro, but strong suppression of tumor growth in vivo

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)

# variables
RANDOM_SEED = 1234

PERTURBATIONS = c(
    "KD_ESRP1_AND_KHDRBS3",
    "KD_EXOSC3_AND_SRRT",
    "KD_ZCCHC8_AND_CBP80",
    "KD_ZCCHC8_AND_SRRT",
    "KO_SCAF4_AND_SCAF8",
    "KD_CLK3_AND_CLK4",
    "KD_CLK1_AND_CLK2_AND_CLK3",
    "KD_CLK1_AND_CLK2_AND_CLK4",
    "KD_CLK1_AND_CLK2_AND_CLK3_AND_CLK4"
)

STUDIES = c("PRJNA223244","PRJNA498529","PRJNA587741","PRJNA321560")


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
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_validation")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","combinatorial_perturbations-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,'metadata','ENASFS.tsv.gz')
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_combinatorial_perturbations")

##### FUNCTIONS #####
plot_activity = function(protein_activity){
    plts = list()
    
    X = protein_activity %>%
        filter(total_avail_sfs>1) %>%
        rowwise() %>%
        mutate(is_regulator_oi = regulator %in% unlist(strsplit(PERT_ENSEMBL, ","))) %>%
        ungroup()
    
    plts[["activity-double_perturbation_rep-ranking-scatter"]] = X %>%
        ggplot(aes(x=activity_ranking, y=activity)) +
        geom_scattermore(data = . %>% filter(!is_regulator_oi), pixels=c(1000,1000), pointsize=4, color=PAL_DARK, alpha=0.5) +
        geom_scattermore(data = . %>% filter(is_regulator_oi), pixels=c(1000,1000), pointsize=8, color=PAL_ACCENT) +
        geom_text_repel(
            aes(label=GENE),
            . %>% filter(is_regulator_oi),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        theme_pubr() +
        facet_wrap(~condition_lab+cell_line_name, ncol=4) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Ranking", y="Protein Activity", color="Is Perturbed")
    
    x = X %>%
        group_by(GENE, is_regulator_oi, condition, cell_line_name, study_accession) %>%
        summarize(
            activity = median(activity, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        group_by(condition) %>%
        arrange(activity) %>%
        mutate(
            activity_ranking = row_number()
        ) %>%
        ungroup() %>%
        mutate(
            condition = factor(condition, levels=PERTURBATIONS),
            study_accession = factor(study_accession, levels=STUDIES)
        ) %>%
        group_by(study_accession, cell_line_name, condition) %>%
        mutate(nudging = 0.25*(activity_ranking - mean(activity_ranking))/sd(activity_ranking)) %>%
        ungroup() %>%
        arrange(condition)
    
    labels_top = x %>% 
        filter(!is_regulator_oi) %>%
        group_by(study_accession, cell_line_name, condition) %>% 
        slice_min(order_by=activity, n=3) %>%
        ungroup()
    
    plts[["activity-double_perturbation_combined-ranking-scatter"]] = x %>%
        ggplot(aes(x=condition, y=activity)) +
        geom_point(color="lightgray", size=1, position=position_nudge(x= x %>% pull(nudging))) +
        geom_point(data = . %>% filter(is_regulator_oi), color="darkred", 
                   size=1, position=position_nudge(x= x %>% filter(is_regulator_oi) %>% pull(nudging))) +
        geom_text_repel(
            aes(label=GENE), 
            . %>% filter(is_regulator_oi),
            color="darkred",
            position = position_nudge(x = (x %>% filter(is_regulator_oi) %>% pull(nudging))),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50, min.segment.length=0
        ) +
        geom_text_repel(
            aes(label=GENE),
            labels_top,
            position = position_nudge(x = (labels_top %>% pull(nudging))),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50, min.segment.length=0
        ) +
        theme_pubr(x.text.angle=70) +
        facet_grid(~study_accession+cell_line_name, scales="free_x", space="free_x") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Perturbation", y="Protein Activity")
    
    return(plts)
}


make_plots = function(protein_activity){
    plts = list(
        plot_activity(protein_activity)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity){
    figdata = list(
        "validation_combinatorial_perturbation" = list(
            "protein_activity" = protein_activity
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
    save_plt(plts, "activity-double_perturbation_rep-ranking-scatter", '.pdf', figs_dir, width=10, height=12)
    save_plt(plts, "activity-double_perturbation_combined-ranking-scatter", '.pdf', figs_dir, width=15, height=10)
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
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    splicing_factors = read_tsv(splicing_factors_file)
    gc()
    
    # prep
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition) %>%
        mutate(
            condition_lab = sprintf("%s_rep%s", condition, replicate)
        ) %>%
        
        # summarize replicates, if available (only 'KO_SCAF4_AND_SCAF8')
        group_by(cell_line_name, condition_lab, condition, PERT_ENSEMBL, PERT_GENE, regulator, study_accession) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        
        # add activity
        group_by(condition_lab) %>%
        arrange(activity) %>%
        mutate(
            activity_ranking = row_number(),
            total_avail_sfs = sum(regulator %in% unlist(strsplit(PERT_ENSEMBL, ",")))
        ) %>%
        ungroup() %>%
        left_join(splicing_factors, by=c("regulator"="ENSEMBL"))
    
    # plot
    plts = make_plots(protein_activity)
    gc()
    
    # make figdata
    figdata = make_figdata(protein_activity)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}