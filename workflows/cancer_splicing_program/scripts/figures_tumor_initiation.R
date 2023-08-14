#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#

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
    "Non-driver"="lightgrey",
    "Tumor suppressor"="#6C98B3",
    "Oncogenic"="#F6AE2D"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","cancer_splicing_program")
# protein_activity_organoids_file = file.path(RESULTS_DIR,"files","protein_activity","Bian2018-EX.tsv.gz")
# metadata_organoids_file = file.path(PREP_DIR,"metadata","Bian2018.tsv.gz")
# protein_activity_mutations_file = file.path(RESULTS_DIR,"files","protein_activity","driver_mutations-EX.tsv.gz")
# metadata_mutations_file = file.path(RESULTS_DIR,"files","metadata","driver_mutations-EX.tsv.gz")
# driver_types_file = file.path("driver_types.tsv")

##### FUNCTIONS #####
plot_brain_organoids = function(protein_activity_organoids){
    plts = list()
    
    X = protein_activity_organoids %>%
        group_by(condition, driver_type, pert_time, replicate) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE)),
               pert_time = factor(pert_time, levels=c(45,130))) %>%
        filter(driver_type!="Non-driver")
    
    plts[["brain_organoids-condition_vs_activity-"]] = X %>%
        ggplot(aes(x=condition, y=activity)) +
        geom_boxplot(aes(fill=driver_type), position=position_dodge(0.9), outlier.shape=NA) +
        geom_jitter(position=position_dodge(0.9)) + 
        facet_wrap(~driver_type+pert_time, ncol=4) +
        stat_compare_means(ref.group="CONTROL", method="t.test", label="p.signif") +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70)
    
    return(plts)
}


plot_driver_mutations = function(protein_activity_mutations){
    plts = list()
    
    X = protein_activity_mutations %>%
        group_by(condition, driver_type, replicate, study_accession, cell_line_name) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(driver_type = factor(driver_type, levels=names(PAL_DRIVER_TYPE)))
    
    order_oi = X %>%
        filter(driver_type == "Tumor suppressor") %>%
        arrange(activity) %>%
        pull(condition)
    
    plts[["driver_mutations-"]] = X %>%
        mutate(condition = factor(condition, levels=order_oi)) %>%
        ggplot(aes(x=condition, y=activity, group=driver_type)) +
        geom_smooth(aes(color=driver_type), fill="lightgray", linetype="dashed", size=LINE_SIZE) +
        geom_point(
            aes(color=driver_type), 
            . %>% group_by(condition, driver_type) %>% summarize(activity=median(activity, na.rm=TRUE)) %>% ungroup(), 
            fill="lightgray", linetype="dashed", size=LINE_SIZE) +
        color_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 25) +
    
    x = X %>% 
        group_by(condition, driver_type, cell_line_name, study_accession) %>% 
        summarize(activity=median(activity, na.rm=TRUE)) %>% 
        ungroup() %>%
        pivot_wider(id_cols = c("study_accession","cell_line_name","condition"), 
                    names_from = "driver_type", values_from = "activity")
        
    x %>%
        ggscatter(x="Tumor suppressor", y="Oncogenic", color="cell_line_name") +
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype="dashed", size=LINE_SIZE) +
        geom_text_repel(
            aes(label=condition),
            . %>% filter(abs(`Tumor suppressor`)>0.15 | abs(`Oncogenic`)>0.15),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        labs(x="median(Tumor suppressor Activity)", y="median(Oncogenic Activity)", color="Cell Line")
    
    return(plts)
}


make_plots = function(
    protein_activity_organoids
){
    plts = list(
        plot_brain_organoids(protein_activity_organoids)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    protein_activity
){
    figdata = list(
        "tumor_initiation" = list(
            "protein_activity_organoids" = protein_activity_organoids
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
    save_plt(plts, "driver_selection-n_signif_vs_driver_type-genexpr-bar", '.pdf', figs_dir, width=12, height=9)
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
        make_option("--diff_activity_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    diff_activity_file = args[["diff_activity_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity_organoids = read_tsv(protein_activity_organoids_file)
    metadata_organoids = read_tsv(metadata_organoids_file)
    protein_activity_mutations = read_tsv(protein_activity_mutations_file)
    metadata_mutations = read_tsv(metadata_mutations_file)
    driver_types = read_tsv(driver_types_file)
    
    # prep
    protein_activity_organoids = protein_activity_organoids %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata_organoids, by="sampleID") %>%
        drop_na(condition) %>%
        mutate(
            condition_lab = sprintf("%s_rep%s (%s%s)", condition, replicate, pert_time, pert_time_units)
        ) %>%
        # summarize technical replicates
        group_by(cell_line_name, condition_lab, condition, PERT_ENSEMBL, PERT_GENE, regulator, study_accession, pert_time, replicate) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        # add driver type
        left_join(
            driver_types %>% distinct(GENE, ENSEMBL, driver_type),
            by=c("regulator"="ENSEMBL")
        ) %>%
        mutate(driver_type = replace_na(driver_type, "Non-driver"))
    
    
    protein_activity_mutations = protein_activity_mutations %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata_mutations, by="sampleID") %>%
        drop_na(condition) %>%
        mutate(
            condition_lab = sprintf("%s_rep%s (%s%s)", condition, replicate, pert_time, pert_time_units)
        ) %>%
        # summarize technical replicates
        group_by(cell_line_name, condition_lab, condition, PERT_ENSEMBL, PERT_GENE, regulator, study_accession, pert_time, replicate) %>%
        summarize(activity = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        # add driver type
        left_join(
            driver_types %>% distinct(GENE, ENSEMBL, driver_type),
            by=c("regulator"="ENSEMBL")
        ) %>%
        mutate(driver_type = replace_na(driver_type, "Non-driver"))
    
    # plot
    plts = make_plots(
        protein_activity_organoids,
        protein_activity_mutations
    )
    
    # make figdata
    figdata = make_figdata(
        protein_activity
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