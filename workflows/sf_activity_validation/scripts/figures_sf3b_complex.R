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
require(scales)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)

# variables
RANDOM_SEED = 1234

MUTATIONS = c(
    "DMSO_AND_PHF5A_Y36C_MUTATION",
    "DMSO_AND_SF3B1_K700E_MUTATION",
    "MUT_SF3B1_K700E",
    "MUT_SUGP1_P636L"
)


DRUGS = c(
    "ISOGINKGETIN",
    "H3B-8800",
    "SPLICEOSTATIN_A",
    "PLADIENOLIDE_B",
    "E7107"
)


SF3b_COMPLEX = c(
    "SF3B1", "SF3B2", "SF3B3", "SF3B4", 
    "SF3B5", "SF3B6", "PHF5A", "DDX42"
)

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
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","sf3b_complex-EX.tsv.gz")
# metadata_file = file.path(RESULTS_DIR,"files","metadata","sf3b_complex-EX.tsv.gz")
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# shortest_paths_file = file.path(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_sf3b_complex.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","validation_sf3b_complex")

##### FUNCTIONS #####
plot_activity_mutations = function(protein_activity){
    plts = list()
    
    # How are spliceosome mutations funneled in the transcriptome?
    X = protein_activity %>%
        filter(condition %in% MUTATIONS) %>%
        # drop PRJNA375102 becasue there is no wt SF3B1
        filter(study_accession != "PRJNA375102") %>%
        rowwise() %>%
        mutate(is_regulator_oi = GENE %in% SF3b_COMPLEX) %>%
        ungroup()
    
    x = X %>%
        mutate(
            condition = factor(condition, levels=MUTATIONS)
        ) %>%
        group_by(study_accession, cell_line_name, condition) %>%
        mutate(nudging = 0.25*(activity_ranking - mean(activity_ranking))/sd(activity_ranking)) %>%
        ungroup() %>%
        arrange(condition)
    
    labels_sf3b_complex = x %>% 
            group_by(study_accession, cell_line_name, condition) %>% 
            slice_max(order_by=condition, n=1) %>%
            ungroup() %>% 
            filter(is_regulator_oi)
    
    labels_top = x %>% 
        filter(!is_regulator_oi) %>%
        group_by(study_accession, cell_line_name, condition) %>% 
        slice_max(order_by=activity, n=3) %>%
        ungroup() %>%
        bind_rows(
            x %>% 
            filter(!is_regulator_oi) %>%
            group_by(study_accession, cell_line_name, condition) %>% 
            slice_min(order_by=activity, n=3) %>%
            ungroup()
        )
    
    plts[["activity_mutations-sf3b_complex-scatter_line"]] = x %>%
        ggplot(aes_string(x="condition", y="activity", group="GENE")) +
        geom_point(color="lightgray", size=1, position=position_nudge(x= x %>% pull(nudging))) +
        geom_point(
            aes(color=GENE), 
            . %>% filter(is_regulator_oi) %>% mutate(),
            size=1, 
            position = position_nudge(x = x %>% filter(is_regulator_oi) %>% pull(nudging))
        ) +
        geom_text_repel(
            aes(label=GENE, color=GENE),
            labels_sf3b_complex,
            position = position_nudge(x = (labels_sf3b_complex %>% pull(nudging))),
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
        labs(x="Mutation", y="Protein Activity", color="SF3b Complex")   
    
    return(plts)
}


plot_activity_drugs = function(protein_activity, shortest_paths){
    plts = list()
    
    # How are spliceosome mutations funneled in the transcriptome?
    X = protein_activity %>%
        filter(condition %in% DRUGS) %>%
        rowwise() %>%
        mutate(is_regulator_oi = GENE %in% SF3b_COMPLEX) %>%
        ungroup()
    
    x = X %>%
        mutate(
            study_accession = factor(
                study_accession, levels = c(
                    # ISOF, SSA
                    "PRJNA292827",
                    # H3B-8800
                    "PRJNA371421",
                    # E7107
                    "PRJNA354957","PRJNA380104",
                    # PLAD B
                    "PRJNA662572","PRJNA685790"
                )
            )
        ) %>%
        group_by(study_accession, cell_line_name, condition) %>%
        mutate(nudging = 0.25*(activity_ranking - mean(activity_ranking))/sd(activity_ranking)) %>%
        ungroup() %>%
        arrange(condition)
    
    labels_sf3b_complex = x %>% 
            group_by(study_accession, cell_line_name, condition_lab) %>% 
            slice_max(order_by=condition_lab, n=1) %>%
            ungroup() %>% 
            filter(is_regulator_oi)
    
    labels_top = x %>% 
        filter(!is_regulator_oi) %>%
        group_by(study_accession, cell_line_name, condition_lab) %>% 
        slice_max(order_by=activity, n=3) %>%
        ungroup() %>%
        bind_rows(
            x %>% 
            filter(!is_regulator_oi) %>%
            group_by(study_accession, cell_line_name, condition_lab) %>% 
            slice_min(order_by=activity, n=3) %>%
            ungroup()
        )
    
    plts[["activity_drugs-sf3b_complex-scatter_line"]] = x %>%
        ggplot(aes_string(x="condition_lab", y="activity", group="GENE")) +
        geom_point(color="lightgray", size=1, position=position_nudge(x= x %>% pull(nudging))) +
        geom_point(
            aes(color=GENE), 
            . %>% filter(is_regulator_oi) %>% mutate(),
            size=1, 
            position = position_nudge(x = x %>% filter(is_regulator_oi) %>% pull(nudging))
        ) +
        geom_text_repel(
            aes(label=GENE, color=GENE),
            labels_sf3b_complex,
            position = position_nudge(x = (labels_sf3b_complex %>% pull(nudging))),
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
        labs(x="Drug", y="Protein Activity", color="SF3b Complex") 
    
    
    # Are genes with altered activities closer to the SF3b complex?
    y = shortest_paths %>%
        group_by(target) %>%
        slice_min(shortest_path_length, n=1, with_ties=FALSE) %>%
        ungroup() %>%
        left_join(x, by=c("target"="GENE")) %>%
        drop_na(study_accession, cell_line_name, condition_lab) %>%
        mutate(shortest_path_length_lab = factor(
            shortest_path_length_lab, levels=c("0","1","2","3",">=4")
        ))
    
    plts[["activity_drugs-shortest_paths-bar"]] = y %>%
        distinct(source, target, shortest_path_length_lab) %>%
        count(shortest_path_length_lab) %>%
        ggbarplot(x="shortest_path_length_lab", y="n", fill=PAL_DARK, color=NA) +
        labs(x="Shortest Path Length", y="Counts")
    
    plts[["activity_drugs-sf3b_complex_vs_shortest_paths-box"]] = y %>%
        ggboxplot(x="shortest_path_length_lab", y="activity", color="condition", outlier.size=0.1) +
        facet_wrap(~study_accession+cell_line_name+condition_lab, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Shortest Path Length", y="Protein Activity", color="Drug")
    
    return(plts)
}


make_plots = function(protein_activity, shortest_paths){
    plts = list(
        plot_activity_mutations(protein_activity),
        plot_activity_drugs(protein_activity, shortest_paths)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity, shortest_paths){
    figdata = list(
        "validation_sf3b_complex" = list(
            "protein_activity" = protein_activity,
            "shortest_paths" = shortest_paths
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
    save_plt(plts, "activity_mutations-sf3b_complex-scatter_line", '.pdf', figs_dir, width=10, height=12)
    save_plt(plts, "activity_drugs-sf3b_complex-scatter_line", '.pdf', figs_dir, width=15, height=12)
    save_plt(plts, "activity_drugs-shortest_paths-bar", '.pdf', figs_dir, width=4, height=2)
    save_plt(plts, "activity_drugs-sf3b_complex_vs_shortest_paths-box", '.pdf', figs_dir, width=10, height=14)
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
        make_option("--shortest_paths_file", type="character"),
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
    shortest_paths_file = args[["shortest_paths_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    splicing_factors = read_tsv(splicing_factors_file)
    shortest_paths = read_tsv(shortest_paths_file)
    gc()
    
    # prep
    metadata = metadata %>% 
        mutate(pert_concentration=as.numeric(gsub(",","\\.",pert_concentration)))
    
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, activity) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, condition, pert_time, pert_time_units, 
                 pert_concentration, pert_concentration_units, cell_line_name, study_accession,
                 PERT_ENSEMBL, PERT_GENE, regulator) %>%
        summarize(
            activity = median(activity, na.rm=TRUE),
            abs_activity = abs(activity),
        ) %>%
        ungroup() %>%
        
        # add activity
        group_by(condition_lab) %>%
        arrange(activity) %>%
        mutate(
            activity_ranking = row_number(),
            total_avail_sfs = sum(regulator %in% unlist(strsplit(PERT_ENSEMBL, ",")))
        ) %>%
        arrange(abs_activity) %>%
        mutate(
            abs_activity_ranking = row_number(),
        ) %>%
        ungroup() %>%
        left_join(splicing_factors, by=c("regulator"="ENSEMBL")) %>%
        mutate(                    
            condition_lab = ifelse(
                condition=="PLADIENOLIDE_B", 
                sprintf("%s (%s%s)", condition, pert_concentration, pert_concentration_units), 
                condition
            )
        )
    
    shortest_paths = shortest_paths %>%
        drop_na(shortest_path_length) %>%
        mutate(
            shortest_path_length_lab = case_when(
                shortest_path_length >= 4 ~ ">=4",
                shortest_path_length < 4 ~ as.character(shortest_path_length)
            )
        )
    
    # plot
    plts = make_plots(protein_activity, shortest_paths)
    gc()
    
    # make figdata
    figdata = make_figdata(protein_activity, shortest_paths)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}