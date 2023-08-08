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
require(ensembldb)
require(EnsDb.Hsapiens.v86)
require(tidyverse)
require(ggpubr)
require(scales)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)

# variables
edb = EnsDb.Hsapiens.v86

RANDOM_SEED = 1234

PTM_METHYLATION = c(
    "GSK3326595",
    "MS023"
)

PTM_PHOSPHORYLATION = c(
    "T3",
    "PALBOCICLIB",
    "T025",
    "KH-CB19",
    "KD_CLK1",
    "KD_CLK2",
    "KD_CLK3",
    "KD_CLK4",
    "KD_CLK3_AND_CLK4",
    "KD_CLK1_AND_CLK2_AND_CLK3",
    "KD_CLK1_AND_CLK2_AND_CLK4",
    "KD_CLK1_AND_CLK2_AND_CLK3_AND_CLK4"
)

PTM_ACETYLATION = c(
    "PHF5A_K29Q_MUTATION"
)

U2_snRNP_COMPLEX = c("PHF5A", "SF3B1", "SF3B3", "SF3B5")

SR_PROTEINS = c(
    "SRSF1",
    "SRSF2",
    "SRSF3",
    "SRSF4",
    "SRSF5",
    "SRSF6",
    "SRSF7",
    "SRSF8",
    "SRSF9",
    "SRSF10",
    "SRSF11",
    "SRSF12",
    "TRA2A",
    "TRA2B"
)

CLK_DRUGS = c("T025","T3","PALBOCICLIB","KH-CB19")

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
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","ptms-EX.tsv.gz")
# metadata_file = file.path(RESULTS_DIR,"files","metadata","ptms-EX.tsv.gz")
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_ptms")

##### FUNCTIONS #####
plot_activity_acetylation = function(protein_activity){
    # PHF5A in the SF3B complex gets acetylated normally, but the mutation prevents that.
    # In the publication (https://doi.org/10.1016/j.molcel.2019.04.009): 
    #     - acetylation increases the expression of KDM3A, which activates Wnt pathway and tumorigenesis
    #     - affects interaction of U2 snRNP complex (PHF5A, SF3B1, SF3B3, SF3B5)
    #     - acetylation of PHF5A at K29Q promotes colorectal cancer cell proliferation and tumor volume
    # (TODO) use this dataset to check activities of cancer-driver splicing factors
    
    plts = list()
    
    X = protein_activity %>%
        filter(condition %in% PTM_ACETYLATION) %>%
        rowwise() %>%
        mutate(is_regulator_oi = GENE %in% U2_snRNP_COMPLEX) %>%
        ungroup()
    
    plts[["activity_acetylation-double_perturbation_rep-ranking-scatter"]] = X %>%
        ggplot(aes(x=activity_ranking, y=activity)) +
        geom_scattermore(data = . %>% filter(!is_regulator_oi), pixels=c(1000,1000), pointsize=4, color=PAL_DARK, alpha=0.5) +
        geom_scattermore(data = . %>% filter(is_regulator_oi), pixels=c(1000,1000), pointsize=8, color=PAL_ACCENT) +
        geom_text_repel(
            aes(label=GENE),
            . %>% filter(is_regulator_oi),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        geom_text_repel(
            aes(label=GENE),
            . %>% group_by(condition_lab) %>% slice_min(activity_ranking, n=5) %>% ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50, min.segment.length=0
        ) +
        theme_pubr() +
        facet_wrap(~condition_lab, ncol=4) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Ranking", y="Protein Activity", color="Is Perturbed")
    
    plts[["activity_acetylation-double_perturbation_combined-ranking-scatter"]] = X %>%
        group_by(GENE, is_regulator_oi, condition, cell_line_name) %>%
        summarize(
            activity = median(activity, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        group_by(condition) %>%
        arrange(activity) %>%
        mutate(
            activity_ranking = row_number(),
            GENE = sprintf("%s (%s)", GENE, activity_ranking),
            GENE = ifelse(is_regulator_oi, sprintf("*%s", GENE), GENE)
        ) %>%
        ungroup() %>%
        ggplot(aes(x=activity_ranking, y=activity)) +
        geom_point(aes(color=condition, shape=cell_line_name), size=1, alpha=0.5) +
        color_palette(get_palette("npg", 15)) + 
        geom_text_repel(
            aes(label=GENE),
            . %>% filter(is_regulator_oi),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50, min.segment.length=0
        ) +
        geom_text_repel(
            aes(label=GENE),
            . %>% slice_min(activity_ranking, n=3),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50, min.segment.length=0
        ) +
        theme_pubr() +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Ranking", y="median(Protein Activity)", color="Perturbation", shape="Cell Line")
    
    return(plts)
}


plot_activity_phosphorylation = function(protein_activity){
    # SR proteins are heavily phosphorylated
    #     - SRSF1 by SRPK1
    #     - CLKs and SRSFs
    #     - DYRKs
    #     - CDK11 phosphorylates
    # CLKs are CDC-like Kinases
    #     - T3 and KH-CB19 are CLK inhibitors, from RNA-seq they define events that monotonically include/exclude
    #       with drug concentration
    #     - T3 high specificity to CLK1-3 proteins
    
    set.seed(RANDOM_SEED)
    
    plts = list()
    
    X = protein_activity %>%
        filter(condition %in% PTM_PHOSPHORYLATION) %>%
        rowwise() %>%
        mutate(is_regulator_oi = GENE %in% SR_PROTEINS) %>%
        ungroup()
    
    # CLKs are known to modulate the activity of SR rich proteins
    for (condition_oi in CLK_DRUGS){
        
        x_var = ifelse(condition_oi=="PALBOCICLIB", "pert_time", "pert_concentration")
        x_lab = ifelse(condition_oi=="PALBOCICLIB", "Time ()", "Concentration ()")
        
        x = X %>%
            filter(condition==condition_oi) %>%
            mutate(
                GENE = factor(GENE, levels=rev(SR_PROTEINS)),
                pert_time = factor(pert_time),
                pert_concentration = factor(pert_concentration, levels=sort(unique(pert_concentration)))
            ) %>%
            group_by(study_accession, cell_line_name, condition) %>%
            mutate(nudging = 0.25*(activity_ranking - mean(activity_ranking))/sd(activity_ranking)) %>%
            ungroup() %>%
            arrange(get(x_var))
        
        labels = x %>% 
            group_by(study_accession, cell_line_name, condition) %>% 
            slice_max(order_by=get(x_var), n=1) %>%
            ungroup()
        
        plts[[sprintf("activity_phosphorylation-sr_proteins-%s-line", condition_oi)]] = x %>%
            ggplot(aes_string(x=x_var, y="activity", group="GENE")) +
            geom_point(color="lightgray", size=1, position=position_nudge(x= x %>% pull(nudging))) +
            geom_path(
                aes(color=GENE), 
                . %>% filter(is_regulator_oi),
                linetype="dashed", size=LINE_SIZE,
                position = position_nudge(x = x %>% filter(is_regulator_oi) %>% pull(nudging))
            ) +
            geom_point(
                aes(color=GENE), 
                . %>% filter(is_regulator_oi),
                size=1, 
                position = position_nudge(x = x %>% filter(is_regulator_oi) %>% pull(nudging))
            ) +
            geom_text_repel(
                aes(label=GENE),
                labels %>% filter(is_regulator_oi),
                position = position_nudge(x = (labels %>% filter(is_regulator_oi) %>% pull(nudging))),
                size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50, min.segment.length=0,
                direction = "y", vjust = .5, hjust = -2
            ) +
            facet_grid(~study_accession+cell_line_name+condition, scales="free_x", space="free_x") +
            theme_pubr() +
            guides(color="none") +
            labs(x=x_lab, y="SR protein", fill="Protein Activity")
    }
    
    ## how ranking changes with concentration and cell line for each SR protein
    for (condition_oi in CLK_DRUGS){
        
        x_var = ifelse(condition_oi=="PALBOCICLIB", "pert_time", "pert_concentration")
        x_lab = ifelse(condition_oi=="PALBOCICLIB", "Time ()", "Concentration ()")
        
        plts[[sprintf("activity_phosphorylation-sr_proteins-%s-heatmap", condition_oi)]] = X %>%
            filter(condition==condition_oi & is_regulator_oi) %>%
            mutate(
                GENE = factor(GENE, levels=rev(SR_PROTEINS)),
                pert_time = factor(pert_time),
                pert_concentration = factor(pert_concentration, levels=sort(unique(pert_concentration)))
            ) %>%
            ggplot(aes_string(x=x_var, y="GENE")) +
            geom_tile(aes(fill=activity)) +
            scale_fill_gradient2(low="blue", mid="white", high="red") +
            facet_grid(~study_accession+cell_line_name+condition, scales="free_x", space="free_x") +
            theme_pubr() +
            labs(x=x_lab, y="SR protein", fill="Protein Activity")
    }
    
    ## compare changes in activity in SR proteins vs a random set of genes
    avail_sr_proteins = X %>% filter(is_regulator_oi) %>% pull(GENE) %>% unique()
    random_genes = sample(X %>% pull(GENE) %>% unique(), length(avail_sr_proteins))
    for (condition_oi in CLK_DRUGS){
        x = X %>%
            filter(condition==condition_oi & is_regulator_oi) %>%
            mutate(gene_set = "real") %>%
            bind_rows(
                X %>%
                    filter(condition==condition_oi) %>%
                    mutate(
                        is_regulator_oi = GENE %in% random_genes,
                        gene_set = "random"
                    ) %>%
                    filter(is_regulator_oi)
            )

        x_var = ifelse(condition_oi=="PALBOCICLIB", "pert_time", "pert_concentration")
        x_lab = ifelse(condition_oi=="PALBOCICLIB", "Time ()", "Concentration ()")

        plts[[sprintf("activity_phosphorylation-sr_proteins-%s-box", condition_oi)]] = x %>%
            mutate(
                GENE = factor(GENE, levels=rev(SR_PROTEINS)),
                pert_time = factor(pert_time),
                pert_concentration = factor(pert_concentration, levels=sort(unique(pert_concentration))),
            ) %>%
            ggplot(aes_string(x=x_var, y="activity", group=sprintf("interaction(%s,gene_set)", x_var))) +
            geom_boxplot(aes(fill=gene_set), outlier.size=0.1) +
            facet_grid(~study_accession+cell_line_name+condition+replicate, scales="free_x", space="free_x") +
            theme_pubr() +
            stat_compare_means(method="wilcox.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
            labs(x=x_lab, y="|Protein Activity|", fill="Gene Set")
    }
    
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
    save_plt(plts, "activity-double_perturbation_combined-ranking-scatter", '.pdf', figs_dir, width=6, height=8)
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
    
    # get protein sequences for each splicing factor
    #     gene_ids = splicing_factors %>% pull(ENSEMBL)
    #     sf_proteins = lapply(gene_ids, function(gene_id){
    #         proteins(edb, filter=GeneIdFilter(gene_id)) %>% as.data.frame()
    #     }) %>% bind_rows()
    #     sf_proteins = sf_proteins %>%
    #         mutate(
    #             sr_n = str_count(protein_sequence, "S|R"),
    #             sr_length = str_length(protein_sequence),
    #             sr_perc = sr_n / sr_length
    #         ) 
    #     sr_content = sf_proteins %>%
    #         group_by(gene_id) %>%
    #         summarize(sr_perc = median(sr_perc)) %>%
    #         ungroup()

    #     splicing_factors = splicing_factors %>%
    #         left_join(sr_content, by=c("ENSEMBL"="gene_id"))
    
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
        left_join(splicing_factors, by=c("regulator"="ENSEMBL"), relationship="many-to-many")
    
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