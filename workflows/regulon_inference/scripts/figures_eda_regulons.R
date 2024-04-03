#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA regulons of interest
#
# Outline
# -------
# - how many targets per regulator?
# - how many regulators per target?
# - which biological processes does each regulator control?
# - which protein impacts does each regulator control?

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(ggrepel)
require(proxy)
require(umap)
require(ggbeeswarm)

# variables
SETS_MAIN = c(
    'aracne_regulons_development',
    'mlr_regulons_development',
    'experimentally_derived_regulons_pruned'
)

SF_CLASS = c("Core", "RBP", "Other")

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_LIGHT = "lightpurple"
PAL_DARK = "darkgreen"
PAL_DUAL = c("grey","orange") # '#1B9E77''#7570B3'

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# SUPPORT_DIR = file.path(ROOT,"support")
# regulons_dir = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX")
# evaluation_ex_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz")
# enrichments_file = file.path(RESULTS_DIR,"files","regulons_eda_gsea","experimentally_derived_regulons_pruned-EX.tsv.gz")
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# spliceosomedb_file = file.path(SUPPORT_DIR,"splicing_factors","literature_suptabs","SpliceosomeDB-human.csv")
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,'figures','eda_regulons-EX')

##### FUNCTIONS #####
plot_regulons = function(regulons, evaluation, splicing_factors){
    plts = list()
    
    X = regulons %>%
        distinct(regulator, target) %>%
        mutate(regulon_id = "Empirical")
    
    plts[["regulons-n_targets_per_regulator-box"]] = X %>%
        count(regulon_id, regulator) %>%
        ggplot(aes(x=regulon_id, y=n)) +
        geom_quasirandom(size=0.1, color="orange", varwidth=0.5) +
        geom_boxplot(width=0.5, outlier.shape=NA, color="black", fill=NA) +
        geom_hline(yintercept = 25, linetype="dashed", color="black", size=LINE_SIZE) +
        theme_pubr(x.text.angle=45) +
        yscale("log10", .format=TRUE) +
        guides(fill="none") +
        labs(x="Regulon ID", y="N. Targets per Regulator")
    
    plts[["regulons-n_regulators_per_target-box"]] = X %>%
        count(regulon_id, target) %>%
        ggplot(aes(x=regulon_id, y=n)) +
        geom_quasirandom(size=0.1, color="orange", varwidth=0.5) +
        geom_boxplot(width=0.5, outlier.shape=NA, color="black", fill=NA) +
        guides(fill="none") +
        theme_pubr(x.text.angle=45) +
        labs(x="Regulon ID", y="N. Regulators per Target")
    
    plts[["regulons-n_targets_per_regulator_vs_evaluation-scatter"]] = evaluation %>%
        left_join(
            X %>%
            count(regulon_id, regulator), 
            by="regulator"
        ) %>%
        ggscatter(x="n", y="ranking_perc", alpha=0.5, size=1, color=PAL_DARK) +
        xscale("log10", .format=TRUE) +
        facet_wrap(~regulon_set+eval_direction) +
        stat_cor()
    
    # RBP vs core spliceosome
    x = X %>%
        count(regulon_id, regulator) %>%
        left_join(splicing_factors, by=c("regulator"="ENSEMBL")) %>%
        mutate(
            sf_class = case_when(
                !is.na(spliceosome_db_complex) ~ "Core",
                in_go_rbp ~ "RBP",
                TRUE ~ "Other"
            )
        )
    plts[["regulons-n_targets_per_regulator_vs_sf_class-box"]] = x %>%
        mutate(sf_class = factor(sf_class, levels=SF_CLASS)) %>%
        ggplot(aes(x=sf_class, y=n)) +
        geom_quasirandom(size=0.1, color="orange", varwidth=0.5) +
        geom_boxplot(width=0.5, outlier.shape=NA, color="black", fill=NA) +
        geom_hline(yintercept = 25, linetype="dashed", color="black", size=LINE_SIZE) +
        theme_pubr(x.text.angle=45) +
        yscale("log10", .format=TRUE) +
        guides(fill="none") +
        geom_text(
            aes(y = 1, label=label), 
            . %>% 
            count(sf_class) %>% 
            mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Splicing Factor Subset", y="N. Targets per Regulator")
    
    return(plts)
}


plot_protein_impact = function(protimp_freqs){
    plts = list()
    
    X = protimp_freqs %>%
        drop_na()
    
    plts[["protein_impact-freqs-violin"]] = X %>%
        mutate(term_clean = fct_reorder(term_clean, n, median)) %>%
        ggplot(aes(x=term_clean, y=n)) +
        geom_violin(trim=TRUE, color=NA, fill="violet") +
        geom_boxplot(width=0.1, fill=NA, outlier.size=0.1) +
        theme_pubr(x.text.angle=70) +
        yscale("log10", .format=TRUE) +
        labs(x="Splicing Impact", y="Targets per Regulator")
    
    return(plts)
}


# plot_enrichments = function(enrichments){
#     plts = list()
    
#     X = enrichments
    
#     terms_oi = X %>%
#         count(Description) %>%
#         slice_max(n, n=10) %>%
#         pull(Description)
    
#     plts[["enrichments-"]] = enrichments %>%
#         filter(Description %in% terms_oi)
    
#     return(plts)
# }


plot_similarities = function(regulons_umap){
    plts = list()
    
    X = regulons_umap
    
    plts[["similarities-umap-scatter"]] = X %>%
        ggscatter(x="UMAP1", y="UMAP2", size=1) +
        theme(aspect.ratio=1)
    
    return(plts)
}


plot_target_lengths = function(regulons, annot){
    plts = list()
    
    X = regulons %>%
        distinct(GENE, target) %>%
        left_join(annot %>% distinct(EVENT, LE_o), by=c("target"="EVENT")) %>%
        group_by(GENE) %>%
        summarize(
            event_length = median(LE_o, na.rm=TRUE),
            n_targets = n()
        ) %>%
        ungroup() %>%
        arrange(event_length) %>%
        mutate(ranking = row_number()) %>%
        drop_na()
    
    plts[["target_lengths-regulators-scatter"]] = X %>%
        ggscatter(x="ranking", y="event_length", size="n_targets", color="brown", alpha=0.5) +
        geom_text_repel(
            aes(label=GENE),
            . %>% slice_min(event_length, n=5),
            size=FONT_SIZE, family=FONT_FAMILY, max.overlaps=50, segment.size=0.1
        ) +
        geom_text_repel(
            aes(label=GENE),
            . %>% slice_max(event_length, n=5),
            size=FONT_SIZE, family=FONT_FAMILY, max.overlaps=50, segment.size=0.1
        ) +
        scale_size(range=c(0.5,3)) +
        theme(aspect.ratio=1) +
        labs(x="Regulator", y="median(Target Length)", size="N. Targets")
    
    return(plts)
}


make_plots = function(regulons, protimp_freqs, regulons_umap, annot, evaluation, splicing_factors){
    plts = list(
        plot_regulons(regulons, evaluation, splicing_factors),
        plot_protein_impact(protimp_freqs),
        plot_similarities(regulons_umap),
        plot_target_lengths(regulons, annot)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(regulons, protimp_freqs, regulons_umap, annot){
    
    figdata = list(
        "eda" = list(
            "regulons" = regulons,
            "protein_impact_frequencies" = protimp_freqs,
            "regulons_umap" = regulons_umap
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
    save_plt(plts, "regulons-n_targets_per_regulator-box", '.pdf', figs_dir, width=2, height=5)   
    save_plt(plts, "regulons-n_regulators_per_target-box", '.pdf', figs_dir, width=2, height=5)
    save_plt(plts, "regulons-n_targets_per_regulator_vs_sf_class-box", '.pdf', figs_dir, width=3, height=4.5)
    save_plt(plts, "protein_impact-freqs-violin", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "similarities-umap-scatter", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "target_lengths-regulators-scatter", '.pdf', figs_dir, width=5, height=6)
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
        make_option("--regulons_dir", type="character"),
        make_option("--protein_impact_file", type="character"),
        make_option("--annotation_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    regulons_dir = args[["regulons_dir"]]
    protein_impact_file = args[["protein_impact_file"]]
    annotation_file = args[["annotation_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    regulons = lapply(list.files(regulons_dir, full.names=TRUE), function(regulons_file){
        regulon_id = basename(regulons_file) %>% gsub("-delta_psi.tsv.gz","",.)
        regulons = read_tsv(regulons_file) %>%
            mutate(regulon_id = regulon_id)
        return(regulons)
    }) %>% bind_rows()
    protein_impact = read_tsv(protein_impact_file) %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT) %>%
            mutate(term_clean=gsub(" \\(.*","",term),
                   term_clean=gsub("ORF disruption upon sequence exclusion",
                                   "ORF disruption (exclusion)",term_clean),
                   term_clean=gsub("ORF disruption upon sequence inclusion",
                                   "ORF disruption (inclusion)",term_clean),
                   term_clean=gsub("In the CDS, with uncertain impact",
                                   "In the CDS (uncertain)",term_clean))
    annot = read_tsv(annotation_file)
    splicing_factors = read_tsv(splicing_factors_file)
    evaluation = read_tsv(evaluation_ex_file)
    
    # prep
    protimp_freqs = regulons %>%
        distinct(regulator, target) %>%
        left_join(protein_impact, by=c("target"="EVENT")) %>%
        count(regulator, term_clean)
    
    evaluation = evaluation %>%
        mutate(regulon_id = gsub("-","_",regulon_id)) %>%
        filter(signature_id!=regulon_id) %>%
        filter(!(str_detect(regulon_id,"ENASFS") & (signature_id=="ENASFS"))) %>%
        # consider only signatures that we know activity 
        # of the splicing factor was altered
        filter(PERT_TYPE %in% c("KNOCKDOWN","KNOCKOUT","OVEREXPRESSION")) %>%
        mutate(
            pert_type_lab = case_when(
                PERT_TYPE=="KNOCKDOWN" ~ "KD",
                PERT_TYPE=="KNOCKOUT" ~ "KO",
                PERT_TYPE=="OVEREXPRESSION" ~ "OE"
            ),
            regulon_set = gsub("-.*","",regulon_set_id)
        ) %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, n_tails, regulon_set_id, pert_type_lab, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MAIN))
    
    # make embedding of common targets
    regulons_mat = regulons %>%
        distinct(regulator, target) %>%
        mutate(value=TRUE) %>%
        pivot_wider(id_cols="target", names_from="regulator", 
                    values_from="value", values_fill=FALSE) %>%
        column_to_rownames("target") %>%
        as.data.frame()
    
    regulons_simil = simil(regulons_mat, method="Jaccard", by_rows=FALSE) %>% as.matrix()
    regulons_dist = 1 - regulons_simil

    regulons_umap = umap(regulons_dist, input="dist", random_state=1234)[["layout"]]
    colnames(regulons_umap) = c("UMAP1", "UMAP2")
    regulons_umap = regulons_umap %>%
        as.data.frame() %>%
        rownames_to_column("regulator") 
    
    # plot
    plts = make_plots(regulons, protimp_freqs, regulons_umap, annot, evaluation, splicing_factors)
    
    # make figdata
    figdata = make_figdata(regulons, protimp_freqs, regulons_umap, annot)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
