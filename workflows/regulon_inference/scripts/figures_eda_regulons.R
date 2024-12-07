require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(ggrepel)
require(proxy)
require(umap)
require(ggbeeswarm)
require(ggvenn)
require(purrr)

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
# regulons_clip_dir = file.path(RESULTS_DIR,"files","postar3_clip_regulons-EX")
# enrichments_file = file.path(RESULTS_DIR,"files","regulons_eda_gsea","experimentally_derived_regulons_pruned-EX.tsv.gz")
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# spliceosomedb_file = file.path(SUPPORT_DIR,"splicing_factors","literature_suptabs","SpliceosomeDB-human.csv")
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,'figures','eda_regulons-EX')

##### FUNCTIONS #####
plot_regulons = function(regulons, splicing_factors){
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
        geom_text(
            aes(y = 1, label=label), 
            . %>% 
            count(regulon_id) %>% 
            mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon ID", y="N. Targets per Regulator")
    
    plts[["regulons-n_regulators_per_target-box"]] = X %>%
        count(regulon_id, target) %>%
        ggplot(aes(x=regulon_id, y=n)) +
        geom_quasirandom(size=0.1, color="orange", varwidth=0.5) +
        geom_boxplot(width=0.5, outlier.shape=NA, color="black", fill=NA) +
        guides(fill="none") +
        theme_pubr(x.text.angle=45) +
        geom_text(
            aes(y = -1, label=label), 
            . %>% 
            count(regulon_id) %>% 
            mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon ID", y="N. Regulators per Target")
    
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


plot_similarities = function(regulons_umap, network_simil){
    plts = list()
    
    X = regulons_umap
    
    plts[["similarities-umap-scatter"]] = X %>%
        ggscatter(x="UMAP1", y="UMAP2", size=1) +
        theme(aspect.ratio=1)
    
    plts[["similarities-jaccard-box"]] = network_simil %>%
        mutate(regulon_id = "Empirical") %>%
        ggplot(aes(x=regulon_id, y=jaccard)) +
        geom_quasirandom(size=0.1, color="orange", varwidth=0.5) +
        geom_boxplot(width=0.5, outlier.shape=NA, color="black", fill=NA) +
        guides(fill="none") +
        theme_pubr(x.text.angle=45) +
        geom_text(
            aes(y = -0.01, label=label), 
            . %>% 
            count(regulon_id) %>% 
            mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon ID", y="Jaccard Similarity")
    
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


plot_emp_vs_clip_regulons = function(emp_vs_clip_interactions, emp_vs_clip_regulators){
    
    plts = list()
    
    plts[["emp_vs_clip-interactions-venn"]] = emp_vs_clip_interactions %>%
        ggvenn(
            c("emp","clip"),
            stroke_color=NA, 
            set_name_size = FONT_SIZE+0.5, text_size = FONT_SIZE) +
        coord_fixed() +
        theme_void()
    
    plts[["emp_vs_clip-regulators-venn"]] = emp_vs_clip_regulators %>%
        ggvenn(
            c("emp","clip"),
            stroke_color=NA, 
            set_name_size = FONT_SIZE+0.5, text_size = FONT_SIZE) +
        coord_fixed() +
        theme_void()
    
    return(plts)
}


plot_clip_eda = function(clip_distances){
    plts = list()
    
    X = clip_distances
    
    plts[["clip_eda-distance_distr-hist"]] = X %>%
        gghistogram(x="rel_distance", fill="rel_position", color=NA, bins=50) +
        facet_wrap(~rel_position, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Distance between CLIP peak to closest exon splice site ", y="N. CLIP Peaks")
    
    return(plts)
}


make_plots = function(regulons, protimp_freqs, regulons_umap, annot, splicing_factors, network_simil,
                      emp_vs_clip_interactions, emp_vs_clip_regulators, clip_distances){
    plts = list(
        plot_regulons(regulons, splicing_factors),
        plot_protein_impact(protimp_freqs),
        plot_similarities(regulons_umap, network_simil),
        plot_target_lengths(regulons, annot),
        plot_emp_vs_clip_regulons(emp_vs_clip_interactions, emp_vs_clip_regulators),
        plot_clip_eda(clip_distances)
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
    save_plt(plts, "similarities-jaccard-box", '.pdf', figs_dir, width=2, height=5)
    save_plt(plts, "target_lengths-regulators-scatter", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "emp_vs_clip-interactions-venn", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "emp_vs_clip-regulators-venn", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "clip_eda-distance_distr-hist", '.pdf', figs_dir, width=7, height=5)
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
        make_option("--regulons_clip_dir", type="character"),
        make_option("--protein_impact_file", type="character"),
        make_option("--annotation_file", type="character"),
        make_option("--splicing_factors_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    regulons_dir = args[["regulons_dir"]]
    regulons_clip_dir = args[["regulons_clip_dir"]]
    protein_impact_file = args[["protein_impact_file"]]
    annotation_file = args[["annotation_file"]]
    splicing_factors_file = args[["splicing_factors_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    ## empirical regulons
    regulons = lapply(list.files(regulons_dir, full.names=TRUE), function(regulons_file){
        regulon_id = basename(regulons_file) %>% gsub("-delta_psi.tsv.gz","",.)
        regulons = read_tsv(regulons_file) %>%
            mutate(regulon_id = regulon_id)
        return(regulons)
    }) %>% bind_rows()
    ## CLIP regulons from POSTAR3
    regulons_clip = lapply(list.files(regulons_clip_dir, full.names=TRUE), function(regulons_file){
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
    
    # prep
    protimp_freqs = regulons %>%
        distinct(regulator, target) %>%
        left_join(protein_impact, by=c("target"="EVENT")) %>%
        count(regulator, term_clean)
    
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
    
    upper_triangle_indices = which(upper.tri(regulons_simil), arr.ind = TRUE)
    network_simil = data.frame(
        regulator_a = rownames(regulons_simil)[upper_triangle_indices[,"row"]],
        regulator_b = colnames(regulons_simil)[upper_triangle_indices[,"col"]],
        jaccard = regulons_simil[upper_triangle_indices]
    )
    
    # CLIP peak distances to closest splice site
    clip_distances = regulons_clip %>%
        distinct(Start, End, peak_start, peak_end) %>%
        mutate(
            rel_position = ifelse(peak_start > End | peak_end < Start, "On Intron", "On Exon"),
            rel_distance = pmin(
              abs(peak_start - Start),
              abs(peak_end - End),
              abs(peak_end - Start),
              abs(peak_start - End)
            )
        )
    
    # empirical vs clip interactions
    interactions_emp = regulons %>% distinct(regulator, target) %>% mutate(edge = paste0(regulator,"_",target))
    interactions_clip = regulons_clip %>% distinct(regulator, target) %>% mutate(edge = paste0(regulator,"_",target))
    ## overlap between interactions?
    emp_vs_clip_interactions = list(
        emp = interactions_emp[["edge"]],
        clip = interactions_clip[["edge"]]
    )
    ## overlap between regulators?
    emp_vs_clip_regulators = list(
        emp = interactions_emp %>% distinct(regulator) %>% pull(regulator),
        clip = interactions_clip %>% distinct(regulator) %>% pull(regulator)
    )
    
    # plot
    plts = make_plots(regulons, protimp_freqs, regulons_umap, annot, splicing_factors, network_simil,
                      emp_vs_clip_interactions, emp_vs_clip_regulators, clip_distances)
    
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
