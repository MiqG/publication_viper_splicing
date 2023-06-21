#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - T025 is a CLK inhibitor
#    - https://www.embopress.org/doi/full/10.15252/emmm.201708289
#    - suppresses phorphorylation of all 14 SR proteins, should be less active
# - PHF5A is a splicing factor acetylated, mutation K29Q reduces acetylation, and should reduce its activity

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)

# variables
RANDOM_SEED = 1234

# formatting
LINE_SIZE = 0.25

PAL_DARK = "darkgrey"

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

GENES_OI = c(
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

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_sf_ptms")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","sf_ptms.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","sf_ptms.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","sf_ptms.tsv.gz")
# annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_sf_ptms")

##### FUNCTIONS #####
plot_validation_phospho = function(genexpr, protein_activity, metadata, annot){
    plts = list()
    
    genes_oi = annot %>%
        filter(GENE %in% GENES_OI)
    studies_oi = c("PRJNA394754")
    conditions_oi = c("DMSO","T025")
    
    X = genes_oi %>%
        left_join(
            genexpr %>%
            filter(ID %in% genes_oi[["ENSEMBL"]]) %>%
            pivot_longer(-ID, names_to="sampleID", values_to="log2_tpm"),
            by=c("ENSEMBL"="ID")
        ) %>%
        left_join(
            protein_activity %>%
            filter(regulator %in% genes_oi[["ENSEMBL"]]) %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity"),
            by=c("ENSEMBL"="regulator","sampleID")
        ) %>%
        left_join(metadata, by="sampleID") %>%
        filter(study_accession%in%studies_oi & condition%in%conditions_oi) %>%
        mutate(GENE = factor(GENE, levels=GENES_OI)) %>%
        drop_na(activity, condition)
    
    
    plts[["validation_phospho-activity-line"]] = X %>%
        group_by(condition, pert_concentration, GENE) %>%
        summarize(
            mean = mean(activity),
            se = sd(activity),
            ymin = mean - se,
            ymax = mean + se
        ) %>%
        ungroup() %>%
        ggplot(aes(x=pert_concentration, y=mean, color=GENE, group=GENE)) +
        geom_line(size=0.1, linetype="dashed") +
        geom_point(size=0.1) +
        geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.01, size=LINE_SIZE) +
        color_palette(get_palette("npg",length(GENES_OI))) +
        theme_pubr() +
        geom_text(
            aes(label=GENE),
            . %>% filter(pert_concentration==max(as.numeric(pert_concentration))),
            size=FONT_SIZE, family=FONT_FAMILY, hjust=0, vjust=0.5, nudge_x=0.1
        ) +
        guides(color="none") +
        theme(aspect.ratio=1) +
        xscale("log10") +
        labs(x="CLK inhibitor T025 Conc. (nM)", y="Protein Activity")   
    
    
    plts[["validation_phospho-genexpr-line"]] = X %>%
        group_by(condition, pert_concentration, GENE) %>%
        summarize(
            mean = mean(log2_tpm),
            se = sd(log2_tpm),
            ymin = mean - se,
            ymax = mean + se
        ) %>%
        ungroup() %>%
        ggplot(aes(x=pert_concentration, y=mean, color=GENE, group=GENE)) +
        geom_line(size=0.1, linetype="dashed") +
        geom_point(size=0.1) +
        geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.01, size=LINE_SIZE) +
        color_palette(get_palette("npg",length(GENES_OI))) +
        theme_pubr() +
        geom_text(
            aes(label=GENE),
            . %>% filter(pert_concentration==max(as.numeric(pert_concentration))),
            size=FONT_SIZE, family=FONT_FAMILY, hjust=0, vjust=0.5, nudge_x=0.1
        ) +
        guides(color="none") +
        theme(aspect.ratio=1) +
        xscale("log10") +
        labs(x="CLK inhibitor T025 Conc. (nM)", y="log2(TPM+1)")        
    
    return(plts)
}


plot_validation_acetylation = function(protein_activity, metadata, annot){
    plts = list()
    
    # we still do not have PHF5A's network
    studies_oi = c("PRJNA506256")
    conditions_oi = c("CONTROL","PHF5A_K29Q_MUTATION")
    X = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        left_join(annot, by=c("regulator"="ENSEMBL")) %>% 
        filter(study_accession%in%studies_oi & condition%in%conditions_oi) %>%
        drop_na(condition) %>%
        group_by(condition, study_accession, cell_line_name, GENE) %>%
        summarize(activity = median(activity)) %>% # take median of activity
        ungroup() %>%
        group_by(condition,study_accession,cell_line_name) %>%
        arrange(activity) %>%
        mutate(activity_ranking = row_number()) %>%
        ungroup()
    
    plts[["validation_acetylation-mut_PHF5A-activity_ranking-scatter"]] = X %>%
            ggplot(aes(x=activity_ranking, y=activity)) +
            geom_scattermore(
                pixels=c(1000,1000), pointsize=4, color=PAL_DARK, alpha=0.5
            ) +
            geom_text_repel(
                aes(label=GENE),
                data= . %>% 
                    group_by(condition,study_accession,cell_line_name) %>% 
                    slice_max(abs(activity), n=10),
                size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
            )+
            color_palette("Dark2") +
            theme_pubr() +
            facet_wrap(~condition+study_accession+cell_line_name) +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(x="Ranking", y="median(Protein Activity)")
        
    return(plts)
}


make_plots = function(genexpr, protein_activity, metadata, annot){
    plts = list(
        plot_validation_phospho(genexpr, protein_activity, metadata, annot),
        plot_validation_acetylation(protein_activity, metadata, annot)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(genexpr, protein_activity, metadata){
    figdata = list(
        "validation_drug_target_activity" = list(
            "genexpr" = genexpr,
            "protein_activity" = protein_activity,
            "metadata" = metadata
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
    save_plt(plts, "validation_phospho-activity-line", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "validation_phospho-genexpr-line", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "validation_acetylation-mut_PHF5A-activity_ranking-scatter", '.pdf', figs_dir, width=10, height=7)
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
        make_option("--genexpr_file", type="character"),
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--annotation_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    genexpr_file = args[["genexpr_file"]]
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    annotation_file = args[["annotation_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    genexpr = read_tsv(genexpr_file)
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    annot = read_tsv(annotation_file)
    gc()
    
    # prep
    ## gene annotations
    annot = annot %>%
        rename(GENE=`Approved symbol`, ENSEMBL=`Ensembl gene ID`) %>%
        distinct(GENE, ENSEMBL) %>%
        drop_na()
    
    # plot
    plts = make_plots(genexpr, protein_activity, metadata, annot)
    gc()
    
    # make figdata
    figdata = make_figdata(genexpr, protein_activity, metadata)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}