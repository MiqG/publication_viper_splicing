#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - Validate that protein activity changes with different mutations in a splicing factor
#     - RBM10
#     - SF3B1
#     - SRSF1
#     - SRSF2
#     - ZRSR2
#     - U2AF1

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)
require(tidytext)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_tcga")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","PANCAN-PrimaryTumor-EX.tsv.gz")
# protein_activity_laml_file = file.path(RESULTS_DIR,"files","protein_activity","PANCAN-PrimaryBloodDerivedCancerPeripheralBlood-EX.tsv.gz")
# snv_file = file.path(RAW_DIR,"UCSCXena","TCGA","snv","mc3.v0.2.8.PUBLIC.xena.gz")
# metadata_file = file.path(RAW_DIR,"TCGA","metadata","PANCAN.tsv.gz")
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","tcga_sf_mutations")

##### FUNCTIONS #####
plot_mut_freq = function(mut_freq){
    plts = list()
    
    X = mut_freq
    
    # which genes and cancers are the top mutated in percentage?
    plts[["mut_freq-top_mutated_sfs-bar"]] = X %>%
        group_by(effect, sample_type, gene) %>%
        mutate(sum_perc = sum(perc_mutated)) %>%
        slice_max(sum_perc, n=10) %>%
        ungroup() %>%
        arrange(perc_mutated) %>%
        ggbarplot(x="gene", y="perc_mutated", fill="cancer_type", color=NA, palette=PALETTE) +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~sample_type+effect, ncol=1, scale="free_x") +
        labs(x="Gene", y="% Mutated Patients")
    
    return(plts)
}


plot_mut_effects = function(protein_activity, snv, gene_oi){
    plts = list()
    
    X = protein_activity %>%
        filter(GENE == gene_oi) %>%
        left_join(snv, by=c("sampleID","cancer_type","sample_type","GENE"="gene")) %>%
        mutate(
            is_mutated = !is.na(effect),
            effect_lab = replace_na(effect, "WT")
        ) %>%
        filter(sample_type == "Primary Tumor")
    
    # which cancers have mutations in that gene?
    plts[["mut_effects-n_mutated-bar"]] = X %>%
        count(GENE, sample_type, cancer_type, is_mutated, effect) %>%
        group_by(cancer_type, sample_type) %>%
        mutate(perc = n / sum(n) * 100) %>%
        ungroup() %>%
        filter(is_mutated) %>%
        mutate(
            cancer_type = factor(cancer_type),
            name = reorder_within(cancer_type, -perc, effect)
        ) %>%
        ggbarplot(x="name", y="perc", fill="cancer_type", color=NA, palette=PALETTE) +
        geom_text(
            aes(label=n), angle=45, vjust=0, hjust=0,
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        scale_x_reordered() +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~effect+sample_type, scale="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Cancer Type", y="% Patients with Mutation", subtitle=gene_oi)
    
    # is the protein activity different when having different types of mutations?
    plts[["mut_effects-activity_vs_mut_effect-boxplot"]] = X %>%
        filter(cancer_type=="BRCA") %>%
        ggplot(aes(x=effect_lab, y=activity)) +
        geom_jitter(aes(color=cancer_type), size=1, alpha=0.5) +
        geom_boxplot(width=0.5, outlier.shape=NA, fill=NA) +
        color_palette(PALETTE) +
        geom_text(
            aes(label=label, y=-4), 
            . %>% count(cancer_type, effect_lab) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle=70) +
        stat_compare_means(method="wilcox.test", label="p.signif", ref.group="WT", 
                           size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~sample_type+cancer_type, scale="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Mutation Effect", y="Protein Activity", subtitle=gene_oi)
    
    # how does the protein activity change having different mutations along the protein?
    plts[["mut_effects-activity_vs_mut_effect-scatter"]] = X %>%
        filter(cancer_type=="BRCA" & effect=="Missense_Mutation") %>%
        arrange(activity) %>%
        mutate(ranking = row_number()) %>%
        ggplot(aes(x=ranking, y=activity)) +
        geom_scattermore(aes(color=cancer_type), pixels=c(1000,1000), pointsize=5, alpha=0.5) +
        color_palette(PALETTE) +
        geom_text_repel(
            aes(label=Amino_Acid_Change), 
            . %>% slice_min(ranking, n=10),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        theme_pubr() +
        facet_wrap(~effect+cancer_type, scale="free") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Ranking", y="Protein Activity", subtitle=gene_oi)
}


make_plots = function(mut_freq){
    plts = list(
        plot_mut_freq(mut_freq)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(mut_freq){
    figdata = list(
        "tcga_sf_mutations" = list(
            "mut_freq" = mut_freq
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
    save_plt(plts, "diff_protein_activity-volcano-scatter", '.pdf', figs_dir, width=5, height=8)
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
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    diff_protein_activity_file = args[["protein_activity_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_activity = read_tsv(protein_activity_file) %>%
        bind_rows(read_tsv(protein_activity_laml_file))
    gene_annotation = read_tsv(gene_annotation_file) %>%
        dplyr::rename(
            GENE = `Approved symbol`,
            ENSEMBL = `Ensembl gene ID`
        )
    snv = read_tsv(snv_file) %>% dplyr::rename(sampleID=sample)
    metadata = read_tsv(metadata_file)
    
    # prep
    ## palette
    cancer_types = unique(metadata[["cancer_type"]])
    PALETTE = setNames(get_palette("Paired", length(cancer_types)), cancer_types)

    ## protein activity
    protein_activity = protein_activity %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("regulator"="ENSEMBL")
        ) %>%
        left_join(
            metadata %>% 
            distinct(sampleID, cancer_type, sample_type),
            by="sampleID"
        ) %>%
        drop_na(cancer_type)
    
    ## mutations
    genes_oi = protein_activity %>% pull(GENE) %>% unique()
    snv = snv %>%
        filter(gene %in% genes_oi) %>%
        left_join(
            metadata %>% 
            distinct(sampleID, cancer_type, sample_type),
            by="sampleID"
        ) %>%
        drop_na(cancer_type)
    
    mut_freq = snv %>%
        distinct(sampleID, gene, effect, cancer_type, sample_type) %>%
        count(gene, effect, cancer_type, sample_type) %>%
        left_join(
            metadata %>% 
            distinct(sampleID, cancer_type, sample_type) %>%
            group_by(cancer_type) %>%
            summarize(n_total=n()) %>%
            ungroup(),
            by="cancer_type"
        ) %>%
        mutate(
            perc_mutated = n/n_total
        )
        
    # plot
    plts = make_plots(mut_freq)
    
    # make figdata
    figdata = make_figdata(mut_freq)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}