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

# protein_activity_tcga_stn_file = file.path(RESULTS_DIR,"files","protein_activity","PANCAN-SolidTissueNormal-EX.tsv.gz")
# genexpr_tcga_stn_file = file.path(PREP_DIR,"genexpr_tpm","PANCAN-SolidTissueNormal.tsv.gz")
# enrichment_scores_tcga_stn_file = file.path(RESULTS_DIR,"files","program_enrichment_scores","PANCAN-SolidTissueNormal-EX.tsv.gz")
# metadata_tcga_file = file.path(PREP_DIR,"metadata","PANCAN.tsv.gz")
# driver_types_file = file.path("driver_types.tsv")

# protein_activity_ccle_file = file.path(RESULTS_DIR,"files","protein_activity","CCLE-EX.tsv.gz")
# genexpr_ccle_file = file.path(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
# metadata_ccle_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# doublings_ccle_file = file.path(PREP_DIR,"doubling_times","CCLE.tsv.gz")
# enrichment_scores_ccle_file = file.path(RESULTS_DIR,"files","program_enrichment_scores","CCLE-EX.tsv.gz")
# metmap_file = file.path(PREP_DIR,"metmap","CCLE.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,"figures","cancer_hallmarks-EX")

##### FUNCTIONS #####
plot_proliferation_stn = function(protein_activity_tcga_stn){
    plts = list()
    
    X = protein_activity_tcga_stn
    
    # median activity of oncogenic and tumor suppressor SFs per sample across cancers
    cancer_order = X %>%
        distinct(cancer_type, sampleID, MKI67) %>%
        group_by(cancer_type) %>%
        summarize(MKI67 = median(MKI67, na.rm=TRUE)) %>%
        ungroup() %>%
        arrange(MKI67) %>%
        pull(cancer_type)
    
    plts[["proliferation_stn-cancer_vs_MKI67-line"]] = X %>%
        distinct(sampleID, MKI67, cancer_type) %>%
        group_by(cancer_type) %>% 
        summarize(
            MKI67_median = median(MKI67, na.rm=TRUE),
            MKI67_q25 = quantile(MKI67, 0.25),
            MKI67_q75 = quantile(MKI67, 0.75),
        ) %>% 
        ungroup() %>%
        mutate(cancer_type = factor(cancer_type, levels=cancer_order)) %>%
        ggplot(aes(x=cancer_type, y=MKI67_median, group=1)) +
        geom_point(color="black", size=1) + 
        geom_line(color="black", size=LINE_SIZE, linetype="dashed") + 
        theme_pubr(x.text.angle = 70) +
        labs(x="STN Tissue Type", y="log2(TPM+1) of MKI67")

    
    plts[["proliferation_stn-cancer_vs_program_activity-line"]] = X %>%
        distinct(sampleID, activity, cancer_type, driver_type, GENE) %>%
        group_by(cancer_type, driver_type, GENE) %>%
        summarize(
            activity = median(activity, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        group_by(cancer_type, driver_type) %>% 
        summarize(
            activity_median = median(activity, na.rm=TRUE),
            activity_q25 = quantile(activity, 0.25),
            activity_q75 = quantile(activity, 0.75),
        ) %>% 
        ungroup() %>%
        mutate(cancer_type = factor(cancer_type, levels=cancer_order)) %>%
        drop_na(driver_type) %>%
        ggplot(aes(x=cancer_type, y=activity_median, group=driver_type)) +
        geom_line(
            aes(color=driver_type), 
            position=position_dodge(0.1),
            size=LINE_SIZE, linetype="dashed") +
        geom_point(
            aes(y=activity_median, color=driver_type), 
            size=1, position=position_dodge(0.1)) + 
        color_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        labs(x="STN Tissue Type", y="median(Protein Activity SFs)", color="Driver Type")
    
    
    plts[["proliferation_stn-cancer_vs_program_enrichment-line"]] = X %>%
        distinct(sampleID, term, enrichment_score, cancer_type) %>%
        group_by(cancer_type, term) %>% 
        summarize( # across samples in each cancer
            enrichment_median = median(enrichment_score, na.rm=TRUE),
            enrichment_q25 = quantile(enrichment_score, 0.25),
            enrichment_q75 = quantile(enrichment_score, 0.75),
        ) %>% 
        ungroup() %>%
        mutate(cancer_type = factor(cancer_type, levels=cancer_order)) %>%
        ggplot(aes(x=cancer_type, y=enrichment_median, group=term)) +
        geom_line(
            aes(color=term), 
            position=position_dodge(0.1),
            size=LINE_SIZE, linetype="dashed") +
        geom_point(
            aes(color=term),
            size=1, position=position_dodge(0.1)) + 
        color_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle = 70) +
        labs(x="STN Tissue Type", y="median(Program Enrichment Score)", color="Driver Type")
        
        
    plts[["proliferation_stn-mki67_vs_activity_median-scatter"]] = X %>%
        distinct(sampleID, activity, cancer_type, driver_type, GENE, MKI67) %>%
        group_by(sampleID, cancer_type, driver_type, MKI67) %>%
        summarize(
            activity_median = median(activity, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        ggscatter(x="MKI67", y="activity_median", color="driver_type", size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) +
        geom_smooth(color="black", method="lm", size=LINE_SIZE, linetype="dashed") +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~driver_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log2(TPM+1) MKI67", y="median(Protein Activity SFs)")
    
    
    plts[["proliferation_stn-corrs_mki67_vs_activity_median-violin"]] = X %>%
        distinct(sampleID, activity, cancer_type, driver_type, GENE, MKI67) %>%
        group_by(sampleID, cancer_type, driver_type, MKI67) %>%
        summarize(
            activity_median = median(activity, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        group_by(cancer_type, driver_type) %>%
        summarize(correlation = cor(activity_median, MKI67, method="spearman")) %>%
        ungroup() %>%
        ggviolin(x="driver_type", y="correlation", fill="driver_type", palette=PAL_DRIVER_TYPE, color=NA, trim=TRUE) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1) +
        geom_text_repel(
            aes(label=cancer_type), 
            segment.size=0.1, max.overlaps=50,
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="Driver Type", y="Spearman Correlation MKI67 vs median(Protein Activity SFs)")
    
    
    plts[["proliferation_stn-mki67_vs_program_enrichment-scatter"]] = X %>%
        distinct(sampleID, MKI67, enrichment_score, term) %>%
        ggscatter(x="MKI67", y="enrichment_score", color="term", size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) +
        geom_smooth(color="black", method="lm", size=LINE_SIZE, linetype="dashed") +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~term) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log2(TPM+1) MKI67", y="Program Enrichment Score")
    
    
    plts[["proliferation_stn-corrs_mki67_vs_program_enrichment-violin"]] = X %>%
        distinct(sampleID, MKI67, enrichment_score, term, cancer_type) %>%
        group_by(cancer_type, term) %>%
        summarize(correlation = cor(enrichment_score, MKI67, method="spearman")) %>%
        ungroup() %>%
        ggviolin(x="term", y="correlation", fill="term", palette=PAL_DRIVER_TYPE, color=NA, trim=TRUE) +
        geom_boxplot(fill=NA, width=0.1, outlier.size=0.1) +
        geom_text_repel(
            aes(label=cancer_type), 
            segment.size=0.1, max.overlaps=50,
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="Driver Type", y="Spearman Correlation MKI67 vs median(Protein Activity SFs)")
    
    
    return(plts)
}


plot_hallmarks_ccle = function(protein_activity_ccle){
    plts = list()
    
    X = protein_activity_ccle
    
    # Does MKI67 correspond to doubling time?
    plts[["hallmarks_ccle-doublings_vs_mki67-scatter"]] = X %>%
        distinct(DepMap_ID, MKI67, doubling_time_ghandi, doubling_time_tsherniak, doubling_time) %>%
        pivot_longer(-c(DepMap_ID, MKI67), values_to="doubling_time", names_to="dataset") %>%
        drop_na() %>%
        ggscatter(x="doubling_time", y="MKI67", size=1, alpha=0.5, color=PAL_DARK) +
        geom_smooth(method="lm", size=LINE_SIZE, color="black", linetype="dashed") +
        facet_wrap(~dataset, ncol=3) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        xscale("log10", .format=TRUE) +
        labs(x="Doubling Time (h)", y="log2(TPM+1) MKI67")
    
    # How do cancer-driver activities relate to proliferation?
    plts[["hallmarks_ccle-mki67_vs_enrichment-scatter"]] = X %>% 
        distinct(DepMap_ID, MKI67, enrichment_score, term) %>%
        ggscatter(x="MKI67", y="enrichment_score", color="term", 
                  size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) + 
        geom_smooth(method="lm", size=LINE_SIZE, color="black", linetype="dashed") +
        stat_cor(method="pearson", size=FONT_SIZE+1, family=FONT_FAMILY) + 
        facet_wrap(~term) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log2(TPM+1) MKI67", y="Program Enrichment Score")
    
    
    plts[["hallmarks_ccle-mki67_vs_activity_median-scatter"]] = X %>% 
        distinct(DepMap_ID, MKI67, activity, driver_type, GENE) %>%
        group_by(DepMap_ID, driver_type, MKI67) %>%
        summarize(activity_median = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        ggscatter(x="MKI67", y="activity_median", color="driver_type", 
                  size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) + 
        geom_smooth(method="lm", size=LINE_SIZE, color="black", linetype="dashed") +
        stat_cor(method="pearson", size=FONT_SIZE+1, family=FONT_FAMILY) + 
        facet_wrap(~driver_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="log2(TPM+1) MKI67", y="median(Protein Activity SFs)")
    
    # How do cancer-driver activities relate to metastasis?
    plts[["hallmarks_ccle-metpotential_vs_enrichment-scatter"]] = X %>% 
        distinct(DepMap_ID, metastatic_tissue, metastatic_potential, enrichment_score, term) %>%
        drop_na(term, metastatic_tissue) %>%
        ggscatter(x="metastatic_potential", y="enrichment_score", color="term", 
                  size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) + 
        geom_smooth(method="lm", size=LINE_SIZE, color="black", linetype="dashed") +
        stat_cor(method="pearson", size=FONT_SIZE+1, family=FONT_FAMILY) + 
        facet_wrap(~term+metastatic_tissue) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Metastatic Potential", y="Program Enrichment Score")
    
    plts[["hallmarks_ccle-metpotential_vs_activity_median-scatter"]] = X %>% 
        distinct(DepMap_ID, metastatic_tissue, metastatic_potential, activity, driver_type, GENE) %>%
        group_by(DepMap_ID, driver_type, metastatic_tissue, metastatic_potential) %>%
        summarize(activity_median = median(activity, na.rm=TRUE)) %>%
        ungroup() %>%
        drop_na(driver_type, metastatic_tissue) %>%
        ggscatter(x="metastatic_potential", y="activity_median", color="driver_type", 
                  size=1, alpha=0.5, palette=PAL_DRIVER_TYPE) + 
        geom_smooth(method="lm", size=LINE_SIZE, color="black", linetype="dashed") +
        stat_cor(method="pearson", size=FONT_SIZE+1, family=FONT_FAMILY) + 
        facet_wrap(~driver_type+metastatic_tissue) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="Metastatic Potential", y="median(Protein Activity SFs)")
    
    return(plts)
}


make_plots = function(protein_activity_tcga_stn, protein_activity_ccle){
    plts = list(
        plot_proliferation_stn(protein_activity_tcga_stn),
        plot_hallmarks_ccle(protein_activity_ccle)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(protein_activity_tcga_stn, protein_activity_ccle){
    figdata = list(
        "cancer_hallmarks" = list(
            "protein_activity_tcga_stn" = protein_activity_tcga_stn,
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
    save_plt(plts, "drivers_activity_stn-cancer_vs_MKI67-line", '.pdf', figs_dir, width=6, height=3)
    save_plt(plts, "drivers_activity_stn-cancer_vs_oncogenic-line", '.pdf', figs_dir, width=6, height=6)
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
    protein_activity_tcga_stn = read_tsv(protein_activity_tcga_stn_file)
    genexpr_tcga_stn = read_tsv(genexpr_tcga_stn_file)
    enrichment_scores_tcga_stn = read_tsv(enrichment_scores_tcga_stn_file)
    metadata_tcga = read_tsv(metadata_tcga_file)
    driver_types = read_tsv(driver_types_file)
    
    protein_activity_ccle = read_tsv(protein_activity_ccle_file)
    genexpr_ccle = read_tsv(genexpr_ccle_file)
    enrichment_scores_ccle = read_tsv(enrichment_scores_ccle_file)
    metadata_ccle = read_tsv(metadata_ccle_file)
    doublings_ccle = read_tsv(doublings_ccle_file)
    metmap = read_tsv(metmap_file)
    
    # prep
    ## TCGA STN
    genexpr_mki67_tcga = genexpr_tcga_stn %>%
        filter(ID == "ENSG00000148773") %>%
        pivot_longer(-ID, names_to="sampleID", values_to="MKI67")
    
    enrichment_scores_tcga_stn = enrichment_scores_tcga_stn %>%
        pivot_longer(-term, names_to="sampleID", values_to="enrichment_score") %>%
        left_join(
            genexpr_mki67_tcga,
            by="sampleID"
        ) %>%
        left_join(
            metadata_tcga,
            by="sampleID"
        )
    
    protein_activity_tcga_stn = protein_activity_tcga_stn %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(
            enrichment_scores_tcga_stn,
            by="sampleID",
            relationship="many-to-many" # each sample, 2 enrichments
        ) %>%
        left_join(
            driver_types,
            by=c("regulator"="ENSEMBL")
        )
    
    ## CCLE
    genexpr_mki67_ccle = genexpr_ccle %>%
        filter(ID == "ENSG00000148773") %>%
        pivot_longer(-ID, names_to="DepMap_ID", values_to="MKI67")
    
    metmap = metmap %>%
        dplyr::select(DepMap_ID, mean_all5, mean_brain, mean_bone, 
                      mean_lung, mean_liver, mean_kidney) %>%
        pivot_longer(-DepMap_ID, names_to="metastatic_tissue", values_to="metastatic_potential")
    
    enrichment_scores_ccle = enrichment_scores_ccle %>%
        pivot_longer(-term, names_to="DepMap_ID", values_to="enrichment_score") %>%
        left_join(
            genexpr_mki67_ccle,
            by="DepMap_ID"
        ) %>%
        left_join(
            metadata_ccle,
            by="DepMap_ID"
        ) %>%
        left_join(
            doublings_ccle,
            by="CCLE_Name"
        ) %>%
        left_join(
            metmap,
            by="DepMap_ID",
            relationship="many-to-many" # each sample, 6 metastatic scores
        )
    
    protein_activity_ccle = protein_activity_ccle %>%
            pivot_longer(-regulator, names_to="DepMap_ID", values_to="activity") %>%
        left_join(
            enrichment_scores_ccle,
            by="DepMap_ID",
            relationship="many-to-many" # each sample, 2 enrichments
        ) %>%
        left_join(
            driver_types,
            by=c("regulator"="ENSEMBL")
        )
    
    # plot
    plts = make_plots(
        protein_activity_tcga_stn,
        protein_activity_ccle
    )
    
    # make figdata
    figdata = make_figdata(
        protein_activity_tcga_stn,
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