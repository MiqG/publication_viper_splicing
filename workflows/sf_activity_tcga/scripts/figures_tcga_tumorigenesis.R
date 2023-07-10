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
THRESH_FDR = 0.05

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_CANCER_DRIVER = setNames(
    c("#6C98B3","#F6AE2D"),
    c("Tumor suppressor","Oncogenic")
)

PAL_CANCER_TYPES = setNames(
    c("#A6CEE3","#3385BB","#84BF96","#6DBD57","#7F9D55","#F57C7C","#E42622",
      "#FBB268","#FE8D19","#DE9E83","#9D7BBA","#977899","#F3E587","#B15928"),
    c("BLCA","BRCA","COAD","HNSC","KICH","KIRC","KIRP",
      "LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_tcga")
# diff_protein_activity_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","tcga_tumorigenesis")
# assocs_gene_dependency_file = file.path(RESULTS_DIR,"..","sf_activity_ccle","files","protein_activity_vs_demeter2","CCLE.tsv.gz")
# survival_analysis_surv_file = file.path(RESULTS_DIR,'files','PANCAN',"survival_analysis-surv.tsv.gz")
# survival_analysis_cat_file = file.path(RESULTS_DIR,'files','PANCAN',"survival_analysis-cat.tsv.gz")

##### FUNCTIONS #####
plot_diff_protein_activity = function(diff_protein_activity){
    plts = list()
    
    X = diff_protein_activity
    cancer_types = unique(X[["cancer_type"]])
    palette = setNames(get_palette("Paired", length(cancer_types)), cancer_types)
    
    plts[["diff_protein_activity-volcano-scatter"]] = X %>%
        ggplot(aes(x=`condition_a-median`, y=log10_padj)) +
        geom_scattermore(aes(color=cancer_type), pixels=c(1000,1000), pointsize=6, alpha=0.5) +
        color_palette(palette) + 
        geom_text_repel(
            aes(label=GENE),
            . %>% slice_max(abs(`condition_a-median`*log10_padj), n=10),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        theme_pubr() +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        labs(x="Protein Activity", 
             y="-log10(FDR)", color="Cancer Type")

    plts[["diff_protein_activity-n_signif-bar"]] = X %>%
        filter(is_significant) %>%
        count(cancer_type) %>%
        arrange(-n) %>%
        ggbarplot(x="cancer_type", y="n", fill="cancer_type", color=NA, palette=palette,
                  label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) +
        guides(fill="none") +
        theme_pubr(x.text.angle=70) +
        labs(x="Cancer Type", y="Count")    
   
   plts[["diff_protein_activity-n_signif_vs_driver_type-bar"]] = X %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup() %>%
        arrange(-n_sum) %>%
        ggbarplot(x="GENE", y="n_sign", fill="driver_type", color=NA) +
        geom_text(
            aes(label=GENE),
            . %>% slice_max(n_sum, n=5) %>% filter(driver_type=="Oncogenic"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            #segment.size=0.1, max.overlaps=50,
            angle=45, hjust=0, vjust=0, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE),
            . %>% slice_min(n_sum, n=5) %>% filter(driver_type=="Tumor suppressor"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            #segment.size=0.1, max.overlaps=50,
            angle=45, hjust=1, vjust=0, nudge_y=-0.25
        ) +
        fill_palette(PAL_CANCER_DRIVER) +
        theme_pubr(x.text.angle=70) +
        labs(x="Splicing Factor", y="Count", fill="Driver Type")
    
    top_doubleagent = X %>%
        filter(is_significant) %>%
        group_by(GENE) %>%
        summarize(std = sd(`condition_a-median`)) %>%
        ungroup() %>%
        slice_max(std, n=10) %>%
        pull(GENE)
    plts[["diff_protein_activity-top_doubleagent-box"]] = X %>%
        filter(GENE %in% top_doubleagent) %>%
        group_by(GENE) %>%
        mutate(std = sd(`condition_a-median`)) %>%
        ungroup() %>%
        mutate(GENE = fct_reorder(GENE, std, max)) %>%
        ggplot(aes(x=GENE, y=`condition_a-median`)) +
        geom_point(aes(color=cancer_type), size=1) +
        geom_boxplot(outlier.shape=NA, fill=NA, width=0.5) +
        color_palette(palette) +
        theme_pubr(x.text.angle=70) +
        labs(x="Splicing Factor", y="Protein Activity", color="Cancer Type")    
        
    top_oncogenes = X %>%
        filter(is_significant & !(GENE%in%top_doubleagent)) %>%
        group_by(GENE) %>%
        summarize(med = median(`condition_a-median`)) %>%
        ungroup() %>%
        slice_max(med, n=10) %>%
        pull(GENE)
    plts[["diff_protein_activity-top_oncogenes-box"]] = X %>%
        filter(GENE %in% top_oncogenes) %>%
        group_by(GENE) %>%
        mutate(med = median(`condition_a-median`)) %>%
        ungroup() %>%
        mutate(GENE = fct_reorder(GENE, med, max)) %>%
        ggplot(aes(x=GENE, y=`condition_a-median`)) +
        geom_point(aes(color=cancer_type), size=1) +
        geom_boxplot(outlier.shape=NA, fill=NA, width=0.5) +
        color_palette(palette) +
        theme_pubr(x.text.angle=70) +
        labs(x="Splicing Factor", y="Protein Activity", color="Cancer Type")    
   
    top_tsgenes = X %>%
        filter(is_significant & !(GENE%in%top_doubleagent)) %>%
        group_by(GENE) %>%
        summarize(med = median(`condition_a-median`)) %>%
        ungroup() %>%
        slice_min(med, n=10) %>%
        pull(GENE)
    plts[["diff_protein_activity-top_tsgenes-box"]] = X %>%
        filter(GENE %in% top_tsgenes) %>%
        group_by(GENE) %>%
        mutate(med = median(`condition_a-median`)) %>%
        ungroup() %>%
        mutate(GENE = fct_reorder(GENE, -med, max)) %>%
        ggplot(aes(x=GENE, y=`condition_a-median`)) +
        geom_point(aes(color=cancer_type), size=1) +
        geom_boxplot(outlier.shape=NA, fill=NA, width=0.5) +
        color_palette(palette) +
        theme_pubr(x.text.angle=70) +
        labs(x="Splicing Factor", y="Protein Activity", color="Cancer Type")
    
    return(plts)
}


plot_prolif_driver = function(diff_protein_activity, assocs_gene_dependency){
    plts = list()
    
    X = diff_protein_activity %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    
    # Do inferred tumor suppressor and oncogenic activities match gene dependencies?
    x = X %>%
        filter(abs(n_sum)>5) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>% 
        left_join(
            assocs_gene_dependency,
            by="GENE"
        )
    
    plts[["prolif_driver-driver_type_vs_demeter2-violin"]] = x %>%
        ggviolin(x="driver_type", y="pearson_coef", color=NA, fill="driver_type", palette=PAL_CANCER_DRIVER, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_text_repel(
            aes(label=GENE),
            x %>% slice_max(pearson_coef, n=5) %>% bind_rows(x %>% slice_min(pearson_coef, n=5)),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1
        ) +
        geom_text(
            aes(y=-0.3, label=label),
            . %>% count(driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="Driver Type", y="Demeter2 Pearson Correlation")
    
    return(plts)
}


plot_survival_analysis = function(diff_protein_activity, survival_analysis_surv, survival_analysis_cat){
    plts = list()
    
    driver_classif = diff_protein_activity %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, regulator, driver_type) %>%
        group_by(GENE, regulator) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        filter(abs(n_sum)>5) %>%
        group_by(GENE, regulator) %>%
        slice_max(n, n=1) %>%
        ungroup()

    # We expect, oncogenic SFs to have positive CoxPH coefficients 
    # (increased hazard with increased protein activity)
    X = driver_classif %>%
        left_join(
            survival_analysis_surv,
            by="regulator"
        ) %>%
        mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")) %>%
        count(GENE, surv_type, driver_type)
    
    plts[["survival_analysis-surv_type_vs_driver_type-violin"]] = X %>%
        ggplot(aes(x=surv_type, y=n, group=interaction(surv_type,driver_type))) +
        geom_violin(aes(fill=driver_type), trim=TRUE, color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, position=position_dodge(0.9)) +
        fill_palette(PAL_CANCER_DRIVER) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY, label="p.signif") +
        geom_text(
            aes(y=3, label=label),
            . %>% count(driver_type, surv_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Risk with High Protein Activity", y="N. Cancers")
        
    return(plts)
}


make_plots = function(diff_protein_activity, assocs_gene_dependency, survival_analysis_surv, survival_analysis_cat){
    plts = list(
        plot_diff_protein_activity(diff_protein_activity),
        plot_prolif_driver(diff_protein_activity, assocs_gene_dependency),
        plot_survival_analysis(diff_protein_activity, survival_analysis_surv, survival_analysis_cat)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_protein_activity, assocs_gene_dependency, survival_analysis_surv, survival_analysis_cat){
    figdata = list(
        "tcga_tumorigenesis" = list(
            "diff_protein_activity" = diff_protein_activity,
            "assocs_gene_dependency" = assocs_gene_dependency,
            "survival_analysis_surv" = survival_analysis_surv,
            "survival_analysis_cat" = survival_analysis_cat
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
    save_plt(plts, "diff_protein_activity-n_signif-bar", '.pdf', figs_dir, width=6, height=5)
    save_plt(plts, "diff_protein_activity-n_signif_vs_driver_type-bar", '.pdf', figs_dir, width=12, height=9)
    save_plt(plts, "diff_protein_activity-top_doubleagent-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "diff_protein_activity-top_oncogenes-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "diff_protein_activity-top_tsgenes-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "diff_protein_activity-top_tsgenes-box", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "prolif_driver-driver_type_vs_demeter2-violin", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "survival_analysis-surv_type_vs_driver_type-violin", '.pdf', figs_dir, width=5, height=7)
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
        make_option("--diff_protein_activity_file", type="character"),
        make_option("--gene_annotation_file", type="character"),
        make_option("--assocs_gene_dependency_file", type="character"),
        make_option("--survival_analysis_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    diff_protein_activity_file = args[["diff_protein_activity_file"]]
    gene_annotation_file = args[["gene_annotation_file"]]
    assocs_gene_dependency_file = args[["assocs_gene_dependency_file"]]
    survival_analysis_file = args[["survival_analysis_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_protein_activity = read_tsv(diff_protein_activity_file)
    gene_annotation = read_tsv(gene_annotation_file) %>%
        dplyr::rename(
            GENE = `Approved symbol`,
            ENSEMBL = `Ensembl gene ID`
        )
    assocs_gene_dependency = read_tsv(assocs_gene_dependency_file)
    survival_analysis_surv = read_tsv(survival_analysis_surv_file)
    survival_analysis_cat = read_tsv(survival_analysis_cat_file)
    
    # prep
    diff_protein_activity = diff_protein_activity %>%
        mutate(
            is_significant = padj < THRESH_FDR
        ) %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("regulator"="ENSEMBL")
        )
    
    # plot
    plts = make_plots(diff_protein_activity, assocs_gene_dependency, survival_analysis_surv, survival_analysis_cat)
    
    # make figdata
    figdata = make_figdata(diff_protein_activity, assocs_gene_dependency, survival_analysis_surv, survival_analysis_cat)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}