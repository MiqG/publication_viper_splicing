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
require(clusterProfiler)
require(org.Hs.eg.db)
library(pROC)

# variables
THRESH_FDR = 0.05
ORGDB = org.Hs.eg.db
THRESH_N_SUM = 5.5

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DRIVER_TYPE = setNames(
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
# sf_cross_regulation_file = file.path(RESULTS_DIR,'files','PANCAN',"sf_cross_regulation.tsv.gz")
# ontology_chea_file = file.path(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz")

##### FUNCTIONS #####
plot_diff_protein_activity = function(diff_protein_activity){
    plts = list()
    
    X = diff_protein_activity

    x = X %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, driver_type) %>%
        group_by(GENE) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        ungroup()
    sf_order = x %>%
        pivot_wider(id_cols="GENE", names_from="driver_type", values_from="n_sign", values_fill=0) %>%
        rowwise() %>%
        mutate(diff=sum(`Tumor suppressor`,`Oncogenic`)) %>%
        ungroup() %>%
        arrange(diff,`Tumor suppressor`,`Oncogenic`) %>%
        pull(GENE)
    
   plts[["diff_protein_activity-n_signif_vs_driver_type-bar"]] = x %>%
        mutate(GENE = factor(GENE, levels=sf_order)) %>%
        ggbarplot(x="GENE", y="n_sign", fill="driver_type", color=NA) +
        geom_text(
            aes(label=GENE),
            . %>% slice_max(n_sum, n=5) %>% filter(driver_type=="Oncogenic"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=1, vjust=1, nudge_y=0.25
        ) +
        geom_text(
            aes(label=GENE),
            . %>% slice_min(n_sum, n=5) %>% filter(driver_type=="Tumor suppressor"),
            size=FONT_SIZE, family=FONT_FAMILY, 
            angle=-45, hjust=0, vjust=1, nudge_y=-0.25
        ) +
        geom_hline(yintercept=c(-5,5), linetype="dashed", color="black", size=LINE_SIZE) +
        geom_text(
            aes(x=x, label=label),
            . %>% 
            filter(abs(n_sum)>THRESH_N_SUM) %>%
            group_by(GENE) %>%
            slice_max(n, n=1) %>%
            ungroup() %>%
            count(driver_type) %>% 
                mutate(
                    label=paste0("n=",n),
                    n_sign=c(10,-10),
                    x=c(120,40)
                ),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        fill_palette(PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=70) +
        labs(x="Splicing Factor", y="Count", fill="Driver Type")
    
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
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE) %>%
        slice_max(n, n=1) %>%
        ungroup() %>% 
        left_join(
            assocs_gene_dependency,
            by="GENE"
        )
    
    plts[["prolif_driver-driver_type_vs_demeter2-violin"]] = x %>%
        ggviolin(x="driver_type", y="pearson_coef", color=NA, fill="driver_type", palette=PAL_DRIVER_TYPE, trim=TRUE) +
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


make_roc_analysis = function(driver_classif, survival_analysis_surv){
    # roc curve based on thresholds of n_sum for each type of cancer
    # We expect to classify Tumor Suppressor as Low Risk and Oncogenic as High Risk
    driver_types = c("Tumor suppressor", "Oncogenic")
    driver_levels = list(
        "Tumor suppressor" = c("High Risk", "Low Risk"),
        "Oncogenic" = c("Low Risk", "High Risk")
    )
    
    driver_classif = driver_classif %>%
        count(GENE, regulator, driver_type) %>%
        group_by(GENE, regulator) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        group_by(GENE, regulator) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    roc_analysis = lapply(driver_types, function(driver_type_oi){
        result = driver_classif %>%
            filter(driver_type==driver_type_oi) %>%
            left_join(
                survival_analysis_surv,
                by="regulator"
            ) %>%
            mutate(surv_type = ifelse(coxph_coef>0, "High Risk", "Low Risk")) %>%
            count(GENE, driver_type, surv_type) %>%
            group_by(driver_type) %>%
            roc(response=surv_type, predictor=n, levels=driver_levels[[driver_type_oi]])
        
        result = result %>%
            coords(transpose=FALSE, re="all") %>% 
            mutate(
                fpr = 1 - specificity,
                auc = result[["auc"]],
                driver_type = driver_type_oi
            ) %>%
            ungroup()
    }) %>% do.call(rbind, .)
    
    return(roc_analysis)
}


plot_survival_analysis = function(survival_roc_analysis){
    plts = list()
    
    X = survival_roc_analysis
    
    plts[["survival_analysis-cancers_all-roc_curves"]] = X %>% 
        filter(cancers_subset=="full") %>%
        mutate(auc_lab = sprintf("AUC(%s)=%s",driver_type,round(auc,2))) %>%
        ggscatter(x="fpr", y="sensitivity", size=1,
                  color="driver_type", palette=PAL_DRIVER_TYPE) +
        geom_path(aes(color=driver_type), linetype="dashed", size=LINE_SIZE) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab, color=driver_type),
            . %>% distinct(auc_lab, driver_type) %>%
                mutate(x=0.12, y=c(0,0.075)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme(aspect.ratio=1) +
        labs(
            x="False Positive Rate (1 - Specificity)", 
            y="True Positive Rate (Sensitivity)",
            color="Driver Type"
        )
    
    plts[["survival_analysis-cancers_differential-roc_curves"]] = X %>% 
        filter(cancers_subset=="differential") %>%
        mutate(auc_lab = sprintf("AUC(%s)=%s",driver_type,round(auc,2))) %>%
        ggscatter(x="fpr", y="sensitivity", size=1,
                  color="driver_type", palette=PAL_DRIVER_TYPE) +
        geom_path(aes(color=driver_type), linetype="dashed", size=LINE_SIZE) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab, color=driver_type),
            . %>% distinct(auc_lab, driver_type) %>%
                mutate(x=0.12, y=c(0,0.075)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text_repel(aes(label=threshold), size=FONT_SIZE, family=FONT_FAMILY,
                        segment.size=0.1, max.overlaps=50) +
        theme(aspect.ratio=1) +
        labs(
            x="False Positive Rate (1 - Specificity)", 
            y="True Positive Rate (Sensitivity)",
            color="Driver Type"
        )
    
    return(plts)
}


plot_sf_cross_regulation = function(diff_protein_activity, sf_cross_regulation){
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
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE, regulator) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    oncogenic = driver_classif %>% filter(driver_type=="Oncogenic") %>% pull(regulator)
    suppressor = driver_classif %>% filter(driver_type=="Tumor suppressor") %>% pull(regulator)
    
    X = sf_cross_regulation %>%
        mutate(
            driver_type_a = case_when(
                regulator_a %in% oncogenic ~ "Oncogenic",
                regulator_a %in% suppressor ~ "Tumor suppressor"
            ),
            driver_type_b = case_when(
                regulator_b %in% oncogenic ~ "Oncogenic",
                regulator_b %in% suppressor ~ "Tumor suppressor"
            ),
            same_driver_type = driver_type_a==driver_type_b
        ) %>%
        drop_na(driver_type_a, driver_type_b) %>%
        rowwise() %>%
        mutate(lab_pair = paste(sort(c(driver_type_a,driver_type_b)), collapse="\n&\n")) %>%
        ungroup() %>%
        group_by(regulator_a,regulator_b,lab_pair) %>%
        summarize(correlation = median(correlation)) %>%
        ungroup()
    
    # are oncogenic vs suppressor corregulated?
    plts[["sf_cross_regulation-correlations-violin"]] = X %>%
        mutate(
            lab_pair=factor(lab_pair, levels=c("Tumor suppressor\n&\nTumor suppressor","Oncogenic\n&\nTumor suppressor","Oncogenic\n&\nOncogenic"))
        ) %>%
        ggviolin(x="lab_pair", y="correlation", fill="lab_pair", color=NA, trim=TRUE,
                 palette=c(PAL_DRIVER_TYPE[[1]],"lightgray",PAL_DRIVER_TYPE[[2]])) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, position=position_dodge(0.9)) +
        stat_compare_means(
            method="wilcox.test", label="p.signif", ref.group="Oncogenic\n&\nTumor suppressor",
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=-0.5, label=label),
            . %>% count(lab_pair) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        guides(fill="none") +
        theme(aspect.ratio=1) +
        labs(x="", y="SF-SF Protein Activity Correlation")
    
    return(plts)
}


make_enrichments = function(driver_classif){
    
    X = driver_classif %>%
        filter(is_significant) %>%
        mutate(driver_type = ifelse(sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")) %>%
        count(GENE, regulator, driver_type) %>%
        group_by(GENE, regulator) %>%
        mutate(
            n_sign = ifelse(driver_type=="Tumor suppressor", -n, n),
            n_sum = sum(n_sign)
        ) %>%
        filter(abs(n_sum)>THRESH_N_SUM) %>%
        group_by(GENE, regulator) %>%
        slice_max(n, n=1) %>%
        ungroup()
    
    oncogenics = X %>% filter(driver_type=="Oncogenic") %>% pull(GENE)
    suppressors = X %>% filter(driver_type=="Tumor suppressor") %>% pull(GENE)
    tf_enrichments = list(
        "oncogenics" = enricher(oncogenics, TERM2GENE=ontology_chea),
        "suppressors" = enricher(suppressors, TERM2GENE=ontology_chea)
    )
    
    return(tf_enrichments)
}


plot_tf_enrichments = function(tf_enrichments){
    plts = list()
    
    X = tf_enrichments
    
    plts[["tf_enrichments-oncogenic-cnet"]] = X[["oncogenics"]] %>% 
        cnetplot(cex_label_category=0.5, cex_label_gene=0.5, 
                 cex_family_category=FONT_FAMILY, cex_family_category=FONT_FAMILY) +
        guides(size="none") +
        theme(aspect.ratio=1)
    
    plts[["tf_enrichments-suppressor-cnet"]] = X[["suppressors"]] %>% 
        cnetplot(cex_label_category=0.5, cex_label_gene=0.5, 
                 cex_family_category=FONT_FAMILY, cex_family_category=FONT_FAMILY) +
        guides(size="none") +
        theme(aspect.ratio=1)
    
    return(plts)
}


make_plots = function(diff_protein_activity, assocs_gene_dependency, survival_roc_analysis, sf_cross_regulation, tf_enrichments){
    plts = list(
        plot_diff_protein_activity(diff_protein_activity),
        plot_prolif_driver(diff_protein_activity, assocs_gene_dependency),
        plot_survival_analysis(survival_roc_analysis),
        plot_sf_cross_regulation(diff_protein_activity, sf_cross_regulation),
        plot_tf_enrichments(tf_enrichments)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_protein_activity, assocs_gene_dependency, survival_roc_analysis, sf_cross_regulation, tf_enrichments){
    figdata = list(
        "tcga_tumorigenesis" = list(
            "diff_protein_activity" = diff_protein_activity,
            "assocs_gene_dependency" = assocs_gene_dependency,
            "survival_analysis_surv" = survival_analysis_surv,
            "survival_analysis_cat" = survival_analysis_cat, 
            "sf_cross_regulation" = sf_cross_regulation
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
    save_plt(plts, "diff_protein_activity-n_signif_vs_driver_type-bar", '.pdf', figs_dir, width=12, height=9)
    save_plt(plts, "prolif_driver-driver_type_vs_demeter2-violin", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "survival_analysis-cancers_all-roc_curves", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "survival_analysis-cancers_differential-roc_curves", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "sf_cross_regulation-correlations-violin", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "sf_cross_regulation-correlations-violin", '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, "tf_enrichments-oncogenic-cnet", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "tf_enrichments-suppressor-cnet", '.pdf', figs_dir, width=4, height=4)
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
        make_option("--survival_analysis_surv_file", type="character"),
        make_option("--survival_analysis_cat_file", type="character"),
        make_option("--sf_cross_regulation_file", type="character"),
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
    survival_analysis_surv_file = args[["survival_analysis_surv_file"]]
    survival_analysis_cat_file = args[["survival_analysis_cat_file"]]
    sf_cross_regulation_file = args[["sf_cross_regulation_file"]]
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
    sf_cross_regulation = read_tsv(sf_cross_regulation_file)
    ontology_chea = read.gmt(ontology_chea_file)
    
    # prep
    diff_protein_activity = diff_protein_activity %>%
        mutate(
            is_significant = padj < THRESH_FDR
        ) %>%
        left_join(
            gene_annotation[,c("ENSEMBL","GENE")],
            by = c("regulator"="ENSEMBL")
        )
    
    driver_classif = diff_protein_activity %>%
        mutate(driver_type = ifelse(
            sign(`condition_a-median`)>0, "Oncogenic", "Tumor suppressor")
        ) %>%
        filter(is_significant)
    
    # enrichment
    tf_enrichments = make_enrichments(driver_classif)
    
    # roc analysis
    survival_roc_analysis = make_roc_analysis(driver_classif, survival_analysis_surv) %>%
        mutate(cancers_subset="full") %>%
        rbind(
            make_roc_analysis(
                driver_classif %>% 
                filter(cancer_type %in% diff_protein_activity[["cancer_type"]]), 
                survival_analysis_surv %>% 
                filter(cancer_type %in% diff_protein_activity[["cancer_type"]])
            ) %>%
            mutate(cancers_subset="differential")
        )
    
    # plot
    plts = make_plots(diff_protein_activity, assocs_gene_dependency, survival_roc_analysis, sf_cross_regulation, tf_enrichments)
    
    # make figdata
    figdata = make_figdata(diff_protein_activity, assocs_gene_dependency, survival_roc_analysis, sf_cross_regulation, tf_enrichments)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}