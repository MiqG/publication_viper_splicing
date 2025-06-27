require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)
require(ComplexHeatmap)
require(ggplotify)
require(gridExtra)
require(circlize)
require(glmnet)

# variables
RANDOM_SEED = 1234
THRESH_FDR = 0.05

MAPPING_DANIELSSON = data.frame(
    sampleID = c('SRR837858','SRR837859','SRR837860','SRR837861','SRR837862','SRR837863','SRR837864','SRR837865'),
    cell_line__replicate = c('BJ_PRIMARY__rep1','BJ_PRIMARY__rep2','BJ_IMMORTALIZED__rep1','BJ_IMMORTALIZED__rep2','BJ_TRANSFORMED__rep1','BJ_TRANSFORMED__rep2','BJ_METASTATIC__rep1','BJ_METASTATIC__rep2')
)

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
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# genexpr_danielsson_file = file.path(PREP_DIR,"genexpr_tpm","tumorigenesis.tsv.gz")
# protein_activity_danielsson_file = file.path(RESULTS_DIR,"files","protein_activity","tumorigenesis-EX.tsv.gz")
# exons_danielsson_file = file.path(PREP_DIR,"event_psi","tumorigenesis-EX.tsv.gz")
# introns_danielsson_file = file.path(PREP_DIR,"event_psi","tumorigenesis-INT.tsv.gz")
# alta_danielsson_file = file.path(PREP_DIR,"event_psi","tumorigenesis-ALTD.tsv.gz")
# altd_danielsson_file = file.path(PREP_DIR,"event_psi","tumorigenesis-ALTA.tsv.gz")
# metadata_danielsson_file = file.path(PREP_DIR,"metadata","tumorigenesis.tsv.gz")
# proteomics_inhouse_file = file.path(PREP_DIR,'inhouse','carcinogenesis','proteomics.tsv.gz')
# phospho_inhouse_file = file.path(PREP_DIR,'inhouse','carcinogenesis','phosphoproteomics.tsv.gz')
# genexpr_matsumoto_file = file.path(PREP_DIR,"genexpr_tpm","Matsumoto2017.tsv.gz")
# protein_activity_matsumoto_file = file.path(RESULTS_DIR,"files","protein_activity","Matsumoto2017-EX.tsv.gz")
# metadata_matsumoto_file = file.path(PREP_DIR,"metadata","Matsumoto2017.tsv.gz")
# driver_types_file = file.path(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# gene_info_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","tumorigenesis-EX")

##### FUNCTIONS #####
plot_tumorigenesis = function(protein_activity, genexpr, proteomics_inhouse, phospho_inhouse){
    plts = list()
    
    X = protein_activity %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE, condition) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        filter(study_accession%in%c("PRJNA193487","PRJDB3148")) %>%
        mutate(cell_line_name=factor(
            cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                     "BJ_TRANSFORMED","BJ_METASTATIC")
        ))
    
    plts[["tumorigenesis-cell_line_vs_activity-danielsson-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_text(
            aes(y=-3, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")
    
    
    plts[["tumorigenesis-cell_line_vs_activity-matsumoto-violin"]] = X %>%
        mutate(cell_line_name = condition) %>%
        filter(cell_line_name!="CONTROL" & study_accession=="PRJDB3148") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_text(
            aes(y=-3, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")

    # genexpr
    X = genexpr %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(genexpr_tpm = median(genexpr_tpm)) %>%
        ungroup() %>%
        filter(study_accession=="PRJNA193487") %>%
        mutate(
            cell_line_name=factor(
                cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                         "BJ_TRANSFORMED","BJ_METASTATIC")
                )
        )
    ctl_genexpr = X %>%
        filter(cell_line_name=="BJ_PRIMARY") %>%
        distinct(genexpr_tpm, GENE) %>%
        dplyr::rename(ctl_tpm = genexpr_tpm)
    X = X %>%
        left_join(ctl_genexpr, by="GENE") %>%
        mutate(genexpr_tpm_fc = genexpr_tpm - ctl_tpm)

    plts[["tumorigenesis-cell_line_vs_genexpr_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=driver_type, y=genexpr_tpm_fc)) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(ref.group="Non-driver", method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=20) + 
        geom_text(
            aes(y=-3.8, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cell Line", y="Gene Expression log2FC", fill="Driver Type")
    
    plts[["tumorigenesis-cell_line_vs_abs_genexpr_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=driver_type, y=abs(genexpr_tpm_fc))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(ref.group="Non-driver", method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=20) + 
        geom_text(
            aes(y=-0.05, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cell Line", y="|Gene Expression log2FC|", fill="Driver Type")
    
    # proteomics
    X = proteomics_inhouse %>%
        group_by(cell_line_name, driver_type, GENE, Protein.Names, ENSEMBL) %>%
        summarize(proteomics_cpm = median(proteomics_cpm, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(
            cell_line_name=factor(
                cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                         "BJ_TRANSFORMED","BJ_METASTATIC")
                )
        )
    ctl_proteomics = X %>%
        filter(cell_line_name=="BJ_PRIMARY") %>%
        distinct(proteomics_cpm, GENE, Protein.Names, ENSEMBL) %>%
        dplyr::rename(ctl_cpm = proteomics_cpm) %>%
        drop_na(ctl_cpm)
    X = X %>%
        left_join(ctl_proteomics, by=c("GENE","Protein.Names","ENSEMBL")) %>%
        mutate(proteomics_cpm_fc = proteomics_cpm - ctl_cpm)

    plts[["tumorigenesis-cell_line_vs_proteomics_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=driver_type, y=proteomics_cpm_fc)) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(ref.group="Non-driver", method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=20) + 
        geom_text(
            aes(y=-3.3, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cell Line", y="Protein Abundance log2FC", fill="Driver Type")
    
    plts[["tumorigenesis-cell_line_vs_abs_proteomics_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=driver_type, y=abs(proteomics_cpm_fc))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(ref.group="Non-driver", method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=20) + 
        geom_text(
            aes(y=-0.05, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cell Line", y="|Protein Abundance log2FC|", fill="Driver Type")
    
    # phosphoproteomics
    X = phospho_inhouse %>%
        mutate(GENE = sprintf("%s_%s_%s", Protein.Names, GENE, Precursor.Id)) %>%
        group_by(cell_line_name, driver_type, GENE) %>%
        summarize(phospho_cpm = median(phospho_cpm, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(
            cell_line_name=factor(
                cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                         "BJ_TRANSFORMED","BJ_METASTATIC")
                )
        )
    ctl_phospho = X %>%
        filter(cell_line_name=="BJ_PRIMARY") %>%
        distinct(phospho_cpm, GENE) %>%
        dplyr::rename(ctl_cpm = phospho_cpm) %>%
        drop_na(ctl_cpm)
    X = X %>%
        left_join(ctl_phospho, by="GENE") %>%
        mutate(phospho_cpm_fc = phospho_cpm - ctl_cpm)

    plts[["tumorigenesis-cell_line_vs_phospho_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=driver_type, y=phospho_cpm_fc)) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(ref.group="Non-driver", method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=20) + 
        geom_text(
            aes(y=-7.4, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cell Line", y="Phosphorylation Abundance log2FC", fill="Driver Type")
    
    plts[["tumorigenesis-cell_line_vs_abs_phospho_fc-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=driver_type, y=abs(phospho_cpm_fc))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(ref.group="Non-driver", method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY, angle=20) + 
        geom_text(
            aes(y=-0.05, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Cell Line", y="|Phosphorylation Abundance log2FC|", fill="Driver Type")
    
    # random control
    X = protein_activity %>%
        group_by(cell_line_name, study_accession) %>%
        mutate(driver_type = sample(driver_type)) %>%
        ungroup() %>%
        drop_na(driver_type) %>%
        group_by(cell_line_name, driver_type, study_accession, GENE) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        filter(study_accession=="PRJNA193487") %>%
        mutate(cell_line_name=factor(
            cell_line_name, levels=c("BJ_PRIMARY","BJ_IMMORTALIZED",
                                     "BJ_TRANSFORMED","BJ_METASTATIC")
        ))
    
    plts[["tumorigenesis-cell_line_vs_activity-random-violin"]] = X %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        ggplot(aes(x=cell_line_name, y=activity, group=interaction(cell_line_name,driver_type))) +
        geom_violin(aes(fill=driver_type), color=NA, position=position_dodge(0.9)) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black", position=position_dodge(0.9)) +
        fill_palette(PAL_DRIVER_TYPE) +
        stat_compare_means(method="wilcox.test", label="p.format", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_text(
            aes(y=-3, label=label, group=driver_type),
            . %>% count(cell_line_name, driver_type) %>% mutate(label=paste0("n=",n)),
            size=FONT_SIZE, family=FONT_FAMILY, position=position_dodge(0.9)
        ) +
        theme_pubr() +
        labs(x="Cell Line", y="Protein Activity", fill="Driver Type")
    
    return(plts)
}


plot_associations = function(associations){
    plts = list()
    
    X = associations %>%
        mutate(
            label = ifelse(str_detect(feature,"ENSG"), sprintf("%sexpr", GENE), feature)
        ) %>%
        filter(feature!="Intercept") # %>%
        # group_by() %>%
        # mutate(n_total = n(), n_zero = coefficient==0) %>%
        # ungroup() %>%
        # filter(n_zero < n_total)
    
    plts[["associations-driver_type_vs_coefficients-bar"]] = X %>%
        mutate(
            abs_coefficient = abs(coefficient),
            regmethod = factor(regmethod, levels=c("lasso","elnet"))
        ) %>%
        arrange(-abs_coefficient) %>%
        filter(coefficient!=0) %>%
        ggbarplot(x="label", y="coefficient", fill="driver_type_feature", color=NA, palette=PAL_DRIVER_TYPE) +
        theme_pubr(x.text.angle=70) +
        facet_grid(~regmethod, scales="free_x", space="free_x") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Associated features", y="Coefficient", fill="SF Type")
    
    return(plts)
}

make_plots = function(
    protein_activity, genexpr, proteomics_inhouse, phospho_inhouse, associations
){
    plts = list(
        plot_tumorigenesis(protein_activity, genexpr, proteomics_inhouse, phospho_inhouse),
        plot_associations(associations)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    protein_activity
){
    figdata = list(
        "tumorigenesis" = list(
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
    save_plt(plts, "tumorigenesis-cell_line_vs_activity-danielsson-violin", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_activity-matsumoto-violin", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_abs_genexpr_fc-violin", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_genexpr_fc-violin", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_abs_proteomics_fc-violin", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_proteomics_fc-violin", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_abs_phospho_fc-violin", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_phospho_fc-violin", '.pdf', figs_dir, width=7, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_activity-random-violin", '.pdf', figs_dir, width=6, height=6)
    
    save_plt(plts, "associations-driver_type_vs_coefficients-bar", '.pdf', figs_dir, width=20, height=16)
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
        make_option("--annotation_file", type="character"),
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--driver_types_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    genexpr_file = args[["genexpr_file"]]
    annotation_file = args[["annotation_file"]]
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    driver_types_file = args[["driver_types_file"]]
    figs_dir = args[["figs_dir"]]
    
    set.seed(RANDOM_SEED)    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    annot = read_tsv(annotation_file)
    driver_types = read_tsv(driver_types_file)
    genexpr_danielsson = read_tsv(genexpr_danielsson_file)
    protein_activity_danielsson = read_tsv(protein_activity_danielsson_file)
    exons_danielsson = read_tsv(exons_danielsson_file)
    introns_danielsson = read_tsv(introns_danielsson_file)
    alta_danielsson = read_tsv(alta_danielsson_file)
    altd_danielsson = read_tsv(altd_danielsson_file)
    metadata_danielsson = read_tsv(metadata_danielsson_file)
    proteomics_inhouse = read_tsv(proteomics_inhouse_file)
    phospho_inhouse = read_tsv(phospho_inhouse_file)
    genexpr_matsumoto = read_tsv(genexpr_matsumoto_file)
    protein_activity_matsumoto = read_tsv(protein_activity_matsumoto_file)
    metadata_matsumoto = read_tsv(metadata_matsumoto_file)
    splicing_factors = read_tsv(splicing_factors_file)
    gene_info = read_tsv(gene_info_file) %>% 
        dplyr::rename(GENE=`Approved symbol`, ENSEMBL=`Ensembl gene ID`) %>%
        distinct(GENE, ENSEMBL)

    # prep
    ## protein activity
    protein_activity_danielsson = protein_activity_danielsson %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata_danielsson, by="sampleID") %>%
        drop_na(condition, activity) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            )
        )
    
    protein_activity_matsumoto = protein_activity_matsumoto %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata_matsumoto, by="sampleID") %>%
        drop_na(condition, activity) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            ),
            PERT_GENE = NA,
            PERT_ENSEMBL = NA
        )
    
    protein_activity = protein_activity_danielsson %>%
        bind_rows(protein_activity_matsumoto) %>%

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
        left_join(driver_types, by=c("regulator"="ENSEMBL"))
    
    ## gene expression
    genexpr_danielsson = genexpr_danielsson %>%
        filter(ID%in%splicing_factors[["ENSEMBL"]]) %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_tpm") %>%
        left_join(metadata_danielsson, by="sampleID") %>%
        drop_na(condition, genexpr_tpm) %>%
        filter(study_accession=="PRJNA193487")
    
    genexpr_matsumoto = genexpr_matsumoto %>%
        filter(ID%in%splicing_factors[["ENSEMBL"]]) %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_tpm") %>%
        left_join(metadata_matsumoto, by="sampleID") %>%
        drop_na(condition, genexpr_tpm) %>%
        mutate(
            PERT_ENSEMBL = NA,
            PERT_GENE = NA
        )
    
    genexpr = genexpr_danielsson %>%
        bind_rows(genexpr_matsumoto) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s%s) (%s%s) | %s | %s", condition, pert_time, pert_time_units, 
                pert_concentration, pert_concentration_units, cell_line_name, study_accession
            )
        ) %>%
        
        # summarize replicates
        group_by(condition_lab, condition, pert_time, pert_time_units, 
                 pert_concentration, pert_concentration_units, cell_line_name, study_accession,
                 PERT_ENSEMBL, PERT_GENE, ID) %>%
        summarize(
            genexpr_tpm = median(genexpr_tpm, na.rm=TRUE),
        ) %>%
        ungroup() %>%
        left_join(gene_info, by=c("ID"="ENSEMBL")) %>%
        left_join(driver_types, by=c("ID"="ENSEMBL", "GENE")) %>%
        mutate(driver_type = replace_na(driver_type, "Non-driver"))

    # splicing of splicing factors
    events_oi = annot %>% filter(GENE%in%driver_types[["GENE"]]) %>% distinct(EVENT,GENE)
    exons_danielsson = exons_danielsson %>%
        filter(EVENT %in% events_oi[["EVENT"]]) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi") %>%
        drop_na(psi) %>%
        left_join(metadata_danielsson, by="sampleID") %>%
        drop_na(condition, psi) %>%
        filter(study_accession=="PRJNA193487")
    introns_danielsson = introns_danielsson %>%
        filter(EVENT %in% events_oi[["EVENT"]]) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi") %>%
        drop_na(psi) %>%
        left_join(metadata_danielsson, by="sampleID") %>%
        drop_na(condition, psi) %>%
        filter(study_accession=="PRJNA193487")
    alta_danielsson = alta_danielsson %>%
        filter(EVENT %in% events_oi[["EVENT"]]) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi") %>%
        drop_na(psi) %>%
        left_join(metadata_danielsson, by="sampleID") %>%
        drop_na(condition, psi) %>%
        filter(study_accession=="PRJNA193487")
    altd_danielsson = altd_danielsson %>%
        filter(EVENT %in% events_oi[["EVENT"]]) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi") %>%
        drop_na(psi) %>%
        left_join(metadata_danielsson, by="sampleID") %>%
        drop_na(condition, psi) %>%
        filter(study_accession=="PRJNA193487")
    
    # proteomics and phosphoproteomics of splicing factors
    duplicated_genes = proteomics_inhouse %>% count(Genes) %>% drop_na() %>% filter(n>1) %>% pull(Genes)
    proteomics_inhouse = proteomics_inhouse %>%
        rename(GENE=Genes) %>%
        drop_na(GENE) %>%
        filter(!GENE%in%duplicated_genes) %>%
        left_join(gene_info, by="GENE") %>%
        filter(ENSEMBL%in%splicing_factors[["ENSEMBL"]]) %>%
        pivot_longer(
            -c(ENSEMBL,GENE,Protein.Group,Protein.Names,First.Protein.Description,N.Sequences,N.Proteotypic.Sequences), 
            names_to="cell_line__replicate", 
            values_to="proteomics_cpm"
        ) %>%
        left_join(MAPPING_DANIELSSON, by="cell_line__replicate") %>%
        separate(cell_line__replicate, sep="__", into=c("cell_line_name","replicate"), remove=FALSE) %>%
        left_join(driver_types, by=c("ENSEMBL", "GENE")) %>%
        mutate(driver_type = replace_na(driver_type, "Non-driver"))
    
    phospho_inhouse = phospho_inhouse %>%
        mutate(
            GENE = Genes,
            Genes = str_c(Genes, Modified.Sequence, sep="_")
        ) %>%
        drop_na(GENE) %>%
        left_join(gene_info, by="GENE") %>%
        filter(ENSEMBL%in%splicing_factors[["ENSEMBL"]]) %>%
        pivot_longer(
            -c(ENSEMBL,GENE,Genes,Protein.Group,Protein.Names,Protein.Ids,First.Protein.Description,Proteotypic,Stripped.Sequence,Modified.Sequence,Precursor.Charge,Precursor.Id), 
            names_to="cell_line__replicate", 
            values_to="phospho_cpm"
        ) %>%
        left_join(MAPPING_DANIELSSON, by="cell_line__replicate") %>%
        separate(cell_line__replicate, sep="__", into=c("cell_line_name","replicate"), remove=FALSE) %>%
        left_join(driver_types, by=c("ENSEMBL", "GENE")) %>%
        mutate(driver_type = replace_na(driver_type, "Non-driver"))
    
    # association protein activity with other omics
    ## Danielsson
    X = protein_activity_danielsson %>%
        distinct(study_accession, sampleID, condition_lab, cell_line_name, regulator, activity) %>%
        filter(study_accession=="PRJNA193487") %>%
        filter(cell_line_name!="BJ_PRIMARY") %>%
        left_join(
            genexpr_danielsson %>% 
                distinct(sampleID, ID, genexpr_tpm),
            by="sampleID", relationship="many-to-many" # SF genexpr vs SF activity
        )
    ### are changes in splicing factor activity explained by the changes in the expression of splicing factors?
    mat_genexpr = genexpr_danielsson %>% 
        distinct(sampleID, ID, genexpr_tpm) %>% 
        pivot_wider(id_cols="sampleID", names_from="ID", values_from="genexpr_tpm") %>%
        column_to_rownames("sampleID")
    
    ### are changes in splicing factor activity explained by changes in splicing?
    mat_splicing = list(
            exons_danielsson,
            introns_danielsson,
            alta_danielsson,
            altd_danielsson
        ) %>% bind_rows() %>%
        left_join(annot, by="EVENT") %>%
        mutate(event_gene = sprintf("%s_%s", EVENT, GENE)) %>%
        distinct(event_gene, sampleID, psi) %>%
        pivot_wider(id_cols="sampleID", names_from="event_gene", values_from="psi") %>%
        column_to_rownames("sampleID")
    #### drop columns with missing values
    mat_splicing = mat_splicing[,colSums(is.na(mat_splicing))==0] # from 4457 to 3320
    #### drop_columns with constant values
    idx_notctt = apply(mat_splicing, 2, function(col){ length(unique(col)) > 1 })
    mat_splicing = mat_splicing[,idx_notctt] # from 3320 to 1481
    
    ### are changes in splicing factor activity explained by changes in proteomics or phosphoproteomics?
    mat_proteomics = proteomics_inhouse %>%
        mutate(GENE = sprintf("%sprot", GENE)) %>%
        distinct(sampleID, GENE, proteomics_cpm) %>%
        pivot_wider(id_cols="sampleID", names_from="GENE", values_from="proteomics_cpm") %>%
        column_to_rownames("sampleID")
    
    mat_phospho = phospho_inhouse %>%
        mutate(GENE = sprintf("%s_%s", GENE, Precursor.Id)) %>%
        distinct(sampleID, GENE, phospho_cpm) %>%
        drop_na() %>%
        pivot_wider(id_cols="sampleID", names_from="GENE", values_from="phospho_cpm") %>%
        column_to_rownames("sampleID")
    #### drop columns with missing values
    mat_proteomics = mat_proteomics[,colSums(is.na(mat_proteomics))==0] # from 4461 to 444
    mat_phospho = mat_phospho[,colSums(is.na(mat_phospho))==0] # from 4511 to 2368
    #### drop_columns with constant values
    idx_notctt = apply(mat_proteomics, 2, function(col){ length(unique(col)) > 1 })
    mat_proteomics = mat_proteomics[,idx_notctt] # from 444 to 444
    idx_notctt = apply(mat_phospho, 2, function(col){ length(unique(col)) > 1 })
    mat_phospho = mat_phospho[,idx_notctt] # from 2368 to 2368
    
    #### subtract control samples
    ctl_samples = metadata_danielsson %>% filter(cell_line_name=="BJ_PRIMARY") %>% pull(sampleID)
    mat_genexpr = sweep(mat_genexpr, 2, colMeans(mat_genexpr[ctl_samples,]), "-")
    mat_splicing = sweep(mat_splicing, 2, colMeans(mat_splicing[ctl_samples,]), "-")
    mat_proteomics = sweep(mat_proteomics, 2, colMeans(mat_proteomics[ctl_samples,]), "-")
    mat_phospho = sweep(mat_phospho, 2, colMeans(mat_phospho[ctl_samples,]), "-")
    
    mat_activity = X %>%
        distinct(sampleID, cell_line_name, regulator, activity) %>% 
        left_join(driver_types %>% distinct(ENSEMBL, driver_type), by=c("regulator"="ENSEMBL")) %>%
        drop_na(driver_type) %>%
        # compute median activity difference between SF programs
        group_by(sampleID, cell_line_name, driver_type) %>%
        summarize(activity = median(activity)) %>%
        ungroup() %>%
        group_by(sampleID, cell_line_name) %>%
        summarize(
            activity = ifelse(driver_type=="Tumor suppressor", -activity, activity),
            activity = sum(activity),
            regulator = "sf_program_diff"
        ) %>% 
        distinct() %>%
        ungroup() %>%
        pivot_wider(id_cols="sampleID", names_from="regulator", values_from="activity") %>%
        column_to_rownames("sampleID")
    
    samples_oi = rownames(mat_activity)
    regulators_oi = colnames(mat_activity)
    alphas = c("lasso"=1, "elnet"=0.5)
    
    associations = lapply(names(alphas), function(regmethod){
        
        assoc = lapply(regulators_oi, function(regulator_oi){

            features = cbind(
                mat_genexpr[samples_oi,,drop=FALSE], 
                mat_splicing[samples_oi,,drop=FALSE],
                mat_proteomics[samples_oi,,drop=FALSE],
                mat_phospho[samples_oi,,drop=FALSE]
            ) %>% scale() %>% as.matrix()
            features = features[, colSums(is.na(features)) == 0]
            features = cbind(Intercept=1, features)
            target = mat_activity[samples_oi, regulator_oi] %>% as.numeric()
            fit = glmnet::glmnet(x=features, y=target, family="gaussian", alpha=alphas[regmethod], nlambda=100) # lasso

            #### get coefficients
            coefs = fit[["beta"]]
            coefs = coefs[,ncol(coefs)]

            #### get predicability metrics
            pred = predict(fit, features)
            pred = pred[,ncol(pred)]
            assoc = data.frame(
                regulator = regulator_oi,
                corr_pearson = cor(pred, target, method="pearson"),
                corr_spearman = cor(pred, target, method="spearman"),
                mse = mean((pred - target)^2),
                coefficient = coefs,
                feature = names(coefs),
                regmethod = regmethod
            )
        
            return(assoc)
        }) %>% bind_rows()
        
        return(assoc)
    }) %>% bind_rows()
    rownames(associations) = NULL
    
    ## add info to features of the models
    annot = annot %>% 
        mutate(event_gene = sprintf("%s_%s", EVENT, GENE)) %>% 
        left_join(gene_info, by="GENE") %>%
        distinct(event_gene, ENSEMBL, GENE) %>%
        drop_na()
    
    associations = associations %>%
        left_join(
            bind_rows(list(
                annot %>% mutate(feature=event_gene),
                gene_info %>% mutate(feature=ENSEMBL) %>% drop_na(ENSEMBL)
            )),
            by="feature"
        ) %>%
        left_join(
            driver_types %>% 
            dplyr::rename(driver_type_feature=driver_type) %>% 
            distinct(ENSEMBL,driver_type_feature), 
            by="ENSEMBL"
        ) %>%
        mutate(
            driver_type_feature = replace_na(driver_type_feature, "Non-driver")
        )
    
    # plot
    plts = make_plots(protein_activity, genexpr, proteomics_inhouse, phospho_inhouse, associations)
    
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