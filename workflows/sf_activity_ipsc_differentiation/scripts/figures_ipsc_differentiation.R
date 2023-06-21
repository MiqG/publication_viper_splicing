#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - Interesting SFs
#    - https://www.nature.com/articles/s41588-021-00851-w :
#        - QKI (opposite brain vs heart)
#    - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6839889/ :
#        - Brain:
#            - PTBP1
#            - PTBP2
#            - SRRM4
#            - RBFOX1
#            - RBFOX3
#            - NOVA2
#            - KHDRBS3
#            - CTCF (brain epigenetics)
#            - TDP43 (neurological disorders)
#        - Muscle:
#            - CELF1
#            - RBFOX1
#            - RBFOX2
#            - RBM24
#            - MBNL1
#            - RBM20
#            - SF3B1
#            - PTBP1
#            - QKI
#        - Pancreas:
#            - NOVA1
#            - RBM4
#            - SRSF3
#            - SRSF10
#            - SLU7
#            - ESRP2
#        - Differentiation
#            - PTBP1 (smooth muscle cells)
#            - HNRNPA
#            - HNRNPB
#            - SNRP70
#            - HNRPLL
#            - MBNL2


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

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

GENES_OI = c(
    "PTBP1",
    "PTBP2",
    "SRRM4",
    "RBFOX1",
    "RBFOX3",
    "NOVA2",
    "KHDRBS3",
    "CTCF",
    "TDP43",
    "CELF1",
    "RBFOX1",
    "RBFOX2",
    "RBM24",
    "MBNL1",
    "RBM20",
    "SF3B1",
    "PTBP1",
    "QKI",
    "NOVA1",
    "RBM4",
    "SRSF3",
    "SRSF10",
    "SLU7",
    "ESRP2",
    "PTBP1",
    "HNRNPA",
    "HNRNPB",
    "SNRP70",
    "HNRPLL",
    "MBNL2"
)


DEV_STAGES = list(
    "PRJEB1195" = c(
        "iPSC","DEFINITIVE_ENDODERM","PRIMITIVE_GUT_TUBE","POSTERIOR_FOREGUT",
        "PANCREATIC_ENDODERM","LATE_PANCREATIC_ENDODERM","ENDOCRINE_CELLS","ADULT_BETA_CELL"
    ),
    "PRJNA379280"=c(
        "iPSC","ALVEOLAR_EPITHELIAL_PROGENITOR","ALVEOLAR_EPITHELIAL_TYPE_II","PRIMARY_ALVEOLAR_EPITHELIAL_TYPE_II"
    ),
    "PRJNA596331"=c(
        "iPSC","ACC_DORSAL",
        "NPC","ROSETTE","NEURONS_ALONE","NEURONS_PLUS_ASTROS"
    ),
    "PRJNA665705"=c(
        "iPSC","MESODERM","CARDIOMYOCYTE"
    )
)


# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_ipsc_differentiation")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","ipsc_differentiation.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","ipsc_differentiation.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","ipsc_differentiation.tsv.gz")
# annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_ipsc_differentiation")

##### FUNCTIONS #####
plot_eda_differentiation = function(genexpr, protein_activity, metadata, annot, study_oi){
    plts = list()
    
    conditions_oi = DEV_STAGES[[study_oi]]
    
    X = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(
            genexpr %>%
            pivot_longer(-ID, names_to="sampleID", values_to="log2_tpm"),
            by=c("regulator"="ID","sampleID")
        ) %>%
        left_join(metadata, by="sampleID") %>%
        filter(study_accession==study_oi & condition%in%conditions_oi) %>%
        left_join(annot, by=c("regulator"="ENSEMBL")) %>%
        mutate(condition = factor(condition, levels=conditions_oi))
    
    x = X %>%
        group_by(study_accession, condition, GENE) %>%
        summarize(activity = mean(activity)) %>%
        ungroup() %>%
        group_by(GENE) %>%
        mutate(corr_activity = cor(as.numeric(condition), activity, method="spearman")*max(abs(activity))) %>%
        ungroup()
    
    sfs_oi = x %>%
        distinct(GENE,corr_activity) %>%
        slice_max(corr_activity, n=5) %>%
        bind_rows(
            x %>%
            distinct(GENE,corr_activity) %>%
            slice_min(corr_activity, n=5)
        ) %>%
        pull(GENE)
    sfs_oi = union(sfs_oi, GENES_OI)
    
    plts[["eda_differentiation-top_assocs-activity-line"]] = X %>%
        filter(GENE %in% sfs_oi) %>%
        mutate(GENE = ifelse(GENE %in% GENES_OI, paste0(GENE,"*"), GENE)) %>%
        group_by(GENE, condition) %>%
        summarize(
            mean = mean(activity),
            se = sd(activity),
            ymin = mean - se,
            ymax = mean + se
        ) %>%
        ungroup() %>%
        ggplot(aes(x=condition, y=mean, color=GENE, group=GENE)) +
        geom_line(size=0.1, linetype="dashed") +
        geom_point(size=0.1) +
        geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.01, size=LINE_SIZE) +
        color_palette(get_palette("lancet", length(sfs_oi))) +
        theme_pubr(x.text.angle = 70) +
        geom_text(
            aes(label=GENE),
            . %>% filter(condition==conditions_oi[length(conditions_oi)]),
            size=FONT_SIZE, family=FONT_FAMILY, hjust=0, vjust=0.5, nudge_x=0.1
        ) +
        guides(color="none") +
        theme(aspect.ratio=1) +
        labs(x="Developmental Stage", y="Protein Activity")
    
    
    plts[["eda_differentiation-top_assocs-genexpr-line"]] = X %>%
        filter(GENE %in% sfs_oi) %>%
        mutate(GENE = ifelse(GENE %in% GENES_OI, paste0(GENE,"*"), GENE)) %>%
        group_by(GENE, condition) %>%
        summarize(
            mean = mean(log2_tpm),
            se = sd(log2_tpm),
            ymin = mean - se,
            ymax = mean + se
        ) %>%
        ungroup() %>%
        ggplot(aes(x=condition, y=mean, color=GENE, group=GENE)) +
        geom_line(size=0.1, linetype="dashed") +
        geom_point(size=0.1) +
        geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.01, size=LINE_SIZE) +
        color_palette(get_palette("lancet", length(sfs_oi))) +
        theme_pubr(x.text.angle = 70) +
        geom_text(
            aes(label=GENE),
            . %>% filter(condition==conditions_oi[length(conditions_oi)]),
            size=FONT_SIZE, family=FONT_FAMILY, hjust=0, vjust=0.5, nudge_x=0.1
        ) +
        guides(color="none") +
        theme(aspect.ratio=1) +
        labs(x="Developmental Stage", y="log2(TPM+1)")
    
    names(plts) = sprintf("%s-%s", names(plts), study_oi)
    
    return(plts)
}


make_plots = function(genexpr, protein_activity, metadata, annot){
    plts = list(
        plot_eda_differentiation(genexpr, protein_activity, metadata, annot, "PRJEB1195"),
        plot_eda_differentiation(genexpr, protein_activity, metadata, annot, "PRJNA379280"),
        plot_eda_differentiation(genexpr, protein_activity, metadata, annot, "PRJNA596331"),
        plot_eda_differentiation(genexpr, protein_activity, metadata, annot, "PRJNA665705")
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
    save_plt(plts, "eda_differentiation-top_assocs-activity-line-PRJEB1195", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "eda_differentiation-top_assocs-genexpr-line-PRJEB1195", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "eda_differentiation-top_assocs-activity-line-PRJNA379280", '.pdf', figs_dir, width=5, height=9)
    save_plt(plts, "eda_differentiation-top_assocs-genexpr-line-PRJNA379280", '.pdf', figs_dir, width=5, height=9)
    save_plt(plts, "eda_differentiation-top_assocs-activity-line-PRJNA596331", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "eda_differentiation-top_assocs-genexpr-line-PRJNA596331", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "eda_differentiation-top_assocs-activity-line-PRJNA665705", '.pdf', figs_dir, width=5, height=8)
    save_plt(plts, "eda_differentiation-top_assocs-genexpr-line-PRJNA665705", '.pdf', figs_dir, width=5, height=8)
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