#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# check drug target protein activity
# - interesting questions:
#   - https://www.nature.com/articles/s41429-021-00450-1
#      - isoginkgetin does not have a clear mechanis of action
#      - Spliceostatin A inhibits SF3B1
# - H3B-8800 targets mutant SF3B1
# - E7107 blocks binding to U2 snRNP to pre-mRNA, essential for SF3B complex formation
# - PLADIENOLIDE_B inhibits SF3B1
# - inhibitors of upstream regulators
#     - GSK3326595 inhibits PRMT5, a methyl-transferase
#     - PALBOCICLIB inhibits CDK4/6 (phosphorylation)

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)

# variables
RANDOM_SEED = 1234

SF3B_COMPLEX = c(
    "SF3B1", "SF3B2", "SF3B3", "SF3B4", 
    "SF3B5", "SF3B6", "PHF5A", "DDX42"
)

SF3B_DRUGS = c(
    "E7107",
    "H3B-8800",
    "SPLICEOSTATIN_A",
    "PLADIENOLIDE_B"
)

OTHER_DRUGS = c(
    "GSK3326595",
    "PALBOCICLIB",
    "ISOGINKGETIN"
)

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "darkgrey"
PAL_ACCENT = "darkred"
PAL_DUAL = c(PAL_DARK, PAL_ACCENT)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_sf_drugs")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","sf_drugs.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","sf_drugs.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","sf_drugs.tsv.gz")
# annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_drug_targets")

##### FUNCTIONS #####
plot_sf3b_complex = function(genexpr, protein_activity, metadata, annot){
    plts = list()
    
    genes_oi = annot %>%
        filter(GENE %in% SF3B_COMPLEX)
    
    # Is the gene expression or protein activity of drug targets different from controls?
    ## combine data
    X = genes_oi %>%
        left_join(
            genexpr %>%
            pivot_longer(-ID, names_to="sampleID", values_to="log2_tpm"),
            by=c("ENSEMBL"="ID")
        ) %>%
        left_join(
            protein_activity %>%
            pivot_longer(-regulator, names_to="sampleID", values_to="activity"),
            by=c("ENSEMBL"="regulator","sampleID")
        ) %>%
        left_join(metadata, by="sampleID") %>%
        mutate(GENE = factor(GENE, levels=SF3B_COMPLEX)) %>%
        drop_na(activity, condition)
    
    ## subset
    x = lapply(
        SF3B_DRUGS, function(drug_oi){
            studies_oi = X %>% filter(condition == drug_oi) %>% pull(study_accession)
            control = ifelse(drug_oi=="SPLICEOSTATIN_A","METHANOL","DMSO")
            x = X %>%
                filter(
                    study_accession%in%studies_oi & 
                    (str_detect(condition, drug_oi) | str_detect(condition, control))
                ) %>%
                mutate(DRUG = drug_oi)
            return(x)
        }) %>% do.call(rbind,.)
    
    ## statistical tests
    stat_tests_genexpr = x %>%
        drop_na(GENE,log2_tpm,activity,condition,study_accession,DRUG,cell_line_name) %>%
        compare_means(
            log2_tpm~condition, data=., method="t.test", label="p.signif", ref.group="DMSO", 
            group.by=c("study_accession","DRUG","GENE","cell_line_name")
        ) %>%
        bind_rows(
            x %>%
            drop_na(GENE,activity,condition,study_accession,DRUG,cell_line_name) %>%
            compare_means(
                log2_tpm~condition, data=., method="t.test", label="p.signif", ref.group="METHANOL", 
                group.by=c("study_accession","DRUG","GENE","cell_line_name")
            )
        )
    
    stat_tests_activity = x %>%
        drop_na(GENE,log2_tpm,activity,condition,study_accession,DRUG,cell_line_name) %>%
        compare_means(
            activity~condition, data=., method="t.test", label="p.signif", ref.group="DMSO", 
            group.by=c("study_accession","DRUG","GENE","cell_line_name")
        ) %>%
        bind_rows(
            x %>%
            drop_na(GENE,activity,condition,study_accession,DRUG,cell_line_name) %>%
            compare_means(
                activity~condition, data=., method="t.test", label="p.signif", ref.group="METHANOL", 
                group.by=c("study_accession","DRUG","GENE","cell_line_name")
            )
        )
    
    plts[["sf3b_complex-drugs_vs_control-genexpr-box"]] = x %>%
        ggplot(aes(x=GENE, y=log2_tpm, group=interaction(GENE,condition))) +
            geom_point(aes(color=condition), position=position_jitterdodge(0.1, dodge.width=0.9), size=1) +
            geom_boxplot(fill=NA, width=0.5, outlier.shape=NA, position=position_dodge(0.9)) +
            color_palette("Paired") +
            theme_pubr(x.text.angle=70) +
            facet_wrap(~study_accession+DRUG+cell_line_name, ncol=1) +
            theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) + 
            geom_text(
                aes(y=2, label=p.signif),
                stat_tests_genexpr %>% mutate(condition=group2), position=position_dodge(0.9),
                size=FONT_SIZE, family=FONT_FAMILY
            ) +
            labs(x="Gene", y="log2(TPM+1)", color="Condition")
    
    plts[["sf3b_complex-drugs_vs_control-activity-box"]] = x %>%
        ggplot(aes(x=GENE, y=activity, group=interaction(GENE,condition))) +
            geom_point(aes(color=condition), position=position_jitterdodge(0.1, dodge.width=0.9), size=1) +
            geom_boxplot(fill=NA, width=0.5, outlier.shape=NA, position=position_dodge(0.9)) +
            color_palette("Paired") +
            theme_pubr(x.text.angle=70) +
            facet_wrap(~study_accession+DRUG+cell_line_name, ncol=1) +
            theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) + 
            geom_text(
                aes(y=2, label=p.signif),
                stat_tests_activity %>% mutate(condition=group2), position=position_dodge(0.9),
                size=FONT_SIZE, family=FONT_FAMILY
            ) +
            labs(x="Gene", y="Protein Activity", color="Condition")
    
    # What is the ranking of our SF3B complex proteins?
    X = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        left_join(annot, by=c("regulator"="ENSEMBL")) %>%
        drop_na(condition) %>%
        group_by(condition, study_accession, cell_line_name, GENE) %>%
        summarize(activity = median(activity)) %>% # take median of activity
        ungroup() %>%
        group_by(condition,study_accession,cell_line_name) %>%
        arrange(activity) %>%
        mutate(activity_ranking = row_number()) %>%
        ungroup() %>%
        mutate(is_regulator_oi = GENE %in% SF3B_COMPLEX)
    
    for (drug_oi in SF3B_DRUGS){
        studies_oi = X %>% filter(condition == drug_oi) %>% pull(study_accession)
        control = ifelse(drug_oi=="SPLICEOSTATIN_A","METHANOL","DMSO")
        x = X %>%
            filter(
                study_accession%in%studies_oi & 
                (str_detect(condition, drug_oi) | str_detect(condition, control))
            )
        
        plts[[paste0("sf3b_complex-drugs_vs_control-activity_ranking-scatter-",drug_oi)]] = x %>%
            ggplot(aes(x=activity_ranking, y=activity)) +
            geom_scattermore(
                data = . %>% filter(!is_regulator_oi), pixels=c(1000,1000), 
                pointsize=4, color=PAL_DARK, alpha=0.5
            ) +
            geom_scattermore(
                aes(color=GENE), 
                data = . %>% filter(is_regulator_oi) %>%
                        mutate(GENE = factor(GENE, levels=SF3B_COMPLEX)), 
                pixels=c(1000,1000), pointsize=10
            ) +
            geom_text_repel(
                aes(label=GENE),
                data= . %>% group_by(condition,study_accession,cell_line_name) %>% slice_min(activity_ranking, n=5),
                size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
            )+
            color_palette("Dark2") +
            theme_pubr() +
            facet_wrap(~condition+study_accession+cell_line_name) +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(x="Ranking", y="median(Protein Activity)", color="SF3B Complex Gene")
    }
    
    return(plts)
}


plot_moa_discovery = function(protein_activity, metadata, annot){
    plts = list()
    
    # Which SFs have the most extreme changes in protein activity?
    X = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        left_join(annot, by=c("regulator"="ENSEMBL")) %>%
        drop_na(condition) %>%
        group_by(condition, study_accession, cell_line_name, GENE) %>%
        summarize(activity = median(activity)) %>% # take median of activity
        ungroup() %>%
        group_by(condition,study_accession,cell_line_name) %>%
        arrange(activity) %>%
        mutate(activity_ranking = row_number()) %>%
        ungroup()
    
    for (drug_oi in OTHER_DRUGS){
        studies_oi = X %>% filter(condition == drug_oi) %>% pull(study_accession)
        control = "DMSO"
        x = X %>%
            filter(
                study_accession%in%studies_oi & 
                (str_detect(condition, drug_oi) | str_detect(condition, control))
            )
        
        plts[[paste0("moa_discovery-drugs_vs_control-activity_ranking-scatter-",drug_oi)]] = x %>%
            ggplot(aes(x=activity_ranking, y=activity)) +
            geom_scattermore(
                pixels=c(1000,1000), pointsize=4, color=PAL_DARK, alpha=0.5
            ) +
            geom_text_repel(
                aes(label=GENE),
                data= . %>% group_by(condition,study_accession,cell_line_name) %>% slice_min(activity_ranking, n=5),
                size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
            )+
            color_palette("Dark2") +
            theme_pubr() +
            facet_wrap(~condition+study_accession+cell_line_name) +
            theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
            labs(x="Ranking", y="median(Protein Activity)")
    }
    
    return(plts)
}


make_plots = function(genexpr, protein_activity, metadata, annot){
    plts = list(
        plot_sf3b_complex(genexpr, protein_activity, metadata, annot),
        plot_moa_discovery(protein_activity, metadata, annot)
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
    save_plt(plts, "sf3b_complex-drugs_vs_control-genexpr-box", '.pdf', figs_dir, width=5, height=30)
    save_plt(plts, "sf3b_complex-drugs_vs_control-activity-box", '.pdf', figs_dir, width=5, height=30)
    save_plt(plts, "sf3b_complex-drugs_vs_control-activity_ranking-scatter-E7107", '.pdf', figs_dir, width=10, height=13)
    save_plt(plts, "sf3b_complex-drugs_vs_control-activity_ranking-scatter-H3B-8800", '.pdf', figs_dir, width=10, height=13)
    save_plt(plts, "sf3b_complex-drugs_vs_control-activity_ranking-scatter-SPLICEOSTATIN_A", '.pdf', figs_dir, width=8, height=13)
    save_plt(plts, "sf3b_complex-drugs_vs_control-activity_ranking-scatter-PLADIENOLIDE_B", '.pdf', figs_dir, width=10, height=13)
    
    save_plt(plts, "moa_discovery-drugs_vs_control-activity_ranking-scatter-GSK3326595", '.pdf', figs_dir, width=15, height=13)
    save_plt(plts, "moa_discovery-drugs_vs_control-activity_ranking-scatter-PALBOCICLIB", '.pdf', figs_dir, width=10, height=13)
    save_plt(plts, "moa_discovery-drugs_vs_control-activity_ranking-scatter-ISOGINKGETIN", '.pdf', figs_dir, width=8, height=13)
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