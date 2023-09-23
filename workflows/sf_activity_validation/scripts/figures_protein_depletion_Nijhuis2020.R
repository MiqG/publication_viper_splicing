#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# See how indisulam affects RBM39 protein and transcriptomic expression and how it compares with inferred protein activities.
# 

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

PAL_DARK = "darkgrey"
PAL_ACCENT = "darkred"
PAL_DUAL = c(PAL_DARK, PAL_ACCENT)
PAL_CONTRAST = c("darkgrey","darkred")
PAL_CELL_LINES = "Dark2"
PAL_TIME = "jco"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_validation")
# proteomics_file = file.path(RAW_DIR,"articles","Nijhuis2020","supplementary_data","proteomics_lfq_intensity.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","Nijhuis2020.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","protein_depletion-Nijhuis2020-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Nijhuis2020.tsv.gz")
# gene_info_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","protein_depletion-Nijhuis2020-EX")
# gene_oi = "ENSG00000131051" # RBM39

##### FUNCTIONS #####
plot_proteomics_gene_oi = function(proteomics, gene_oi){
    plts = list()
    
    X = proteomics %>%
        filter(GENE == gene_oi)
    
    plts[["proteomics-indisulam_vs_dmso-box"]] = X %>%
        ggplot(aes(x=condition_lab, y=lfq)) +
        geom_point(aes(color=cell_line_name), position=position_jitter(0.1), size=0.5) +
        geom_boxplot(fill=NA, width=0.25, outlier.shape=NA) +
        color_palette(PAL_CELL_LINES) +
        theme_pubr() +
        facet_wrap(~pert_time_lab, nrow=1) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_compare_means(method="t.test", label="p.signif", size=FONT_SIZE, family=FONT_FAMILY, ref.group="DMSO") + 
        labs(x="Condition", y="log2(LFQ+1)", color="Cell Line", subtitle=gene_oi)

    names(plts) = sprintf("%s-%s",names(plts),gene_oi)
    
    return(plts)
}


plot_diff_proteomics_gene_oi = function(diff_proteomics, gene_oi){
    plts = list()
    
    X = diff_proteomics %>%
        mutate(
            log10_pvalue = -log10(p),
            is_gene_oi = GENE == gene_oi
        )
    
    plts[["diff_proteomics-indisulam_vs_dmso-volcano"]] = X %>%
        ggplot(aes(x=diff_lfq, y=log10_pvalue)) +
        geom_scattermore(data = . %>% filter(!is_gene_oi), pixels=c(1000,1000), pointsize=4, color=PAL_DARK, alpha=0.5) +
        geom_scattermore(data = . %>% filter(is_gene_oi), pixels=c(1000,1000), pointsize=8, color=PAL_ACCENT) +
        theme_pubr() +
        facet_wrap(~pert_time_lab, ncol=2, scales="free_y") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text_repel(
            aes(label=GENE),
            . %>% group_by(pert_time_lab) %>% slice_max(abs(diff_lfq) * log10_pvalue, n=10) %>% ungroup(),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1, max.overlaps=50
        ) +
        labs(x="Log2(Fold Change LFQ)", y="-log10(p-value)", color=sprintf("Is %s", gene_oi))
    
    names(plts) = sprintf("%s-%s",names(plts),gene_oi)
    
    return(plts)
}


plot_genexpr_gene_oi = function(genexpr, gene_oi){
    plts = list()
    
    X = genexpr %>%
        filter(ID == gene_oi)
    
    plts[["genexpr-indisulam_vs_dmso-box"]] = X %>%
        ggplot(aes(x=condition_lab, y=genexpr_tpm, group=pert_time_lab)) +
        geom_point(aes(color=pert_time_lab), position=position_jitter(0.1), size=0.5) +
        geom_line(aes(color=pert_time_lab), size=LINE_SIZE, linetype="dashed") +
        color_palette(PAL_TIME) +
        theme_pubr() +
        facet_wrap(~cell_line_name, nrow=1) +
        theme(aspect.ratio=NULL, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Condition", y="log2(TPM+1)", color="Pert. Time", subtitle=gene_oi)
    
    plts[["genexpr-indisulam_vs_dmso-bar"]] = X %>%
        filter(condition_lab!="INDISULAM\n(1micromolar)") %>%
        ggbarplot(x="condition_lab", y="genexpr_tpm", fill=PAL_ACCENT, color=NA) +
        geom_text(
            aes(label=round(genexpr_tpm,2)),
            vjust=0, hjust=0.5, nudge_y=0.1,
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr() +
        facet_wrap(~cell_line_name+pert_time_lab, nrow=1, scales="free_y") +
        theme(aspect.ratio=NULL, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Condition", y="log2(TPM+1)", color="Pert. Time", subtitle=gene_oi)
    
    names(plts) = sprintf("%s-%s",names(plts),gene_oi)
    
    return(plts)
}


plot_activity_gene_oi = function(protein_activity, gene_oi){
    plts = list()
    
    X = protein_activity %>%
        mutate(is_regulator_oi = regulator == gene_oi)
    
    plts[["activity-indisulam_vs_dmso-ranking-scatter"]] = X %>%
        ggplot(aes(x=activity_ranking, y=activity)) +
        geom_scattermore(data = . %>% filter(!is_regulator_oi), pixels=c(1000,1000), pointsize=4, color=PAL_DARK, alpha=0.5) +
        geom_scattermore(data = . %>% filter(is_regulator_oi), pixels=c(1000,1000), pointsize=8, color=PAL_ACCENT) +
        theme_pubr() +
        facet_wrap(~condition_lab+pert_time_lab, ncol=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Ranking", y="Protein Activity", color=sprintf("Is %s", gene_oi))
    
    names(plts) = sprintf("%s-%s",names(plts),gene_oi)
    
    return(plts)
}


make_plots = function(proteomics, diff_proteomics, genexpr, protein_activity){
    plts = list(
        plot_proteomics_gene_oi(proteomics, "RBM39"),
        plot_diff_proteomics_gene_oi(diff_proteomics, "RBM39"),
        plot_genexpr_gene_oi(genexpr, "ENSG00000131051"),
        plot_activity_gene_oi(protein_activity, "ENSG00000131051")
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(proteomics, diff_proteomics, genexpr, protein_activity){
    figdata = list(
        "validation_drug_target_activity" = list(
            "proteomics" = proteomics,
            "diff_proteomics" = diff_proteomics,
            "genexpr" = genexpr,
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
    save_plt(plts, "proteomics-indisulam_vs_dmso-box-RBM39", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "diff_proteomics-indisulam_vs_dmso-volcano-RBM39", '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, "genexpr-indisulam_vs_dmso-box-ENSG00000131051", '.pdf', figs_dir, width=5.5, height=6)
    save_plt(plts, "genexpr-indisulam_vs_dmso-bar-ENSG00000131051", '.pdf', figs_dir, width=4.5, height=6)
    save_plt(plts, "activity-indisulam_vs_dmso-ranking-scatter-ENSG00000131051", '.pdf', figs_dir, width=9, height=15)
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
        make_option("--proteomics_file", type="character"),
        make_option("--genexpr_file", type="character"),
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    print(args)
    
    proteomics_file = args[["proteomics_file"]]
    genexpr_file = args[["genexpr_file"]]
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    proteomics = read_tsv(proteomics_file)
    genexpr = read_tsv(genexpr_file)
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    gene_info = read_tsv(gene_info_file)
    gc()
    
    # prep
    proteomics = proteomics %>%
        filter(GENE %in% gene_info[["Approved symbol"]]) %>%
        pivot_longer(-GENE, values_to="lfq", names_to="sampleID_rep") %>%
        mutate(lfq = log2(lfq+1)) %>%
        separate(sampleID_rep, sep="_", into=c("sampleID","technical_replicate")) %>%
        mutate(technical_replicate = sprintf("rep%s", technical_replicate)) %>%
        drop_na() %>%
        group_by(sampleID, GENE) %>%
        summarize(lfq = mean(lfq, na.rm=TRUE)) %>%
        ungroup() %>%
        left_join(metadata, by="sampleID") %>%
        mutate(
            condition_lab = ifelse(
                is.na(pert_concentration), condition, 
                sprintf("%s\n(%s%s)", condition, pert_concentration, pert_concentration_units)
            ),
            condition_lab = factor(condition_lab, levels=c("DMSO","INDISULAM\n(5micromolar)")),
            pert_time_lab = sprintf("%s%s", pert_time, pert_time_units),
            pert_time_lab = factor(pert_time_lab, levels=c("6hours","16hours"))
        )
    
    # differential protein levels
    diff_proteomics = proteomics %>% 
        filter(cell_line_name=="IMR32_AUTONOMIC_GANGLIA") %>%
        compare_means(
            lfq ~ condition_lab, 
            data=., 
            group.by=c("pert_time_lab","cell_line_name","GENE"), 
            method="t.test"
        )
    fold_changes = proteomics %>% 
        group_by(pert_time_lab,cell_line_name,GENE,condition_lab) %>% 
        summarize(
            lfq = mean(lfq),
            lfq = ifelse(condition_lab=="DMSO", -lfq, lfq)
        ) %>% 
        ungroup() %>%
        group_by(pert_time_lab,cell_line_name,GENE) %>%
        summarize(
            diff_lfq = sum(lfq),
            diff_lfq = ifelse(diff_lfq<0, -log2(abs(diff_lfq)+1), log2(diff_lfq+1))
        ) %>%
        ungroup() 
    diff_proteomics = diff_proteomics %>%
        left_join(
            fold_changes, 
            by=c("pert_time_lab","cell_line_name","GENE")
        ) %>%
        group_by(pert_time_lab,cell_line_name) %>%
        mutate(
            fdr = p.adjust(p, method="fdr"),
            is_target = GENE=="RBM39"
        ) %>%
        ungroup()
    
    genexpr = genexpr %>%
        pivot_longer(-ID, names_to="sampleID", values_to="genexpr_tpm") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition) %>%
        mutate(
            condition_lab = ifelse(
                is.na(pert_concentration), condition, 
                sprintf("%s\n(%s%s)", condition, pert_concentration, pert_concentration_units)
            ),
            condition_lab = factor(
                condition_lab, levels=c("DMSO","INDISULAM\n(1micromolar)","INDISULAM\n(5micromolar)")
            ),
            pert_time_lab = sprintf("%s%s", pert_time, pert_time_units),
            pert_time_lab = factor(pert_time_lab, levels=c("6hours","16hours"))
        )
    
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition) %>%
        mutate(
            condition_lab = ifelse(
                is.na(pert_concentration), condition, 
                sprintf("%s\n(%s%s)", condition, pert_concentration, pert_concentration_units)
            ),
            condition_lab = factor(
                condition_lab, levels=c("DMSO","INDISULAM\n(1micromolar)","INDISULAM\n(5micromolar)")
            ),
            pert_time_lab = sprintf("%s%s", pert_time, pert_time_units),
            pert_time_lab = factor(pert_time_lab, levels=c("6hours","16hours"))
        ) %>%
        group_by(condition_lab, pert_time_lab) %>%
        arrange(activity) %>%
        mutate(activity_ranking = row_number()) %>%
        ungroup()
    
    
    # plot
    plts = make_plots(proteomics, diff_proteomics, genexpr, protein_activity)
    gc()
    
    # make figdata
    figdata = make_figdata(proteomics, diff_proteomics, genexpr, protein_activity)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}