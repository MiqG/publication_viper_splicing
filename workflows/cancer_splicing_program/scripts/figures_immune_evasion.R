#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Outline
# ----
# 1. Differential inclusion in responders (CR/PR) vs non-responders (PD)
# 2. Which differentially included exons affect functionally important genes?
# 3. Does inclusion of prioritized exon(s) correlate with patient survival?

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)
require(survival)
require(survminer)
require(readxl)

# variables
RANDOM_SEED = 1234
THRESH_FDR = 0.05

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
# RESULTS_DIR = file.path(ROOT,"results","cancer_splicing_program")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","Riaz2017-PRE.tsv.gz")
# splicing_file = file.path(PREP_DIR,"event_psi","Riaz2017-PRE-EX.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","Riaz2017-PRE-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Riaz2017.tsv.gz")
# driver_types_file = file.path(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","immune_evasion")
# enrichments_reactome_file = file.path(RESULTS_DIR,"figures","cancer_program","figdata","cancer_program","enrichments_reactome.tsv.gz")
# immune_screen_file = file.path(SUPPORT_DIR,"supplementary_tables_literature","Dubrot2022-suptabs-41590_2022_1315_MOESM2_ESM.xlsx") # Sup. Tab. 13
# human2mouse_file = file.path(RAW_DIR,"BIOMART","human2mouse.tsv")
# survival_analysis_file = file.path(RESULTS_DIR,'files',"survival_analysis",'splicing-EX-Riaz2017-PRE-surv.tsv.gz')
# annotation_file = file.path(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# REGINF_DIR = file.path(ROOT,"results","regulon_inference")
# regulons_path = file.path(REGINF_DIR,"files","experimentally_derived_regulons_pruned-EX")


##### FUNCTIONS #####
load_regulons = function(regulons_path, patt=NULL){
    if (file.exists(regulons_path) && !dir.exists(regulons_path)){
        # regulons_path is a file, we load only that regulon (we'll tun regular VIPER)
        regulon_files = list(regulons_path)
    }else if (dir.exists(regulons_path)){
        # regulons_path is a directory, we load all regulons contained (we'll run metaVIPER)
        regulon_files = list.files(regulons_path, pattern=patt, full.names=TRUE)
    }else {
        stop("Invalid regulons_path.")
    }
    
    regulons = sapply(regulon_files, function(regulon_file){
        regulon = read_tsv(regulon_file)
        return(regulon)
    }, simplify=FALSE)
    
    return(regulons)
}


plot_diff_response = function(diff_response, immune_screen, splicing){
    plts = list()
    
    X = diff_response
    
    plts[["diff_response-median_diff_vs_pvalue-scatter"]] = X %>%
        ggplot(aes(x=median_diff, y=-log10(p))) +
        geom_scattermore(pixels=c(1000,1000), pointsize=15, alpha=0.5) +
        geom_hline(yintercept=-log10(0.05), size=LINE_SIZE, linetype="dashed", color="black") +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Median PSI difference Responder (CR/PR) vs Non-responder (PD)", y="-log10(p-value)")
        
    X = diff_response %>%
        left_join(immune_screen, by=c("GENE"="human_symbol")) %>%
        mutate(
            label = sprintf("%s_%s (%s)", EVENT, GENE, Gene),
            term_clean = replace_na(term_clean,"Unknown")
        ) %>%
        filter(Comparison=="ICB vs NSG" & p<0.05)
    
    events_oi = X %>% 
        group_by(Comparison,Dataset) %>% 
        slice_max(abs(score*median_diff), n=6) %>%
        ungroup()
    
    plts[["diff_response-median_diff_vs_immune_screen_score-scatter"]] = X %>%
        ggplot(aes(x=median_diff, y=score)) +
        geom_scattermore(aes(color=term_clean), pixels=c(1000,1000), pointsize=15, alpha=0.8) +
        geom_hline(yintercept=0, size=LINE_SIZE, linetype="dashed", color="black") +
        geom_vline(xintercept=0, size=LINE_SIZE, linetype="dashed", color="black") +
        color_palette("jco") +
        geom_text_repel(
            aes(label=label), events_oi ,
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1
        ) +
        theme_pubr() +
        facet_wrap(~Comparison+Dataset) +
        theme(aspect.ratio=1) +
        labs(
            x="Median PSI difference Responder (CR/PR) vs Non-responder (PD)", 
            y="Fitness Score Gene KO",
            color="Exon Impact"
        )
        
    EXONS_OI = events_oi %>% slice_max(abs(score), n=1) %>% pull(EVENT)
    for (event_oi in EXONS_OI){
        x = splicing %>%
            filter(EVENT %in% event_oi) %>%
            surv_cutpoint(time="OS_time", event="OS_event", variables="psi") 
        
        cutpoint = x[["cutpoint"]][["cutpoint"]]
        x = x %>% surv_categorize()

        fit = survfit(Surv(OS_time, OS_event) ~psi, data=x)
        
        gene_oi = splicing %>% filter(EVENT %in% event_oi) %>% pull(GENE) %>% unique()
        event_gene = sprintf("%s_%s",event_oi,gene_oi)
        plts[[sprintf("diff_response-psi_vs_survival-km-%s",event_gene)]] = x %>%
            ggsurvplot(
                fit, data=., risk.table=TRUE, conf.int=TRUE, pval=TRUE, pval.size=FONT_SIZE+2,
                risk.table.fontsize=FONT_SIZE+2, risk.table.font.family=FONT_FAMILY,
                palette = get_palette("Dark2", 2)
            ) + labs(title=sprintf("Best Cutpoint: PSI=%s",cutpoint), subtitle=event_gene)
        plts[[sprintf("diff_response-psi_vs_survival-km-%s",event_gene)]] = plts[[sprintf("diff_response-psi_vs_survival-km-%s",event_gene)]][["plot"]]
    }
    
    return(plts)
}


make_plots = function(
    diff_response, immune_screen, splicing
){
    plts = list(
        plot_diff_response(diff_response, immune_screen, splicing)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    diff_response, immune_screen, splicing
){
    figdata = list(
        "tumorigenesis" = list(
            "splicing" = splicing,
            "diff_response" = diff_response, 
            "immune_screen" = immune_screen
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
    save_plt(plts, "diff_response-median_diff_vs_pvalue-scatter", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "diff_response-median_diff_vs_immune_screen_score-scatter", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "diff_response-psi_vs_survival-km-HsaEX1036341_SEC22B", '.pdf', figs_dir, width=5, height=7)
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
        make_option("--metadata_file", type="character"),
        make_option("--driver_types_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    driver_types_file = args[["driver_types_file"]]
    figs_dir = args[["figs_dir"]]
    
    set.seed(RANDOM_SEED)    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    annot = read_tsv(annotation_file)
    regulons = load_regulons(regulons_path)
    splicing = read_tsv(splicing_file)
    metadata = read_tsv(metadata_file)
    enrichments_reactome = read_tsv(enrichments_reactome_file)
    immune_screen = read_excel(immune_screen_file, sheet="Supplementary Table 13")
    human2mouse = read_tsv(
        human2mouse_file, col_names=c("human_ensembl","human_symbol","mouse_ensembl","mouse_symbol")
    )
    survival_analysis = read_tsv(survival_analysis_file)
    protein_impact = read_tsv(protein_impact_file) %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT) %>%
            mutate(term_clean=gsub(" \\(.*","",term),
                   term_clean=gsub("ORF disruption upon sequence exclusion",
                                   "ORF disruption (exclusion)",term_clean),
                   term_clean=gsub("ORF disruption upon sequence inclusion",
                                   "ORF disruption (inclusion)",term_clean),
                   term_clean=gsub("In the CDS, with uncertain impact",
                                   "In the CDS (uncertain)",term_clean),
                   term_clean=replace_na(term_clean, "Unknown")
            )
    # prep
    regulons = regulons %>% 
        bind_rows() %>%
        distinct(regulator, target)
    
    enrichments_reactome = enrichments_reactome %>%
                filter(str_detect(ID, "ANTIGEN")) %>%
                separate_rows(geneID)
    
    annot = annot %>%
        mutate(
            in_regulons = EVENT %in% regulons[["target"]],
            in_enrichment = GENE %in% enrichments_reactome[["geneID"]]
        ) %>%
        left_join(protein_impact, by="EVENT")
    
    metadata = metadata %>%
        mutate(
            is_responder = NA,
            is_responder = case_when(
                treatment_response %in% c("PR","CR") ~ "Responder",
                treatment_response == "PD" ~ "Non-responder"
            )
        )
    
    splicing = splicing %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi") %>%
        drop_na(psi) %>%
        left_join(annot, by="EVENT") %>%
        left_join(metadata, by="sampleID") %>%
        filter(in_enrichment & in_regulons)
    
    immune_screen = immune_screen %>%
        # translate mouse genes to human
        left_join(
            human2mouse %>% drop_na(human_symbol, mouse_symbol), 
            by=c("Gene"="mouse_symbol"),
            relationship = "many-to-many" # WARNING! There are duplicates!
        ) %>%
        mutate(score = Sign*`Average Score`) %>% 
        distinct(Gene, human_symbol, score, Dataset, Comparison) %>%
        drop_na(human_symbol)    
    
    # plot
    plts = make_plots(diff_response, immune_screen, splicing)
    
    # make figdata
    figdata = make_figdata(diff_response, immune_screen, splicing)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}