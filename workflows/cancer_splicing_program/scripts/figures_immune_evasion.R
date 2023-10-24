#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# TODO
# ----
# - check which exons in the selected genes correlate with patient survival

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)
require(ggrepel)
require(survival)
require(survminer)

# variables
RANDOM_SEED = 1234
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

GENES_CALR = c('ENSG00000084463','ENSG00000021776','ENSG00000189091','ENSG00000215301')
GENES_TAPBP = c(
    'ENSG00000100056','ENSG00000125352',
    'ENSG00000139793',
    'ENSG00000174231','ENSG00000143368','ENSG00000160201','ENSG00000165119',
    'ENSG00000183431','ENSG00000197111','ENSG00000115524',
    'ENSG00000131013','ENSG00000074201','ENSG00000104897','ENSG00000189091','ENSG00000179950'
)

EXONS_OI = c(
    "HsaEX0063618", # TAPBP
    "HsaEX6093277" # CALR
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","cancer_splicing_program")
# splicing_file = file.path(PREP_DIR,"event_psi","Riaz2017-PRE_CombPD1_CTLA4-EX.tsv.gz")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","Riaz2017-PRE-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Riaz2017.tsv.gz")
# driver_types_file = file.path(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","immune_evasion")

##### FUNCTIONS #####
# plot kaplan meier of median protein activity of tumor suppressor, oncogenic, and regulators of strong immune evasion genes
plot_survival = function(protein_activity){
    plts = list()
    
    X = 
    
    
    return(plts)
}


make_plots = function(
    protein_activiy
){
    plts = list(
        plot_tumorigenesis(protein_activiy)
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
    save_plt(plts, "tumorigenesis-cell_line_vs_activity-violin", '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, "tumorigenesis-cell_line_vs_activity-random-violin", '.pdf', figs_dir, width=6, height=6)
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
    splicing = read_tsv(splicing_file)
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    driver_types = read_tsv(driver_types_file)
    
    # prep
    splicing = splicing %>%
        filter(EVENT %in% EXONS_OI) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, psi)        
    
    protein_activity = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(metadata, by="sampleID") %>%
        drop_na(condition, activity) %>%
        mutate(
            condition_lab = sprintf(
                "%s (%s | %s) | %s", condition, sampleID, patientID, study_accession
            )
        ) %>%
        # add activity
        group_by(condition_lab) %>%
        arrange(activity) %>%
        mutate(
            abs_activity = abs(activity),
            activity_ranking = row_number()
        ) %>%
        arrange(abs_activity) %>%
        mutate(
            abs_activity_ranking = row_number(),
        ) %>%
        ungroup() %>%
        left_join(driver_types, by=c("regulator"="ENSEMBL"))
    
    # survival analysis
    for (event_oi in EXONS_OI){
        x = splicing %>%
            filter(EVENT %in% event_oi) %>%
            surv_cutpoint(time="OS_time", event="OS_event", variables="psi", minprop=0.05) %>%
            surv_categorize()

        fit = survfit(Surv(OS_time, OS_event) ~psi, data=x)
        plt = x %>%
            ggsurvplot(
                fit, data=., risk.table=TRUE, conf.int=TRUE, pval=TRUE, pval.size=FONT_SIZE+2,
                risk.table.fontsize=FONT_SIZE+2, risk.table.font.family=FONT_FAMILY,
                palette = get_palette("Dark2", 2)
            ) + labs(subtitle=event_oi)
        print(plt)
    }
    
    
    for (gene_oi in GENES_TAPBP){
        x = protein_activity %>%
            filter(regulator %in% gene_oi) %>%
            drop_na(driver_type) %>%
            group_by(driver_type, sampleID, patientID, OS_time, OS_event) %>%
            summarize(activity = median(activity)) %>%
            ungroup() %>%
            filter(driver_type == "Tumor suppressor") %>%
            surv_cutpoint(time="OS_time", event="OS_event", variables="activity") %>%
            surv_categorize()

        fit = survfit(Surv(OS_time, OS_event) ~activity, data=x)
        plt = x %>%
            ggsurvplot(
                fit, data=., risk.table=TRUE, conf.int=TRUE, pval=TRUE, pval.size=FONT_SIZE+2,
                risk.table.fontsize=FONT_SIZE+2, risk.table.font.family=FONT_FAMILY,
                palette = get_palette("Dark2", 2)
            ) + labs(subtitle=gene_oi)
        print(plt)
    }
    
    
    # plot
    plts = make_plots(protein_activity)
    
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