require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_EVAL = setNames(c("#EB9486","#7E7F9A"), c(TRUE, FALSE))

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# experimental_pruned_path = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX")
# aracne_and_experimental_path = file.path(RESULTS_DIR,"files","aracne_and_experimental_regulons-EX")
# mlr_and_experimental_path = file.path(RESULTS_DIR,"files","mlr_and_experimental_regulons-EX")
# figs_dir = file.path(RESULTS_DIR,"figures","inference_troubleshooting")

##### FUNCTIONS #####
load_networks = function(network_path, patt=NULL){
    if (file.exists(network_path) && !dir.exists(network_path)){
        # network_path is a file, we load only that network (we'll tun regular VIPER)
        network_files = list(network_path)
    }else if (dir.exists(network_path)){
        # network_path is a directory, we load all networks contained (we'll run metaVIPER)
        network_files = list.files(network_path, pattern=patt, full.names=TRUE)
    }else {
        stop("Invalid network_path.")
    }
    
    networks = sapply(network_files, function(network_file){
        network = read_tsv(network_file) %>%
            mutate(network_id = network_file)
        return(network)
    }, simplify=FALSE) %>% bind_rows()
    
    return(networks)
}


plot_inference_evaluation = function(eval_likelihood, eval_tfmode){
    plts = list()
    
    # Are likelihoods predictive of true interactions?
    X = eval_likelihood
    
    plts[["inference_eval-likelihood-box"]] = X %>% 
        ggboxplot(x="in_ground_truth", y="likelihood", fill="in_ground_truth", 
                  palette=PAL_EVAL, width=0.5, outlier.size=0.1) + 
        geom_text(
            aes(y=0.1, label=label),
            . %>% count(method, in_ground_truth) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~method, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="In Ground Truth", y="Interaction Likelihood")
    
    # Do MoR correlate?
    X = eval_tfmode
    
    plts[["inference_eval-tfmode-bar"]] = X %>% 
        ggbarplot(x="tfmode_experimental", y="perc", fill="label", color=NA, palette=PAL_EVAL) +
        facet_wrap(~method, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="MoR Ground Truth", y="Percentage", fill="Correct MoR Prediction")
    
    return(plts)
}


make_plots = function(eval_likelihood, eval_tfmode){
    plts = list(
        plot_inference_evaluation(eval_likelihood, eval_tfmode)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(eval_likelihood, eval_tfmode){
    figdata = list(
        "inference_evaluation" = list(
            "eval_likelihood" = eval_likelihood,
            "eval_tfmode" = eval_tfmode
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
    save_plt(plts, "inference_eval-likelihood-box", '.pdf', figs_dir, width=5, height=5.5)
    save_plt(plts, "inference_eval-tfmode-bar", '.pdf', figs_dir, width=5, height=7)
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
        make_option("--experimental_pruned_path", type="character"),
        make_option("--aracne_and_experimental_path", type="character"),
        make_option("--mlr_and_experimental_path", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    experimental_pruned_path = args[["experimental_pruned_path"]]
    aracne_and_experimental_path = args[["aracne_and_experimental_path"]]
    mlr_and_experimental_path = args[["mlr_and_experimental_path"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    experimental_pruned = load_networks(experimental_pruned_path)
    aracne_and_experimental = load_networks(aracne_and_experimental_path)
    mlr_and_experimental = load_networks(mlr_and_experimental_path)
    
    # prep
    experimental_pruned = experimental_pruned %>%
        distinct(regulator, target) %>%
        mutate(in_ground_truth = TRUE)
        # 121071 total interactions
    
    aracne_and_experimental = aracne_and_experimental %>%
        left_join(experimental_pruned, by=c("regulator","target")) %>%
        mutate(
            in_ground_truth = replace_na(in_ground_truth, FALSE),
            tfmode_aracne = sign(tfmode_aracne)
        )
    
    mlr_and_experimental = mlr_and_experimental %>%
        left_join(experimental_pruned, by=c("regulator","target")) %>%
        mutate(in_ground_truth = replace_na(in_ground_truth, FALSE))
    
    # evaluate
    ## Are likelihoods predictive of true interactions?
    eval_likelihood = aracne_and_experimental %>%
        distinct(regulator, target, likelihood, in_ground_truth) %>%
        mutate(method = "ARACNe") %>%
        bind_rows(
            mlr_and_experimental %>%
            distinct(regulator, target, likelihood, in_ground_truth) %>%
            mutate(method = "MLR")
        )
    
    ## Can we predict correct MoR?
    eval_tfmode = aracne_and_experimental %>%
        dplyr::rename(tfmode_method = tfmode_aracne) %>%
        distinct(regulator, target, tfmode_experimental, tfmode_method) %>%
        count(tfmode_experimental, tfmode_method) %>%
        group_by(tfmode_experimental) %>%
        mutate(perc = n / sum(n)) %>%
        ungroup() %>%
        mutate(method = "ARACNe") %>%
        bind_rows(
            mlr_and_experimental %>%
            dplyr::rename(tfmode_method = tfmode_mlr) %>%
            distinct(regulator, target, tfmode_experimental, tfmode_method) %>%
            count(tfmode_experimental, tfmode_method) %>%
            group_by(tfmode_experimental) %>%
            mutate(perc = n / sum(n)) %>%
            ungroup() %>%
            mutate(method = "MLR")
        ) %>%
        mutate(
            label = tfmode_experimental==tfmode_method
        )
    
    # plot
    plts = make_plots(eval_likelihood, eval_tfmode)
    
    # make figdata
    figdata = make_figdata(eval_likelihood, eval_tfmode)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}