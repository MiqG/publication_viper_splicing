# Note
# ----
# WE ARE RUNNING METAVIPER!

require(optparse)
require(tidyverse)
require(viper)

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# signature_file = file.path(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-EX.tsv.gz')
# regulons_path = file.path(RESULTS_DIR,"files","mlr_regulons_development-EX")
# eval_labels_file = file.path(RESULTS_DIR,"files","regulon_evaluation_labels","ENCOREKO_HepG2.tsv.gz")

##### FUNCTIONS #####
as_regulon_network = function(regulons){

    uregs = regulons[['regulator']] %>% unique()
    reg = sapply(uregs, function(ureg){
            X = regulons %>%
                filter(regulator %in% ureg)
            X = list(
                tfmode = setNames(X[['tfmode']], X[['target']]),
                likelihood = X[['likelihood']]
            )
            return(X)
        }, simplify=FALSE)
    
    return(reg)
}


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
        network = read_tsv(network_file)
        network = as_regulon_network(network)
        return(network)
    }, simplify=FALSE)
    
    return(networks)
}


run_viper = function(signature, regulons){
    # runs VIPER or metaVIPER depending on whether there are multiple regulons
    # in `regulons`
    
    result = viper(
        signature, 
        regulons,
        verbose=FALSE
    ) 
    
    return(result)
}


evaluate_protein_activities = function(protein_activities, eval_labels){
    # prep
    X = protein_activities %>%
        as.data.frame() %>%
        rownames_to_column("regulator") %>%
        pivot_longer(-regulator, names_to="PERT_ID", values_to="activity") %>%
        left_join(eval_labels, by="PERT_ID") %>%
        # change sign of protein activities from OVEREXPRESSION signatures
        mutate(activity = ifelse(PERT_TYPE=="OVEREXPRESSION",-activity,activity)) %>%
        drop_na()
    
    # Can we gess in which perturbations was a regulator perturbed?
    evaluation_between_real = X %>%
        group_by(regulator) %>%
        arrange(-activity) %>% # high score better
        mutate(
            ranking_raw = row_number(),
            ranking_perc = ranking_raw / n()
        ) %>%
        filter(regulator == PERT_ENSEMBL) %>%
        mutate(
            eval_direction = "between",
            eval_type = "real"
        )
    
    evaluation_between_random = X %>%
        group_by(regulator) %>%
        mutate(activity = sample(activity)) %>%
        arrange(-activity) %>% # high score better
        mutate(
            ranking_raw = row_number(),
            ranking_perc = ranking_raw / n()
        ) %>%
        filter(regulator == PERT_ENSEMBL) %>%
        mutate(
            eval_direction = "between",
            eval_type = "random"
        )
    
    # Can we gess which regulator was perturbed in each perturbation?
    evaluation_within_real = X %>%
        group_by(PERT_ID) %>%
        arrange(-activity) %>% # high score better
        mutate(
            ranking_raw = row_number(),
            ranking_perc = ranking_raw / n()
        ) %>%
        filter(regulator == PERT_ENSEMBL) %>%
        mutate(
            eval_direction = "within",
            eval_type = "real"
        )
    
    evaluation_within_random = X %>%
        group_by(PERT_ID) %>%
        mutate(activity = sample(activity)) %>%
        arrange(-activity) %>% # high score better
        mutate(
            ranking_raw = row_number(),
            ranking_perc = ranking_raw / n()
        ) %>%
        filter(regulator == PERT_ENSEMBL) %>%
        mutate(
            eval_direction = "within",
            eval_type = "random"
        )
    
    # merge
    evaluation = evaluation_between_real %>%
        bind_rows(evaluation_between_random) %>%
        bind_rows(evaluation_within_real) %>%
        bind_rows(evaluation_within_random)
    
    return(evaluation)
}


run_viper_and_evaluate = function(signature, regulons, eval_labels){
    # run viper for each regulon set
    result = lapply(names(regulons), function(regulons_oi){
        # compute protein activities
        protein_activities = viper(signature, regulons[[regulons_oi]], verbose=FALSE)
        
        # evaluate protein activities
        evaluation = evaluate_protein_activities(protein_activities, eval_labels)
        
        # add info
        evaluation[["regulon_id"]] = basename(regulons_oi) %>% gsub(".tsv.gz","",.)
        
        return(evaluation)
        
    }) %>% do.call(rbind, .)
    
    return(result)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character"),
        make_option("--regulons_path", type="character"),
        make_option("--eval_labels_file", type="character", default=NULL),
        make_option("--output_file", type="character"),
        make_option("--random_seed", type="integer", default=1234)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    regulons_path = args[["regulons_path"]]
    eval_labels_file = args[["eval_labels_file"]]
    output_file = args[["output_file"]]
    
    set.seed(args[["random_seed"]])
    
    # load
    signature = read_tsv(signature_file)
    regulons = load_networks(regulons_path)
    if (!is.null(eval_labels_file)){
        eval_labels = read_tsv(eval_labels_file)
    }else{
        eval_labels = eval_labels_file
    }
    
    # prep
    ## signature
    signature = signature %>% as.data.frame()
    rownames(signature) = signature[,1]
    signature = signature[,2:ncol(signature)]
    signature = signature %>% 
        dplyr::select(where(is.numeric))
    
    if (is.null(eval_labels)){
        
        # run regular viper
        result = run_viper(signature, regulons)
        result = result %>% as.data.frame() %>% rownames_to_column('regulator')
    
    }else{
        print("Evaluation mode...")
        
        # run viper and evaluate predicted protein activities
        result = run_viper_and_evaluate(signature, regulons, eval_labels)
        result[["regulon_set_id"]] = basename(regulons_path)
        result[["signature_id"]] = basename(eval_labels_file) %>% gsub(".tsv.gz","",.)
    }
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
