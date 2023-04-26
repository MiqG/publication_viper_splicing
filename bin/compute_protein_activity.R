# Note
# ----
# WE ARE RUNNING METAVIPER!

require(optparse)
require(tidyverse)
require(viper)

# Development
# -----------
# ROOT = here::here()
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_ccle")
# ACTVAL_DIR = file.path(ROOT,"results","validation_activity")
# signature_file = file.path(RESULTS_DIR,"files","signatures","CCLE-EX.tsv.gz")
# regulons_dir = file.path(ACTVAL_DIR,"files","subsetted_regulons","regulons_selected")

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

load_networks = function(networks_dir, patt=NULL){
    
    network_files = list.files(networks_dir, pattern=patt, full.names=TRUE)
    networks = sapply(network_files, function(network_file){
        network = read_tsv(network_file)
        network = as_regulon_network(network)
        return(network)
    }, simplify=FALSE)
    
    return(networks)
}

run_viper = function(signature, regulons){
    
    result = viper(
        signature, 
        regulons,
        verbose=FALSE
    ) 
    
    return(result)
}

parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character"),
        make_option("--regulons_dir", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    regulons_dir = args[["regulons_dir"]]
    output_file = args[["output_file"]]
    
    # load
    signature = read_tsv(signature_file)
    regulons = load_networks(regulons_dir)
    
    # prep
    ## signature
    signature = signature %>% as.data.frame()
    rownames(signature) = signature[,1]
    signature = signature[,2:ncol(signature)]
    signature = signature %>% 
        dplyr::select(where(is.numeric))
    
    # run
    result = run_viper(signature, regulons)
    
    # save
    result = result %>% as.data.frame() %>% rownames_to_column('regulator')
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
