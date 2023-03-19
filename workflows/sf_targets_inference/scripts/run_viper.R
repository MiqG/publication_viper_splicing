#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Run VIPER using splicing factor regulons.

require(optparse)
require(viper)
require(tidyverse)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","sf_targets_inference")

# signature_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"K562",'log2_fold_change_tpm.tsv.gz')
# regulons_file = file.path(RESULTS_DIR,'files','target_inference','aracne','genexpr_tpm','CCLE.tsv.gz')

##### FUNCTIONS #####
prep_regulons = function(regulons){
    
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

parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character"),
        make_option("--regulons_file", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    regulons_file = args[["regulons_file"]]
    output_file = args[["output_file"]]
    
    # load
    signature = read_tsv(signature_file)
    regulons = read_tsv(regulons_file)
    
    # prep
    signature = signature %>% as.data.frame()
    rownames(signature) = signature[,1]
    signature = signature[,2:ncol(signature)]
    signature = signature %>% 
        dplyr::select(where(is.numeric))
    
    regulons = prep_regulons(regulons)
    
    # run viper
    result = viper(
        signature, 
        regulons,
        verbose=FALSE
    ) %>% 
    as.data.frame() %>%
    rownames_to_column('regulator')
    
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
