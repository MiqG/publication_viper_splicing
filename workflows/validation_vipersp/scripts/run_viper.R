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
# RESULTS_DIR = file.path(ROOT,"results","validation_vipersp")

# delta_psi_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"K562",'delta_psi-EX.tsv.gz')
# regulon_file = file.path(RESULTS_DIR,'files','ground_truth_regulon','ENCORE','K562.tsv.gz')

##### FUNCTIONS #####
prep_regulons = function(regulons){
    
    sfs = regulons[['splicing_factor']] %>% unique()
    reg = sapply(sfs, function(sf){
            X = regulons %>%
                filter(splicing_factor %in% sf)
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
        make_option("--delta_psi_file", type="character"),
        make_option("--regulons_file", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    delta_psi_file = args[["delta_psi_file"]]
    regulons_file = args[["regulons_file"]]
    output_file = args[["output_file"]]
    
    # load
    delta_psi = read_tsv(delta_psi_file)
    regulons = read_tsv(regulons_file)
    
    # prep
    delta_psi = delta_psi %>%
        column_to_rownames('EVENT') %>%
        as.matrix()
    
    regulons = prep_regulons(regulons)
    
    # run viper
    result = viper(
        delta_psi, 
        regulons,
        verbose=FALSE
    ) %>% 
    as.data.frame() %>%
    rownames_to_column('splicing_factor')
    
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}