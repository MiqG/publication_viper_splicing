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
require(aracne.networks)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","tfs_splicing")

# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","PAAD.tsv.gz")
# cancer_type = "PAAD"
# gene_annotation_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")

##### FUNCTIONS #####
parseargs = function(){
    
    option_list = list( 
        make_option("--genexpr_file", type="character"),
        make_option("--gene_annotation_file", type="character"),
        make_option("--cancer_type", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    genexpr_file = args[["genexpr_file"]]
    gene_annotation_file = args[["gene_annotation_file"]]
    cancer_type = args[["cancer_type"]]
    output_file = args[["output_file"]]
    
    # load
    genexpr = read_tsv(genexpr_file)
    regulons_name = paste0("regulon",tolower(cancer_type))
    gene_annotation = read_tsv(gene_annotation_file) %>%
        rename(ensembl=`Ensembl gene ID`, entrez=`NCBI Gene ID`, symbol=`Approved symbol`) %>%
        dplyr::select(ensembl, entrez, symbol) %>%
        mutate(entrez = as.character(entrez))
    
    # prep
    genexpr = genexpr %>% 
        left_join(gene_annotation %>% dplyr::select(ensembl,entrez), by=c("ID"="ensembl")) %>%
        dplyr::select(-ID) %>%
        drop_na(entrez) %>%
        as.data.frame()
    rownames(genexpr) = genexpr[,ncol(genexpr)]
    genexpr = genexpr %>% 
        dplyr::select(where(is.numeric))
    
    regulons = get(regulons_name)
    
    print(regulons_name)
    
    # run viper
    result = viper(
        genexpr, 
        regulons,
        verbose=FALSE,
        method="scale"
    ) %>% 
    as.data.frame() %>%
    rownames_to_column('regulator')
    
    # annotate genes
    result = result %>%
        left_join(gene_annotation, by=c("regulator"="entrez"))
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
