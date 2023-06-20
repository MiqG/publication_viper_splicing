#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# make figures of evaluate target inference algorithms.
# median accuracy of 0.973 and 0.975 in HepG2 and K562 with threshold of 0.2

require(optparse)
require(tidyverse)
require(qep)

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,'data','prep')
# matrix_file = file.path(PREP_DIR,'genexpr_tpm','LGG.tsv.gz')

##### FUNCTIONS #####
quantize <- function(data, qscheme,
                     decreasing=T, verbose=T)
{
    
    if(decreasing) {
        data <- apply(-data,2,rank,ties.method="random")
        } else data <- apply(data,2,rank,ties.method="random")

    data <- (data-1)/(nrow(data)-1)


    if(!is.null(dim(qscheme))) {
        if(abs(sum(qscheme)-ncol(qscheme))>(10^-15)*ncol(qscheme))
            stop("All schemes must sum to 1")
        if(ncol(qscheme)!=ncol(data))
            stop("data and scheme must have same number of columns")
        binned <- matrix(NA,nrow(data),ncol(data))
        for(i in 1:ncol(data)) {
          cs <- cumsum(qscheme[,i])
          ## make sure 1 is 1 and not 1-10^16...    
          cs[length(cs)] <- 1
          binned[,i] <- cut(data[,i], c(0,cs), F, T)
        }
    } else {
      cs <- cumsum(qscheme)
      ## make sure 1 is 1 and not 1-10^16...    
      cs[length(cs)] <- 1
      binned <- apply(data, 2, cut, c(0,cs), F, T)
    }

    rownames(binned) <- rownames(data)
    colnames(binned) <- colnames(data)
    
    return(binned) # simplified output
}


discretize = function(X){
    # get binning scheme (defaults from paper)
    # m = 36, s = 14.31, w = −9.39, v = 0.04.
    qscheme = logis.qs(nbins=36, steep=14.31, width=-9.39, baseline=0.04)
    
    # bin each sample
    bins = quantize(X, qscheme)
    
    return(bins)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--matrix_file", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    matrix_file = args[["matrix_file"]]
    output_file = args[["output_file"]]
    
    # load
    matrix = read_tsv(matrix_file)
    matrix = matrix %>%
        column_to_rownames(colnames(matrix)[1]) %>%
        as.matrix()
    
    # discretize
    result = discretize(matrix)
    result = result %>%
        as.data.frame() %>%
        rownames_to_column("index")
    
    # save
    write_tsv(result, output_file)
}

##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}