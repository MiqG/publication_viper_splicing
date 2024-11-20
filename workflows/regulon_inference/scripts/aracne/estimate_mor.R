Sys.setenv(VROOM_CONNECTION_SIZE = 500000)
require(optparse)
require(tidyverse)
require(viper)

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# splicing_file = file.path(PREP_DIR,'genexpr_tpm','CardosoMoreira2020.tsv.gz')
# genexpr_file = file.path(PREP_DIR,'event_psi_imputed','CardosoMoreira2020-EX.tsv.gz')
# aracne_network_file = file.path(RESULTS_DIR,"files","aracne_regulons","CardosoMoreira2020-EX","pruned","network.txt")

##### FUNCTIONS #####
readAracneAdj = function (afile) {
    aracne = read_tsv(afile)
    aracne = aracne %>%
        dplyr::rename(
            tf=Regulator,
            target=Target,
            mi=MI
        ) %>%
        dplyr::select(tf,target,mi) %>%
        as.data.frame() 
    return(aracne)
}


aracne2regulon = function (afile, eset, gene = FALSE, format = c("adj", "3col"), 
    verbose = TRUE) 
{
    format <- match.arg(format)
    if (is(eset, "ExpressionSet")) 
        eset <- exprs(eset)
    if (verbose) 
        message("\nLoading the dataset...")
    if (length(eset) == 1) {
        tmp <- strsplit(readLines(eset), "\t")
        dset <- t(sapply(tmp[-1], function(x) as.numeric(x[-(1:2)])))
        colnames(dset) <- tmp[[1]][-(1:2)]
        rownames(dset) <- sapply(tmp[-1], function(x) x[1])
        annot <- t(sapply(tmp[-1], function(x) x[1:2]))
    } else {
        dset <- eset
        annot <- rownames(eset)
        names(annot) <- rownames(eset)
        rm(eset)
    }
    switch(format, adj = {
        aracne <- readAracneAdj(afile)
    }, `3col` = {
        tmp <- t(sapply(strsplit(readLines(afile), "\t"), function(x) x[1:3]))
        aracne <- data.frame(tf = tmp[, 1], target = tmp[, 2], 
            mi = as.numeric(tmp[, 3])/max(as.numeric(tmp[, 3])))
    })
    if (gene) {
        if (verbose) 
            message("Collapsing the interactomes to the gene level...")
        tmp <- aracne[order(aracne$mi, decreasing = TRUE), ]
        tmp$tf <- annot[match(tmp$tf, annot[, 1]), 2]
        tmp$target <- annot[match(tmp$target, annot[, 1]), 2]
        aracne <- tmp[!duplicated(paste(tmp$tf, tmp$target, sep = "_")), 
            ]
        rownames(dset) <- annot[match(rownames(dset), annot[, 
            1]), 2]
        dset <- filterCV(dset)
    }
    if (verbose) 
        message("Generating the regulon objects...")
    tmp <- aracne[!is.na(aracne$mi), ]
    tmp <- tmp[rowSums(matrix(as.matrix(tmp[, 1:2]) %in% rownames(dset), 
        nrow(tmp), 2)) == 2, ]
    aracne <- tapply(1:nrow(tmp), as.vector(tmp$tf), function(pos, 
        tmp) {
        tfmode <- rep(0, length(pos))
        names(tfmode) <- tmp$target[pos]
        list(tfmode = tfmode, likelihood = tmp$mi[pos])
    }, tmp = tmp)
    names(aracne) <- levels(factor(as.vector(tmp$tf)))
    aracne <- viper:::TFmode1(aracne, dset)
    rm(dset)
    aracne <- aracne[names(aracne) != "NA"]
    aracne <- lapply(aracne, function(x) {
        filtro <- !(names(x$tfmode) == "NA" | is.na(x$tfmode) | 
            is.na(x$likelihood))
        x$tfmode <- x$tfmode[filtro]
        x$likelihood <- x$likelihood[filtro]
        return(x)
    })
    aracne <- aracne[sapply(aracne, function(x) length(names(x$tfmode))) > 
        0]
    regul <- viper:::TFscore(aracne, verbose = verbose)
    class(regul) <- "regulon"
    return(regul)
}
                            
regulons_to_edgelist = function(regulons){
    edgelist = lapply(names(regulons), function(regulator_oi){
        
        tfmode = regulons[[regulator_oi]][["tfmode"]]
        likelihood = regulons[[regulator_oi]][["likelihood"]]
        
        edges = data.frame(
            "regulator" = regulator_oi,
            "target" = names(tfmode),
            "likelihood" = likelihood,
            "tfmode" = tfmode
        )
        return(edges)
    }) %>% do.call(rbind, .)
    
    return(edgelist)
}
                            
                            
parseargs = function(){
    
    option_list = list( 
        make_option("--splicing_file", type="character"),
        make_option("--genexpr_file", type="integer", default=10),
        make_option("--aracne_network_file", type="numeric", default=0.5),
        make_option("--output_as_edgelist", type="character", default="FALSE"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    splicing_file = args[["splicing_file"]]
    genexpr_file = args[["genexpr_file"]]
    aracne_network_file = args[["aracne_network_file"]]
    output_as_edgelist = args[["output_as_edgelist"]]
    output_file = args[["output_file"]]
    
    # load
    genexpr = read_tsv(genexpr_file)
    splicing = read_tsv(splicing_file)
    
    # prep
    genexpr = genexpr %>% column_to_rownames(colnames(genexpr)[1])
    splicing = splicing %>% column_to_rownames(colnames(splicing)[1])
    X = genexpr %>% bind_rows(splicing) %>% distinct() %>% as.matrix()
    
    # make regulons
    regulons = aracne2regulon(aracne_network_file, X, format="adj")
    
    # save
    if (output_as_edgelist=="TRUE"){
        
        regulons = regulons_to_edgelist(regulons)
        write_tsv(regulons, output_file)
        
    }else{
        saveRDS(regulons, output_file)    
    }
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}

