require(optparse)
require(tidyverse)
require(DESeq2)

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,"data","prep")
# data_file = file.path("~/databases/data/GDC/gene_counts","LUAD.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","LUAD.tsv.gz")
# sample_col = "sampleID"
# comparison_col = "sample_type_clean"
# condition_a = "PrimaryTumor"
# condition_b = "SolidTissueNormal"
# padj_method = "fdr"

##### FUNCTIONS #####
run_deseq2 = function(count_data, metadata, sample_col, comparison_col, condition_a, condition_b, padj_method, n_jobs) {
    # Use BiocParallel to specify the number of cores
    bpparam = BiocParallel::MulticoreParam(workers = n_jobs)
    
    # Modify column names to keep only the first 15 characters
    colnames(count_data) = substr(colnames(count_data), 1, 15)
    
    # Make sure rownames are the gene IDs and the column names match sample names
    count_data = count_data %>% as.data.frame()
    rownames(count_data) = count_data[,1]
    count_data = count_data[,2:ncol(count_data)]
    count_data = count_data %>% 
        dplyr::select(where(is.numeric)) %>% 
        mutate_if(is.numeric, replace_na, replace = 0) %>%
        mutate_if(is.numeric, as.integer)
    
    # Match metadata and count data by sample column
    common_samples = intersect(metadata[["sampleID"]],colnames(count_data))
    col_data = metadata %>%
        filter(!!sym(comparison_col) %in% c(condition_a, condition_b)) %>%  # filter relevant conditions
        filter(!!sym(sample_col) %in% common_samples) %>%  # filter relevant samples
        column_to_rownames(sample_col)  # Make sample column rownames
    
    count_data = count_data[,rownames(col_data)]  # Match sample names in count_data to metadata

    # Create DESeq2 dataset
    dds = DESeqDataSetFromMatrix(
        countData = count_data, 
        colData = col_data, 
        design = as.formula(paste("~", comparison_col))
    )
    
    # Filter low counts
    dds = dds[rowSums(counts(dds)) > 1, ]
    
    # Run DESeq2 analysis
    dds = DESeq(dds, parallel=TRUE, BPPARAM=bpparam)
    
    # Extract results for condition_a vs condition_b
    result = results(dds, contrast = c(comparison_col, condition_a, condition_b), pAdjustMethod = padj_method)
    
    # Turn result into a tibble
    result_df = as_tibble(result, rownames = "gene")
    
    # add info
    ## normalized counts
    normalized_counts = counts(dds, normalized=TRUE) %>% as.data.frame()
    mean_counts_a = rowMeans(normalized_counts[, col_data[[comparison_col]] == condition_a], na.rm = TRUE)
    mean_counts_b = rowMeans(normalized_counts[, col_data[[comparison_col]] == condition_b], na.rm = TRUE)
    median_counts_a = matrixStats::rowMedians(as.matrix(normalized_counts[, col_data[[comparison_col]] == condition_a]), na.rm = TRUE)
    median_counts_b = matrixStats::rowMedians(as.matrix(normalized_counts[, col_data[[comparison_col]] == condition_b]), na.rm = TRUE)
    
    ## add
    result_df = result_df %>%
        mutate(
            comparison_col = comparison_col,
            condition_a = condition_a,
            condition_b = condition_b,
            `condition_a-mean` = mean_counts_a[gene],
            `condition_b-mean` = mean_counts_b[gene],
            `condition_a-median` = median_counts_a[gene],
            `condition_b-median` = median_counts_b[gene]
        )
    
    
    return(result_df)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--data_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--sample_col", type="character"),
        make_option("--comparison_col", type="character"),
        make_option("--condition_a", type="character"),
        make_option("--condition_b", type="character"),
        make_option("--padj_method", type="character"),
        make_option("--n_jobs", type="integer", default=1),
        make_option("--output_file", type="character"),
        make_option("--random_seed", type="integer", default=1234)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    data_file = args[["data_file"]]
    metadata_file = args[["metadata_file"]]
    sample_col = args[["sample_col"]]
    comparison_col = args[["comparison_col"]]
    condition_a = args[["condition_a"]]
    condition_b = args[["condition_b"]]
    padj_method = args[["padj_method"]]
    n_jobs = args[["n_jobs"]]
    output_file = args[["output_file"]]
    
    set.seed(args[["random_seed"]])
    
    # load
    data = read_tsv(data_file)
    metadata = read_tsv(metadata_file)
    
    # differential gene expression analysis
    result = run_deseq2(data, metadata, sample_col, comparison_col, condition_a, condition_b, padj_method, n_jobs)
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
