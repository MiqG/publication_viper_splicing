require(optparse)
require(tidyverse)
require(clusterProfiler)

# Development
# -----------
# ROOT = here::here()
# RESULTS_DIR = file.path(ROOT,"results","cancer_splicing_program")
# omic_table_file = file.path(RESULTS_DIR,"files","protein_activity","CCLE-EX.tsv.gz")
# gene_sets_file = file.path(ROOT,"workflows","cancer_splicing_program","scripts","driver_types-gene_sets.tsv")

##### FUNCTIONS #####
add_diff_probs = function(df){
    
    # real gene set
    ecdf_in_gene_set_real = ecdf(df %>% filter(in_gene_set) %>% pull(ranking))
    ecdf_notin_gene_set = ecdf(df %>% filter(!in_gene_set) %>% pull(ranking))
    df = df %>% mutate(
        diff_probs = ecdf_in_gene_set_real(ranking) - ecdf_notin_gene_set(ranking)
    )

    return(df)
}


compute_enrichment_score_single = function(x, gene_sets){
    
    result = x %>% 
        # prepare data
        enframe(name="gene") %>%
        drop_na() %>%
        left_join(gene_sets, by="gene") %>%
        mutate(in_gene_set = ifelse(is.na(term), -1, 1)) %>%
        pivot_wider(id_cols=c("gene","value"), names_from = "term", values_from="in_gene_set", values_fill=-1) %>%
        pivot_longer(-c(gene,value), names_to="term", values_to="score") %>%
        filter(term!="NA") %>%
        # compute running scores
        group_by(term) %>%
        arrange(desc(value)) %>%
        mutate(
            in_gene_set = score==1,
            ranking = row_number(),
        ) %>%
        ungroup(term) %>%
        group_by(in_gene_set) %>%
        mutate(
            score = 1,
            running_score = cumsum(score) / n()
        ) %>%
        ungroup(in_gene_set)
    
    # fit ecdf on the rankings of each term ad
    result = lapply(result %>% pull(term) %>% unique(), function(term_oi){ 
        df = result %>% filter(term==term_oi)
        df = add_diff_probs(df)
        return(df)
    }) %>% bind_rows()
    
    # compute enrichment score
    enrichment_score = result %>%
        group_by(term) %>%
        slice_max(abs(diff_probs), n=1, with_ties=FALSE) %>%
        ungroup() %>%
        column_to_rownames("term") %>%
        distinct(diff_probs)
    
    return(enrichment_score)
}


# compute_enrichment_score_single = function(x, gene_sets){
    
#     enrichment_score = x %>%
#         sort(decreasing=TRUE) %>%
#         GSEA(TERM2GENE=gene_sets, pvalueCutoff=1, verbose=FALSE) %>% 
#         as.data.frame() %>% 
#         mutate(score = sign(enrichmentScore)*(1-pvalue)) %>%
#         distinct(score)
    
#     return(enrichment_score)
# }


compute_enrichment_score = function(omic_table, gene_sets){
    
    result = apply(omic_table, 2, function(x){
        enrichment_score = compute_enrichment_score_single(x, gene_sets)
        return(enrichment_score)
    }, simplify=FALSE) 
    sample_names = names(result)
    result = result %>% bind_cols()
    colnames(result) = sample_names
    
    return(result)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--omic_table_file", type="character"),
        make_option("--gene_sets_file", type="character"),
        make_option("--output_file", type="character"),
        make_option("--random_seed", type="integer", default=1234)
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    omic_table_file = args[["omic_table_file"]]
    gene_sets_file = args[["gene_sets_file"]]
    output_file = args[["output_file"]]
    
    set.seed(args[["random_seed"]])
    
    # load
    omic_table = read_tsv(omic_table_file)
    gene_sets = read_tsv(gene_sets_file)
    
    # prep
    omic_table = omic_table %>% as.data.frame()
    rownames(omic_table) = omic_table[,1]
    omic_table = omic_table[,2:ncol(omic_table)]
    omic_table = omic_table %>% dplyr::select(where(is.numeric))
    
    colnames(gene_sets) = c("term","gene")
    
    # compute enrichment_score for each sample/column
    result = compute_enrichment_score(omic_table, gene_sets)
    
    # prepare outputs
    result = result %>%
        rownames_to_column("term")
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
