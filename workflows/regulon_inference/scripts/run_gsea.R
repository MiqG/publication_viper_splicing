#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA of gene dependencies regressed on event PSI and gene TPMs.
# 
# Outline
# -------

require(optparse)
require(tidyverse)
require(clusterProfiler)
require(org.Hs.eg.db)

# variables
THRESH_FDR = 0.05
ORGDB = org.Hs.eg.db

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','regulon_inference')
# regulons_path = file.path(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX")
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# msigdb_dir = file.path(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')

##### FUNCTIONS #####
load_regulons = function(regulons_path, patt=NULL){
    if (file.exists(regulons_path) && !dir.exists(regulons_path)){
        # regulons_path is a file, we load only that regulon (we'll tun regular VIPER)
        regulon_files = list(regulons_path)
    }else if (dir.exists(regulons_path)){
        # regulons_path is a directory, we load all regulons contained (we'll run metaVIPER)
        regulon_files = list.files(regulons_path, pattern=patt, full.names=TRUE)
    }else {
        stop("Invalid regulons_path.")
    }
    
    regulons = sapply(regulon_files, function(regulon_file){
        regulon = read_tsv(regulon_file)
        return(regulon)
    }, simplify=FALSE)
    
    return(regulons)
}


get_genes_lists = function(regulons){
    df = regulons
    
    genes_lists = regulons %>%
        drop_na() %>%
        with(., split(GENE, regulator))
        
    return(genes_lists)
}


get_events_lists = function(regulons){
    df = regulons 
    
    events_lists = regulons %>%
        drop_na() %>%
        with(., split(target, regulator))
        
    return(events_lists)
}


get_universe = function(regulons){
    df = regulons
    universe = df[,c('target','GENE')] %>% apply(., 2, unique)
    names(universe) = c('events','genes')
    return(universe)
}


run_enrichment = function(genes, events, universe, ontologies){
    enrichments = list()
    if(length(genes)>0){
        enrichments[['hallmarks']] = enricher(genes, TERM2GENE=ontologies[['hallmarks']], universe=universe[['genes']])
        enrichments[['oncogenic_signatures']] = enricher(genes, TERM2GENE=ontologies[['oncogenic_signatures']], universe=universe[['genes']])
        enrichments[['GO_BP']] = enrichGO(genes, OrgDb=ORGDB, keyType='SYMBOL', universe=universe[['genes']])
    }
    
    if(length(events)>0){
        enrichments[['protein_impact']] = enricher(events, TERM2GENE=ontologies[['protein_impact']], universe=universe[['events']])
    }
    return(enrichments)
}


run_enrichments = function(genes_lists, events_lists, universe, ontologies){
    enrichments = sapply(names(genes_lists), function(cond){
            genes = genes_lists[[cond]]
            events = events_lists[[cond]]
            run_enrichment(genes, events, universe, ontologies)
        }, simplify=FALSE)
    return(enrichments)
}


get_enrichment_result = function(enrich_list, thresh){
    ## regulators are extracted from names
    if (length(enrich_list)>0){
        ontos = names(enrich_list[[1]])
        regulators = names(enrich_list)
        results = sapply(
            ontos, function(onto){
            result = lapply(regulators, function(regulator){
                result = enrich_list[[regulator]][[onto]]
                if(!is.null(result)){
                    result = result@result
                    result[['regulator']] = regulator
                    result[['ontology']] = onto
                }
                return(result)
            })
            result[sapply(result, is.null)] = NULL
            result = do.call(rbind,result)
            ## filter by p.adjusted
            result = result %>% filter(p.adjust<thresh)
        }, simplify=FALSE)
        results = do.call(rbind,results)
    }else{
        results = data.frame(
            ID=NA, Description=NA, GeneRatio=NA, BgRatio=NA, pvalue=NA,
            p.adjust=NA, qvalue=NA, geneID=NA, Count=NA, regulator=NA, ontology=NA
        ) %>% drop_na()
    }
    
    return(results)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--regulons_path", type="character"),
        make_option("--annotation_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--protein_impact_file", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    regulons_path = args[["regulons_path"]]
    annotation_file = args[["annotation_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    protein_impact_file = args[["protein_impact_file"]]
    output_file = args[["output_file"]]
    
    # load
    regulons = load_regulons(regulons_path)
    annot = read_tsv(annotation_file)
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt")),
        "protein_impact" = read_tsv(protein_impact_file) %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT)
    )
    
    # prep
    regulons = regulons %>% 
        bind_rows() %>%
        distinct(regulator, target) %>%
        left_join(annot %>% distinct(EVENT, GENE), by=c("target"="EVENT"))
    
    # run enrichments
    genes_lists = get_genes_lists(regulons)
    events_lists = get_events_lists(regulons)
    universe = get_universe(regulons)
    enrichments = run_enrichments(genes_lists, events_lists, universe, ontologies)
    results_enrich = get_enrichment_result(enrichments, THRESH_FDR)
    
    # save
    write_tsv(results_enrich, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}