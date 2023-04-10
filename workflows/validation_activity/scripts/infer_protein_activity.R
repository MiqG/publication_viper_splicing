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
require(clusterProfiler)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","regulon_evaluation")

# signature_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"HepG2",'delta_psi_rel-EX-masked.tsv.gz')
# regulons_file = file.path(RESULTS_DIR,"files","associations","LIHC","genexpr_vs_psi_imputed","correlation_spearman.tsv.gz")
# assoc_method = "correlation_spearman"
# actinf_method = "viper"

##### FUNCTIONS #####
prep_regulons = function(regulons, assoc_method){
    
    if (assoc_method=="correlation_spearman"){
        
        regulons = regulons %>% drop_na(padj, association)
        regulons[["likelihood"]] = abs(regulons[["association"]])
        regulons[["tfmode"]] = sign(regulons[["association"]]) * (-log10(regulons[["padj"]]))
        #regulons[["tfmode"]] = sign(regulons[["association"]])
        regulons = regulons %>%
            distinct(regulator, target, likelihood, tfmode)
        
    } else if (assoc_method=="aracne_w_spearman"){
        
        regulons = regulons %>% drop_na(association_aracne, association_spearman)
        regulons[["likelihood"]] = regulons[["association_aracne"]]
        regulons[["tfmode"]] = regulons[["association_spearman"]]
        regulons = regulons %>%
            distinct(regulator, target, likelihood, tfmode)
        
    } else if (assoc_method=="linear_model"){
        
        regulons = regulons %>% drop_na(lr_padj, association)
        regulons[["likelihood"]] = abs(regulons[["association"]])
        regulons[["tfmode"]] = sign(regulons[["association"]]) * (-log10(regulons[["lr_padj"]]))
        #regulons[["tfmode"]] = sign(regulons[["association"]])
        regulons = regulons %>%
            distinct(regulator, target, likelihood, tfmode)
        
    }
    
    return(regulons)
}

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

run_viper = function(signature, regulons){
    
    regulon_network = as_regulon_network(regulons)
    
    # run viper
    result = viper(
        signature, 
        regulon_network,
        verbose=FALSE
    ) 
    
    return(result)
}


run_correlation = function(signature, regulons, method="spearman"){
    
    # correlate each regulon with each signatures
    regulators = regulons %>% pull(regulator) %>% unique()
    result = lapply(
        regulators, function(regulator_oi){
            
            # get associations for the regulator
            assoc_signature = regulons %>%
                filter(regulator == regulator_oi) %>% 
                mutate(association = likelihood * tfmode) %>%
                distinct(target, association) %>%
                deframe() %>%
                as.matrix()
            colnames(assoc_signature) = regulator_oi
            
            # compute correlation
            common_targets = intersect(rownames(assoc_signature), rownames(signature))
            corrs = cor(
                signature[common_targets,], 
                assoc_signature[common_targets,], 
                use="pairwise.complete.obs",
                method=method
            )
            colnames(corrs) = regulator_oi
            corrs = t(corrs)
            return(corrs)
    })
    result = do.call(rbind, result)
    
    return(result)
}


run_enrichment_test = function(signature, regulons){
    
    # are the changing features enriched in some regulon?
    # overlap enrichment analysis for each sample
    # save enrichment scores
    term2gene = regulons %>% distinct(regulator, target)
    universe = rownames(signature)
    samples = colnames(signature)
    result = lapply(
        samples, function(sample_oi){
            
            idx_changes = !is.na(signature[,sample_oi])
            changes_oi = rownames(signature[idx_changes,sample_oi,drop=FALSE])
                
            enrichment = enricher(
                changes_oi,
                universe=universe,
                pvalueCutoff=1,
                qvalueCutoff=1,
                minGSSize=0,
                maxGSSize=Inf,
                TERM2GENE=term2gene
            )
            
            enrich_ratios = data.frame(row.names = enrichment@result[["ID"]])
            enrich_ratios[[sample_oi]] = sapply(
                enrichment@result[["GeneRatio"]], function(x){eval(parse(text=x))}
            )
            return(enrich_ratios)  
        })
    result = do.call(cbind,result)
    
    return(result)
}


run_fisher_test = function(signature, regulons){
    
    # do the signs of changes match tf mode?
    regulators = regulons %>% pull(regulator) %>% unique()
    result = lapply(
        regulators, function(regulator_oi){
            
            # get associations for the regulator
            assoc_signature = regulons %>%
                filter(regulator == regulator_oi) %>% 
                mutate(association = sign(tfmode)) %>%
                distinct(target, association) %>%
                deframe() %>%
                as.matrix()
            colnames(assoc_signature) = regulator_oi
            
            # get fisher test p-values
            common_targets = intersect(rownames(assoc_signature), rownames(signature))
            assoc_signature = assoc_signature[common_targets,,drop=FALSE]
            
            fisher_odds = lapply(
                colnames(signature), function(sample_oi){
                    
                    column = signature[common_targets,sample_oi,drop=FALSE]
                    
                    X = cbind(sign(column), assoc_signature)
                    X = X[rowSums(abs(X))==2,]
                    X = table(X[,1], X[,2])
                    
                    test = fisher.test(X)
                    odds = data.frame(row.names = regulator_oi)
                    odds[[sample_oi]] = test[["estimate"]]
                    
                    return(odds)
                }
            )
            fisher_odds = do.call(cbind, fisher_odds)
            return(fisher_odds)
    })
    result = do.call(rbind, result)

    return(result)
}


infer_protein_activity = function(signature, regulons, actinf_method){
    
    if (actinf_method=="viper"){
        
        result = run_viper(signature, regulons)
    
    } else if (actinf_method=="correlation"){
        
        result = run_correlation(signature, regulons)
        
    } else if (actinf_method=="enrichment_test"){
        
        result = run_enrichment_test(signature, regulons)
        
    } else if (actinf_method=="fisher_test"){
        
        result = run_fisher_test(signature, regulons)
        
    }
    
    return(result)
}

parseargs = function(){
    
    option_list = list( 
        make_option("--signature_file", type="character"),
        make_option("--regulons_file", type="character"),
        make_option("--output_file", type="character"),
        make_option("--assoc_method", type="character"),
        make_option("--actinf_method", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    signature_file = args[["signature_file"]]
    regulons_file = args[["regulons_file"]]
    output_file = args[["output_file"]]
    assoc_method = args[["assoc_method"]]
    actinf_method = args[["actinf_method"]]
    
    # load
    signature = read_tsv(signature_file)
    regulons = read_tsv(regulons_file)
    
    # prep
    ## signature
    signature = signature %>% as.data.frame()
    rownames(signature) = signature[,1]
    signature = signature[,2:ncol(signature)]
    signature = signature %>% 
        dplyr::select(where(is.numeric))
    ## regulons
    regulons = prep_regulons(regulons, assoc_method)
    
    # run
    result = infer_protein_activity(signature, regulons, actinf_method)
    
    # save
    result = result %>% as.data.frame() %>% rownames_to_column('regulator')
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
