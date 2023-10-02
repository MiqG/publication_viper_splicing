# VIPER algorithm applied to splicing cancer-driver program

## Overview
1. [ ] Inference and validation of splicing factor activities
    - [X] download GO splicing factors (https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=GOBP_RNA_SPLICING)
    - [X] download GO RNA binding (https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOMF_RNA_BINDING.html)
    - [X] download CCLE splicing factor omics
    - [X] download SpliceAid RBP exon interactions (http://srv00.recas.ba.infn.it/SpliceAidF/)
    - [X] download eCLIP seq from ENCODE
    - [X] validate concordance CLIP seq vs RBP vs SpliceAid interactions
        - eCLIP vs RBP KD --> strong changes in splicing do not occur where RBP binds
    - [ ] network inference: save raw outputs and process them to make inference of regulation sign and magnitude
        - [ ] with correlations: mRNA levels SF vs exon inclusion
        - [ ] with WGCNA
        - [ ] with ARACNE
        - [ ] with robustica
        - [ ] with multiomic signatures
    - [ ] network evaluation:
        - we want to find those exons sensitive to changes in splicing factor activity
        - "correlation" between association scores and observed inclusion changes upon KD of that splicing factor
        - or precision and recall generating different sets of exons
            - from RBP KDs, consider different delta PSI thresholds (consider also possible % change thresholds) to define targets and non-targets
            - from eCLIP, we can only have true positives
            - from SpliceAid, we can only have true positives
        - evaluate considering or not the directionality of the SF-exon interactions
2. [ ] Chart and validate splicing factor activities in cancer
3. [ ] tissue specific vs pan-cancer regulatory networks

## Structure
1. `obtain_data`
2. `preprocess_data`
3. `regulon_inference`
    - experimentally derived regulons
    - aracne regulons
    - regulon evaluation
    - regulon eda
4. `sf_activity_validation`: validation of splicing factor activity estimation
    - combinatorial perturbations
    - SF3b complex
    - PTMs
5. `cancer_splicing_program`: cancer-driver splicing program
    - definition using TCGA
    - hallmarks
        - proliferation
        - (?) immune evasion
        - (?) metastasis
    - tumorigenesis
        - from fibroblasts to cancer
    - program regulators
        - perturbations that reverse oncogenic and tumor suppressor activities
   

## Requirements

## Usage

## Supporting Data

## Authors

## License

## References
- Giulietti M, Piva F, D'Antonio M, D'Onorio De Meo P, Paoletti D, Castrignanò T, D'Erchia AM, Picardi E, Zambelli F, Principato G, Pavesi G, Pesole G. SpliceAid-F: a database of human splicing factors and their RNA-binding sites. Nucleic Acids Res. 2013 Jan;41(Database issue):D125-31. doi: 10.1093/nar/gks997. Epub 2012 Oct 30. PMID: 23118479; PMCID: PMC3531144.