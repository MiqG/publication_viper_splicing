# Regulon Inference - Splicing factor network inference workflows

This set of Snakemake workflows make up a pipeline for constructing and evaluating splicing factor networks generated through different approaches. These workflows use raw and preprocessed files to populate the `results/regulon_inference` directory.

## Outline (run in the following order)
1. Regulon inference
    - `aracne_regulons.smk`
        - infers computational splicing factor networks using the ARACNe algorithm
        - executed for three observational datasets: CardosoMoreira2020, TCGA-STN, and TCGA-PT
        - outputs in `results/regulon_inference/files/arance_regulons`
        
    - `mlr_regulons.smk`
        - infers computational splicing factor networks using multivariate linear regression (MLR) models
        - executed for three observational datasets: CardosoMoreira2020, TCGA-STN, and TCGA-PT
        - outputs in `results/regulon_inference/files/mlr_regulons`
    
    - `experimentally_derived_regulons.smk`
        - infers empirical splicing factor networks from RNA-seq experiments perturbing individual splicing factors
        - outputs in `results/regulon_inference/files/experimentally_derived_regulons_pruned-EX`
    
    - `postar3_clip_regulons.smk`
        - infers splicing factor networks from CLIP experiments for RBP splicing factors
        - outputs in `results/regulon_inference/files/postar3_clip_regulons-EX`
    
    - `splicinglore_regulons.smk`
        - maps SplicingLore database to VastDB exon identifiers
        - infers splicing factor networks from SplicingLore data
        - outputs in `results/regulon_inference/files/splicinglore_regulons-EX`
        
    - `inferred_and_experimental_regulons.smk`
        - infers splicing factor networks combining features from other splicing factor networks:
            - ARACNe MI + experimental MoR
            - MLR beta + experimental MoR
            - ARACNe MI + ARACNe MoR (negative) + MLR MoR (positive)
        - outputs in 
            - `results/regulon_inference/files/aracne_and_experimental_regulons-EX`        
            - `results/regulon_inference/files/mlr_and_experimental_regulons-EX`        
            - `results/regulon_inference/files/aracne_and_mlr_regulons-EX`        
            
    - `postar3_clip_regulons.smk`
        - infers splicing factor networks from CLIP experiments for RBP splicing factors
        - outputs in 
            - `results/regulon_inference/files/postar3_and_experimental_regulons-EX`        
            - `results/regulon_inference/files/experimental_without_postar3_regulons-EX`        
            
2. Benchmarks of empirical splicing factor networks
    - `regulon_robustness.smk`: 
        - generates empirical networks considering top N edges
        - outputs in 
            - `results/regulon_inference/files/top{N}_experimentally_derived_regulons_pruned-EX`
            - `results/regulon_inference/files/regulon_properties`
        
    - `regulon_thresholds.smk`: 
        - generates empirical networks considering PSI thresholds
        - outputs in 
            - `results/regulon_inference/files/dPSIthresh{N}_experimentally_derived_regulons_pruned-EX`
            - `results/regulon_inference/files/regulon_properties`
            
    - `aracne_regulons_benchmark_boostrap_pruning.smk`
        - generates computational networks using different "max_targets" parameter in ARACNe
        - outputs in `results/regulon_inference/files/aracne_regulons_CardosoMoreira2020_{max_targets}-EX`

        
3. Regulon Evaluation workflow (`regulon_evaluation.smk`): computes performance evaluation metrics for each type of network
    
4. Regulon EDA (`regulon_eda.smk`): 
    - EDA of inferred empirical splicing factor networks
    - outputs in `results/regulon_inference/files`
    
## Notes
To run `aracne_regulons.smk`, adapted ARACNe-AP was obtained from https://github.com/chaolinzhanglab/ARACNe-AP and `scripts/aracne/filter_arachne_bootstraps.pl` and `scripts/aracne/estimate_mor.R` were obtained from https://github.com/chaolinzhanglab/mras.

## Expected outputs
```{shell}
results/regulon_inference/
├── figures
│   ├── eda_regulons-EX
│   │   ├── clip_eda-distance_distr-hist.pdf
│   │   ├── emp_vs_clip-interactions-venn.pdf
│   │   ├── emp_vs_clip-regulators-venn.pdf
│   │   ├── figdata
│   │   │   └── eda
│   │   │       ├── protein_impact_frequencies.tsv.gz
│   │   │       ├── regulons.tsv.gz
│   │   │       └── regulons_umap.tsv.gz
│   │   ├── protein_impact-freqs-violin.pdf
│   │   ├── regulons-n_regulators_per_target-box.pdf
│   │   ├── regulons-n_targets_per_regulator-box.pdf
│   │   ├── regulons-n_targets_per_regulator_vs_sf_class-box.pdf
│   │   ├── similarities-jaccard-box.pdf
│   │   ├── similarities-umap-scatter.pdf
│   │   └── target_lengths-regulators-scatter.pdf
│   ├── inference_troubleshooting
│   │   ├── inference_eval-likelihood-box.pdf
│   │   └── inference_eval-tfmode-bar.pdf
│   ├── regulon_evaluation
│   │   ├── eval_splicinglore-held_out_ds-raw_auc_roc-box.pdf
│   │   ├── evaluation-clip-median_auc_roc-box.pdf
│   │   ├── evaluation-comp_best_aracne_mlr-median_auc_roc-box.pdf
│   │   ├── evaluation-comp_best_vs_in_empirical-median_auc_roc-box.pdf
│   │   ├── evaluation-comp_best_w_empirical_likelihood-median_auc_roc-box.pdf
│   │   ├── evaluation-comp_split-median_auc_roc-box.pdf
│   │   ├── evaluation-general-median_auc_roc-box.pdf
│   │   ├── evaluation-held_out_ds-raw_auc_roc-box.pdf
│   │   ├── evaluation-max_targets_aracne-median_auc_roc-box.pdf
│   │   ├── evaluation-robustness-median_auc_roc-box.pdf
│   │   ├── evaluation-sf_class-raw_auc_roc-box.pdf
│   │   ├── evaluation-sf_class_vs_held_out_ds-raw_auc_roc-box.pdf
│   │   ├── evaluation-threshold-median_auc_roc-box.pdf
│   │   └── figdata
│   │       └── regulon_evaluation
│   │           └── evaluation.tsv.gz
│   └── regulon_inference
│       ├── figdata
│       │   └── inference_evaluation
│       │       ├── eval_likelihood.tsv.gz
│       │       └── eval_tfmode.tsv.gz
│       ├── inference_eval-likelihood-box.pdf
│       └── inference_eval-tfmode-bar.pdf
└── files
    ├── aracne_and_experimental_regulons-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── aracne_and_mlr_regulons-EX
    │   └── PANCAN_PT.tsv.gz
    ├── aracne_regulons
    │   ├── CardosoMoreira2020-EX
    │   │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │   ├── genexpr.tsv
    │   │   ├── miThreshold_p1E-8_samples313.txt
    │   │   ├── pruned
    │   │   │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │   │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │   │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │   │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │   │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │   │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │   │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │   │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │   │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │   │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │   │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │   │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │   │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │   │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │   │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │   │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │   │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │   │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │   │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │   │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │   │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │   │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │   │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │   │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │   │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │   │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │   │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │   │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │   │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │   │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │   │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │   │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │   │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │   │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │   │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │   │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │   │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │   │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │   │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │   │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │   │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │   │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │   │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │   │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │   │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │   │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │   │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │   │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │   │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │   │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │   │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │   │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │   │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │   │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │   │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │   │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │   │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │   │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │   │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │   │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │   │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │   │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │   │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │   │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │   │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │   │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │   │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │   │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │   │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │   │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │   │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │   │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │   │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │   │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │   │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │   │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │   │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │   │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │   │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │   │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │   │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │   │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │   │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │   │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │   │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │   │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │   │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │   │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │   │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │   │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │   │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │   │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │   │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │   │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │   │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │   │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │   │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │   │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │   │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │   │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │   │   ├── network.txt
    │   │   │   └── regulons.tsv.gz
    │   │   └── splicing.tsv
    │   ├── CardosoMoreira2020-EX-100
    │   │   └── pruned
    │   │       ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │       ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │       ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │       ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │       ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │       ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │       ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │       ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │       ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │       ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │       ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │       ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │       ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │       ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │       ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │       ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │       ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │       ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │       ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │       ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │       ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │       ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │       ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │       ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │       ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │       ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │       ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │       ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │       ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │       ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │       ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │       ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │       ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │       ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │       ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │       ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │       ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │       ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │       ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │       ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │       ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │       ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │       ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │       ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │       ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │       ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │       ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │       ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │       ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │       ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │       ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │       ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │       ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │       ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │       ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │       ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │       ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │       ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │       ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │       ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │       ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │       ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │       ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │       ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │       ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │       ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │       ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │       ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │       ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │       ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │       ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │       ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │       ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │       ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │       ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │       ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │       ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │       ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │       ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │       ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │       ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │       ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │       ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │       ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │       ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │       ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │       ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │       ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │       ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │       ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │       ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │       ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │       ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │       ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │       ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │       ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │       ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │       ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │       ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │       ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │       ├── network.txt
    │   │       └── regulons.tsv.gz
    │   ├── CardosoMoreira2020-EX-1000
    │   │   └── pruned
    │   │       ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │       ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │       ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │       ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │       ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │       ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │       ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │       ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │       ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │       ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │       ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │       ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │       ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │       ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │       ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │       ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │       ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │       ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │       ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │       ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │       ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │       ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │       ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │       ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │       ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │       ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │       ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │       ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │       ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │       ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │       ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │       ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │       ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │       ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │       ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │       ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │       ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │       ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │       ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │       ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │       ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │       ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │       ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │       ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │       ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │       ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │       ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │       ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │       ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │       ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │       ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │       ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │       ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │       ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │       ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │       ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │       ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │       ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │       ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │       ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │       ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │       ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │       ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │       ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │       ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │       ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │       ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │       ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │       ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │       ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │       ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │       ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │       ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │       ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │       ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │       ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │       ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │       ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │       ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │       ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │       ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │       ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │       ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │       ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │       ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │       ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │       ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │       ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │       ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │       ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │       ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │       ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │       ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │       ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │       ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │       ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │       ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │       ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │       ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │       ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │       ├── network.txt
    │   │       └── regulons.tsv.gz
    │   ├── CardosoMoreira2020-EX-2000
    │   │   └── pruned
    │   │       ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │       ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │       ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │       ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │       ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │       ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │       ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │       ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │       ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │       ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │       ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │       ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │       ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │       ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │       ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │       ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │       ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │       ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │       ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │       ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │       ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │       ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │       ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │       ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │       ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │       ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │       ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │       ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │       ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │       ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │       ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │       ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │       ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │       ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │       ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │       ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │       ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │       ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │       ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │       ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │       ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │       ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │       ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │       ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │       ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │       ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │       ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │       ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │       ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │       ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │       ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │       ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │       ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │       ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │       ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │       ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │       ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │       ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │       ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │       ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │       ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │       ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │       ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │       ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │       ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │       ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │       ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │       ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │       ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │       ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │       ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │       ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │       ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │       ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │       ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │       ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │       ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │       ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │       ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │       ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │       ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │       ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │       ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │       ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │       ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │       ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │       ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │       ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │       ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │       ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │       ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │       ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │       ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │       ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │       ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │       ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │       ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │       ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │       ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │       ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │       ├── network.txt
    │   │       └── regulons.tsv.gz
    │   ├── CardosoMoreira2020-EX-500
    │   │   └── pruned
    │   │       ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │       ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │       ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │       ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │       ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │       ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │       ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │       ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │       ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │       ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │       ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │       ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │       ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │       ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │       ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │       ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │       ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │       ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │       ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │       ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │       ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │       ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │       ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │       ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │       ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │       ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │       ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │       ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │       ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │       ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │       ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │       ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │       ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │       ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │       ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │       ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │       ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │       ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │       ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │       ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │       ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │       ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │       ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │       ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │       ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │       ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │       ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │       ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │       ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │       ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │       ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │       ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │       ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │       ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │       ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │       ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │       ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │       ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │       ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │       ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │       ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │       ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │       ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │       ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │       ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │       ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │       ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │       ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │       ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │       ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │       ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │       ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │       ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │       ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │       ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │       ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │       ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │       ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │       ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │       ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │       ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │       ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │       ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │       ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │       ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │       ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │       ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │       ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │       ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │       ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │       ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │       ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │       ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │       ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │       ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │       ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │       ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │       ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │       ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │       ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │       ├── network.txt
    │   │       └── regulons.tsv.gz
    │   ├── CardosoMoreira2020-EX-5000
    │   │   └── pruned
    │   │       ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │       ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │       ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │       ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │       ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │       ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │       ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │       ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │       ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │       ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │       ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │       ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │       ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │       ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │       ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │       ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │       ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │       ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │       ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │       ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │       ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │       ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │       ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │       ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │       ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │       ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │       ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │       ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │       ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │       ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │       ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │       ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │       ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │       ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │       ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │       ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │       ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │       ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │       ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │       ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │       ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │       ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │       ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │       ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │       ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │       ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │       ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │       ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │       ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │       ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │       ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │       ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │       ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │       ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │       ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │       ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │       ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │       ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │       ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │       ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │       ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │       ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │       ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │       ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │       ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │       ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │       ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │       ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │       ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │       ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │       ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │       ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │       ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │       ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │       ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │       ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │       ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │       ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │       ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │       ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │       ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │       ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │       ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │       ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │       ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │       ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │       ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │       ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │       ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │       ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │       ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │       ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │       ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │       ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │       ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │       ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │       ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │       ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │       ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │       ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │       ├── network.txt
    │   │       └── regulons.tsv.gz
    │   ├── CardosoMoreira2020-genexpr
    │   │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │   ├── miThreshold_p1E-8_samples313.txt
    │   │   ├── pruned
    │   │   │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │   │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │   │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │   │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │   │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │   │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │   │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │   │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │   │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │   │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │   │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │   │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │   │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │   │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │   │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │   │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │   │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │   │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │   │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │   │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │   │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │   │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │   │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │   │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │   │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │   │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │   │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │   │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │   │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │   │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │   │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │   │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │   │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │   │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │   │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │   │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │   │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │   │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │   │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │   │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │   │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │   │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │   │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │   │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │   │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │   │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │   │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │   │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │   │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │   │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │   │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │   │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │   │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │   │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │   │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │   │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │   │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │   │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │   │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │   │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │   │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │   │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │   │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │   │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │   │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │   │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │   │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │   │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │   │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │   │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │   │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │   │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │   │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │   │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │   │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │   │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │   │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │   │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │   │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │   │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │   │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │   │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │   │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │   │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │   │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │   │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │   │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │   │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │   │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │   │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │   │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │   │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │   │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │   │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │   │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │   │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │   │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │   │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │   │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │   │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │   │   ├── network.txt
    │   │   │   └── regulons.tsv.gz
    │   │   ├── regulators.tsv
    │   │   └── targets.tsv
    │   ├── PANCAN_PT-EX
    │   │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │   ├── miThreshold_p1E-8_samples9854.txt
    │   │   ├── pruned
    │   │   │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │   │   │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │   │   │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │   │   │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │   │   │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │   │   │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │   │   │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │   │   │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │   │   │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │   │   │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │   │   │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │   │   │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │   │   │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │   │   │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │   │   │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │   │   │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │   │   │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │   │   │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │   │   │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │   │   │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │   │   │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │   │   │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │   │   │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │   │   │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │   │   │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │   │   │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │   │   │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │   │   │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │   │   │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │   │   │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │   │   │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │   │   │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │   │   │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │   │   │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │   │   │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │   │   │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │   │   │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │   │   │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │   │   │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │   │   │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │   │   │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │   │   │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │   │   │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │   │   │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │   │   │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │   │   │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │   │   │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │   │   │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │   │   │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │   │   │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │   │   │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │   │   │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │   │   │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │   │   │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │   │   │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │   │   │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │   │   │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │   │   │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │   │   │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │   │   │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │   │   │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │   │   │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │   │   │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │   │   │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │   │   │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │   │   │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │   │   │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │   │   │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │   │   │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │   │   │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │   │   │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │   │   │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │   │   │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │   │   │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │   │   │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │   │   │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │   │   │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │   │   │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │   │   │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │   │   │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │   │   │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │   │   │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │   │   │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │   │   │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │   │   │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │   │   │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │   │   │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │   │   │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │   │   │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │   │   │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │   │   │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │   │   │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │   │   │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │   │   │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │   │   │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │   │   │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │   │   │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │   │   │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │   │   │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │   │   │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │   │   │   ├── network.txt
    │   │   │   └── regulons.tsv.gz
    │   │   ├── regulators.tsv
    │   │   └── targets.tsv
    │   └── PANCAN_STN-EX
    │       ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │       ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │       ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │       ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │       ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │       ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │       ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │       ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │       ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │       ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │       ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │       ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │       ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │       ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │       ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │       ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │       ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │       ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │       ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │       ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │       ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │       ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │       ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │       ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │       ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │       ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │       ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │       ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │       ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │       ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │       ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │       ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │       ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │       ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │       ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │       ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │       ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │       ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │       ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │       ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │       ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │       ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │       ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │       ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │       ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │       ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │       ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │       ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │       ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │       ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │       ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │       ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │       ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │       ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │       ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │       ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │       ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │       ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │       ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │       ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │       ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │       ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │       ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │       ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │       ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │       ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │       ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │       ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │       ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │       ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │       ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │       ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │       ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │       ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │       ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │       ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │       ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │       ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │       ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │       ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │       ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │       ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │       ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │       ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │       ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │       ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │       ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │       ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │       ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │       ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │       ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │       ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │       ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │       ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │       ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │       ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │       ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │       ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │       ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │       ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │       ├── miThreshold_p1E-8_samples736.txt
    │       ├── pruned
    │       │   ├── bootstrapNetwork_1o16trgvqmt17tec7b1ma8njof.txt
    │       │   ├── bootstrapNetwork_1paet1fa307s0v7c9clkigh9eh.txt
    │       │   ├── bootstrapNetwork_1rjmsnbk3aiqq0scre9qro8uki.txt
    │       │   ├── bootstrapNetwork_1ssusd9u3ktlj2ld5ftp402k2k.txt
    │       │   ├── bootstrapNetwork_31hatvjvp5vku1ck9grn8sebs6.txt
    │       │   ├── bootstrapNetwork_32qit5ga1gajn33kjhflhk81a7.txt
    │       │   ├── bootstrapNetwork_333qtbck1qleg4qktj3rqs1mg9.txt
    │       │   ├── bootstrapNetwork_34d2shaua4099mjl7knq33rc6b.txt
    │       │   ├── bootstrapNetwork_4aamtpja7vd3dn1slm9mgo0pdu.txt
    │       │   ├── bootstrapNetwork_4cjutffk09o26oosvotopfqek0.txt
    │       │   ├── bootstrapNetwork_4dt6sldu8j2svqht1phr27i421.txt
    │       │   ├── bootstrapNetwork_4j4qtpg1fqhu2ogb0itbji9jng.txt
    │       │   ├── bootstrapNetwork_5k42tjgkepqhtcn4ptnpo3j6vm.txt
    │       │   ├── bootstrapNetwork_5ldat9euf35cmee53vbs1bcs5o.txt
    │       │   ├── bootstrapNetwork_5nmisvd8fdgbfg55m0vu9j4hjp.txt
    │       │   ├── bootstrapNetwork_5su6t3hbmjvcie7jcobertq0p8.txt
    │       │   ├── bootstrapNetwork_5u7et9flmtabbfsjmqvh45lmfa.txt
    │       │   ├── bootstrapNetwork_6v6mt3e8ltiv5k3dg5pv8mt9ng.txt
    │       │   ├── bootstrapNetwork_70fuspaim7tpv5qdq7e1humuti.txt
    │       │   ├── bootstrapNetwork_75eatnibl3208i5rmt5fqhip4v.txt
    │       │   ├── bootstrapNetwork_76nitdelldcv1jqrovpi39ceb0.txt
    │       │   ├── bootstrapNetwork_770qt3cvtnnpqljs30dgc163p2.txt
    │       │   ├── bootstrapNetwork_8n1uthgd3nkdd99idv745ne0g9.txt
    │       │   ├── bootstrapNetwork_8ob6t7en41v86b0ig1r6ef7lub.txt
    │       │   ├── bootstrapNetwork_8qkestd14ba6vcniq2f8n71bcc.txt
    │       │   ├── bootstrapNetwork_8rtms39bcll1oeej443b0ep0ie.txt
    │       │   ├── bootstrapNetwork_a1ratbhn2g1vsusqi6lbd2ue21.txt
    │       │   ├── bootstrapNetwork_a24it1c1arcqlglr4799mao383.txt
    │       │   ├── bootstrapNetwork_a3dqsnabb5nleicre9tfuihou4.txt
    │       │   ├── bootstrapNetwork_a8c6tli4a1rvouk8r0kq7lfith.txt
    │       │   ├── bootstrapNetwork_b9betvin104jj2r2sbf8c6p65o.txt
    │       │   ├── bootstrapNetwork_bakmt5f19afec4i36d3alegrbp.txt
    │       │   ├── bootstrapNetwork_bctusrdb9kq956b38enctmagpr.txt
    │       │   ├── bootstrapNetwork_bi5itvhe8q9e84bhf62tfh00fa.txt
    │       │   ├── bootstrapNetwork_bjeqt5foh4k9162hh8mvo8pllb.txt
    │       │   ├── bootstrapNetwork_cke2tfeb84ssrq9bajhdsq38ti.txt
    │       │   ├── bootstrapNetwork_clnaslclge7rkrubkl5c5hsubj.txt
    │       │   ├── bootstrapNetwork_cruut9gonkmonq0pjdgsncidh2.txt
    │       │   ├── bootstrapNetwork_ct86tff2nu1ngrnpte5304a374.txt
    │       │   ├── bootstrapNetwork_cuheslbco8cm9teq7gp18s5od5.txt
    │       │   ├── bootstrapNetwork_dvgmsvbvn7l641lk0rjjdddblb.txt
    │       │   ├── bootstrapNetwork_e4f2ttjom3pcedv1lib1m0b5so.txt
    │       │   ├── bootstrapNetwork_e5oat3g2ue4b7fk1vkv3uo2r2q.txt
    │       │   ├── bootstrapNetwork_e61it9ccuof601d29lj67vsg8s.txt
    │       │   ├── bootstrapNetwork_f7qthhlrsfio7p3r4jj2d56en.txt
    │       │   ├── bootstrapNetwork_fd8etnh2st6utjia1pp0trtj6h.txt
    │       │   ├── bootstrapNetwork_ffhmtdfct7hpmlbabqd76jl8ki.txt
    │       │   ├── bootstrapNetwork_fgqusjdn5hskfn2als15fretqk.txt
    │       │   ├── bootstrapNetwork_gh2t7fvs6qdh9g4567hb4ursp.txt
    │       │   ├── bootstrapNetwork_h0s2t1h4bhp82ao0orqp8hmqpr.txt
    │       │   ├── bootstrapNetwork_h15at7debr46rcf12serhpeg7t.txt
    │       │   ├── bootstrapNetwork_h2eistboc5f1ke61cu2tqh85lu.txt
    │       │   ├── bootstrapNetwork_h7cutrjhb1jbtqff1kqc347vlb.txt
    │       │   ├── bootstrapNetwork_h8m6t1frjbu6nc6fbmeecbvkrd.txt
    │       │   ├── bootstrapNetwork_hqastc9sh58ab94f7rnjsmh2q.txt
    │       │   ├── bootstrapNetwork_i9letbeeab6qhgd9d18sgt983j.txt
    │       │   ├── bootstrapNetwork_ibumshcoilhlai49n3sup50thl.txt
    │       │   ├── bootstrapNetwork_ih6atlgrhr0qdg4nlr8fbfmcv3.txt
    │       │   ├── bootstrapNetwork_iifitbf5i5bl6htnvsshjng2d5.txt
    │       │   ├── bootstrapNetwork_ijoqshdfqfmjvjio1ugjsv9nj7.txt
    │       │   ├── bootstrapNetwork_jko2src2hfv7q7phr9b210jb3d.txt
    │       │   ├── bootstrapNetwork_jqvmtfg5gle8t5rvq1miir8q8s.txt
    │       │   ├── bootstrapNetwork_js8ut5efovp7m7h043agrj2fut.txt
    │       │   ├── bootstrapNetwork_jti6srapp942f980e5un4aq54v.txt
    │       │   ├── bootstrapNetwork_l3fqtjj5v4gojpo7s7gjhv1ici.txt
    │       │   ├── bootstrapNetwork_l4p2t9hfverncrf8684hqmr7qk.txt
    │       │   ├── bootstrapNetwork_l52atfdpvo6i5t688aok3ekt0l.txt
    │       │   ├── bootstrapNetwork_l7bisla403hguuv8ibcmc6cimn.txt
    │       │   ├── bootstrapNetwork_mc96ttgftuub3fdg0dumpajvua.txt
    │       │   ├── bootstrapNetwork_meiet3eq6895s14gifil22dlcc.txt
    │       │   ├── bootstrapNetwork_mfrmspd46ik0l2tgsg6rba5aie.txt
    │       │   ├── bootstrapNetwork_nm2itnhq4nbpil2okkcm164d83.txt
    │       │   ├── bootstrapNetwork_nnbqtde452mkbmroul0oadu2m4.txt
    │       │   ├── bootstrapNetwork_npl2sjcedc1j4ogp0nkqilno46.txt
    │       │   ├── bootstrapNetwork_nquas98odmcdtq9pao8srthda8.txt
    │       │   ├── bootstrapNetwork_o062tdercsrj0o879hkdd86svm.txt
    │       │   ├── bootstrapNetwork_o1fasjb5l66dpq17ji8fmfui5o.txt
    │       │   ├── bootstrapNetwork_o6dmthiuk2ak36al89vtuiscd5.txt
    │       │   ├── bootstrapNetwork_o7mut7f8kclis81liajs7qm1j7.txt
    │       │   ├── bootstrapNetwork_o906tddikm0dlpolsc82gifn18.txt
    │       │   ├── bootstrapNetwork_pavesnc5jm91gdvfln2cl3pa9e.txt
    │       │   ├── bootstrapNetwork_pg72trg8iso6irvtkfe16uepmt.txt
    │       │   ├── bootstrapNetwork_phgat1eir631btmtuh23fm6f4v.txt
    │       │   ├── bootstrapNetwork_pipisncsrgds5ffu0im5oe04b0.txt
    │       │   ├── bootstrapNetwork_qp0etlhipl5l21l60ms0ea170l.txt
    │       │   ├── bootstrapNetwork_qr9mtbfspvgjrja6aog6n1osmn.txt
    │       │   ├── bootstrapNetwork_qsiusha72arekl36kp4509ihsp.txt
    │       │   ├── bootstrapNetwork_s2gitpiio588olje2rm1ddpv4c.txt
    │       │   ├── bootstrapNetwork_s3pqtfgsofj7hn8ecta7m5hkie.txt
    │       │   ├── bootstrapNetwork_s432t5d70pu2b91emuu5utb9of.txt
    │       │   ├── bootstrapNetwork_s6casrbh138t4aoep0i87l4veh.txt
    │       │   ├── bootstrapNetwork_tb9utjht6uln8b6mf244l9acm4.txt
    │       │   ├── bootstrapNetwork_tdj6t9e7790i1svmh3oath2246.txt
    │       │   ├── bootstrapNetwork_tesesvch7jbcqummr5c96otna7.txt
    │       │   ├── bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt
    │       │   ├── bootstrapNetwork_umcit3fhe2e0h2kv5a6a5skfdu.txt
    │       │   ├── bootstrapNetwork_uolqspdrecova4bvfbqce4e4s0.txt
    │       │   ├── bootstrapNetwork_upv2sf85en3q362vhdeenc7q21.txt
    │       │   ├── bootstrapNetwork_uutetdhudi80cicde45t0f5k1e.txt
    │       │   ├── bootstrapNetwork_vvsmt7ghcigk7mj77f0b4gd7hl.txt
    │       │   ├── network.txt
    │       │   └── regulons.tsv.gz
    │       ├── regulators.tsv
    │       └── targets.tsv
    ├── aracne_regulons_CardosoMoreira2020_1000-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── aracne_regulons_CardosoMoreira2020_100-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── aracne_regulons_CardosoMoreira2020_2000-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── aracne_regulons_CardosoMoreira2020_5000-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── aracne_regulons_CardosoMoreira2020_500-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── aracne_regulons_CardosoMoreira2020-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── aracne_regulons_combined-EX
    │   ├── CardosoMoreira2020.tsv.gz
    │   ├── PANCAN_PT.tsv.gz
    │   └── PANCAN_STN.tsv.gz
    ├── aracne_regulons_PANCAN_PT-EX
    │   └── PANCAN_PT.tsv.gz
    ├── aracne_regulons_PANCAN_STN-EX
    │   └── PANCAN_STN.tsv.gz
    ├── dPSIthresh10_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh15_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh20_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh25_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh30_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh35_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh40_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh45_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh50_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── dPSIthresh5_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── experimentally_derived_regulons_raw-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── experimental_without_postar3_regulons-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── mlr_and_experimental_regulons-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── mlr_regulons
    │   ├── CardosoMoreira2020-EX
    │   │   └── pruned
    │   │       └── regulons.tsv.gz
    │   ├── CardosoMoreira2020-genexpr
    │   │   └── pruned
    │   │       └── regulons.tsv.gz
    │   ├── PANCAN_PT-EX
    │   │   └── pruned
    │   │       └── regulons.tsv.gz
    │   └── PANCAN_STN-EX
    │       └── pruned
    │           └── regulons.tsv.gz
    ├── mlr_regulons_CardosoMoreira2020-EX
    │   └── CardosoMoreira2020.tsv.gz
    ├── mlr_regulons_combined-EX
    │   ├── CardosoMoreira2020.tsv.gz
    │   ├── PANCAN_PT.tsv.gz
    │   └── PANCAN_STN.tsv.gz
    ├── mlr_regulons_PANCAN_PT-EX
    │   └── PANCAN_PT.tsv.gz
    ├── mlr_regulons_PANCAN_STN-EX
    │   └── PANCAN_STN.tsv.gz
    ├── postar3_and_experimental_regulons-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── postar3_clip_regulons-EX
    │   ├── POSTAR3-metaexperiment0-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment10-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment11-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment12-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment13-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment14-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment15-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment16-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment17-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment18-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment19-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment1-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment2-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment3-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment4-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment5-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment6-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment7-clip_peak.tsv.gz
    │   ├── POSTAR3-metaexperiment8-clip_peak.tsv.gz
    │   └── POSTAR3-metaexperiment9-clip_peak.tsv.gz
    ├── regulon_evaluation_labels
    │   ├── ENASFS.tsv.gz
    │   ├── ENCOREKD_HepG2.tsv.gz
    │   ├── ENCOREKD_K562.tsv.gz
    │   ├── ENCOREKO_HepG2.tsv.gz
    │   └── ENCOREKO_K562.tsv.gz
    ├── regulon_evaluation_scores
    │   ├── correlation_pearson
    │   │   ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   └── top90_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   ├── correlation_spearman
    │   │   ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   └── top90_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   ├── gsea
    │   │   ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── aracne_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── experimental_without_postar3_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── mlr_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_clip_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_clip_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_clip_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_clip_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── postar3_clip_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── splicinglore_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │   │   ├── top90_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │   │   └── top90_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │   ├── merged-EX.tsv.gz
    │   ├── merged-genexpr.tsv.gz
    │   └── viper
    │       ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_and_mlr_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_1000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_1000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_100-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_100-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_2000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_2000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_5000-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_5000-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_500-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020_500-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── aracne_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── aracne_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh10_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh15_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh20_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh25_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh30_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh35_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh45_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── dPSIthresh5_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimental_without_postar3_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimental_without_postar3_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimental_without_postar3_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimental_without_postar3_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── experimental_without_postar3_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_CardosoMoreira2020-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_CardosoMoreira2020-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_CardosoMoreira2020-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_combined-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_combined-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_combined-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_combined-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_combined-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENASFS-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-one_tailed.tsv.gz
    │       ├── mlr_regulons_development-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_PT-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_PT-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_PT-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_PT-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_PT-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_STN-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_STN-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_STN-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_STN-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── mlr_regulons_PANCAN_STN-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_and_experimental_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_and_experimental_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_and_experimental_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_and_experimental_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_and_experimental_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_clip_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_clip_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_clip_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_clip_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── postar3_clip_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── splicinglore_regulons-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── splicinglore_regulons-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── splicinglore_regulons-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── splicinglore_regulons-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── splicinglore_regulons-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top100_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top100_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top100_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top40_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top40_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top40_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top50_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top50_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top50_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top60_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top60_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top60_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top70_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top70_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top70_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top80_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top80_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top80_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top90_experimentally_derived_regulons_pruned-ENASFS-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top90_experimentally_derived_regulons_pruned-ENCOREKD_K562-EX-shadow_no-two_tailed.tsv.gz
    │       ├── top90_experimentally_derived_regulons_pruned-ENCOREKO_HepG2-EX-shadow_no-two_tailed.tsv.gz
    │       └── top90_experimentally_derived_regulons_pruned-ENCOREKO_K562-EX-shadow_no-two_tailed.tsv.gz
    ├── regulon_properties
    │   ├── dPSIthresh-regulators_per_target-EX.tsv.gz
    │   ├── dPSIthresh-targets_per_regulator-EX.tsv.gz
    │   ├── regulators_per_target-EX.tsv.gz
    │   ├── regulators_per_target-genexpr.tsv.gz
    │   ├── targets_per_regulator-EX.tsv.gz
    │   └── targets_per_regulator-genexpr.tsv.gz
    ├── regulons_eda_gsea
    │   └── experimentally_derived_regulons_pruned-EX.tsv.gz
    ├── regulons_eda_jaccard
    │   └── experimentally_derived_regulons_pruned-EX.tsv.gz
    ├── splicinglore_regulons-benchmarkable.tsv.gz
    ├── splicinglore_regulons-EX
    │   ├── SplicingLore-metaexperiment0-delta_psi.tsv.gz
    │   ├── SplicingLore-metaexperiment1-delta_psi.tsv.gz
    │   ├── SplicingLore-metaexperiment2-delta_psi.tsv.gz
    │   ├── SplicingLore-metaexperiment3-delta_psi.tsv.gz
    │   └── SplicingLore-metaexperiment4-delta_psi.tsv.gz
    ├── top100_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── top40_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── top50_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── top60_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── top70_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    ├── top80_experimentally_derived_regulons_pruned-EX
    │   ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
    │   ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
    │   ├── ENCOREKD-HepG2-delta_psi.tsv.gz
    │   ├── ENCOREKD-K562-delta_psi.tsv.gz
    │   ├── ENCOREKO-HepG2-delta_psi.tsv.gz
    │   └── ENCOREKO-K562-delta_psi.tsv.gz
    └── top90_experimentally_derived_regulons_pruned-EX
        ├── ENASFS-metaexperiment0-delta_psi.tsv.gz
        ├── ENASFS-metaexperiment1-delta_psi.tsv.gz
        ├── ENASFS-metaexperiment2-delta_psi.tsv.gz
        ├── ENASFS-metaexperiment3-delta_psi.tsv.gz
        ├── ENCOREKD-HepG2-delta_psi.tsv.gz
        ├── ENCOREKD-K562-delta_psi.tsv.gz
        ├── ENCOREKO-HepG2-delta_psi.tsv.gz
        └── ENCOREKO-K562-delta_psi.tsv.gz
```