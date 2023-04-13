import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### RULES #####
rule all:
    input:
        # RNA seq perturbation profiles
        ## create ENA query to find samples with splicing factors
        ".done_query_sfs",
        ## download ENA metadata table to find samples
        expand(os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments-{field}.tsv"), field=["alias","title"]),
        ## merge tables
        os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments_merged.tsv.gz"),
        ## clean SF ENA metadata
        ## create ENA query to get corresponding controls
        
        
rule create_ena_query_sfs:
    input:
        sfs = os.path.join(SUPPORT_DIR,"splicing_factors.tsv")
    output:
        touch(".done_query_sfs")
    run:
        import pandas as pd
        
        sfs = pd.read_table(input.sfs)
        
        # alias
        sfs_query = (
            " OR ".join((
                'sample_alias="*'+sfs["GENE"]+'*" OR ' +\
                'run_alias="*'+sfs["GENE"]+'*"'
            ).values)
        )
        query = 'tax_eq(9606) AND library_strategy="RNA-Seq" AND library_source="TRANSCRIPTOMIC" AND ('+sfs_query + ")"
        print(query)        
        
        # title
        sfs_query = (
            " OR ".join((
                #'study_title="*'+sfs["GENE"]+'*" OR ' +\
                'sample_title="*'+sfs["GENE"]+'*" OR ' +\
                'experiment_title="*'+sfs["GENE"]+'*"'
            ).values)
        )
        query = 'tax_eq(9606) AND library_strategy="RNA-Seq" AND library_source="TRANSCRIPTOMIC" AND ('+sfs_query + ")"
        print(query)
                
        print("Done!")
        
        
rule download_metadata_sfs:
    input:
        ".done_query_sfs"
    output:
        os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments-{field}.tsv")
    params:
        field = "{field}"
    shell:
        """
        bash scripts/ENA-download_sf_exps_metadata-{params.field}.sh > {output}
        
        echo "Done!"
        """
    
rule merge_ena_sfs_metadata:
    input:
        metadata_alias = os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments-alias.tsv"),
        metadata_title = os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments-title.tsv"),
        sfs = os.path.join(SUPPORT_DIR,"splicing_factors.tsv")
    output:
        metadata = os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments_merged.tsv.gz")
    threads: 24
    run:
        import pandas as pd
        from joblib import Parallel, delayed
        from tqdm import tqdm
        
        metadata_alias = pd.read_table(input.metadata_alias, low_memory=False)
        metadata_title = pd.read_table(input.metadata_title, low_memory=False)
        sfs = pd.read_table(input.sfs)
        threads = threads
        
        # combine metadata
        metadata = pd.concat([metadata_title, metadata_alias]).drop_duplicates()
        
        # filter by instrument platform
        metadata = metadata.loc[metadata["instrument_platform"]=="ILLUMINA"]
        
        # filter by whether we can download the data
        metadata = metadata.loc[~metadata["fastq_ftp"].isnull()].copy()

        # filter out studies we know we don't want
        studies_todrop = [
            "PRJNA379992", # rheumatoid patients
            "PRJNA486023", # triple negative breast cancer patients
            "PRJNA325938", # cancer cell lines
            "PRJNA495476", # HBV-induced hepatocellular carcinoma (single cell)
            "PRJEB38744", # timecourse COVID
            "PRJNA412363", # cancer cell lines in different extracellular conditions
            "PRJEB32263", # modulation innate immune response during infection
            "PRJNA523380", # the CCLE
            "PRJEB38015", # organoids with mutations
            "PRJEB33498", # KDELR1 KO
            "PRJEB38070", # KO of INPP5D (leukemia)
            "PRJEB42732", # cells treated with drugs at different timepoints
            "PRJEB13892", # drugged xenografs patient derived
            "PRJEB47484", # KLF17 KO timecourse
            "PRJNA818797", # PD1/CTLA4 immunotherapy
            "PRJNA842677", # PDX + patient tumor + normal tissue
            "PRJNA417221", # lisozyme treatment
            "PRJNA398556", # pregnancy
            "PRJEB19655", # AML
            "PRJNA701374", # papilary cancer paediatric
            "PRJNA339968", # monocytes
            "PRJNA722536", # pregnancy
            "PRJEB41156", # time course MCF7 grouwth
            "PRJEB36353", # EMT
            "PRJEB48549", # treatment with GTPP
            "PRJNA884818", # Covid
            "PRJNA542922", # viral infection
            "PRJNA672723", # infection
            "PRJNA807823", # infection
            "PRJNA881214", # covid
            "PRJNA677887", # Epstein Barr Virus
            "PRJEB35113", # MEF2D-HNRNPUL1 fusion
            "PRJNA172563", # lymphoma samples
            "PRJNA356760", # Zika virus
            "PRJNA698960", # Zika virus
            "PRJNA638035", # Zika virus
            "PRJEB43894", # endometrial cancer cell line Ishikawa treated with BMP2
            "PRJNA888457", # MEN1-KO NCI-H460 cells
            "PRJNA782615", # schizophrenia
            "PRJEB32530", # influenza
            "PRJNA559719", # screen in HAP cells
            "PRJNA507472", # virus
            "PRJNA34535", # Human Reference Epigenome Mapping Project
            "PRJNA574055", # SMARCB1 epigenetics
            "PRJNA30709", # ENCODE
            "PRJNA504037", # RNA seq different cells
            "PRJNA883202", # ZIKA
            "PRJNA315289", # human preadipocytes
            "PRJNA667475", # virus
            "PRJNA883202", # ZIKA
            "PRJNA726637", # drug
            "PRJNA606334", # zika
            "PRJNA662490", # zika
            "PRJEB13167", # parkinson's
            "PRJNA374549", # zika
            "PRJEB4337", # HPA RNA-seq normal tissues
            "PRJNA669277", # leukemia treatment
            "PRJNA642003", # innate immunity
            "PRJEB36921", # Tau protein
            "PRJNA665692", # PTBP1 with HALO tag, covid
            "PRJNA665581", # ribodeplete
            "PRJNA260224", # not interesting
            "PRJNA379358", # virus
            "PRJNA662366", # zika
            "PRJEB33642", # neurodegeneration
            "PRJNA513099", # CDKN1A
            "PRJNA859517", # BHLHE40
            "PRJNA418335", # NF95
            "PRJNA517339", # circulating free mRNA
            "PRJNA321028", # endometrial carcinoma
            "PRJNA522052", # immune cells
            "PRJNA344703", # zika
            "PRJNA396095", # io torrent
            "PRJNA645773", # anemia
            "PRJEB55015", # PBMCs
            "PRJNA811414", # metabolic labeling methylation
            "PRJNA659884", # shikonin
            "PRJNA335350", # translation KD IBTK
            "PRJNA764684", # healthy human tissues
            "PRJNA932451", # RBFOX1 cardiomyocites
            "PRJNA633293", # PIKIII treatment in RKO cells
            "PRJNA748615", # parkinson
            "PRJNA830988", # KD KDKHNYN
            "PRJNA751732", # IKZF1
            "PRJNA764368", # human keloid pathogenesis
            "PRJNA275642", # virus
            "PRJNA795888", # VISTA
            "PRJNA782183", # autosomy
            "PRJNA522627", # osteonecrosis
            "PRJNA412440", # fusion cells
            "PRJNA477380", # skin
            "PRJEB31881", # CD4 differentiaiton
            "PRJNA594493", # NK cells
            "PRJNA382346", # alzheimer
            "PRJNA209481", # ovarian cancer cell lines
            "PRJNA314266", # co culture with bacteria
            "PRJNA415678", # stress granules
            "PRJNA414471", # microcephaly
            "PRJEB25956", # IDH2 mutated acute myeloid
            "SRP006674", # T reg cells
            "PRJNA709641", # neural crest
            "PRJNA679773", # solid tumors
            "PRJEB40462", # cancer placenta antigen
            "PRJEB47471", # lateral sclerosis
            "PRJEB49820", # KD JMJD6
            "PRJNA396985", # Resistance to tyrosine kinase inhibitors (TKIs)
            "PRJNA596245", # BCR sequencing
            "PRJNA636928", # DDX53
            "PRJEB57753", # hydrocortisone treatment
            "PRJNA564062", # placental metastasis
            "PRJNA591091", # neural crest
            "PRJNA612416", # cell line HEK293
            "PRJNA662434", # blood
            "PRJNA681741", # patient tumor sample
            "PRJEB25108", # enteroids
            "PRJNA722770", # microglia
            "PRJNA731139", # drug combinations, low read count
            "PRJNA63443", # ENCODE
            "PRJNA799044", # HSF1,2 in different cell lines
            "PRJNA532934", # metabolism
            "PRJNA555321", # ES cells
            "PRJEB49074", # enteroids drug toxicity
            "PRJNA605264", # clonal human renal
            "PRJNA612305", # tumor patient samples
            "PRJNA614494", # EMT timecourse
            "PRJNA729731", # nuclear export combination therapies
            "PRJNA310000", # breast cancer
            "PRJDB5361", # TF OE in hPSCs
            "PRJNA589697", # lung epithelial
            "PRJNA263687", # unperturbed human sample
            "PRJNA765665", # sammson gene
            "PRJEB54896", # patients
        ]
        metadata = metadata.loc[~metadata["study_accession"].isin(studies_todrop)].copy()
        
        # filter by whether run_alias, sample_alias, sample_title, experiment_title, or study_title
        # contain any of the splicing factors
        cols_oi = ["run_alias","sample_alias","sample_title","experiment_title","study_title"]
        metadata[cols_oi] = metadata[cols_oi].apply(lambda col: col.str.upper())
        
        query = "|".join(sfs["GENE"]).upper()
        new_cols = ["match_in_"+col for col in cols_oi]
        metadata[new_cols] = metadata[cols_oi].apply(lambda col: col.str.contains(query))
        metadata = metadata.loc[metadata[new_cols].any(axis=1)].copy()
        
        # inform on the match
        def get_gene_matches(row, cols_oi, genes_oi):
            result = {}
            for col_oi in cols_oi:
                result["found_sfs_in_"+col_oi] = ",".join(
                    [gene for gene in genes_oi if row[[col_oi]].str.contains(gene).any()]
                )
            result = pd.Series(result)
            row = pd.concat([row,result])
            return row
            
        found_sfs = Parallel(threads)(
                delayed(get_gene_matches)(
                    row,
                    cols_oi,
                    sfs["GENE"]
                )
                for idx, row in tqdm(metadata.iterrows(), total=metadata.shape[0])
            )
        metadata = pd.DataFrame(found_sfs).copy()
        
        # save
        metadata.to_csv(output.metadata, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule clean_ena_metadata:
    input:
        metadata = os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments_merged.tsv.gz"),
        sfs = os.path.join(SUPPORT_DIR,"splicing_factors.tsv")
    output:
        metadata = os.path.join(RAW_DIR,"ENA","splicing_factors","metadata","sf_experiments_clean.tsv.gz")
    run:
        import pandas as pd
        
        # load
        metadata = pd.read_table(input.metadata, low_memory=False)
        sfs = pd.read_table(input.sfs)
        
        # combine iformation of matches
        cols_oi = [col for col in metadata.columns if "found_sfs" in col]
        metadata["found_sfs"] = metadata[cols_oi].apply(
            lambda row: list(set(",".join(row.dropna().unique()).split(","))), axis=1
        )
        metadata = metadata.explode("found_sfs").drop_duplicates().copy()

        # consider only those SFs hand curated, in Rogalska, in Hegele
        idx = sfs["source"].isin(["Rogalska2022","hand-curated","Hegele2012"]) & sfs["in_hegele2012"]
        sfs = sfs.loc[idx].copy()
        metadata = metadata.loc[metadata["found_sfs"].isin(sfs["GENE"])].copy()

        # filter out bad matches
        cols_oi = ["run_alias","sample_alias","sample_title","experiment_title","study_title"]
        is_bad_match = metadata[cols_oi].apply(
            lambda col: col.str.contains(
                "PREGNA|IMMUN|CD4|CD8|VIRUS|ZIKA|SCHIZO|ALZHEI|COVID|INFECTION|CREST|SINGLE-CELL|SINGLE CELL|T CELL|SKIN|NEUTROPH|FETAL|PLACEN|MONONU|EOSINOPHILS|PARKINSO|DIFFERENTIAT|PAEDIATRIC|INTERFERON|IPS|BUDESONIDE|TRANSPOSON|CORTISON|CITOKIN|ASTHMA|EPIDERMIS|INTERLEUKINS|FUSION GENES|MONOCYTE|TUBERCULOSIS|STEM CELL|ANGIOGE|EOSINOPHIL|ZIKV"
            )
        ).any(axis=1)
        metadata = metadata.loc[~is_bad_match].copy()
        
        # studies to drop
        studies_todrop = [
            "PRJNA809587", # Riboseq
            "PRJNA597343", # 
        ]
        metadata = metadata.loc[~metadata["study_accession"].isin(studies_todrop)].copy()        
        
        # how many SFs we found data for?
        found_sfs = metadata.loc[metadata["found_sfs"].isin(sfs["GENE"]),"found_sfs"].unique()
        new_sfs = set(found_sfs) - set(sfs.loc[~sfs["in_encore_kd"],"GENE"])
        missing_sfs = set(sfs["GENE"]) - set(found_sfs)
        print("Found %s SFs in total, %s of them not in ENCORE KDs. We are missing %s SFs" %\
              (len(found_sfs), len(new_sfs), len(missing_sfs)))
        
        # consider only studies with high diversity in SF matches, how many/which ones are we missing? Are they covered by ENCORE KDs?
        # consider larger studies, how much are we covering?
        n_sfs = metadata[["study_accession","found_sfs"]].drop_duplicates()["study_accession"].value_counts()
        studies_oi = set()
        prev_found_missing = set()
        prev_found_total = set()
        for top_n in range(len(n_sfs)):
            studies = set([n_sfs.index[top_n]])
            found_total = metadata.loc[
                metadata["study_accession"].isin(studies), "found_sfs"
            ].unique()
            found_missing = set(new_sfs).intersection(set(found_total))
            
            if len(found_missing - prev_found_missing) > 0:  
                print(top_n)
                # save the study if we found something new
                studies_oi = studies_oi.union(studies)
                prev_found_missing = prev_found_missing.union(found_missing)
                prev_found_total = prev_found_total.union(found_total)

        print("""
            Taking %s studies we are would be covering a total of %s different SFs, 
            %s out of %s for which we don't have data now.
        """ % (len(studies_oi), len(prev_found_total), len(prev_found_missing), len(new_sfs)))
        
        print("Selected studies:\n", n_sfs[studies_oi].sort_values())
        
        metadata.loc[metadata["study_accession"].isin(studies_oi), ["study_title","study_accession"]].drop_duplicates().values
        
        print(metadata.loc[metadata["study_accession"].isin(studies_oi),"found_sfs"].value_counts().to_string())
        
        
        # problematics: IK, KIN, SON      
        
        studies_tokeep = [
            "PRJDB6952", # drug perturbation
            "PRJNA642291", # KD RBPs
            "PRJNA494172", # KD SMNDC1
            "PRJNA674890", # OE SNRPA1
            "PRJNA420970", # KD SFs
            "PRJNA290340", # KD SFs
        ] + [
            "PRJEB22885", # PRPF31 mut https://doi.org/10.1038%2Fs41467-018-06448-y
            "PRJNA321560", # CLK inhibitor, KDs https://doi.org/10.1038/s41467-016-0008-7
            "PRJEB46899", # siRNA screen https://doi.org/10.1161/CIRCRESAHA.120.318141
            "PRJEB33667", # iCLIP + KDs
            "PRJNA384899", # TRRAP knockdowns https://doi.org/10.1083%2Fjcb.201706106
            "PRJEB36354", # WT1 knockdown
            "PRJNA300429", # CPSF6 KO
            "PRJNA388736", # CLIP + KD EIF4A3
            "PRJNA551123", # SRRM4 OE
            "PRJEB23439", # KO PTBP1/2 MATR3(fractionated)
            "PRJNA733583", # KD DDX42
            "PRJEB6651", # U2AF1 depletion
            "PRJNA394049", # HNRNPK depletion (fractionated)
            "PRJNA394053", # SRSF6 KD
            "PRJEB39343", # PTBP1, ESRP2, and MBNL1 KDs
            "PRJEB41023", # LUC7L KD
            "PRJNA530774", # drug CDKs
            "PRJEB57225", # KO RBM15
            "PRJNA400256", # KO TIA1
            "PRJEB41023", # KD LUC7L
            "PRJNA738336", # KD SFs
            "PRJNA223244", # KD SFs
            "PRJNA561411", # KD SFs
            "PRJNA738336", # KD SFs
            "PRJNA521683", # KD RBFOX2
            "PRJNA678286", # KD HNRNPF
            "PRJNA624911", # KO SMN2
            "PRJEB41929", # KD ERN1
            "PRJEB44104", # KD PPARGC1A
            "PRJNA283786", # OE RBM8A
            "PRJEB41850", # KD KEAP1 ZRANB2
            "PRJNA789180", # DDX27
            "PRJNA645345", # EIF4A3
            "PRJEB41777", # ZRANB2
            "PRJEB41872", # KD NCBP1
            "PRJNA523954", # KD SRSF3
            "PRJNA782063", # KD RNF213
            "PRJNA779469", # KD RNF213 (weird)
            "PRJEB53909", # KD RBM42
            "PRJNA622794", # PRPF31 KO
            "PRJNA886829", # KD SART3
            "PRJNA474911", # OE SRRM4
            "PRJNA789153", # KO CDK12
            "PRJNA797682", # KO SNIP1
            "PRJNA809499", # KD USP22
            "PRJEB4007", # KD SON
            "PRJNA669300", # KO PPIL1
            "PRJEB23554", # KD DDX54
            "PRJEB28561", # KD RNF40
            "PRJEB25085", # KD MSI1 MSI2
            "PRJNA540831", # KO PRPF40B
            "PRJNA689247", # KO SRSF7
            "PRJNA633677", # KD UPF1
            "PRJNA520804", # KD CPSF6
            "PRJEB36921", # KD HNRNPH1
            "PRJNA655345", # KD HNRNPA2B1
            "PRJNA353381", # KD RPUSD4
            "PRJNA610182", # KD HNRNPL
            "PRJNA361121", # mutations RBM39
            "PRJNA258436", # KD PRMT1 and PRMT5
            "PRJNA722728", # KD YTHDC1
            "PRJNA915416", # KD METTL3
            "PRJNA781783", # KD METTL3
            "PRJNA428150", # KD SFs
            "PRJNA241095", # KD RBPMS
            "PRJNA293234", # SRSF1
            "PRJNA505781", # KD SFPQ
            "PRJNA471771", # KD AHNAK2 EVPL
            "PRJEB7544", # KD PRPF8
            "PRJNA850080", # KD PTBP1
            "PRJNA314922", # KD PTBP1 HNRNPK
            "PRJEB33860", # KD USP22
            "PRJEB33859", # KD USP22
            "PRJNA943438", # KD WTAP
            "PRJNA551046", # KD MOV1 TRA1
            "PRJNA797585", # OE RBFOX2
            "PRJEB3048", # KD HNRNPC
            "PRJEB29040", # KD EIF4A3
            "PRJNA523477", # KD NOVA2
            "PRJNA674660", # KD SF3B1
            "PRJNA682912", # KD CPSF1
            "PRJNA639929", # KD PABPC1
            "PRJNA278702", # KD CWC22
            "PRJNA722701", # KD USP22
            "PRJNA720522", # RBM10
            "PRJNA764507", # KD NCL
            "PRJNA309732", # PTBP1 and 2 KD
            "PRJNA525954", # KD WTAP
            "PRJEB34009", # KD USP22
            "PRJNA510082", # KD CELF2
            "PRJNA579088", # KD EIF4A3
            "PRJNA780251", # KD KDM1A
            "PRJNA401938", # KD EIF4A3
            "PRJEB48601", # KD DHX38
            "PRJNA218550", # RBFOX1,2,3
            "PRJEB40550", # KD MSI1
            "PRJNA610157", # KO KAT2B
            "PRJEB34918", # KD KLF5
            "PRJEB34048", # KD USP22
            "PRJEB32567", # KD MBNL1
            "PRJNA164615", # KD MSI1
            "PRJEB38971", # KD RNF40
            "PRJEB10298", # KD SNRPB
            "PRJEB9546", # KD HNRNP + iCLIP
            "PRJEB14614", # KD SF1
            "PRJNA301134", # KD HNRNP
            "PRJEB29794", # KD USP22
            "PRJNA666753", # KD METTL3
            "PRJNA294972", # KD UPF1
            "PRJNA914900", # KD RBBP6
            "PRJNA774953", # KD SRPK1
            "PRJNA484809", # KD RBMS3 + CLIP
            "PRJNA752633", # KD XRN2
            "PRJNA494617", # KD METTL3
            "PRJNA706906", # KD SRSF6 + Bleomycin
            "PRJNA482875", # KD METTL3
            "PRJNA660570", # KD U2AF1
            "PRJEB30370", # KD SMU1 RED IK MFAP1
            "PRJEB32489", # KD YBX3
            "PRJEB49865", # KD HNRNPL
            "PRJNA234443", # KD RBFOX2
            "PRJNA386251", # KD RBM10
            "PRJNA413916", # KD SRSF5
            "PRJNA414346", # KD SRSF3
            "PRJNA422160", # (KD) CLIP RBM10
            "PRJNA185008", # KD FUS
            "PRJEB8225", # KD QKI
            "PRJNA591294", # KD UPF1 and G3BP1
        ] + [
            "PRJNA683080", # Lu2021, indisulam treatment
            "PRJNA140779", # KD ELAVL1, https://doi.org/10.1016/j.molcel.2011.06.008 (listed in Desai2020)
            "PRJNA192838", # KD and OE RBM10, https://doi.org/10.1093/nar/gkx508 (listed in Desai2020)
        ]
        
        # filter by library selection?
        
        print("Done!")

