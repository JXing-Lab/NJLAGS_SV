"""
CNV_Summarize.py
"""

# FIND MEDIAN LENGTH OF CNVs
def medianfinder(input_tsv):
    "Running medianfinder() method..."
    
    # import pandas
    import pandas as pd
    
    # load in tsv
    df_cnv = pd.read_csv(input_tsv, sep = "\t")
    # calculate length
    df_cnv["length"] = df_cnv["stop"] - df_cnv["start"]
    # calculate median
    cnv_median = df_cnv["length"].median()
    
    "medianfinder() method finished running."
    
    return(cnv_median)

 
# FIND # CALLED PATHOGENIC FOR StrVCTVRE AND SvAnna
def pathocounter(input_tsv):
    "Running pathocounter() method..."
    
    # import pandas
    import pandas as pd
    import numpy as np
    
    # declare counting variables
    StrVCTVRE_count = 0
    SvAnna_count = 0
    
    # load in tsv
    df_cnv = pd.read_csv(input_tsv, sep = "\t")
    
    # count pathogenic for StrVCTVRE
    df_cnv["StrVCTVRE"] = df_cnv["StrVCTVRE"].str.replace("liftOver_to_GRCh38_failed", "0")
    df_cnv["StrVCTVRE"] = df_cnv["StrVCTVRE"].str.replace("not_exonic", "0")
    StrVCTVRE_column = df_cnv["StrVCTVRE"].astype("float")
    StrVCTVRE_count = StrVCTVRE_column[StrVCTVRE_column > 0].count()
    
    # count pathogenic for SvAnna
    SvAnna_column = df_cnv["SvAnna"].astype("float")
    SvAnna_count = SvAnna_column[SvAnna_column > 0].count()
    
    # output
    counts_list = [StrVCTVRE_count, SvAnna_count]
    return(counts_list)
    
    "pathocounter() method finished running."

# BUILD TABLE METHOD
def build_varsum_table():
    print("Building table...")
    
    import pandas as pd
    
    # open output file
    sum_out = open("../Results/CNV_summary.tsv", "w")
    
    phenotypes = ["ASD_only", "ASD_LI", "ASD_RI"]
    
    # print out header
    header_list = ["Phenotype", "Raw_CNVs", "Prioritized_CNVs", "Median_CNV_length", "#duplications", "#deletions", "#patho_StrVCTVRE", "#patho_SvAnna"]
    header_string = "\t".join(header_list)
    sum_out.write(header_string + "\n")

    
    # iterate
    for a in range(0, len(phenotypes)):
        # import results file
        results_file = "../Results/" + phenotypes[a] + "_results.tsv"
        # import as dataframe
        dfr = pd.read_csv(results_file, sep = "\t")
        
        
        CNV_median = medianfinder(results_file)
        CNV_StrVCTVRE = pathocounter(results_file)[0]
        CNV_SvAnna = pathocounter(results_file)[1]
        results_in = open(results_file, "r")
        results_string = results_in.read()
        results_list = results_string.split("\n")
        raw_CNV_count = 2538
        CNV_count = len(results_list)-2
        del_count = len(dfr[dfr.sv_type == "<CN0>"]) + len(dfr[dfr.sv_type == "<CN1>"])
        dup_count = len(dfr[dfr.sv_type == "<CN3>"]) + len(dfr[dfr.sv_type == "<CN4>"])
        
        # output
        output_list = [phenotypes[a], str(raw_CNV_count), str(CNV_count), str(CNV_median), str(del_count), str(dup_count), str(CNV_StrVCTVRE), str(CNV_SvAnna)]
        output_string = "\t".join(output_list) + "\n"
        sum_out.write(output_string)
        
    # close open file
    sum_out.close()
    
    
    print("Table built.")
    

# create counts
def create_counts(phenotype1, phenotype2):
    
    # import libraries
    import pandas as pd
    
    # read in table as dataframe
    df = pd.read_csv("../Data/Starting_Manifest/alternative_table.tsv", sep = "\t")
    
    # collate counts
    patient_count = len(df[(df["New PARTID"].notnull()) & 
                           ((df[phenotype1] == "2") | 
                           (df[phenotype2] == "2"))])
    male_count =    len(df[(df["SEX"] == "1") & 
                           ((df[phenotype1] == "2") | 
                           (df[phenotype2] == "2"))])
    female_count =    len(df[(df["SEX"] == "2") & 
                           ((df[phenotype1] == "2") |
                           (df[phenotype2] == "2"))])
    family_count = len(df["New Family ID"].unique())
    counts = [patient_count, male_count, female_count, family_count]
    
    # print counts out
    print("For the combination of phenotypes " + phenotype1 + " & " + phenotype2)
    print("Patient count: " + str(patient_count))
    print("Male count: " + str(male_count))
    print("Female count: " + str(female_count))
    print("Family count: " + str(family_count))
    

    
# gmn counts
def gmn_counts(phenotype1, phenotype2):
    
    # import libraries
    import pandas as pd
    
    # read Geminesque output as dataframe
    df = pd.read_csv("../Data/" + phenotype1 + "_" + phenotype2 + "/GMN_out_" + phenotype1 + "_" + phenotype2 + ".tsv", sep = "\t")
    
    
    # collate counts
    CNV_count = df["SVsN"].sum()
    denovo_count = len(df[df["de_novo__single_affected"].notnull()])
    rec_count = len(df[df["autosomal_recessive__single_affected"].notnull()])
    dom_count = len(df[df["autosomal_dominant__single_affected"].notnull()])
    
    # print counts out
    print("For the combination of phenotypes " + phenotype1 + " & " + phenotype2)
    print("CNV count: " + str(CNV_count))
    print("denovo_count: " + str(denovo_count))
    print("rec_count: " + str(rec_count))
    print("dom_count: " + str(dom_count))
    
#gmn_counts("ASD", "RI")

# method to extract gene list from results
def gene_extract(results_tsv, annotsv_file, output1_tsv, output2_tsv, output3_tsv):
    print("Running gene_extract() method...")
    
    # import libraries
    import pandas as pd
    import numpy as np
    
    # read in results file as dataframe
    dfr = pd.read_csv(results_tsv, sep = "\t")
    
    # extract list of genes
    gene_list = dfr["genes"].tolist()
    gene_list = [gene.split("|") for gene in gene_list] # extract genes from multiple genes in same variant
    gene_list = [gene for sublist in gene_list for gene in sublist] # collapse list
    
    # extract list of transcripts
    tx_list = dfr["transcript_id"].tolist()
    tx_list = [tx.split("|") for tx in tx_list]
    tx_list = [tx for sublist in tx_list for tx in sublist]
    
    # extract list of overlap CDS percentages
    OCDS_list = dfr["overlap_CDS_percent"].tolist()
    OCDS_list = [OCDS.split("|") for OCDS in OCDS_list]
    OCDS_list = [OCDS for sublist in OCDS_list for OCDS in sublist]
    
    
    # build new dataframe
    dfo = pd.DataFrame(
                        {"gene": gene_list,
                         "transcript_id": tx_list,
                         "overlap_CDS_percent": OCDS_list}
                        )
    
    # drop duplicates
    dfo = dfo.drop_duplicates(subset = ["gene"])
    
    # output1 prior to CDS filtering 
    dfo.to_csv(output1_tsv, sep = "\t")
    
    # filter on CDS > 0
    dfo = dfo[dfo["overlap_CDS_percent"].astype("int") > 0]
    dfo.to_csv(output2_tsv, sep = "\t")

    # select protein coding
    dfo["protein_coding"] = np.where(dfo["transcript_id"].str.contains("NM"), True, False)
    dfo = dfo[dfo["protein_coding"] == True]    
    
    # output2 for genes for network analysis
    dfo.to_csv(output3_tsv, sep = "\t")
    
    
    print("gene_extract() method finished running.")
 


# RESULTS SORTER METHOD
def sort_results():
    print("Running sort_results() method...")
    
    import pandas as pd

    # import the files
    ASD_only_in = "../Results/intermediary_files/result_genes_3_ASD_only.txt"
    ASD_LI_in = "../Results/intermediary_files/result_genes_3_ASD_LI.txt"
    ASD_RI_in = "../Results/intermediary_files/result_genes_3_ASD_RI.txt"

    # read in as pandas dataframes
    dfao = pd.read_csv(ASD_only_in, sep = "\t")
    dfal = pd.read_csv(ASD_LI_in, sep = "\t")
    dfar = pd.read_csv(ASD_RI_in, sep = "\t")

    # sort alphabetically
    dfao = dfao.sort_values("gene")
    dfal = dfal.sort_values("gene")
    dfar = dfar.sort_values("gene")

    # output all
    ao_genes = dfao["gene"].unique().tolist()
    al_genes = dfal["gene"].unique().tolist()
    ar_genes = dfar["gene"].unique().tolist()

    ao_out = open("../Results/ASD_only_genes_TPM_UNFILTERED.txt", "w")
    al_out = open("../Results/ASD_LI_genes_TPM_UNFILTERED.txt", "w")
    ar_out = open("../Results/ASD_RI_genes_TPM_UNFILTERED.txt", "w")

    ao_out.write("\n".join(ao_genes))
    al_out.write("\n".join(al_genes))
    ar_out.write("\n".join(ar_genes))
    
    # get unique genes
    all_genes = ao_genes + al_genes + ar_genes
    df_all = pd.DataFrame(all_genes, columns = ["genes"])
    df_all = df_all.sort_values("genes")
    uniq_gene_count = str(len(df_all["genes"].unique())) 
    uniq_genes_out = open("../Results/cnv_uniq_genes.txt", "w")
    uniq_genes_list = df_all["genes"].unique()
    uniq_genes_out.write("\n".join(uniq_genes_list))
    
    
    # print out gene counts
    print("ASD_only gene count: " + str(len(ao_genes)))
    print("LI* gene count: " + str(len(al_genes)))
    print("RI* gene count: " + str(len(ar_genes)))
    print("Unique gene count: " + uniq_gene_count)

    print("sort_results() method finished running.")


# FAMILY COUNTER METHOD
def count_families(results_tsv, gene_in, fam_out):
    import pandas as pd
     
    # read in results
    dfr = pd.read_csv(results_tsv, sep = "\t")
    dfg = pd.read_csv(gene_in).values.tolist()
    dfg = [''.join(ele) for ele in dfg]
    
    match_list = []
    # require genes be expressed
    for row in dfr.iterrows():
        row_list = row[1].tolist()
        row_genes = row_list[6].split("|")
        gene_match = any(elem in row_genes for elem in dfg)
        match_list.append([row_list[2], gene_match])
    dfm = pd.DataFrame(match_list, columns = ["stop", "expressed"])
    dfr = pd.merge(dfr, dfm)
    dfr = dfr[dfr["expressed"] == True]
    
    # extract columns as lists
    dn_list = dfr["de_novo_anno"].dropna().tolist()
    ar_list = dfr["auto_rec_anno"].dropna().tolist()
    ad_list = dfr["auto_dom_anno"].dropna().tolist()
    
    # extract families
    fam_list = []
    
    for a in range(0, len(dn_list)):
        fam_id = dn_list[a].split("family")[1].split("[")[0].strip(":")
        if fam_id not in fam_list:
            fam_list.append(fam_id)
    for b in range(0, len(ar_list)):
        fam_id = ar_list[b].split("family")[1].split("[")[0].strip(":")
        if fam_id not in fam_list:
            fam_list.append(fam_id)
    for a in range(0, len(ad_list)):
        fam_id = ad_list[a].split("family")[1].split("[")[0].strip(":")
        if fam_id not in fam_list:
            fam_list.append(fam_id)

    fam_output = open(fam_out, "w")
    fam_output.write("\n".join(fam_list))

    fam_count = str(len(fam_list))
        
    return fam_count



# ACQUIRE UNIQUE FAMILIES METHOD
def get_uniq_fam(ao_file, al_file, ar_file):
    import pandas as pd
    
    ao_in = open(ao_file, "r")
    al_in = open(al_file, "r")
    ar_in = open(ar_file, "r")
    
    ao_list = ao_in.read().split("\n")
    al_list = al_in.read().split("\n")
    ar_list = ar_in.read().split("\n")
    
    fam_list = ao_list + al_list + ar_list
    dff = pd.DataFrame(fam_list, columns = ["families"])
    dff = dff.sort_values("families")
    uniq_fam_count = str(len(dff["families"].unique())) 
    
    uniq_fams_out = open("../Results/cnv_uniq_fams.txt", "w")
    uniq_fams_list = dff["families"].unique()
    uniq_fams_out.write("\n".join(uniq_fams_list))
    
    
    return uniq_fam_count
    
    
    # families from CNV and SV
    
    
    
    # RI* unique genes and families
    
    # genes from CNV and SV
    
    # families from CNV and SV
    
    
    

# FULL SAMPLE UNIQUES METHOD (genes + families)
def full_sample_uniq(cnv_genes_file, cnv_fams_file, sv_file):
    import pandas as pd
    
    # import sv data
    sv_in = open(sv_file, "r")
    sv_rows = sv_in.read().split("\n")
    
    dfsv = pd.read_csv(sv_file, sep = "\t")
    
    # ASD only unique genes and families
    # genes from CNV and SV
    cnv_ao_genes_in = open("../Results/AO_genes.txt", "r")
    cnv_ao_genes = cnv_ao_genes_in.read().split("\n")
    sv_ao_genes = dfsv["Unique_Genes"][1].strip("[").strip("]").split("', '")
    sv_ao_genes = [x.strip("'") for x in sv_ao_genes]
    
    # families from CNV and SV
    cnv_ao_fams_in = open("../Results/ao_fams.txt", "r")
    cnv_ao_fams = cnv_ao_fams_in.read().split("\n")
    sv_ao_fams = dfsv["Unique_Fams"][1].strip("[").strip("]").split(", ")
    
    # LI* unique genes and families
    # genes from CNV and SV
    cnv_al_genes_in = open("../Results/AL_genes.txt", "r")
    cnv_al_genes = cnv_al_genes_in.read().split("\n")
    sv_al_genes = dfsv["Unique_Genes"][2].strip("[").strip("]").split("', '")
    sv_al_genes = [x.strip("'") for x in sv_al_genes]
    
    # families from CNV and SV
    cnv_al_fams_in = open("../Results/al_fams.txt", "r")
    cnv_al_fams = cnv_al_fams_in.read().split("\n")
    sv_al_fams = dfsv["Unique_Fams"][2].strip("[").strip("]").split(", ")
    
    # RI* unique genes and families
    # genes from CNV and SV
    cnv_ar_genes_in = open("../Results/AR_genes.txt", "r")
    cnv_ar_genes = cnv_ar_genes_in.read().split("\n")
    sv_ar_genes = dfsv["Unique_Genes"][3].strip("[").strip("]").split("', '")
    sv_ar_genes = [x.strip("'") for x in sv_ar_genes]
    
    # families from CNV and SV
    cnv_ar_fams_in = open("../Results/ar_fams.txt", "r")
    cnv_ar_fams = cnv_ar_fams_in.read().split("\n")
    sv_ar_fams = dfsv["Unique_Fams"][3].strip("[").strip("]").split(", ")
    
    # unique genes for the whole sample
    # ASD ONLY
    all_ao_genes = cnv_ao_genes + sv_ao_genes
    dfao = pd.DataFrame(all_ao_genes, columns = ["ao_genes"])
    all_ao_genes_uniq = dfao["ao_genes"].unique().tolist()
    ao_genes_count = len(all_ao_genes_uniq)
    # LI*
    all_al_genes = cnv_al_genes + sv_al_genes
    dfal = pd.DataFrame(all_al_genes, columns = ["al_genes"])
    all_al_genes_uniq = dfal["al_genes"].unique().tolist()
    al_genes_count = len(all_al_genes_uniq)
    # RI*
    all_ar_genes = cnv_ar_genes + sv_ar_genes
    dfar = pd.DataFrame(all_ar_genes, columns = ["ar_genes"])
    all_ar_genes_uniq = dfar["ar_genes"].unique().tolist()
    ar_genes_count = len(all_ar_genes_uniq)
    
    # unique families for the whole sample
    all_ao_fams = cnv_ao_fams + sv_ao_fams
    dfao = pd.DataFrame(all_ao_fams, columns = ["ao_fams"])
    all_ao_fams_uniq = dfao["ao_fams"].unique().tolist()
    ao_fams_count = len(all_ao_fams_uniq)
    # LI*
    all_al_fams = cnv_al_fams + sv_al_fams
    dfal = pd.DataFrame(all_al_fams, columns = ["al_fams"])
    all_al_fams_uniq = dfal["al_fams"].unique().tolist()
    al_fams_count = len(all_al_fams_uniq)
    # RI*
    all_ar_fams = cnv_ar_fams + sv_ar_fams
    dfar = pd.DataFrame(all_ar_fams, columns = ["ar_fams"])
    all_ar_fams_uniq = dfar["ar_fams"].unique().tolist()
    ar_fams_count = len(all_ar_fams_uniq)
    
    # calculate totals
    # total genes
    total_genes_list = all_ao_genes_uniq + all_al_genes_uniq + all_ar_genes_uniq
    dftg = pd.DataFrame(total_genes_list, columns = ["total genes"])
    total_genes_uniq = dftg["total genes"].sort_values().unique().tolist()
    total_genes_count = len(total_genes_uniq)
    
    # total fams
    total_fams_list = all_ao_fams_uniq + all_al_fams_uniq + all_ar_fams_uniq
    dftf = pd.DataFrame(total_fams_list, columns = ["total fams"])
    total_fams_uniq = dftf["total fams"].unique().tolist()
    total_fams_count = len(total_fams_uniq)
    total_fams_out = open("../Results/total_uniq_fams_strict.txt", "w")
    total_fams_out.write("\n".join(total_fams_uniq))
    
    # output matrix
    output_matrix = [
                    ["uniq_genes", "uniq_fams"],
                    [str(ao_genes_count), str(ao_fams_count)],
                    [str(al_genes_count), str(al_fams_count)],
                    [str(ar_genes_count), str(ar_fams_count)],
                    [str(total_genes_count), str(total_fams_count)]
                    ]
    out_file = open("../Results/uniq_matrix.txt", "w")
    for g in range(0, len(output_matrix)):
        line_string = "\t".join(output_matrix[g])
        out_file.write(line_string + "\n")

    # output total unique genes list
    total_genes_out = open("../Results/total_uniq_genes_strict.txt", "w")
    total_genes_out.write("\n".join(total_genes_uniq))
    
    
    print("ao_genes_count: " + str(ao_genes_count))
    print("al_genes_count: " + str(al_genes_count))
    print("ar_genes_count: " + str(ar_genes_count))
    print("ao_fams_count: " + str(ao_fams_count))
    print("al_fams_count: " + str(al_fams_count))
    print("ar_fams_count: " + str(ar_fams_count))
    print("total_fams_uniq: " + str(len(total_fams_uniq)))
    print("total_genes_uniq: " + str(len(total_genes_uniq)))
    
    # close all files
    sv_in.close()
    cnv_ao_genes_in.close()
    cnv_ao_fams_in.close()
    cnv_al_genes_in.close()
    cnv_al_fams_in.close()
    cnv_ar_genes_in.close()
    cnv_ar_fams_in.close()
    out_file.close()
    total_genes_out.close()
    
# MERGE CPDB OUTPUT DATA
def CPDB_merge(path_file, GO_file, merged_path_out, merged_GO_out):
    import pandas as pd
    
    # read in pathway file as dataframe
    dfp = pd.read_csv(path_file, sep = "\t")
    
    # merge on basis of matching column
    agg_functions = {"pathway":" | ".join, "q-value":"first", "effective_size":"first"}
    dfp = dfp.groupby(dfp["members_input_overlap"]).aggregate(agg_functions)
    dfp.sort_values(by = "q-value", inplace = True)

    # output pathway dataframe as tab delimited table
    dfp.to_csv(merged_path_out, sep = "\t")
    
    # read in GO term file as dataframe
    dfg = pd.read_csv(GO_file, sep = "\t")

    # merge on basis of matching column
    agg_functions = {"term_name":" | ".join, "q-value":"first", "effective_size":"first"}
    dfg = dfg.groupby(dfg["members_input_overlap"]).aggregate(agg_functions)
    dfg.sort_values(by = "q-value", inplace = True)

    # output pathway dataframe as tab delimited table
    dfg.to_csv(merged_GO_out, sep = "\t")

# filter out genes TPM < 5
def gene_expression_filter():
    import pandas as pd 
    import pickle
    import math
    
    # read in gene lists
    aog = pd.read_csv("../Results/ASD_only_genes_TPM_UNFILTERED.txt", sep = "\t").values.tolist()
    alg = pd.read_csv("../Results/ASD_LI_genes_TPM_UNFILTERED.txt", sep = "\t").values.tolist()
    arg = pd.read_csv("../Results/ASD_RI_genes_TPM_UNFILTERED.txt", sep = "\t").values.tolist()
    
    # process tpm1
    tpm1 = pd.read_csv("../Data/tpm1.txt", sep = "\t")
    columns_GTEx = list(tpm1.columns[1:-1]) # select GTEx columns from dataframe
    columns_GTEx_brain = [e for e in columns_GTEx if "Brain" in e] # retain only the ones from brain 
    tpm1['tpm1_max'] = tpm1[columns_GTEx_brain].max(axis=1) # declare new column that collects maximum tpm value from all brain columns
    # reduce df_tpm1 to only the relevant columns
    tpm1.rename(columns={'Description': 'GENE'}, inplace=True)
    tpm1 = tpm1[['GENE', 'tpm1_max']]
    
    # process tpm2    
    tpm2 = pd.read_csv("../Data/tpm2.txt", sep = "\t")
    # treat all columns as necessary for now
    columns_brainspan = list(tpm2.columns[2:]) 
    tpm2['tpm2_max'] = tpm2[columns_brainspan].max(axis=1)
    
    # process tpm3    
    tpm3 = pickle.load(open("../Data/tpm3.txt",'rb')) 
    tpm3['tpm3_max'] = tpm3[tpm3.columns[2:]].max(axis=1) 
    tpm3 = tpm3[['Gene_Name', 'tpm3_max'] + [e for e in tpm3.columns if "cerebral_cortex" in e]]
    # The second column is what we care about: TPM3
    tpm3.rename(columns={'Gene_Name': 'GENE'}, inplace = True)
    tpm3 = tpm3[['GENE', 'tpm3_max']]
    
    # merge tpm files into a single dataframe referencing every gene
    tpm_1and2 = tpm1.merge(tpm2, how = 'left', left_on = "GENE", right_on ='geneSymbol')
    tpm_all = tpm_1and2.merge(tpm3, how = 'left', left_on = "GENE", right_on ='GENE') 
    tpm_all = tpm_all[["GENE", "tpm1_max", "tpm2_max", "tpm3_max"]]
    tpm_all.to_csv("../Data/tpm_all.txt", sep = "\t")
    
    # read in tmp_all and compare 
    tpm_in = open("../Data/tpm_all.txt", "r")
    tpm_list = tpm_in.read().split("\n")
    
    # compile all genes
    genes_all = aog + alg + arg
    allg = pd.DataFrame(genes_all, columns = ["genes"])
    allg_uniq = allg["genes"].unique().tolist()
    print("Unique genes in CNV, TPM unfiltered: " + str(len(allg_uniq)))
    all_passing = filter_tpm(allg_uniq, tpm_all, "all")
    
    
    # phenotype specific results
    aog = [''.join(ele) for ele in aog]
    aog_passing = filter_tpm(aog, tpm_all, "ASD_only")
    alg = [''.join(ele) for ele in alg]
    alg_passing = filter_tpm(alg, tpm_all, "LI*")
    arg = [''.join(ele) for ele in arg]
    arg_passing = filter_tpm(arg, tpm_all, "RI*")
    
# =============================================================================
#     passing_genes = []
#     # filter by tpm > 5
#     for a in range(0, len(allg_uniq)):
#         gene = allg_uniq[a]
#         gene_tpm = tpm_all.loc[tpm_all['GENE'] == gene].values.tolist()
#         if len(gene_tpm) < 1:
#             continue
#         else:
#             gene_tpm = gene_tpm[0]
#             for b in range(1, len(gene_tpm)):
#                 if math.isnan(gene_tpm[b]):
#                     continue
#                 else:
#                     if int(gene_tpm[b]) > 5:
#                         passing_genes.append(gene)
#                         break
#                     else:
#                         continue        
#     print("Unique genes in CNV expressed in brain: " + str(len(passing_genes)))
# =============================================================================
    
    # export
    all_genes_out = open("../Results/CNV_candidate_genes.txt", "w") # all three phenotypes, unique, brain expression
    all_genes_out.write("\n".join(all_passing))
    
    ao_out = open("../Results/AO_genes.txt", "w")
    al_out = open("../Results/AL_genes.txt", "w")
    ar_out = open("../Results/AR_genes.txt", "w")

    ao_out.write("\n".join(aog_passing))
    al_out.write("\n".join(alg_passing))
    ar_out.write("\n".join(arg_passing))
    
    ao_out.close()
    al_out.close()
    ar_out.close()
    
    
    # close files
    tpm_in.close()
    
def filter_tpm(gene_input_list, tpm_all, phenotype):
    import pandas as pd
    import math
    
    passing_genes = []
    # filter by tpm > 5
    for a in range(0, len(gene_input_list)):
        gene = gene_input_list[a]
        gene_tpm = tpm_all.loc[tpm_all['GENE'] == gene].values.tolist()
        if len(gene_tpm) < 1:
            continue
        else:
            gene_tpm = gene_tpm[0]
            for b in range(1, len(gene_tpm)):
                if math.isnan(gene_tpm[b]):
                    continue
                else:
                    if int(gene_tpm[b]) > 5:
                        passing_genes.append(gene)
                        break
                    else:
                        continue        
    print("Unique genes in CNV expressed in brain for " + phenotype + ": " + str(len(passing_genes)) + "")
    return(passing_genes)

def variant_table():
    print("Building variant table for all three phenotypes' results...")
    import pandas as pd
    
    # import variant results
    dfao = pd.read_csv("../Results/ASD_only_results.tsv", sep = "\t")
    dfal = pd.read_csv("../Results/ASD_LI_results.tsv", sep = "\t")
    dfar = pd.read_csv("../Results/ASD_RI_results.tsv", sep = "\t")
    
    # add phenotype tag
    dfao["phenotype"] = "ASD"
    dfal["phenotype"] = "LI*"
    dfar["phenotype"] = "RI*"
    
    # combined dataframe
    dfm = pd.DataFrame()
    dfao["family"] = dfao["individuals"]
    dfal["family"] = dfal["individuals"]
    dfar["family"] = dfar["individuals"]
    
    # extract columns to retain
    dfao = dfao[["chrom", "start", "stop", "sv_type", "family", "overlap_CDS_percent", "transcript_id", "genes", "patho_gain_list[phen|hpo|source|coord]", "patho_loss_list[phen|hpo|source|coord]", "cohort_freq", "StrVCTVRE", "SvAnna", "phenotype"]]
    dfal = dfal[["chrom", "start", "stop", "sv_type", "family", "overlap_CDS_percent", "transcript_id", "genes", "patho_gain_list[phen|hpo|source|coord]", "patho_loss_list[phen|hpo|source|coord]", "cohort_freq", "StrVCTVRE", "SvAnna", "phenotype"]]
    dfar = dfar[["chrom", "start", "stop", "sv_type", "family", "overlap_CDS_percent", "transcript_id", "genes", "patho_gain_list[phen|hpo|source|coord]", "patho_loss_list[phen|hpo|source|coord]", "cohort_freq", "StrVCTVRE", "SvAnna", "phenotype"]]
    
    # concatenate dataframes and sort
    dfm = pd.concat([dfao, dfal, dfar], axis=0)
    dfm.sort_values(by=["chrom", "start", "stop"], inplace=True)
    
    # output dataframe
    dfm.to_csv("../Results/CNV_variant_table_unmerged.tsv", sep = "\t", index = None)
    
    
    # MERGE MATCHING VARIANTS 
    
    # open tables
    table_input = open("../Results/CNV_variant_table_unmerged.tsv", "r")
    input_list = table_input.read().split("\n")
    table_output = open("../Results/CNV_variant_table.tsv", "w")
    
    # acquire indices
    chrom_index = dfm.columns.get_loc("chrom")
    SvAnna_index = dfm.columns.get_loc("SvAnna")
    pheno_index = dfm.columns.get_loc("phenotype")
    fam_index = dfm.columns.get_loc("family")
    
    # declare holding list
    holding_list = []
    output_list = []
    output_list.append(input_list[0])
    
    # standalone check
    standalones = 0
    
    # iterate and merge
    for a in range(1, len(input_list)-3):
        first_line = input_list[a].split("\t")
        second_line = input_list[a + 1].split("\t")
        third_line = input_list[a + 2].split("\t")
        
        # determine first positions
        first_pos = "_".join(first_line[chrom_index : SvAnna_index + 1])
        first_pheno = first_line[pheno_index]
        
        # determine second positions
        second_pos = "_".join(second_line[chrom_index : SvAnna_index + 1])
        second_pheno = second_line[pheno_index]
        
        # determine third positions
        third_pos = "_".join(third_line[chrom_index : SvAnna_index + 1])
        third_pheno = third_line[pheno_index]
        
        # if match
        if first_pos == second_pos and first_pos not in holding_list:
            if second_pos == third_pos:
                merge_pheno = first_pheno + ", " + second_pheno + ", " + third_pheno
                merge_line = first_line
                merge_line[pheno_index] = merge_pheno
                output_list.append("\t".join(merge_line))
                holding_list.append(first_pos)    
                continue
            else:
                merge_pheno = first_pheno + ", " + second_pheno
                merge_line = first_line
                merge_line[pheno_index] = merge_pheno
                output_list.append("\t".join(merge_line))
                holding_list.append(first_pos)    
                continue
        elif first_pos not in holding_list:
            output_list.append("\t".join(first_line))
            holding_list.append(first_pos)
            standalones +=1
    output_list2 = []
    output_list2.append(input_list[0])
    
    print("# of standalones: " + str(standalones))
    
    # correct fam IDs
    for b in range(1, len(output_list)):
        line_list = output_list[b].split("\t")
        fam_field = line_list[fam_index]
        if "|" in fam_field:
            fam1 = fam_field.split("|")[0][:4]
            fam2 = fam_field.split("|")[1][:4]
            if fam1 == fam2:
                fam = convert_fam(fam1)
                
            else:
                fam = convert_fam(fam1) + ", " + convert_fam(fam2)
        else:
            fam = convert_fam(fam_field[:4])
        line_list[fam_index] = fam
        output_list2.append("\t".join(line_list))
            
    # output merged table
    table_output.write("\n".join(output_list2))
    
    # close open files
    table_input.close()
    table_output.close()
    
    print("Variant table built.")
        

def convert_fam(input_fam):
    import pandas as pd
    
    dfy = pd.read_csv("../../CNV_Pipeline/Data/Starting_Manifest/alternative_table.tsv", sep = "\t")
    dfy = dfy[["Family ID", "New Family ID"]]
    dfy["Family ID"] = dfy["Family ID"].astype(str)
    dfy = dfy.drop_duplicates()
    fam_dict = dfy.set_index('Family ID')['New Family ID'].to_dict()
    new_ID = fam_dict[input_fam]
    
    return(new_ID)

    
def count_CNV2_indv():
    import pandas as pd
    
    df = pd.read_csv("../Data/CNV2.vcf", skiprows = 49, sep = "\t")
    
    vcf_in = open("../Data/CNV2.vcf", "r")
    vcf_list = vcf_in.read().split("\n")
    vcf_out = open("../Data/CNV3.tsv", "w")
    
    invariable_out = open("../Data/CNV_invariable.tsv", "w")
    variable_count = 0
    
    for a in range(49, len(vcf_list)-1):
        var_list = vcf_list[a].split("\t")
        invariable = True
        for b in range(9, len(var_list)):
            if var_list[b] != "0/0:0:0:0:0:0:null:0":
                invariable = False
                break
        if invariable == True:
            invariable_out.write(vcf_list[a] + "\n")
        else:
            vcf_out.write(vcf_list[a] + "\n")
            variable_count += 1
    
    
            
    print("# of variants with invariable sites removed: " + str(variable_count))
    vcf_in.close()
    vcf_out.close()
    invariable_out.close()

# find median of QC filtered sample
def median_finder():
    import pandas as pd
    import numpy as np
    
    dfv = pd.read_csv("../Data/CNV3.vcf", skiprows = 49, sep = "\t")
    indv_list = dfv.columns
    print(len(indv_list[9:]))
    print(len(dfv))
    
    
    dfv["start"] = dfv["POS"]
    dfv[["chr", "start", "stop"]] = dfv["ID"].str.split("_", expand = True)
    dfv["length"] = dfv["stop"].astype(int) - dfv["start"].astype(int)
    
    print(dfv["ALT"].value_counts())
    
    CNV_median = dfv["length"].median()
    CNV_min = dfv["length"].min()
    CNV_max = dfv["length"].max()
    
    indv_list = dfv.columns
    print(indv_list[9:])
    
    print("Median of QC filtered CNVs: " + str(CNV_median))
    print("Minimum of QC filtered CNVs: " + str(CNV_min))
    print("Maximum of QC filtered CNVs: " + str(CNV_max))
 
def all_fams():
    import pandas as pd
    
    dfm = pd.read_csv("../Data/Starting_Manifest/SampleSummary_NewID_2023_05_24.tsv", sep = "\t")
    dfm = dfm[(dfm["caller"] != "No variants called") | (dfm["VCF"] != 0)]
    fams_count = len(dfm["Family ID"].unique())
    print(fams_count)
    
all_fams()
 
    
# COMMANDS
build_varsum_table()
 



phenotype = "ASD_only"
gene_extract("../Results/" + phenotype + "_results.tsv", "../Data/CNV2_anno_corrected.tsv", "../Results/intermediary_files/result_genes_1_" + phenotype +".txt", "../Results/intermediary_files/result_genes_2_" + phenotype +".txt", "../Results/intermediary_files/result_genes_3_" + phenotype +".txt")
     
phenotype = "ASD_LI"
gene_extract("../Results/" + phenotype + "_results.tsv", "../Data/CNV2_anno_corrected.tsv", "../Results/intermediary_files/result_genes_1_" + phenotype +".txt", "../Results/intermediary_files/result_genes_2_" + phenotype +".txt", "../Results/intermediary_files/result_genes_3_" + phenotype +".txt")# 

phenotype = "ASD_RI"
gene_extract("../Results/" + phenotype + "_results.tsv", "../Data/CNV2_anno_corrected.tsv", "../Results/intermediary_files/result_genes_1_" + phenotype +".txt", "../Results/intermediary_files/result_genes_2_" + phenotype +".txt", "../Results/intermediary_files/result_genes_3_" + phenotype +".txt")

sort_results()

gene_expression_filter()
 
variant_table()


print("ASD_only family count: " + count_families("../Results/ASD_only_results.tsv", "../Results/AO_genes.txt", "../Results/ao_fams.txt"))
print("LI* family count: " + count_families("../Results/ASD_LI_results.tsv", "../Results/AL_genes.txt", "../Results/al_fams.txt"))
print("RI* family count: " + count_families("../Results/ASD_RI_results.tsv","../Results/AR_genes.txt", "../Results/ar_fams.txt"))
print("Unique families: " + get_uniq_fam("../Results/ao_fams.txt", "../Results/al_fams.txt", "../Results/ar_fams.txt")) 
print("\n\nFULL SAMPLE:")
full_sample_uniq("../Results/cnv_uniq_genes.txt", "../Results/cnv_uniq_fams.txt", "../Results/counts_by_phenotype.tsv")

CPDB_merge("../Results/pathways_20230609.tab", "../Results/GOterms_20230609.tab", "../Results/merged_pathways.tsv", "../Results/merged_GOterms.tsv")



#median_finder()

    

