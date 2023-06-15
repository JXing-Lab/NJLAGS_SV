# table build method
def enrichment_table(outliers_input, manifest_input, ASD_proband_output, LI_proband_output, RI_proband_output, unaffected_output):
    print("Building enrichment table...")
    
    # import libraries
    import pandas as pd
    import numpy as np
    
    # import file and remove outliers
    dfo = pd.read_csv(outliers_input, sep = "\t")
    outliers_list = dfo["PARTID"].tolist()
    
    print(len(outliers_list))
    
    dfm = pd.read_csv(manifest_input, sep = "\t")
    
    dfp_LI_only = dfm[(dfm["LI"] == 2) & (dfm["ASD"] == 1) & (dfm["RI"] == 1)]

    print(len(dfm))
    dfm = dfm[~dfm["PARTID"].isin(outliers_list)] # remove outliers  
    dfm = dfm[dfm["caller"] != "No variants called"] # remove variants with no CNV calling data
    print(len(dfm))
    # define unaffected
    dfu = dfm[(dfm["ASD"] == "1") & ((dfm["LI"] == "1") | (dfm["LI"] == "x")) & ((dfm["RI"] == "1") | (dfm["RI"] == "x"))]
    
    # define affected
    dfp_ASD = dfm[dfm["ASD"] == "2"]
    dfp_LI = dfm[dfm["LI"] == "2"]
    dfp_RI = dfm[dfm["RI"] == "2"]


    # output 
    dfu.to_csv(unaffected_output, sep = "\t", index=False)
    dfp_ASD.to_csv(ASD_proband_output, sep = "\t", index=False)
    dfp_LI.to_csv(LI_proband_output, sep = "\t", index=False)
    dfp_RI.to_csv(RI_proband_output, sep = "\t", index=False)
    
    print("Enrichment table built.")




# length based enrichment
def length_enrichment(vcf_file, pedigree_file, table_out,  ASD_proband_output, LI_proband_output, RI_proband_output, unaffected_output):
    
    # import packages
    import pandas as pd
    import numpy as np
    
    # import vcf as list
    vcf_in = open(vcf_file, "r")
    vcf_list = vcf_in.read().split("\n")
    
    # import vcf as dataframe
    dfv = pd.read_csv(vcf_file, skiprows = 49, sep = "\t") # skip the vcf header
    
    # build new dataframe to hold CN counts
    partid_list = dfv.columns.tolist()[9:]
    dfi = pd.DataFrame(partid_list, columns = ["PARTID"])
    dfi["length_impact_1"] = 0
    dfi["length_impact_2"] = 0
    dfi["CN0"] = 0
    dfi["CN1"] = 0
    dfi["CN3"] = 0
    dfi["CN4"] = 0
    dfi.set_index("PARTID")
    
    # VCF header list
    header_list = vcf_list[49].split("\t")

    # iterate through vcf list, incrementing corresponding individual as variant is found
    for b in range(50, len(vcf_list) - 1):
        var_list = vcf_list[b].split("\t")
        CN = var_list[4].strip(">").strip("<") # assign CN number for comparison    
        
        # acquire length of CNV
        CNV_length = int(var_list[2].split("_")[2]) - int(var_list[2].split("_")[1])
        
        # iterate through columns until a nonzero genotype is found
        
        for c in range(9, len(var_list)):
            geno_list = var_list[c].split(":")
            genotype = geno_list[0]
            if genotype != "0/0": 
                #print(str(b) + "| " + header_list[c] + ": " + str(CN))
                indv_index = int(dfi[dfi["PARTID"] == header_list[c]].index[0])
                
                # increment individual's length impact
                dfi.at[indv_index, "length_impact_1"] += CNV_length
                
                # increment individual's length impact by CN number
                if CN == "CN1" or CN == "CN3":
                    dfi.at[indv_index, "length_impact_2"] += CNV_length
                elif CN == "CN0" or CN == "CN4":
                    dfi.at[indv_index, "length_impact_2"] += CNV_length * 2
                else:
                    print("ERROR: " + str(CN))
                
                dfi.at[indv_index, CN] += 1

    # calculate CN sum
    CN_columns = ["CN0","CN1", "CN3", "CN4"]
    dfi["CN_sum"] = dfi[CN_columns].sum(axis = 1)
    
    # merge tables
    dfp = pd.read_csv(pedigree_file, sep = "\t")
    dfi["PARTID"] = dfi["PARTID"].astype(int)
    dfp = dfp.dropna(subset = ["PARTID"])
    dfp["PARTID"] = dfp["PARTID"].astype(int)
    dfm = pd.merge(dfp, dfi, on = "PARTID")
    dfm = dfm.astype(str)
  
    # define unaffected
    dfu = dfm[(dfm["ASD"] == "1") & ((dfm["LI"] == "1") | (dfm["LI"] == "x")) & ((dfm["RI"] == "1") | (dfm["RI"] == "x"))]
    # define affected
    dfp_ASD = dfm[dfm["ASD"] == "2"]
    dfp_LI = dfm[dfm["LI"] == "2"]
    dfp_RI = dfm[dfm["RI"] == "2"]

    print("length_impact_1")
    unaff_median = dfu["length_impact_1"].astype(int).median()
    print("unaff median: " + str(unaff_median))
    ASD_median = dfp_ASD["length_impact_1"].astype(int).median()
    print("ASD median: " + str(ASD_median))
    LI_median = dfp_LI["length_impact_1"].astype(int).median()
    print("LI median: " + str(LI_median))
    RI_median = dfp_RI["length_impact_1"].astype(int).median()
    print("RI median: " + str(RI_median))
    
    unaff_mean = dfu["length_impact_1"].astype(int).mean()
    print("unaff mean: " + str(unaff_mean))
    ASD_mean = dfp_ASD["length_impact_1"].astype(int).mean()
    print("ASD mean: " + str(ASD_mean))
    LI_mean = dfp_LI["length_impact_1"].astype(int).mean()
    print("LI mean: " + str(LI_mean))
    RI_mean = dfp_RI["length_impact_1"].astype(int).mean()
    print("RI mean: " + str(RI_mean))


    print("\n\nlength_impact_2")
    unaff_median = dfu["length_impact_2"].astype(int).median()
    print("unaff median: " + str(unaff_median))
    ASD_median = dfp_ASD["length_impact_2"].astype(int).median()
    print("ASD median: " + str(ASD_median))
    LI_median = dfp_LI["length_impact_2"].astype(int).median()
    print("LI median: " + str(LI_median))
    RI_median = dfp_RI["length_impact_2"].astype(int).median()
    print("RI median: " + str(RI_median))
    
    unaff_mean = dfu["length_impact_2"].astype(int).mean()
    print("unaff mean: " + str(unaff_mean))
    ASD_mean = dfp_ASD["length_impact_2"].astype(int).mean()
    print("ASD mean: " + str(ASD_mean))
    LI_mean = dfp_LI["length_impact_2"].astype(int).mean()
    print("LI mean: " + str(LI_mean))
    RI_mean = dfp_RI["length_impact_2"].astype(int).mean()
    print("RI mean: " + str(RI_mean))


    # output 
    dfu.to_csv(unaffected_output, sep = "\t", index=False)
    dfp_ASD.to_csv(ASD_proband_output, sep = "\t", index=False)
    dfp_LI.to_csv(LI_proband_output, sep = "\t", index=False)
    dfp_RI.to_csv(RI_proband_output, sep = "\t", index=False)

    
    # export to table
    dfm.to_csv(table_out, sep = "\t", index = False)
    
    




# CDS overlap based enrichment
def CDS_enrichment(anno_file, pedigree_file, table_out,  ASD_proband_output, LI_proband_output, RI_proband_output, unaffected_output):
    
    # import packages
    import pandas as pd
    import numpy as np
    
    # import vcf as list
    anno_in = open(anno_file, "r")
    anno_list = anno_in.read().split("\n")
    
    # import vcf as dataframe
    dfa = pd.read_csv(anno_file, sep = "\t") # skip the vcf header
    
    # set indices
    anno_mode_index = dfa.columns.get_loc("Annotation_mode")
    CDS_index = dfa.columns.get_loc("Overlapped_CDS_length")
    CN_index = dfa.columns.get_loc("SV_type")
    indv_start_index = dfa.columns.get_loc("1002001")
    
    # build new dataframe to hold CN counts
    partid_list = dfa.columns.tolist()[indv_start_index:anno_mode_index]
    dfi = pd.DataFrame(partid_list, columns = ["PARTID"])
    dfi["CDS_impact"] = 0
    dfi["CN0"] = 0
    dfi["CN1"] = 0
    dfi["CN3"] = 0
    dfi["CN4"] = 0
    dfi.set_index("PARTID")
    
    # VCF header list
    header_list = dfa.columns.tolist()
    
    # iterate through vcf list, incrementing corresponding individual as variant is found
    for b in range(1, len(anno_list) - 1):
    
        var_list = anno_list[b].split("\t")
        CN = var_list[CN_index].strip(">").strip("<") # assign CN number for comparison    
        
        # only run on split columns
        if var_list[anno_mode_index] == "full":
            continue
        elif var_list[anno_mode_index] == "split":        
            # acquire length of CNV
            CDS_length = int(var_list[CDS_index])
            
            # iterate through columns until a nonzero genotype is found
            for c in range(indv_start_index, anno_mode_index):     
                geno_list = var_list[c].split(":")
                genotype = geno_list[0]
        
                if genotype != "0/0":                
                    # match to individual
                    indv_index = int(dfi[dfi["PARTID"] == header_list[c]].index[0])
                    
                    # increment individual's length impact
                    dfi.at[indv_index, "CDS_impact"] += CDS_length
                    dfi.at[indv_index, CN] += 1

    # calculate CN sum
    CN_columns = ["CN0","CN1", "CN3", "CN4"]
    dfi["CN_sum"] = dfi[CN_columns].sum(axis = 1)
    
    # merge tables
    dfp = pd.read_csv(pedigree_file, sep = "\t")
    dfi["PARTID"] = dfi["PARTID"].astype(int)
    dfp = dfp.dropna(subset = ["PARTID"])
    dfp["PARTID"] = dfp["PARTID"].astype(int)
    dfm = pd.merge(dfp, dfi, on = "PARTID")
    print(dfm)
    dfm = dfm.astype(str)
  
    # define unaffected
    dfu = dfm[(dfm["ASD"] == "1") & ((dfm["LI"] == "1") | (dfm["LI"] == "x")) & ((dfm["RI"] == "1") | (dfm["RI"] == "x"))]
    # define affected
    dfp_ASD = dfm[dfm["ASD"] == "2"]
    dfp_LI = dfm[dfm["LI"] == "2"]
    dfp_RI = dfm[dfm["RI"] == "2"]

    print("CDS_impact")
    unaff_median = dfu["CDS_impact"].astype(int).median()
    print("unaff median: " + str(unaff_median))
    ASD_median = dfp_ASD["CDS_impact"].astype(int).median()
    print("ASD median: " + str(ASD_median))
    LI_median = dfp_LI["CDS_impact"].astype(int).median()
    print("LI median: " + str(LI_median))
    RI_median = dfp_RI["CDS_impact"].astype(int).median()
    print("RI median: " + str(RI_median))
    
    unaff_mean = dfu["CDS_impact"].astype(int).mean()
    print("unaff mean: " + str(unaff_mean))
    ASD_mean = dfp_ASD["CDS_impact"].astype(int).mean()
    print("ASD mean: " + str(ASD_mean))
    LI_mean = dfp_LI["CDS_impact"].astype(int).mean()
    print("LI mean: " + str(LI_mean))
    RI_mean = dfp_RI["CDS_impact"].astype(int).mean()
    print("RI mean: " + str(RI_mean))

    # output 
    dfu.to_csv(unaffected_output, sep = "\t", index=False)
    dfp_ASD.to_csv(ASD_proband_output, sep = "\t", index=False)
    dfp_LI.to_csv(LI_proband_output, sep = "\t", index=False)
    dfp_RI.to_csv(RI_proband_output, sep = "\t", index=False)

    
    # export to table
    dfm.to_csv(table_out, sep = "\t", index = False)


# COMMANDS    
enrichment_table("../../../Data/outlier_IDs.txt", "SampleSummary_NewID_2023_05_25.tsv",  "ASD_proband_counts.tsv", "LI_proband_counts.tsv", "RI_proband_counts.tsv", "unaffected_counts.tsv")    
#length_enrichment("../../../Data/CNV3.vcf", "../../../Data/NJLAGS_CNV.ped", "./length_enrichment_table.tsv", "ASD_proband_lengths.tsv", "LI_proband_lengths.tsv", "RI_proband_lengths.tsv", "unaffected_lengths.tsv")    
#CDS_enrichment("../../../Data/CNV3_anno_corrected.tsv", "../../../Data/NJLAGS_CNV.ped", "./CDS_enrichment_table.tsv", "ASD_proband_CDS_lengths.tsv", "LI_proband_CDS_lengths.tsv", "RI_proband_CDS_lengths.tsv", "unaffected_CDS_lengths.tsv")