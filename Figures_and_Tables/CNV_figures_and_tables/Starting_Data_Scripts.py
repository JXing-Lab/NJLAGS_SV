"""
Starting_Data_Scripts.py
> Builds a table that summarizes all starting data by individual
> Assigns each individual the corresponding anonymized ID and pedigree information
> Determines which callers and which batch the individual's data was covered by
> Assigns the count of CNVs according to caller
> Assings the count of CNVs after merging and QC filtering (but BEFORE annotation)

Created 2023.03.29
Updated 2023.06.16
"""

# import modules
import pandas as pd
import numpy as np
import os

# MAIN METHOD
def build_table(input_manifest, input_pedigree, output_table):
    print("Running build_table() method...")
    
    # load in manifest file and pedigree
    dfm = pd.read_csv(input_manifest, sep = "\t", skipfooter = 1, engine = "python")
    dfp = pd.read_csv(input_pedigree, sep = "\t")
    # declare output dataframe
    dfo = pd.DataFrame()
    
    # assign columns from manifest dataframe to output dataframe
    """
    dfo["PARTID"] = dfm["PARTID"].astype(int) # to be replaced with anonymized IDs
    dfo["FAMID"] = dfm["New PARTID"]
    dfo["WGSID"] = dfm["WGS ID"]
    dfo["SEX"] = dfm["SEX"]
    # match phenotypes from pedigree file
    dfm["ASD"] = dfm["PARTID"].map(dfp.set_index(["PARTID"])["ASD"])
    dfm["LI"] = dfm["PARTID"].map(dfp.set_index(["PARTID"])["LI"])
    dfm["RI"] = dfm["PARTID"].map(dfp.set_index(["PARTID"])["RI"])
    dfm["SRS"] = dfm["PARTID"].map(dfp.set_index(["PARTID"])["SRS"])
    """
    # remove duplicates
    dfm = dfm.drop_duplicates(subset = "PARTID", keep = "first")

    # assign counts of CNVs
    partid_list = dfm["PARTID"].values.tolist()
    CNV_counter_raw("../starting_call_files/in_active_use", partid_list, "./batch_caller_counts.tsv")
    CNV_counter_filtered("../CNV3.vcf", partid_list, "./filtered_counts.tsv")
    
    # import raw counts table
    dfr = pd.read_csv("./batch_caller_counts.tsv", sep = "\t")
    dfm["Axiom_CNVs"] = dfm["PARTID"].map(dfr.set_index(["PARTID"])["Axiom_CNVs"])
    dfm["Illumina_CNVs"] = dfm["PARTID"].map(dfr.set_index(["PARTID"])["Illumina_CNVs"])
    dfm["quantisnp_calls"] = dfm["PARTID"].map(dfr.set_index(["PARTID"])["quantisnp_calls"])
    dfm["penncnv_calls"] = dfm["PARTID"].map(dfr.set_index(["PARTID"])["penncnv_calls"])
    dfm["caller"] = dfm["PARTID"].map(dfr.set_index(["PARTID"])["caller"])
    
    # import filtered counts table
    dff = pd.read_csv("./filtered_counts.tsv", sep = "\t")
    dfm["CN0"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN0"])
    dfm["CN1"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN1"])
    dfm["CN2"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN2"])
    dfm["CN3"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN3"])
    dfm["CN4"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN4"])
    dfm["CN5"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN5"])
    dfm["CN6"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN6"])
    dfm["CN7"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN7"])
    dfm["CN8"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN8"])
    dfm["CN9"] = dfm["PARTID"].map(dff.set_index(["partid"])["CN9"])
    
    # print out count of CN1s
    print("Most common genotype is CN1 at " + str(dfm["CN1"].sum()))
    
    # create new column that sums up all other CN counts
    CN_columns = ["CN0","CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9"]
    dfm["CN_sum"] = dfm[CN_columns].sum(axis = 1)
    
    # add LI* and RI* columns
    dfm["LI*"] = np.where((dfm["ASD"] == "2") | (dfm["LI"] == "2"), True, False)
    dfm["RI*"] = np.where((dfm["ASD"] == "2") | (dfm["RI"] == "2"), True, False)
    
    # acquire individual lists from results
    ASD_only_ids = find_indv_results("ASD_only")
    ASD_LI_ids = find_indv_results("ASD_LI")
    ASD_RI_ids = find_indv_results("ASD_RI")
    all_results_ids = ASD_LI_ids + ASD_RI_ids
    
    # reference individual lists from results
    dfm.loc[dfm["PARTID"].isin(ASD_only_ids), "results"] = "ASD_only"
    dfm.loc[dfm["PARTID"].isin(ASD_LI_ids), "results"] = "ASD_LI"
    dfm.loc[dfm["PARTID"].isin(ASD_RI_ids), "results"] = "ASD_RI"
    dfm.loc[dfm["PARTID"].isin(all_results_ids), "results"] = "ASD_LI, ASD_RI"
    
    # output dataframe to file
    dfm.to_csv(output_table, sep = "\t", index = False)
    
    print("build_table() method finished running.")
    
# CNV COUNTS PREANNO METHOD
def CNV_counter_raw(data_folder, partid_list, output_file):
    print("Running CNV_counter_raw() method...")
    
    # import all starting data files
    file_list = [f for f in os.listdir(data_folder) if os.path.isfile(os.path.join(data_folder, f))]
    file_out = open(output_file, "w")
    
    # for each individual, go through all starting call files and increment CNV counts by caller
    header_list = ["PARTID", "Axiom_CNVs", "Illumina_CNVs", "quantisnp_calls", "penncnv_calls", "caller"]
    header_string = "\t".join(header_list) + "\n"
    file_out.write(header_string)
    output_list = []
    # iterate through individuals
    for a in range(0, len(partid_list)):
        output_list = [partid_list[a], 0, 0, 0, 0, "caller"] # build holding list
        # batch booleans
        axiom_batch = False
        illm_batch = False
    
        
        # iterate through starting call files
        for b in range(0, len(file_list)):
            call_file = "../starting_call_files/in_active_use/" + file_list[b]
            print("--Now reading in call file: " + call_file + "--")
            call_in = open(call_file, "r")
            call_list = call_in.read().split("\n")
            # iterate through lines in current call file
            for c in range(1, len(call_list)):
                # if file is Axiom, increment Axiom_CNVs count on any match
                if str(partid_list[a]) in call_list[c] and ("axiom" in call_file or "psych_run4" in call_file):
                    output_list[1] += 1
                    axiom_batch = True
                # if file is Illumina, increment Illumina_CNVs on any match
                if str(partid_list[a]) in call_list[c] and ("illm" in call_file or "may" in call_file or "jan" in call_file):
                    output_list[2] += 1
                    illm_batch = True
                # if file was called by quantisnp, increment quantisnp_calls on any match
                if str(partid_list[a]) in call_list[c] and ("axiom_run3" in call_file or "illm_jan2015" in call_file or "illm_may2017" in call_file):
                    output_list[3] += 1
                # if file was called by penncnv, increment penncnv_calls on any match
                if str(partid_list[a]) in call_list[c] and ("axiom_run7" in call_file or "psych_run4" in call_file or "psych_may2017" in call_file):
                    output_list[4] += 1
        
        # assign batch
        if axiom_batch and illm_batch:
            output_list[5] = "Axiom and Illumina"
        elif axiom_batch:
            output_list[5] = "Axiom"
        elif illm_batch:
            output_list[5] = "Illumina"
        else:
            output_list[5] = "No variants called"
        
        for d in range(0, len(output_list)):
            output_list[d] = str(output_list[d])
        output_string = "\t".join(output_list) + "\n"
        file_out.write(output_string)
        
    file_out.close()    
    
    print("CNV_counter_raw() method finished running.")
    
def CNV_counter_filtered(input_vcf, partid_list, output_file):
    print("Running CNV_counter_filtered() method...")
    
    # import vcf
    vcf_in = open(input_vcf, "r")
    vcf_list = vcf_in.read().split("\n")
    
    # declare header_list
    header_list = ["partid", "CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9"]
    
    # open export file
    out_file = open(output_file, "w")
    out_file.write("\t".join(header_list) + "\n") # output header
   
    # sort through partid_list, increment counts upon finding matching variants
    for a in range(0, len(partid_list)):
        summ_indv = str(partid_list[a]) # individual ID of summary file
        output_list = [summ_indv, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        # set up partID list of individuals from vcf
        vcf_header = vcf_list[49].split("\t")
        indv_list = vcf_header[9:len(vcf_header)]
        
        # iterate through vcf
        for b in range(50, len(vcf_list) - 1):
            var_list = vcf_list[b].split("\t")
            CN = int(var_list[4].split("N")[1].strip(">")) # assign CN number for comparison    
            
            # iterate through columns until a nonzero genotype is found
            for c in range(9, len(var_list)):
                geno_list = var_list[c].split(":")
                genotype = geno_list[0]
                if genotype != "0/0": 
                    # if in individual of interest, increment corresponding CN
                    vcf_indv = indv_list[c - 9] # individual ID of vcf
                    if vcf_indv == summ_indv:
                        output_list[CN + 1] += 1 # increment the corresponding output list CN column
                    else:
                        continue
                else:
                    continue
        
        # write out output
        out_file.write("\t".join(str(x) for x in output_list) + "\n")
                    
    print("CNV_counter_filtered() method finished running.")
    
# method to determine if an individual is included in CNV results
def find_indv_results(phenotype):
    print("Collating individuals found in results for " + phenotype + "...")
    
    # import libraries
    import pandas as pd
    
    # declare final_list
    final_list = []
    
    # open each results file for each phenotype    
    dfr = pd.read_csv("../../Results/" + phenotype + "_results.tsv", sep = "\t")
    indv_list = dfr["individuals"].unique().tolist()
    for a in range(0, len(indv_list)):
        if "|" in indv_list[a]:
            id_list = indv_list[a].split("|")
            for b in range(0, len(id_list)):
                if id_list[b] not in final_list:
                    final_list.append(int(id_list[b]))
        elif indv_list[a] not in final_list:
            final_list.append(int(indv_list[a]))
        else:
            continue
    
    return final_list

# method to identify outliers
def find_outliers(input_file):
    print("Running find_outliers() method...")
    
    # import pandas
    import pandas as pd
    import numpy as np
    
    # load in counts file as dataframe of counts (dfc)
    dfc = pd.read_csv(input_file, sep = "\t")
    # create new column that sums up all other CN counts
    CN_columns = ["CN0","CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9"]
    dfc["CN_sum"] = dfc[CN_columns].sum(axis = 1)
    
    # find median
    count_median = dfc["CN_sum"].median()
    print("cohort median :" + str(count_median))
    # calculate interquartile range
    Q3 = np.quantile(dfc["CN_sum"], 0.75)
    Q1 = np.quantile(dfc["CN_sum"], 0.25)
    IQR = Q3 - Q1
    # calculate outlier threshold
    out_thr = count_median + (IQR * 1.5)
    
    print(str(out_thr))
    
    # extract outliers
    dfo = dfc[dfc["CN_sum"] > out_thr]
    dfo.to_csv("./outliers.tsv", sep = "\t")
    
    outlier_IDs = dfo["partid"].tolist()
    CN_index = dfo.columns.get_loc("CN_sum")
    print(type(CN_index))
    
    print("Finished running find_outliers() method.")
    
# FIND SIBLING PAIRS METHOD 
def find_sibs(input_file):
    """
    Parameters
    ----------
    input_file : TABLE
        alternative_table.tsv

    Returns
    -------
    sibling pairs that have sequencing data, microarray data, and one proband and one unaffected
    """
    
    import pandas as pd
    
    # build table for relevant families
    dft = pd.read_csv(input_file, sep = "\t") # import table
    dfd = dft[dft["VCF"] == 1] # extract with sequence data
    dfd = dfd[dfd["caller"] != "No variants called"] # extract indvs with microarray data
    dfd = dfd[dfd["Father ID"] != 0] # remove parents
    dfa = dfd[dfd["ASD"] == "2"] # acquire list of affected
    dfu = dfd[dfd["ASD"] == "1"] # acquire list of unaffected
    
    pair_list = [] # declare holding list
    # iterate through aff list, find matching unaffected siblings
    for a in dfa.index:
        aff_partid = dfa["PARTID"][a]
        aff_father = dfa["Father ID"][a]
        aff_mother = dfa["Mother ID"][a]
        
        # iterate through unaff list
        for b in dfu.index:
            unaff_partid = dfu["PARTID"][b]
            unaff_father = dfu["Father ID"][b]
            unaff_mother = dfu["Mother ID"][b]
            
            # test for match
            if aff_father == unaff_father and aff_mother == unaff_mother:
                sib_pair = [str(aff_partid), str(unaff_partid)]
                pair_list.append(sib_pair)
    
    # write out pairs
    out_file = open("./sib_pairs.tsv", "w")
    header_list = ["aff", "unaff"]
    out_file.write("\t".join(header_list) + "\n")
    for c in range(0, len(pair_list)):
        sib_string = "\t".join(pair_list[c]) + "\n"
        out_file.write(sib_string)
    out_file.close()


# build sample family summary table
def build_sample_fam(input_file, ASD_GMN, LI_GMN, RI_GMN, output_file):
    import pandas as pd
    import math

    # build table for relevant families
    dft = pd.read_csv(input_file, sep = "\t") # import table
    
    # populate CNV section
    df_CNV = dft[dft["caller"] != "No variants called"] # remove individuals with no microarray data
    df_SV = dft[dft["VCF"] == 1]
    candidate_fams = len(df_CNV["Family ID"].unique())
    candidate_CNVs = df_CNV["CN_sum"].sum()
    df_ASD = df_CNV[df_CNV["ASD"] == "2"]
    df_ASD_LI = df_CNV[(df_CNV["ASD"] == "2")|(df_CNV["LI"] == "2")]
    df_ASD_RI = df_CNV[(df_CNV["ASD"] == "2")|(df_CNV["RI"] == "2")]
    
   
    
    # count families
    
    # ASD
    ASD_probands = df_ASD["PARTID"].values.tolist()
    # holding list for dominant segregating families
    ASD_dom_fams = []
    # holding list for denovo/recessive segregating families
    ASD_denrec_fams = []
    # iterate through list of affected individuals
    for a in range(0, len(ASD_probands)):
        # holding booleans
        dominant = False
        
        partid = ASD_probands[a]
        famid = str(df_CNV["Family ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".")
        # check both parents for affected status
        father_id = str(df_CNV["Father ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".") # find father
        mother_id = str(df_CNV["Mother ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".") # find father
        # find parent affected status
        father_status = dft.loc[dft["PARTID"] == float(father_id)]["ASD"].tolist()
        mother_status = dft.loc[dft["PARTID"] == float(mother_id)]["ASD"].tolist()
        #print("Son is " + str(partid) + ", father is " + str(father_id) + ", mother is " + str(mother_id))
      
        # if either parent is affected, assign family ID field to list of dominant (run this check on dft)
        if len(father_status) > 0:
            father_status = father_status[0]
            if father_status == "2":
                dominant = True
        if len(mother_status) > 0:
            if mother_status == "2":
                dominant = True
        
        # assign family to category
        if dominant == True:
            ASD_dom_fams.append(famid)
        else:
            ASD_denrec_fams.append(famid)
    
    df_ASD_dom_fams = pd.DataFrame(ASD_dom_fams, columns = ["families"])
    ASD_dom_fams_count = str(len(df_ASD_dom_fams["families"].unique()))
    df_ASD_denrec_fams = pd.DataFrame(ASD_denrec_fams, columns = ["families"])
    ASD_denrec_fams_count = str(len(df_ASD_denrec_fams["families"].unique()))
            
    # LI*
    ASD_LI_probands = df_ASD_LI["PARTID"].values.tolist()
    # holding list for dominant segregating families
    ASD_LI_dom_fams = []
    # holding list for denovo/recessive segregating families
    ASD_LI_denrec_fams = []
    # iterate through list of affected individuals
    for a in range(0, len(ASD_LI_probands)):
        # holding booleans
        dominant = False
        
        partid = ASD_LI_probands[a]
        famid = str(df_CNV["Family ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".")
        # check both parents for affected status
        father_id = str(df_CNV["Father ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".") # find father
        mother_id = str(df_CNV["Mother ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".") # find father
        # find parent affected status
        father_status = dft.loc[dft["PARTID"] == float(father_id)]["LI"].tolist()
        mother_status = dft.loc[dft["PARTID"] == float(mother_id)]["LI"].tolist()
        #print("Son is " + str(partid) + ", father is " + str(father_id) + ", mother is " + str(mother_id))
      
        # if either parent is affected, assign family ID field to list of dominant (run this check on dft)
        if len(father_status) > 0:
            father_status = father_status[0]
            if father_status == "2":
                dominant = True
        if len(mother_status) > 0:
            if mother_status == "2":
                dominant = True
        
        # assign family to category
        if dominant == True:
            ASD_LI_dom_fams.append(famid)
        else:
            ASD_LI_denrec_fams.append(famid)
    
    df_ASD_LI_dom_fams = pd.DataFrame(ASD_LI_dom_fams, columns = ["families"])
    ASD_LI_dom_fams_count = str(len(df_ASD_LI_dom_fams["families"].unique()))
    df_ASD_LI_denrec_fams = pd.DataFrame(ASD_LI_denrec_fams, columns = ["families"])
    ASD_LI_denrec_fams_count = str(len(df_ASD_LI_denrec_fams["families"].unique()))
    
    # RI*
    ASD_RI_probands = df_ASD_RI["PARTID"].values.tolist()
    # holding list for dominant segregating families
    ASD_RI_dom_fams = []
    # holding list for denovo/recessive segregating families
    ASD_RI_denrec_fams = []
    # iterate through list of affected individuals
    for a in range(0, len(ASD_RI_probands)):
        # holding booleans
        dominant = False
        
        partid = ASD_RI_probands[a]
        famid = str(df_CNV["Family ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".")
        # check both parents for affected status
        father_id = str(df_CNV["Father ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".") # find father
        mother_id = str(df_CNV["Mother ID"].loc[df_CNV["PARTID"] == partid].values).strip("[").strip("]").strip(".") # find father
        # find parent affected status
        father_status = dft.loc[dft["PARTID"] == float(father_id)]["LI"].tolist()
        mother_status = dft.loc[dft["PARTID"] == float(mother_id)]["LI"].tolist()
        #print("Son is " + str(partid) + ", father is " + str(father_id) + ", mother is " + str(mother_id))
      
        # if either parent is affected, assign family ID field to list of dominant (run this check on dft)
        if len(father_status) > 0:
            father_status = father_status[0]
            if father_status == "2":
                dominant = True
        if len(mother_status) > 0:
            if mother_status == "2":
                dominant = True
        
        # assign family to category
        if dominant == True:
            ASD_RI_dom_fams.append(famid)
        else:
            ASD_RI_denrec_fams.append(famid)
    
    #
    df_ASD_RI_dom_fams = pd.DataFrame(ASD_RI_dom_fams, columns = ["families"])
    ASD_RI_dom_fams_count = str(len(df_ASD_RI_dom_fams["families"].unique()))
    df_ASD_RI_denrec_fams = pd.DataFrame(ASD_RI_denrec_fams, columns = ["families"])
    ASD_RI_denrec_fams_count = str(len(df_ASD_RI_denrec_fams["families"].unique()))

    
    
    # Count sex
    ASD_patients = len(df_ASD["PARTID"])
    ASD_male = len(df_ASD[df_ASD["SEX"] == "1"])
    ASD_female = len(df_ASD[df_ASD["SEX"] == "2"])
    
    ASD_LI_patients = len(df_ASD_LI["PARTID"])
    ASD_LI_male = len(df_ASD_LI[df_ASD_LI["SEX"] == "1"])
    ASD_LI_female = len(df_ASD_LI[df_ASD_LI["SEX"] == "2"])
    
    ASD_RI_patients = len(df_ASD_RI["PARTID"])
    ASD_RI_male = len(df_ASD_RI[df_ASD_RI["SEX"] == "1"])
    ASD_RI_female = len(df_ASD_RI[df_ASD_RI["SEX"] == "2"])
    
    # count families with segregation patterns
    # import GMN
    AO_gmn = pd.read_csv(ASD_GMN, sep = "\t")
    AL_gmn = pd.read_csv(LI_GMN, sep = "\t")
    AR_gmn = pd.read_csv(RI_GMN, sep = "\t")
    
    
    # SORT THROUGH EACH SEG LIST AND EXTRACT FAM IDS
    
    # ASD only
    AO_seg = AO_gmn["de_novo__single_affected"].values.tolist()
    AO_fam = []
    for a in range(0, len(AO_seg)):
        if(type(AO_seg[a]) != float): # skip NaN
            seg_list = AO_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AO_fam.append(fam_ID)
        else:
           continue
    dfao_fam = pd.DataFrame(AO_fam, columns = ["families"])
    dfao_fam_count = str(len(dfao_fam["families"].unique()))
    dfao_dom_count = "0"
    dfao_rec_count = "0"
    dfao_denovo_count = dfao_fam_count
      
    # LI* all
    AL_seg = AL_gmn["de_novo__single_affected"].values.tolist() + AL_gmn["autosomal_recessive__single_affected"].values.tolist() + AL_gmn["autosomal_dominant__single_affected"].values.tolist()
    AL_fam = []
    for a in range(0, len(AL_seg)):
        if(type(AL_seg[a]) != float): # skip NaN
            seg_list = AL_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AL_fam.append(fam_ID)
        else:
           continue
    dfal_fam = pd.DataFrame(AL_fam, columns = ["families"])
    dfal_fam_count = str(len(dfal_fam["families"].unique()))
    
    # LI* de novo
    AL_denovo_seg = AL_gmn["de_novo__single_affected"].values.tolist()
    AL_denovo_fam = []
    for a in range(0, len(AL_denovo_seg)):
        if(type(AL_denovo_seg[a]) != float): # skip NaN
            seg_list = AL_denovo_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AL_denovo_fam.append(fam_ID)
        else:
           continue
    dfal_denovo_fam = pd.DataFrame(AL_denovo_fam, columns = ["families"])
    dfal_denovo_count = str(len(dfal_denovo_fam["families"].unique()))
    
    # LI* dom only
    AL_dom_seg = AL_gmn["autosomal_dominant__single_affected"].values.tolist()
    AL_dom_fam = []
    for a in range(0, len(AL_dom_seg)):
        if(type(AL_dom_seg[a]) != float): # skip NaN
            seg_list = AL_dom_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AL_dom_fam.append(fam_ID)
        else:
           continue
    dfal_dom_fam = pd.DataFrame(AL_dom_fam, columns = ["families"])
    dfal_dom_count = str(len(dfal_dom_fam["families"].unique()))
    
    # LI* rec only
    AL_rec_seg = AL_gmn["autosomal_recessive__single_affected"].values.tolist()
    AL_rec_fam = []
    for a in range(0, len(AL_rec_seg)):
        if(type(AL_rec_seg[a]) != float): # skip NaN
            seg_list = AL_rec_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AL_rec_fam.append(fam_ID)
        else:
           continue
    dfal_rec_fam = pd.DataFrame(AL_rec_fam, columns = ["families"])
    dfal_rec_count = str(len(dfal_rec_fam["families"].unique()))
       
    # RI* all
    AR_seg = AR_gmn["de_novo__single_affected"].values.tolist() + AR_gmn["autosomal_recessive__single_affected"].values.tolist() + AR_gmn["autosomal_dominant__single_affected"].values.tolist()
    AR_fam = []
    for a in range(0, len(AR_seg)):
        if(type(AR_seg[a]) != float): # skip NaN
            seg_list = AR_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AR_fam.append(fam_ID)
        else:
           continue
    dfar_fam = pd.DataFrame(AR_fam, columns = ["families"])
    dfar_fam_count = str(len(dfar_fam["families"].unique()))
    
    # RI* de novo
    AR_denovo_seg = AR_gmn["de_novo__single_affected"].values.tolist()
    AR_denovo_fam = []
    for a in range(0, len(AR_denovo_seg)):
        if(type(AR_denovo_seg[a]) != float): # skip NaN
            seg_list = AR_denovo_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AR_denovo_fam.append(fam_ID)
        else:
           continue
    dfar_denovo_fam = pd.DataFrame(AR_denovo_fam, columns = ["families"])
    dfar_denovo_count = str(len(dfar_denovo_fam["families"].unique()))
    
    # RI* dom only
    AR_dom_seg = AR_gmn["autosomal_dominant__single_affected"].values.tolist()
    AR_dom_fam = []
    for a in range(0, len(AR_dom_seg)):
        if(type(AR_dom_seg[a]) != float): # skip NaN
            seg_list = AR_dom_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AR_dom_fam.append(fam_ID)
        else:
           continue
    dfar_dom_fam = pd.DataFrame(AR_dom_fam, columns = ["families"])
    dfar_dom_count = str(len(dfar_dom_fam["families"].unique()))
    
    # RI* rec only
    AR_rec_seg = AR_gmn["autosomal_recessive__single_affected"].values.tolist()
    AR_rec_fam = []
    for a in range(0, len(AR_rec_seg)):
        if(type(AR_rec_seg[a]) != float): # skip NaN
            seg_list = AR_rec_seg[a].split("family:")
            fam_ID = seg_list[1].split("[")[0]
            AR_rec_fam.append(fam_ID)
        else:
           continue
    dfar_rec_fam = pd.DataFrame(AR_rec_fam, columns = ["families"])
    dfar_rec_count = str(len(dfar_rec_fam["families"].unique()))
       
    # calculate cohort counts
    cohort_patients = len(df_CNV["PARTID"])
    cohort_male = len(df_CNV[df_CNV["SEX"] == "1"])
    cohort_female = len(df_CNV[df_CNV["SEX"] == "2"])
    SV_cohort_patients = len(df_SV["PARTID"])
    SV_cohort_male = len(df_SV[df_SV["SEX"] == "1"])
    SV_cohort_female = len(df_SV[df_SV["SEX"] == "2"])
    SV_candidate_fams = 73
    
    # table
    summary_table = [
                    ["Phenotype", "Patients", "Male", "Female", "Families", "Dominant", "Recessive/de novo"],
                    ["ASD", ASD_patients, ASD_male, ASD_female, candidate_fams, ASD_dom_fams_count, ASD_denrec_fams_count],
                    ["LI*", ASD_LI_patients, ASD_LI_male, ASD_LI_female, candidate_fams, ASD_LI_dom_fams_count, ASD_LI_denrec_fams_count],
                    ["RI*", ASD_RI_patients, ASD_RI_male, ASD_RI_female, candidate_fams, ASD_RI_dom_fams_count, ASD_RI_denrec_fams_count],
                    ["Cohort", cohort_patients, cohort_male, cohort_female, candidate_fams, "--", "--"],
                    ["SV_cohort", SV_cohort_patients, SV_cohort_male, SV_cohort_female, SV_candidate_fams, "--", "--"]
                    ]
    #print(candidate_fams)
    
    output = open(output_file, "w")
    for a in range(0, len(summary_table)):
        output.write("\t".join([str(elem) for elem in summary_table[a]]) + "\n")
  
    

    
# zero out caller counts in excluded individuals
def revise_S1():
    import pandas as pd
    
    # import table
    dfm = pd.read_csv("./SampleSummary_NewID_2023_05_25.tsv", sep = "\t")
    
    # individuals and families to drop
    """
    - 2018 family - father is missing
    - 2108 family - not in pedigree file
    - 2117001 - not in pedigree file
    - 2078 family - don’t have microarray data for the ASD proband (2078003)
    - 2094 family - no variants called on ASD proband (209400CNV_vcf3), plate failed
    - 2111 family - no one in family has ASD (2111005 has ADHD and RI)
    - 2119003 - no data
    """
    
    # EXTRACT ALL INDIVIDUALS FROM FAMILIES
    exclusion_list = []
    
    # 2018 family
    df_2018 = dfm[dfm["Family ID"] == 2018]
    fam_2018 = df_2018["PARTID"].tolist()
    exclusion_list = exclusion_list + fam_2018
    
    # 1029 family
    df_1029 = dfm[dfm["Family ID"] == 1029]
    fam_1029 = df_1029["PARTID"].tolist()
    exclusion_list = exclusion_list + fam_1029
    
    # 2108 family
    df_2108 = dfm[dfm["Family ID"] == 2108]
    fam_2108 = df_2108["PARTID"].tolist()
    exclusion_list = exclusion_list + fam_2108
    
    # 2117001 individual
    exclusion_list.append(2117001)
    
    # 2078 family
    df_2078 = dfm[dfm["Family ID"] == 2078]
    fam_2078 = df_2078["PARTID"].tolist()
    exclusion_list = exclusion_list + fam_2078

    # 2094 family
    df_2094 = dfm[dfm["Family ID"] == 2094]
    fam_2094 = df_2094["PARTID"].tolist()
    exclusion_list = exclusion_list + fam_2094
    
    # 2111 family
    df_2111 = dfm[dfm["Family ID"] == 2111]
    fam_2111 = df_2111["PARTID"].tolist()
    exclusion_list = exclusion_list + fam_2111    
    
    # 2119003 individual
    exclusion_list.append(2119003)
    
    
    # zero out rows in xclusion list
    dfm.loc[dfm["PARTID"].isin(exclusion_list), "Axiom_CNVs"] = 0
    dfm.loc[dfm["PARTID"].isin(exclusion_list), "Illumina_CNVs"] = 0
    dfm.loc[dfm["PARTID"].isin(exclusion_list), "quantisnp_calls"] = 0
    dfm.loc[dfm["PARTID"].isin(exclusion_list), "penncnv_calls"] = 0
    dfm.loc[dfm["PARTID"].isin(exclusion_list), "caller"] = "No variants called"
    
    
    
    
    
    
    # mark ones in the included 482
    dfx = pd.read_csv("../CNV3_anno_corrected.tsv", sep = "\t") # converters={'first_column': convert_dtype,'second_column': convert_dtype}
    
    first_indv = dfx.columns.get_loc("1002001")
    last_indv = dfx.columns.get_loc("2118017") + 1
    
    retained = dfx.columns.tolist()[first_indv:last_indv]
    retained = [int(y) for y in retained]
    
    # assign field for CNV_vcf
    dfm["CNV_vcf"] = 0
    dfm.loc[dfm["PARTID"].isin(retained), "CNV_vcf"] = 1
    dfm.loc[dfm["caller"] != "No variants called", "CNV_vcf"] = 1
    
    # assign field VCF/CNV_VCF
    dfm["CNV/SV/MEI_union"] = 0
    dfm["CNV/SV/MEI_intersection"] = 0
    dfm.loc[(dfm["CNV_vcf"] == 1) | (dfm["VCF"] == 1), "CNV/SV/MEI_union"] = 1    
    dfm.loc[(dfm["CNV_vcf"] == 1) & (dfm["VCF"] == 1), "CNV/SV/MEI_intersection"] = 1    
    
    # print out counts
    #print(dfm['CNV/SV/MEI_data'].value_counts()[1])

    dff = pd.DataFrame()
    dff["Family ID"] = dfm.loc[dfm["CNV/SV/MEI_union"] == 1, "Family ID"]
    dff_list = dff["Family ID"].unique().tolist()
    
    dffc = pd.DataFrame()
    dffc["Family ID"] = dfm.loc[dfm["CNV_vcf"] == 1, "Family ID"]
    dffc_list  = dffc["Family ID"].unique().tolist()
    
    dffs = pd.DataFrame()
    dffs["Family ID"] = dfm.loc[dfm["VCF"] == 1, "Family ID"]
    dffs_list  = dffs["Family ID"].unique().tolist()
    
    df1 = dfm.loc[dfm["PARTID"].isin(retained)]
    df2 = dfm.loc[dfm["caller"] != "No variants called"]
  
    retained_families = df1["Family ID"].unique().tolist()
    microarray_families = df2["Family ID"].unique().tolist()
    
    #for j in range(0, len(retained_families)):
     #   if retained_families[j] not in microarray_families:
      #      print(retained_families[j])
            
    for k in range(0, len(microarray_families)):
        if microarray_families[k] not in retained_families:
            print(microarray_families[k])
    
    #print(len(retained_families))
    #print(len(microarray_families))
    
    
# =============================================================================
#     for a in range(0, len(dff_list)):
#         #print(dff_list[a])
#         if dff_list[a] not in dffc_list:
#             print(dff_list[a])
# =============================================================================
        
    dfm.to_csv("./SampleSummary_NewID_revised.tsv", sep = "\t")
    
"""[_______COMMANDS_______]"""
#build_table("./sample_summary_anonymized.tsv","./NJLAGS_CNV.ped", "./alternative_table.tsv")


#find_sibs("./alternative_table.tsv")



#build_sample_fam("./SampleSummary_NewID_2023_05_25.tsv", "../ASD_only/GMN_out_ASD_only_final.tsv",  "../ASD_LI/GMN_out_ASD_LI_final.tsv", "../ASD_RI/GMN_out_ASD_RI_final.tsv", "./sample_family_summary.tsv")
revise_S1()



