"""
CNV_Builder.py
> prepares a .vcf for AnnotSV annotation
> calls AnnotSV to annotate it
> prepares a .vcf for Geminesque annotation
> calls Geminesque to annotate it

Controlled by CNV_Master.py
Encompasses methods previously corresponding to Manipulator1.py, AnnotSV, Manipulator2.py, Geminesque.py

Created: 2023.02.24
Updated: 2023.06.07
"""

# CALCULATE OUTLIER THRESHOLD
def find_outliers(strictness, vcf_file, outliers_file):
    print("Running find_outliers() at strictness [" + strictness + "]...")
    
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
        # iterate through columns until a nonzero genotype is found
        for c in range(9, len(var_list)):
            geno_list = var_list[c].split(":")
            genotype = geno_list[0]
            if genotype != "0/0": 
                #print(str(b) + "| " + header_list[c] + ": " + str(CN))
                indv_index = int(dfi[dfi["PARTID"] == header_list[c]].index[0])
                dfi.at[indv_index, CN] += 1

    # calculate CN sum
    CN_columns = ["CN0","CN1", "CN3", "CN4"]
    dfi["CN_sum"] = dfi[CN_columns].sum(axis = 1)
    #dfi.to_csv("CNV_counts_prefilter.tsv", sep = "\t", index = False)

    # find median
    count_median = dfi["CN_sum"].median()
    count_mean = dfi["CN_sum"].mean()
    
    # calculate interquartile range
    Q3 = np.quantile(dfi["CN_sum"], 0.75)
    Q1 = np.quantile(dfi["CN_sum"], 0.25)
    IQR = Q3 - Q1
    # calculate outlier threshold
    if strictness == "standard":
        out_thr = Q3 + (IQR * 1.5)
    elif strictness == "extreme":
        out_thr = Q3 + (IQR * 3)
    else:
        print("Unknown strictness specified, aborting.")
    
    # output list of outliers
    dfo = dfi[dfi["CN_sum"] > out_thr]
    outlier_list = dfo["PARTID"].tolist()
    dfo.to_csv(outliers_file, sep = "\t", index = False)
    
    # print results
    print("cohort median: " + str(count_median))
    print("cohort mean: " + str(count_mean))
    print("outlier_threshold: " + str(out_thr))
    
    print("Finished running find_outliers().")

projects_folder = "/home/rohan/CNV_Project/2023/Submission_Pipeline"
find_outliers("standard", projects_folder + "/Data/CNV.vcf", projects_folder + "/Data/outlier_IDs.txt") # strictness, vcf input, outlier ID output

# SANITIZE VCF METHOD
def vcf_sanitize(input_file, output_file, outliers_file):
    print("Running vcf_sanitize() on file " + input_file)
    """
    

    Parameters
    ----------
    input_file : .vcf
        A VCF (variant call format) file with complete .vcf specification, in its original form designed to be annotated
        using ANNOVAR
    output_file : .vcf
        A modified .vcf that can be imported into AnnotSV, and is prepared to be used for structural variant (SV) analysis
        rather than single nucleotide variation (SNV) analysis

    Returns
    -------
    None.

    """
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_file, "w")
    
    # declare list for positions to be removed
    to_remove = []
    
    # iterate through every line of the file
    for line in file_in:
        # skip past specification info
        if(line[1] == "#"):
            file_out.write(line)
            continue
        
        # identify index position of family to remove
        if(line[0] == "#"):
            # divide header into list by tab
            header_list = line.strip("\n").split("\t")
        
            # iterate through header_list to find desired column
            for a in range(0, len(header_list)): # remove individuals that should be ignored due to missing data
                if("2018" in header_list[a]):
                    to_remove.append(a)
                if("2117001" in header_list[a]):
                    to_remove.append(a)
                if("2108" in header_list[a]):
                    to_remove.append(a)
                if("2078" in header_list[a]):
                    to_remove.append(a)
                if("2094" in header_list[a]):
                    to_remove.append(a)
                if("2111" in header_list[a]):
                    to_remove.append(a)
                if("2119003" in header_list[a]):
                    to_remove.append(a)
                
                # remove outliers
                outliers_in = open(outliers_file, "r")
                outliers_list = outliers_in.read().split("\n")
                if header_list[a] in outliers_list:
                    to_remove.append(a)
            
            header_removed = []
            # remove target column in header
            for b in range(0, len(to_remove)):
                header_removed.append(header_list.pop(to_remove[b]-b))
            
            #print(header_removed)
                
            
            # write out header
            header_string = "\t".join(header_list)
            file_out.write(header_string + "\n")
        else:
            # divide line into list by tab
            line_list = line.strip("\n").split("\t")
            
            # remove target columns in body
            for c in range(0, len(to_remove)):
                line_list.pop(to_remove[c]-c)
        
            # write out edited line
            line_string = "\t".join(line_list)
            file_out.write(line_string + "\n")
    
    # close files
    file_in.close()
    file_out.close()
    
    print("vcf_sanitize() on file " + input_file + " has completed running.")

# remove invariable sites
def remove_invariable(VCF2_input, VCF3_output):
    import pandas as pd
    
    df = pd.read_csv(VCF2_input, skiprows = 49, sep = "\t")
    
    vcf_in = open(VCF2_input, "r")
    vcf_list = vcf_in.read().split("\n")
    vcf_out = open(VCF3_output, "w")
    
    # write out VCF header
    header = "\n".join(vcf_list[0:49])
    vcf_out.write(header + "\n")
    
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


# CALL ANNOTSV METHOD
def call_AnnotSV():
    print("Calling AnnotSV from annotator script...")
    
    # import subprocess package
    import subprocess 
    
    # call AnnotSV
    subprocess.run('./AnnotSV_call.sh')
    
    print("Annotator script has finished calling AnnotSV.")
    
# MANIPULATE AN ANNOTATED VCF METHOD
def tsv_manipulate(input_file, output_file):
    print("Running tsv_manipulate from annotator script...")
    
    # import
    file_in = open(input_file, "r")
    file_out = open(output_file, "w")
    
    # iterate through file line by line
    for line in file_in:
        line_list = line.split("\t")
        
        
        #iterate through line
        for a in range(9, len(line_list)):
            # replace all nulls with homozygous reference
            geno_list = line_list[a].split(":")
            if(geno_list[0] == "./."):
                geno_list[0] = "0/0"
            # and reassign genotypes according to indel counts
            else:
                # 0.5/. should be called as 0/1 for heterozygous deletion
                if(geno_list[0] == "0.5/."):
                    geno_list[0] = "0/1"
                # 1.0/. should be called as 0/1 for heterozygous deletion
                elif(geno_list[0] == "1.0/."):
                    geno_list[0] = "0/1"
                # 3.0/. should be called 0/1 for heterozygous insertion
                elif(geno_list[0] == "3.0/."):
                    geno_list[0] = "0/1"
                # 3.5/. and 4.0 should both be called 1/1 for homozygous insertion
                elif(geno_list[0] == "3.5/." or geno_list[0] == "4.0/."):
                    geno_list[0] = "1/1"
            line_list[a] = ":".join(geno_list)
                
        # prepare output
        output_string = "\t".join(line_list)
        output_string = output_string.strip("\n")
        output_string = output_string + "\n"
        file_out.write(output_string)
            
    
    file_in.close()
    file_out.close()
    print("Annotator script has finished running tsv_manipulate()")

# CALL GEMINESQUE METHOD
def call_Geminesque(phenotype, AnnotSV_input):
    print("Calling Geminesque from annotator script...")
    
    # import Geminesque
    import Geminesque_2023
    
    Geminesque_2023.Geminesque_main(phenotype, AnnotSV_input)
    Geminesque_2023.add_Overlapped(phenotype, AnnotSV_input)
    
    print("Annotator script has finished calling Geminesque.")
    
    
    
    
    
    
    
    
    