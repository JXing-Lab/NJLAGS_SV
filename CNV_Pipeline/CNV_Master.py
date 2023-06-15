"""
CNV_Master.py
> calls other files and scripts to run CNV pipeline

Created: 2023.01.12
Updated: 2023.06.08

"""

# begin
print("CNV_Master.py script is running...\n")

# DECLARE GLOBAL VARIABLES
projects_folder = "/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline"
phenotype = "ASD_RI" # set phenotype
print("Projects folder located at: " + projects_folder)



"""
╒◖════════════════════◗╕
    VCF BUILD STEP
╘◖════════════════════◗╛
"""


# Call pre-annotation assembly and formatting scripts
import CNV_Builder as builder

# PHASE 1 COMMANDS
# merge starting call files by caller
builder.CNV_combine(projects_folder + "/Data/starting_call_files/")
# sort and merge by batch
builder.CNV_sort(projects_folder + "/Data/merge_outputs/axiom_unsorted.tsv", projects_folder)
builder.caller_merge(projects_folder + "/Data/merge_outputs/axiom_sorted.tsv", projects_folder)
builder.CNV_sort(projects_folder + "/Data/merge_outputs/jan_unsorted.tsv", projects_folder)
builder.caller_merge(projects_folder + "/Data/merge_outputs/jan_sorted.tsv", projects_folder)
builder.CNV_sort(projects_folder + "/Data/merge_outputs/may_unsorted.tsv", projects_folder)
builder.caller_merge(projects_folder + "/Data/merge_outputs/may_sorted.tsv", projects_folder)

# merge to write and form all_sorted
axiom_in = open(projects_folder + "/Data/merge_outputs/axiom_unsorted.tsv", "r")
jan_in = open(projects_folder + "/Data/merge_outputs/jan_unsorted.tsv", "r")
may_in = open(projects_folder + "/Data/merge_outputs/may_unsorted.tsv", "r")
all_out = open(projects_folder + "/Data/merge_outputs/all_unsorted.tsv", "w")
# read in all files as stri
axiom_string = axiom_in.read()
jan_string = jan_in.read()
may_string = may_in.read()
all_string = axiom_string + "\n" + jan_string + "\n" + may_string
# write out combined file, merge, and sort
all_out.write(all_string)
builder.CNV_sort(projects_folder + "/Data/merge_outputs/all_unsorted.tsv", projects_folder)
builder.caller_merge(projects_folder + "/Data/merge_outputs/all_sorted.tsv", projects_folder) # this produces the all_merged.tsv file
# close open files
axiom_in.close()
jan_in.close()
may_in.close()
all_out.close()


# PHASE 2 COMMANDS
# compile list of individuals
builder.indv_find(projects_folder + "/Data/merge_outputs/", projects_folder + "/Data/individuals.txt")

# batch merge
builder.batch_merge(projects_folder + "/Data/merge_outputs/axiom_sorted.tsv", projects_folder + "/Data/individuals.txt", "_merge_intermediary_1", projects_folder)
builder.batch_merge(projects_folder + "/Data/merge_outputs/jan_sorted.tsv", projects_folder + "/Data/individuals.txt", "_merge_intermediary_1", projects_folder)
builder.batch_merge(projects_folder + "/Data/merge_outputs/may_sorted.tsv", projects_folder + "/Data/individuals.txt", "_merge_intermediary_1", projects_folder)

# TO DO: Fix merging error
builder.batch_merge(projects_folder + "/Data/merge_outputs/all_sorted.tsv", projects_folder + "/Data/individuals.txt", "_merge_intermediary_1", projects_folder)

# remerge consecutive positions
builder.CNV_remerge(projects_folder + "/Data/merge_outputs/all_merge_intermediary_1.tsv", projects_folder + "/Data/merge_outputs/all_merge_intermediary_2.tsv", 1, projects_folder + "/Data/error_reports/error1.txt") # returns 16 remerges
builder.CNV_remerge(projects_folder + "/Data/merge_outputs/all_merge_intermediary_2.tsv", projects_folder + "/Data/merge_outputs/all_merge_intermediary_3.tsv", 2, projects_folder + "/Data/error_reports/error2.txt") # returns 0 remerges

# adjacency-agnostic merge
remerges_required = 99999
file_counter = 2
cycle_count = 1

# remerger loop
while remerges_required > 0:
    remerge_input = projects_folder + "/Data/merge_outputs/all_merge_intermediary_" + str(file_counter) + ".tsv"
    remerge_output = projects_folder + "/Data/merge_outputs/all_merge_intermediary_" + str(file_counter + 1) + ".tsv"
    error_output = projects_folder + "/Data/error_reports/error_" + str(file_counter + 1) + ".txt"
    remerges_required = builder.AJA_remerge(remerge_input, remerge_output, cycle_count, error_output) # returns 120 remerges
    file_counter += 1
    cycle_count += 1
final_remerge = remerge_output

    
# PHASE 3 COMMANDS
# convert to vcf
builder.convert_to_vcf(final_remerge, projects_folder + "/Data/CNV.vcf", projects_folder + "/Data/human_g1k_v37.fasta")


"""
╒◖════════════════════◗╕
    ANNOTATION STEP
╘◖════════════════════◗╛
"""


import CNV_Annotator as annotator


# # sanitize .vcf into structural variant vcf 
annotator.find_outliers("standard", projects_folder + "/Data/CNV.vcf", projects_folder + "/Data/outlier_IDs.txt") # strictness, vcf input, outlier ID output
annotator.vcf_sanitize(projects_folder + "/Data/CNV.vcf", projects_folder + "/Data/CNV2.vcf", projects_folder + "/Data/outlier_IDs.txt")
annotator.remove_invariable(projects_folder + "/Data/CNV2.vcf", projects_folder + "/Data/CNV3.vcf")
annotator.call_AnnotSV()
annotator.tsv_manipulate(projects_folder + "/Data/CNV3_anno.tsv", projects_folder + "/Data/CNV3_anno_corrected.tsv")
# call post-annotation formatting and segregation analysis
annotator.call_Geminesque(phenotype, projects_folder + "/Data/CNV3_anno_corrected.tsv")



"""
╒◖═══════════════════════◗╕
    PRIORITIZATION STEP
╘◖═══════════════════════◗╛
"""

import CNV_Prioritizer as prioritizer

# GMNP COMMANDS
filepath = projects_folder + "/Data/"
leniency = "strict"
prioritizer.var_process(filepath + phenotype + "/GMN_out_" + phenotype + "_final.tsv", filepath + "CNV3_anno_corrected.tsv", filepath + phenotype + "/" + phenotype + ".ped","../Data/" + phenotype + "/GMNP_intermediate1.tsv")
prioritizer.var_sort(filepath + phenotype + "/GMNP_intermediate1.tsv", filepath + phenotype + "/GMNP_intermediate2.tsv")
prioritizer.var_prioritize(filepath + phenotype + "/GMNP_intermediate2.tsv", filepath + phenotype + "/GMNP_intermediate3.tsv", filepath + "dispensable_genes_and_muc.txt", filepath + "NDD_genes.tsv", leniency)
prioritizer.add_overlap_CDS(filepath + phenotype + "/GMNP_intermediate3.tsv" + "_" + leniency, filepath + "CNV3_anno_corrected.tsv", filepath + phenotype + "/GMNP_intermediate4.tsv")
prioritizer.add_gtex(filepath + phenotype + "/GMNP_intermediate4.tsv", filepath + phenotype + "/prioritized_variants.tsv", filepath + "GTEx_2019_12_12.txt")
prioritizer.add_expression(filepath + phenotype + "/GMNP_intermediate4.tsv", filepath + phenotype + "/GMNP_intermediate5.tsv", filepath + "tpm1.txt", filepath + "tpm2.txt", filepath + "tpm3.txt")

# call StrVCTVRE scripts
prioritizer.build_bed(filepath + phenotype + "/prioritized_variants.tsv", filepath + phenotype + "/GMN.bed")
prioritizer.call_StrVCTVRE(phenotype, projects_folder)
prioritizer.add_STR(filepath + phenotype + "/GMN_StrVCTVRE.tsv", filepath + phenotype + "/GMNP_intermediate5.tsv", filepath + phenotype + "/GMNP_out_" + phenotype + ".tsv")

# Call SvAnna pathogenicity ranking script
prioritizer.VCF_to_BED(filepath + "CNV3.vcf", filepath + "CNV_b37.bed")
prioritizer.b37_to_hg38(filepath + "CNV_b37.bed", filepath + "CNV_hg19.bed", filepath + "CNV_hg38.bed")
prioritizer.BED_to_VCF(filepath + "CNV_hg38.bed", filepath + "SvAnna_input_" + phenotype + ".vcf", filepath + "CNV3.vcf")
prioritizer.call_SvAnna(filepath + "SvAnna_input_" + phenotype + ".vcf", filepath + phenotype, phenotype)
prioritizer.SvAnna_assign(filepath + phenotype + "/SvAnna_input_" + phenotype + ".SVANNA.csv", filepath + phenotype + "/GMNP_out_" + phenotype + ".tsv", filepath + phenotype + "/" + phenotype + "_results.tsv")

# final bash command to output results
import subprocess
subprocess.run(["cp", filepath + phenotype + "/" + phenotype + "_results.tsv", projects_folder + "/Results/"])
    
    
"""
╒◖═══════════════════════◗╕
       SUMMARY STEP
╘◖═══════════════════════◗╛
"""

import CNV_Summarizer as summarizer
