import os
import re

data_dir="/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data"
MEI_dir = os.fsencode(f"{data_dir}/merge_prep-MEI_callset/step-1_merge_ready_MEI_vcfs")
SV_dir = os.fsencode(f"{data_dir}/starting_callsets/SV_callset")
sample_sets_dir = f"{data_dir}/merge_step_one/sample_sets/"

for MEI_file in os.listdir(MEI_dir):
    MEI_filename = os.fsdecode(MEI_file)
    if MEI_filename == "MEI_merged_duplicate_cleaned.vcf":
        continue
    sample = re.search('MEI_merged_duplicate_cleaned.(.*).vcf', MEI_filename).group(1)
    for SV_file in os.listdir(SV_dir):
        SV_filename = os.fsdecode(SV_file)
        if SV_filename.startswith(sample):
            # Create file
            f = open(sample_sets_dir + "/" + sample, "a")
            f.write(SV_filename + "\n")
            f.write(MEI_filename + "\n")
            f.close()