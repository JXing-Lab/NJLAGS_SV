## Two steps
DATA_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data
# A. Make the sample name file, containing the names of all 259 individuals to be merged into one big vcf
cd ${DATA_DIR}/merge_step_one/merged_sample_vcfs
ls > ${DATA_DIR}/merge_step_two/sample_set.txt

# B. Inter-individual merge
cd ${DATA_DIR}/merge_step_one/merged_sample_vcfs
/lab01/Tools/SURVIVOR/Debug/SURVIVOR merge ${DATA_DIR}/merge_step_two/sample_set.txt 100 0 1 1 0 0 ${DATA_DIR}/merge_step_two/SURVIVOR_merged_callset.vcf
