DATA_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data


mkdir -p ${DATA_DIR}/merge_step_one/pre-merge_sample_vcfs

# Getting SVs
cp ${DATA_DIR}/starting_callsets/SV_callset/* ${DATA_DIR}/merge_step_one/pre-merge_sample_vcfs/

# Getting MEIs
cp ${DATA_DIR}/merge_prep-MEI_callset/step-1_merge_ready_MEI_vcfs/* ${DATA_DIR}/merge_step_one/pre-merge_sample_vcfs/
rm ${DATA_DIR}/merge_step_one/pre-merge_sample_vcfs/MEI_merged_duplicate_cleaned.vcf # This is the leftover file that was split into the 259 MEI VCFs

# Creating sample_sets dir for next script in pipeline
mkdir -p ${DATA_DIR}/merge_step_one/sample_sets