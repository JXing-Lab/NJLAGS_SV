# For reference, starting callsets used for this project are located at 
# /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/starting_callsets
SRC_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/src

## 1. Generating final individual level MEI files for merging
# 1.a. Clean up the MEI callset by removing duplicate position variants on same chromosome
python ${SRC_DIR}/CleanupDuplicateMEIsInVCF.py
# Before cleaning duplicates: 12,587 sites.
# After cleaning duplicates: 12,466 sites.

# 1.b. Split up the MEI callset into 1 VCF file per individual, and decompress the gzips after
bash ${SRC_DIR}/SplitMEIsIntoMultipleVCFs.sh



## 2. Generating final individual level SV files for merging
# No data cleaning needed because individual SV VCFs do not have multiple variants at the same position.
# Checked by taking amount of unique variants and comparing with amount of variants, using uniq and wc -l



## 3. First step of SURVIVOR merging: Intra-individual merging
# 3.a. Make folder containing all individuals' VCFs (259 MEI VCFs and 272 SV VCFs)
bash ${SRC_DIR}/SURVIVOR/CompileStepOneVCFs.sh

# 3.b. Make 259 sample names files (Other 13, Knome Batch, will not be merged since low quality in MEI callset)
#      used as an input for SURVIVOR merge
python ${SRC_DIR}/SURVIVOR/MakeStepOneSampleNames.py

# 3.c. Intra-individual merge
bash ${SRC_DIR}/SURVIVOR/MakeStepOneMergedVCFs.sh

# 3.d. Copy un-merged Knome batch individuals from SV into the set of merged VCFs. 
DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_one
cp ${DIR}/pre-merge_sample_vcfs/LP600* ${DIR}/merged_sample_vcfs/


## 4. Second step of SURVIVOR merging: Inter-individual merging
# Two steps: Make the sample name file, containing the names of all 272 individuals to be merged into one big vcf
bash ${SRC_DIR}/SURVIVOR/StepTwoFull.sh
