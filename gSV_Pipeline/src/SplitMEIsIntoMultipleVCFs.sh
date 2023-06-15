MAIN_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_prep-MEI_callset/
mkdir -p ${MAIN_DIR}/step-1_merge_ready_MEI_vcfs
cp ${MAIN_DIR}/MEI_merged_duplicate_cleaned.vcf ${MAIN_DIR}/step-1_merge_ready_MEI_vcfs/MEI_merged_duplicate_cleaned.vcf
cd ${MAIN_DIR}/step-1_merge_ready_MEI_vcfs

# Splitting
FILE=MEI_merged_duplicate_cleaned.vcf
for sample in `bcftools query -l $FILE`; do
    bcftools view -c0 -Oz -s $sample -o ${FILE/.vcf*/.$sample.vcf.gz} $FILE
done

# Decompression
for file in *.gz; do
    gzip -d $file
done
