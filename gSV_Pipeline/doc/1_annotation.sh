### Modify the SURVIVOR_merged_callset.vcf, so that it works with AnnotSV and GEMINI.
END_OF_HEADER=127
MERGED_VCF_PATH=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/SURVIVOR_merged_callset.vcf
MOD_VCF_PATH=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/modified_SURVIVOR_merged_callset.vcf


# if [ ! -f "$MOD_VCF_PATH" ]; then
#     # Get rid of CIPOS and CIEND fields
#     awk -v END_OF_HEADER=$END_OF_HEADER 'BEGIN {FS=OFS="\t"} NR<=END_OF_HEADER {print} NR>END_OF_HEADER {gsub(/;CIPOS=[^;]+/,"",$8); gsub(/;CIEND=[^;]+/,"",$8); print}' $MERGED_VCF_PATH > $MOD_VCF_PATH
#     # Add SVTYPE of IDP to the fields that need it
#     awk 'BEGIN {FS=OFS="\t"} $8 ~ /SVTYPE=NA/ {sub(/SVTYPE=NA/, "SVTYPE=IDP", $8)} {print}' $MOD_VCF_PATH > tmp && mv tmp $MOD_VCF_PATH
#     # Add header tag for IDP
#     sed -i '101a\##ALT=<ID=IDP,Description="Interspersed Duplication">' $MOD_VCF_PATH

# else
#     echo "Modified SURVIVOR merged vcf already exists."
# fi

# Get rid of CIPOS and CIEND fields
awk -v END_OF_HEADER=$END_OF_HEADER 'BEGIN {FS=OFS="\t"} NR<=END_OF_HEADER {print} NR>END_OF_HEADER {gsub(/;CIPOS=[^;]+/,"",$8); gsub(/;CIEND=[^;]+/,"",$8); print}' $MERGED_VCF_PATH > $MOD_VCF_PATH
# Add SVTYPE of IDP to the fields that need it
awk 'BEGIN {FS=OFS="\t"} $8 ~ /SVTYPE=NA/ {sub(/SVTYPE=NA/, "SVTYPE=IDP", $8)} {print}' $MOD_VCF_PATH > tmp && mv tmp $MOD_VCF_PATH
# Add header tag for IDP
sed -i '101a\##ALT=<ID=IDP,Description="Interspersed Duplication">' $MOD_VCF_PATH

### Create annotation file
ANNOTATION_PATH=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/annotations/annotated.tsv
export ANNOTSV=/lab01/Tools/AnnotSV/

# if [ ! -f "$ANNOTATION_PATH" ]; then
#     $ANNOTSV/bin/AnnotSV -SVinputFile $MOD_VCF_PATH -outputFile $ANNOTATION_PATH -SVinputInfo 1
# else
#     echo "Annotation file already exists."
# fi

$ANNOTSV/bin/AnnotSV -SVinputFile $MOD_VCF_PATH -outputFile $ANNOTATION_PATH -SVinputInfo 1