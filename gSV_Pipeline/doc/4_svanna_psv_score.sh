
#1. Make bed file
SVANNA_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/SvAnna
SRC_DIR=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/src

mkdir -p $SVANNA_DIR
python ${SRC_DIR}/MakeSvannaBed.py


#2. Copy the bed file and make a concatenated column
cp ${SVANNA_DIR}/liftOver_ready.bed ${SVANNA_DIR}/liftOver_ready_concat.bed
awk '{$5=$1";"$2";"$3";"$4; print}' ${SVANNA_DIR}/liftOver_ready_concat.bed > tmp && mv tmp ${SVANNA_DIR}/liftOver_ready_concat.bed


#3. Run liftOver on both bed files
CHAIN=/lab01/Tools/liftOver/hg19ToHg38.over.chain

mkdir -p ${SVANNA_DIR}/waste
/lab01/Tools/liftOver/liftOver ${SVANNA_DIR}/liftOver_ready.bed $CHAIN ${SVANNA_DIR}/lifted.bed ${SVANNA_DIR}/waste/unlifted.bed
/lab01/Tools/liftOver/liftOver ${SVANNA_DIR}/liftOver_ready_concat.bed $CHAIN ${SVANNA_DIR}/lifted_concat.bed ${SVANNA_DIR}/waste/unlifted_concat.bed


# 4. Take the bed file without a 5th column, and turn it into a VCF file using my own script
python ${SRC_DIR}/MakeSvannaVCF.py


# 5. Run SvAnna on the VCF file and decompress the gzipped output
java -jar /lab01/Tools/SvAnna/svanna-cli-1.0.2-distribution/svanna-cli-1.0.2/svanna-cli-1.0.2.jar \
    prioritize -d /lab01/Tools/SvAnna/svanna-cli-1.0.2-distribution/svanna-cli-1.0.2/svanna-data \
    --output-format csv \
    --vcf ${SVANNA_DIR}/SvAnna_ready.vcf \
    --out-dir ${SVANNA_DIR}/Output
cp ${SVANNA_DIR}/Output/SvAnna_ready.SVANNA.csv.gz ${SVANNA_DIR}/Output/copy_SvAnna_ready.SVANNA.csv.gz
gzip -d ${SVANNA_DIR}/Output/SvAnna_ready.SVANNA.csv.gz


# 6. Data formatting and merging to make svanna candidates
# 7. Load in the bed file with 5 columns as a dataframe, and the csv as another dataframes.
# 8. Merge the two on the first 4 columns of the bed file.
# 9. Now you have a dataframe with the desired PSV values, and the corresponding hg19 positions. Merge with
python ${SRC_DIR}/FormatAndMergeSvanna.py



