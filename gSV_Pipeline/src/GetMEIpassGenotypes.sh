main_dir='/lab01/Projects/Sammy_Projects/'
# Insertions
grep -v "#" ${main_dir}/NJLAGS_MEIs/data/20200623_RawData/20200623ALU.final_comp.vcf | awk -F'\t' '{ for(i=10; i<=NF; i++) printf "%s ", $i; printf "\n" }' > ${main_dir}/NJLAGS_SVs/data/MEI_pass_genotypes/ALU.csv
grep -v "#" ${main_dir}/NJLAGS_MEIs/data/20200623_RawData/20200623LINE1.final_comp.vcf | awk -F'\t' '{ for(i=10; i<=NF; i++) printf "%s ", $i; printf "\n" }' > ${main_dir}/NJLAGS_SVs/data/MEI_pass_genotypes/LINE1.csv
grep -v "#" ${main_dir}/NJLAGS_MEIs/data/20200623_RawData/20200623SVA.final_comp.vcf | awk -F'\t' '{ for(i=10; i<=NF; i++) printf "%s ", $i; printf "\n" }' > ${main_dir}/NJLAGS_SVs/data/MEI_pass_genotypes/SVA.csv
# Deletions
grep -v "#" ${main_dir}/NJLAGS_MEIs/data/20200623_RawData/20200623dALU.final_comp.vcf | awk -F'\t' '{ for(i=10; i<=NF; i++) printf "%s ", $i; printf "\n" }' > ${main_dir}/NJLAGS_SVs/data/MEI_pass_genotypes/dALU.csv
grep -v "#" ${main_dir}/NJLAGS_MEIs/data/20200623_RawData/20200623dLINE1.final_comp.vcf | awk -F'\t' '{ for(i=10; i<=NF; i++) printf "%s ", $i; printf "\n" }' > ${main_dir}/NJLAGS_SVs/data/MEI_pass_genotypes/dLINE1.csv
grep -v "#" ${main_dir}/NJLAGS_MEIs/data/20200623_RawData/20200623dSVA.final_comp.vcf | awk -F'\t' '{ for(i=10; i<=NF; i++) printf "%s ", $i; printf "\n" }' > ${main_dir}/NJLAGS_SVs/data/MEI_pass_genotypes/dSVA.csv