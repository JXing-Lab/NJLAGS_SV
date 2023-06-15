# Part 3 - Segregation

### Set variables
GEMINI=/lab01/Tools/gemini/bin/gemini
GEM_DATA=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI
MOD_VCF_PATH=/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/modified_SURVIVOR_merged_callset.vcf    ### This was created in annotation.sh

### Remove suffixes from the copied entire_pedigree.ped file, if they still have 'em
# entire_pedigree.py is copied from /lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20210715_Gemini
awk 'BEGIN{FS=OFS="\t"} {for (i=1; i<=NF; i++) if ($i ~ /[a-zA-Z]/) sub(/\..*/, "", $i)} 1' $GEM_DATA/entire_pedigree.ped > tmp && mv tmp $GEM_DATA/entire_pedigree.ped

### Make ASD_only, ASD_LI, and ASD_RI pedigrees, if hasn't been done already
python /lab01/Projects/Sammy_Projects/NJLAGS_SVs/src/MakeCustomPedigrees.py

### Create databases (in parallel) with GEMINI load, if they don't exist already
for PHENOTYPE in ASD_only ASD_LI ASD_RI; do
    $GEMINI load -v $MOD_VCF_PATH -p $GEM_DATA/${PHENOTYPE}.ped $GEM_DATA/${PHENOTYPE}.db &
done
wait

### Adjust the databases to account for GEMINI weirdness
for DB in ASD_only ASD_LI ASD_RI; do
    MOD_DB_PATH=$GEM_DATA/modified_${DB}.db
    # # Modify a test database so don't corrupt the original
    cp "$GEM_DATA/${DB}.db" "$MOD_DB_PATH"
    sqlite3 -line "$MOD_DB_PATH" 'update variants set gene="test"'
    sqlite3 -line "$MOD_DB_PATH" 'update variants set start=start+1'
done

### Calculate inheritance patterns (Note: No mendel errors calculations because [1] calculating them here is giving me errors [2] i remove those inh patterns downstream anyways. Will not affect my results)
for PHENOTYPE in ASD_only ASD_LI ASD_RI; do
    MOD_DB_PATH=$GEM_DATA/modified_${PHENOTYPE}.db
    for MODE in autosomal_recessive autosomal_dominant x_linked_recessive x_linked_dominant x_linked_de_novo de_novo; do
        $GEMINI $MODE $MOD_DB_PATH --columns "chrom,start,end,sub_type" > $GEM_DATA/Inheritance/${PHENOTYPE}_${MODE}.tsv
    done
done