import pandas as pd
results_dir = "/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/recyclebin/Results"

missingindv_alu = pd.read_csv(f"{results_dir}/HQ_missingindv_alu.imiss", sep="\t")
missingindv_aluDel = pd.read_csv(f"{results_dir}/HQ_missingindv_aluDel.imiss", sep="\t")
missingindv_line1 = pd.read_csv(f"{results_dir}/HQ_missingindv_line1.imiss", sep="\t")
missingindv_line1Del = pd.read_csv(f"{results_dir}/HQ_missingindv_line1Del.imiss", sep="\t")
missingindv_sva = pd.read_csv(f"{results_dir}/HQ_missingindv_sva.imiss", sep="\t")
missingindv_svaDel = pd.read_csv(f"{results_dir}/HQ_missingindv_svaDel.imiss", sep="\t")

missing_indv = missingindv_alu + missingindv_aluDel + missingindv_line1 + missingindv_line1Del + missingindv_sva + missingindv_svaDel
missing_indv["F_MISS"] = missing_indv["N_MISS"] / missing_indv["N_DATA"]

# Saving table
missing_indv.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/MEI_missing_rate_by_individual.csv", index=None)

# Missing rates for knome versus other batches
knome_sites = missing_indv['INDV'].str.contains('LP6005170')
print("Avg. Knome batch missing rate: ", missing_indv[knome_sites]['F_MISS'].mean()*100)
print("Avg. Individual (excluding Knome batch) missing rate: ", missing_indv[~knome_sites]['F_MISS'].mean()*100)

