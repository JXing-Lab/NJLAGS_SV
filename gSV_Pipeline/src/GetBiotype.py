import pandas as pd

# Importing dataframes
candidates_sfari = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/sfari_candidates.csv")
#gene_biotypes = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211228_Biotypes/gene_biotypes.csv')
gene_biotypes = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/gencode_biotypes.csv')
# Saving gene list
candidates_sfari['Gene_name'].to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/candidates_sfari_geneList', sep='\t', index=None)

# Merging
#candidates_biotypes = pd.merge(candidates_sfari, gene_biotypes, how='left', left_on='Gene_name', right_on='hgnc_symbol')
candidates_biotypes = pd.merge(candidates_sfari, gene_biotypes, how='left', left_on='Gene_name', right_on='Gene Name')
# Save csv
candidates_biotypes_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/biotypes_candidates.csv"
candidates_biotypes.to_csv(candidates_biotypes_path, index=False)
