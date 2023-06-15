import numpy as np
import pandas as pd
import openpyxl

''' 
Note that tourettes list is used as a filter in the intronic AF-focused pipeline.
SFARI annotation added here is not used later on, since I add another SFARI annotation when creating the variant level table:
See the GetVariantLevelSummaryTable.py script.
'''

candidates_eqtl = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/eqtl_candidates.csv")

# Tourette Related Genes list comes from Xiaolong Cao's paper: https://www.nature.com/articles/s41380-021-01094-1
ndds = pd.read_excel('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211201_NDDs/20190319TouretteRelatedGenes.all.xlsx', engine='openpyxl')

# Formatting messed up gene names
ndds.at[678,'Gene'] = 'SEPTIN10'
ndds.at[1061,'Gene'] = 'SEPTIN5'
ndds.at[1319,'Gene'] = 'MARCH11'

ndd_list = ndds['Gene'].tolist()
candidates_ndd=pd.merge(candidates_eqtl, ndds, how='left', left_on='Gene_name', right_on='Gene')

sfari = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211201_NDDs/SFARI.csv')
candidates_sfari = pd.merge(candidates_ndd, sfari, how='left', left_on='Gene_name', right_on='gene-symbol')
candidates_sfari.loc[pd.isna(candidates_sfari['Gene'])==False, 'ndd_tourettes'] = True
candidates_sfari.loc[pd.isna(candidates_sfari['gene-symbol'])==False, 'ndd_sfari'] = True

# Save csv
candidates_sfari_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/sfari_candidates.csv"
candidates_sfari.to_csv(candidates_sfari_path, index=False)