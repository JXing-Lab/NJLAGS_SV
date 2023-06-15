import pandas as pd
from biomart import BiomartServer
import io

# ## Keep commented out since reading the excel file takes forever
# For potential brain eQTLs, let’s use the Table S2 from the GTEx SV paper:
# https://genome.cshlp.org/content/31/12/2249.long.
# Table S2 direct link:
# https://genome.cshlp.org/content/suppl/2021/11/12/gr.275488.121.DC1/Supplemental_Table_S2.xlsx 
supp_table_2 = pd.read_excel("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GTEx_SV_data/Supplemental_Table_S2.xlsx", sheet_name='eQTLs')

# Formatting raw excel file to be usable with pandas
supp_table_2 = supp_table_2.drop(0)
supp_table_2.columns = supp_table_2.iloc[0]
supp_table_2 = supp_table_2.drop(1).reset_index()
supp_table_2 = supp_table_2.drop(columns=["index"])
# Save to CSV for ease of use
supp_table_2.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GTEx_SV_data/Supplemental_Dataframe_S2.csv", index=None)