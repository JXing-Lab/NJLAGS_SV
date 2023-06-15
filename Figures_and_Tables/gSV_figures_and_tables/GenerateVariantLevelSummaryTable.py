import pandas as pd
import numpy as np
from ast import literal_eval

def replace_ids(row):
    if isinstance(row, str):
        numerical_ids = row.split(',')
        string_ids = [fam_mapping.get(int(id_), id_) for id_ in numerical_ids]
        return ','.join(string_ids)
    else:
        return fam_mapping.get(row, row)


data_dir = '/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data'
eaf = pd.read_csv(f"{data_dir}/results/full_exons_af-focused.csv")
iaf = pd.read_csv(f"{data_dir}/results/full_introns_af-focused.csv")
ief = pd.read_csv(f"{data_dir}/results/full_introns_eqtl-focused.csv")
imf = pd.read_csv(f"{data_dir}/results/full_intergenics_eqtl-focused.csv")

eaf['Pipeline'] = "ExonicAF"
iaf['Pipeline'] = "IntronicAF"
ief['Pipeline'] = "IntronicEQTL"
imf['Pipeline'] = "IntergenicEQTL"

# Changing gene name column to only look at eQTLs
ief['Gene_name'] = ief['brain_eQTL_genes']
imf['Gene_name'] = imf['brain_eQTL_genes']

df = pd.concat([eaf, iaf, ief, imf]).reset_index(drop=True)
df['ASD_only'] = ~pd.isna(df['ASD_only_unique_fams'])
df['ASD_LI'] = ~pd.isna(df['ASD_LI_unique_fams'])
df['ASD_RI'] = ~pd.isna(df['ASD_RI_unique_fams'])

variant_level_summary = df[[
    'SV_chrom', 'SV_start', 'SV_end', 'SV_length', 'SV_type', 'REF', 'ALT',
    'Gene_name', 'ACMG_class', 'psv', 'AF', 'popAF', 'tissue_count_ME_isCausal','causal_tissues',
    'ASD_only', 'ASD_only_unique_fams', 'ASD_LI', 'ASD_LI_unique_fams', 'ASD_RI', 'ASD_RI_unique_fams', 
    'Inheritance_Patterns', 'Pipeline'
]].drop_duplicates()

# Saving table
variant_level_summary.rename(columns={"AF": "Sample_AF", "popAF": "Pop_AF", "SV_chrom": "Chrom", "SV_start": "Start", "SV_end": "End", "SV_length": "Length", "SV_type": "Type", "ACMG_class": "ACMG_Class", "psv": "SvAnna_PSV"}, inplace=True)
variant_level_summary['ACMG_Class'] = variant_level_summary['ACMG_Class'].replace('full=3','3').replace('full=NA', np.nan)

# Read ID mapping as a df and make it a dictionary
samp_summ_by_id = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/2023_05_09_sample_summary_by_id.csv')[['Family ID','New Family ID']]
fam_mapping = {k: v for k,v in zip(samp_summ_by_id['Family ID'].tolist(), samp_summ_by_id['New Family ID'].tolist())}

columns_to_use_mapping = [x for x in variant_level_summary.columns if 'unique_fams' in x]
# Replace numerical IDs with string IDs using the dictionary
for col in columns_to_use_mapping:
    variant_level_summary[col] = variant_level_summary[col].map(replace_ids)


variant_level_summary['Gene_name'] = variant_level_summary['Gene_name'].apply(lambda x: x if '[' not in x else ';'.join(literal_eval(x)))
variant_level_summary.to_csv(f"{data_dir}/tables/variant_level_summary_table.csv", index=None)
