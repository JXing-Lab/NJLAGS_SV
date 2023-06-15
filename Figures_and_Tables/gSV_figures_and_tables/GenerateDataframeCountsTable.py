import pandas as pd

# ik = knome_included, ek = knome_excluded
ik_dir = '/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/'
ik_seg = pd.read_csv(f'{ik_dir}/seg_candidates.csv')
ik_impc = pd.read_csv(f'{ik_dir}/impc_candidates.csv')
ik_af = pd.read_csv(f'{ik_dir}/af_candidates.csv')
ik_eqtl = pd.read_csv(f'{ik_dir}/eqtl_candidates.csv')
ik_sfari = pd.read_csv(f'{ik_dir}/sfari_candidates.csv')
ik_biotypes = pd.read_csv(f'{ik_dir}/biotypes_candidates.csv')
ik_svanna = pd.read_csv(f'{ik_dir}/svanna_candidates.csv')
ik_inh = pd.read_csv(f'{ik_dir}/inh_candidates.csv')
ik_nonbenign = pd.read_csv(f'{ik_dir}/nonbenign_candidates.csv')

name_col = ['seg', 'impc', 'af', 'eqtl', 'sfari', 'biotypes', 'svanna', 'inh', 'nonbenign']
ik_sites_col = []
ik_genes_col = []

for name in name_col:
    ik_df = eval(f'ik_{name}')

    ik_sites_col.append(len(ik_df['AnnotSV_ID'].unique()))
    ik_genes_col.append(len(ik_df['Gene_name'].unique()))

df = pd.DataFrame({'df_name': name_col, 'ik_sites': ik_sites_col, 'ik_genes': ik_genes_col})
df.to_csv(f'{ik_dir}/tables/dataframe_counts_table.csv', index=None)