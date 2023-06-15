import pandas as pd

ea = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/exons_af-focused.csv")
eaf = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_exons_af-focused.csv")

ia = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/introns_af-focused.csv")
iaf = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_introns_af-focused.csv")

ie = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/introns_eqtl-focused.csv")
ief = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_introns_eqtl-focused.csv")
iea = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/introns_eqtl-focused_annotations.csv")

im = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/intergenics_eqtl-focused.csv")
imf = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_intergenics_eqtl-focused.csv")
ima = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/intergenics_eqtl-focused_annotations.csv")

'''
###
eQTL-focused candidates being added
###
'''

fin_df = pd.DataFrame()
dfs = {
    'ExonicAF': eaf,
    'IntronicAF': iaf
}

cols_to_copy = [
    'Gene_name', 'Name', 'Ensembl', 'MGD', 'Gene Type', 
    'Synonym','GnomAD_pLI', 'IMPC',
    'GTEx_brain_max', 'brainPrenatal_maxTPM', 'brainspan_maxTPM',
    'ndd_tourettes', 'mp_term_name', 'ndd_sfari', 'DDD_disease',
    'OMIM_phenotype', 'ExAC_misZ', 
    'ASD_only', 'ASD_LI', 'ASD_RI', 'Pipeline'
]

result_df = pd.DataFrame()

for name, df in dfs.items():
    df['ASD_only'] = ~pd.isna(df['ASD_only_unique_fams'])
    df['ASD_LI'] = ~pd.isna(df['ASD_LI_unique_fams'])
    df['ASD_RI'] = ~pd.isna(df['ASD_RI_unique_fams'])
    df['Pipeline'] = name
    
    df['Name'] = ''
    df['Ensembl'] = ''
    df['MGD'] = ''
    df['Strand'] = 'NA'
    df['Synonym'] = ''

    temp = df[cols_to_copy]
    result_df = pd.concat([result_df, temp])

result_df.reset_index(drop=True, inplace=True)


'''
###
eQTL-focused candidates being added
###
'''

def process_eqtl_dataframes(gene_df, full_annot, pipeline):
    rows = []
    for g, a in zip(full_annot['Gene_name'].tolist(), full_annot['AnnotSV_ID'].tolist()):
        df = gene_df[(gene_df['brain_eQTL_genes'].str.contains(g)) & (gene_df['AnnotSV_ID'] == a)]
        annot_df = full_annot[(full_annot['Gene_name'].str.contains(g)) & (full_annot['AnnotSV_ID'] == a)]

        ASD_only = ~pd.isna(df['ASD_only_unique_fams'])
        ASD_LI = ~pd.isna(df['ASD_LI_unique_fams'])
        ASD_RI = ~pd.isna(df['ASD_RI_unique_fams'])

        vals = [
            g, '', '', '', annot_df.iloc[0]['Gene Type'], '',
            df.iloc[0]['GnomAD_pLI'], annot_df.iloc[0]['IMPC'], annot_df.iloc[0]['GTEx_brain_max'],
            annot_df.iloc[0]['brainPrenatal_maxTPM'], annot_df.iloc[0]['brainspan_maxTPM'],
            annot_df.iloc[0]['ndd_tourettes'], annot_df.iloc[0]['mp_term_name'],
            annot_df.iloc[0]['ndd_sfari'], df.iloc[0]['DDD_disease'], df.iloc[0]['OMIM_phenotype'],
            df.iloc[0]['ExAC_misZ'], ASD_only.iloc[0], ASD_LI.iloc[0], ASD_RI.iloc[0], pipeline
        ]
        rows.append(vals)

    return pd.DataFrame(rows, columns=cols)

cols = [
    'Gene_name', 'Name', 'Ensembl', 'MGD', 'Gene Type', 
    'Synonym', 'GnomAD_pLI', 'IMPC', 'GTEx_brain_max', 
    'brainPrenatal_maxTPM', 'brainspan_maxTPM', 'ndd_tourettes', 
    'mp_term_name', 'ndd_sfari', 'DDD_disease', 'OMIM_phenotype', 
    'ExAC_misZ', 'ASD_only', 'ASD_LI', 'ASD_RI', 'Pipeline'
]
# Process iea and ief dataframes
result_ie = process_eqtl_dataframes(ief, iea, 'IntronicEQTL')

# Process ima and imf dataframes
result_im = process_eqtl_dataframes(imf, ima, 'IntergenicEQTL')

# Concatenate the results
gene_level_summ_MEISV = pd.concat([result_df, result_ie, result_im]).reset_index(drop=True)

# Save results
gene_level_summ_MEISV.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/gene_level_summary_table_MEISV.csv", index=None)

