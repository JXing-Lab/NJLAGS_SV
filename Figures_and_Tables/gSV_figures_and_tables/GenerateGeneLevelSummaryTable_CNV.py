import os
import numpy as np
import pandas as pd
import os
import pickle
import openpyxl


'''
#####
Adding Brain Expression Data
#####
'''

# Copied from: /lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Data/CNV3_anno_corrected.tsv
rohan_annot = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/rohan_data/CNV3_anno_corrected.tsv", sep='\t')
seg_candidates = rohan_annot.copy()

###TPM 1
file_expression = '/lab02/Data_Raw/Xiaolong/HumanCommon/GTEx_v8/20191212GeneExpressionMedianRecalculateMinCoding.3.GTEx_v8.codingNoMT.txt' 
df_expression = pd.read_csv(file_expression,sep='\t')#54 tissues 
df_expression = df_expression[~df_expression['Name'].apply(lambda x:x.endswith('_PAR_Y'))]
del df_expression['Name']

columns_GTEx = list(df_expression.columns[1:-1]) 
columns_GTEx_brain = [e for e in columns_GTEx if "Brain" in e]
df_expression['GTEx_brain_max'] = df_expression[columns_GTEx_brain].max(axis=1) 
df_expression.rename(columns={'Description': 'GENE'}, inplace=True)

#The last column is what we care about: TPM1

df_expression = df_expression[['GENE', 'GTEx_brain_max']]
df_expression
candidates_tpm1 = seg_candidates.merge(df_expression, how = 'left', left_on = "Gene_name", right_on ='GENE') 
del candidates_tpm1['GENE']


###TPM 2
file_brainspan = '/lab02/Data_Raw/Xiaolong/HumanCommon/Expression/BrainSpan/20191216BrainSpanGeneExpression.codingNomalized.csv' 
df_brainspan = pd.read_csv(file_brainspan, sep='\t')#379 tissues 

# only keep primary_motor_cortex 

columns_brainspan = list(df_brainspan.columns[2:]) 
columns_primarymotor = [e for e in columns_brainspan if 'primary_motor_cortex' in e] 
df_brainspan['brainspan_maxTPM'] = df_brainspan[columns_brainspan].max(axis=1) 
df_primarymotor = df_brainspan[['geneSymbol', 'brainspan_maxTPM'] + columns_primarymotor] 

# The second column is what we care about: TPM2

df_primarymotor = df_primarymotor[['geneSymbol', 'brainspan_maxTPM']]
df_primarymotor
candidates_tpm2 = candidates_tpm1.merge(df_primarymotor, how = 'left', left_on = "Gene_name", right_on ='geneSymbol') 
del candidates_tpm2['geneSymbol']


###TPM3
## add gene expression in prenatal human brain development. only include "cerebral cortex" 

file_prenatal = '/lab02/Data_Raw/Xiaolong/HumanCommon/Expression/brainDevelopment/20190415brainDevelopment.binary' 
df_prenatal = pickle.load(open(file_prenatal,'rb'))#165 tissues 

df_prenatal['brainPrenatal_maxTPM'] = df_prenatal[df_prenatal.columns[2:]].max(axis=1) 
df_cerebral = df_prenatal[['Gene_Name', 'brainPrenatal_maxTPM'] + [e for e in df_prenatal.columns if "cerebral_cortex" in e]]

# The second column is what we care about: TPM3

df_cerebral = df_cerebral[['Gene_Name', 'brainPrenatal_maxTPM']]
df_cerebral
candidates_tpm3 = candidates_tpm2.merge(df_cerebral, how = 'left', left_on = "Gene_name", right_on ='Gene_Name') 
del candidates_tpm3['Gene_Name']


###IMPC
df = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20210808_IMPC/genotype-phenotype-assertions-ALL.csv') 
df1 = df[df['top_level_mp_term_name'].apply(lambda x:'behavior/neurological phenotype' in str(x) or "nervous system phenotype" in str(x))]
df2 = df1[~df1['p_value'].isna()]
df2

# get human-mouse gene ortholog relationships 
# df_homolog = pd.read_csv('/lab01/Projects/Sammy_Projects/VCF/IMPC/HOM_AllOrganism.rpt', sep='\t')
df_homolog = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20210808_IMPC/HOM_AllOrganism.rpt', sep='\t')



df2[['marker_symbol', 'top_level_mp_term_name', 'mp_term_name']]
df2['IMPC'] = True
df3 = df2[['marker_symbol', 'top_level_mp_term_name', 'mp_term_name', 'IMPC']]
df3
df4 = df_homolog.merge(df3, how='left', left_on='Symbol', right_on='marker_symbol')
#df4[(df4['Boolean']==True) & (df4['Common Organism Name']=='human')]
keylist = df4[df4['IMPC']==True]['DB Class Key'].tolist()
df4.loc[df4['DB Class Key'].isin(keylist), 'IMPC'] = True
df4[(df4['IMPC']==True) & (df4['Common Organism Name']=='human')]

#Adding mouse knockout phenotype data to matching genes in humans
impc_final = df4[(df4['IMPC']==True) & (df4['Common Organism Name']=='human')][['Symbol','IMPC']]
df5 = df4[df4['Common Organism Name']=='mouse, laboratory'][['Symbol', 'top_level_mp_term_name', 'mp_term_name']]
df5['Symbol'] = df5['Symbol'].str.upper()
df6 = df5.groupby('Symbol').agg({'top_level_mp_term_name': lambda x: list(x),'mp_term_name': lambda x: list(x)}).reset_index()
impc_final = impc_final.merge(df6, how='left', left_on='Symbol', right_on='Symbol')

#impc_final.drop_duplicates(inplace=True)
impc_final
candidates_impc = candidates_tpm3.merge(impc_final, how='left', left_on='Gene_name', right_on='Symbol')
del candidates_impc['Symbol']
candidates_impc.rename(columns={'violation': 'violation_mendel_errors', 'violation_prob': 'violation_prob_mendel_errors'}, inplace=True)
candidates_impc = candidates_impc.astype({'SV_chrom': str,'SV_start': str, 'SV_end': str})



'''
#####
Adding NDD annotations
#####
'''


ndds = pd.read_excel('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211201_NDDs/20190319TouretteRelatedGenes.all.xlsx', engine='openpyxl')

# Formatting messed up gene names
ndds.at[678,'Gene'] = 'SEPTIN10'
ndds.at[1061,'Gene'] = 'SEPTIN5'
ndds.at[1319,'Gene'] = 'MARCH11'

ndd_list = ndds['Gene'].tolist()
candidates_ndd=pd.merge(candidates_impc, ndds, how='left', left_on='Gene_name', right_on='Gene')

sfari = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211201_NDDs/SFARI.csv')
candidates_sfari = pd.merge(candidates_ndd, sfari, how='left', left_on='Gene_name', right_on='gene-symbol')
candidates_sfari.loc[pd.isna(candidates_sfari['Gene'])==False, 'ndd_tourettes'] = True
candidates_sfari.loc[pd.isna(candidates_sfari['gene-symbol'])==False, 'ndd_sfari'] = True



'''
#####
Adding Biotype annotations
#####
'''


# Importing dataframes
gene_biotypes = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/gencode_biotypes.csv')
# Merging
candidates_biotypes = pd.merge(candidates_sfari, gene_biotypes, how='left', left_on='Gene_name', right_on='Gene Name')



'''
#####
Creating the Gene-level table
#####
'''



# Getting gene_lists for each phenotype
# These files were originally copied from: /lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/A*_genes.txt
file_paths = '''
/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/rohan_data/ASD_only_genes.txt
/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/rohan_data/ASD_LI_genes.txt
/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/rohan_data/ASD_RI_genes.txt
'''.split()
rohans_candidates = {}
for file in file_paths:
    with open(file) as f:
        v = [x.strip() for x in f.readlines()]
        rohans_candidates[os.path.basename(file).split('.')[0]] = v
rohans_candidates['All'] = [y for x in rohans_candidates.values() for y in x]

# Getting the three phenotype columns
ASD_only = []
ASD_LI = []
ASD_RI = []
df = candidates_biotypes
df = df[(~pd.isna(df['Gene_name'])) & (df['Annotation_mode']=='split')]
genes = df['Gene_name'].tolist()
for g in genes:
    if g in rohans_candidates['ASD_only_genes']:
        ASD_only.append(True)
    else:
        ASD_only.append(False)

    if g in rohans_candidates['ASD_LI_genes']:
        ASD_LI.append(True)
    else:
        ASD_LI.append(False)

    if g in rohans_candidates['ASD_RI_genes']:
        ASD_RI.append(True)
    else:
        ASD_RI.append(False)

# Adding the three phenotype columns
df['ASD_only'] = ASD_only
df['ASD_LI'] = ASD_LI
df['ASD_RI'] = ASD_RI
filt = (df['ASD_only']==False) & (df['ASD_LI']==False) & (df['ASD_RI']==False)
df = df[~filt]

# Other added columns
df['Pipeline'] = 'CNV'
df['Name'] = ''
df['Ensembl'] = ''
df['MGD'] = ''
df['Strand'] = 'NA'
df['Synonym'] = ''

# Displaying the final product
gene_level_summ_CNV = df[[
    'Gene_name', 'Name', 'Ensembl', 'MGD', 'Gene Type', 'Synonym',
    'GnomAD_pLI', 'IMPC', 'GTEx_brain_max', 'brainPrenatal_maxTPM',
    'brainspan_maxTPM', 'ndd_tourettes', 'mp_term_name', 'ndd_sfari',
    'DDD_disease', 'OMIM_phenotype', 'ExAC_misZ', 'ASD_only', 'ASD_LI',
    'ASD_RI', 'Pipeline'
]]

gene_level_summ_CNV.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/gene_level_summary_table_CNV.csv", index=None)
