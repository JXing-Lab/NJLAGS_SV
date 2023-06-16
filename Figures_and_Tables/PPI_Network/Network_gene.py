import pandas as pd 

# get genes 

file_ASD = '/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230521xiaolong/20230524SFARI.genes.csv' 

file_other = "/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230521xiaolong/Table S2 Previously implicated NDD Genes.xlsx" 

file_target = '/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/total_uniq_genes_strict.txt' 

 

df_ASD = pd.read_csv(file_ASD) 

df_other = pd.read_excel(file_other, skiprows=12) 

df_target = pd.read_csv(file_target,header=None) 

df_target.columns = ['Gene_name'] 

 

genes_ADHD = set(df_other[df_other['ADHD_evidence'] > 0]['gene_symbol'])#354 

genes_other_Neuro = set(df_other[df_other['other_Neuro'] > 0]['gene_symbol'])#233 

genes_SFARI = set(df_ASD['gene-symbol'])#1128 

df_target['target'] = 1 

df_target['ADHD_evidence'] = df_target['Gene_name'].isin(genes_ADHD).astype(int) 

df_target['other_Neuro'] = df_target['Gene_name'].isin(genes_other_Neuro).astype(int) 

df_target['SFARI'] = df_target['Gene_name'].isin(genes_SFARI).astype(int) 

 

columns_list = ['target','ADHD_evidence','other_Neuro','SFARI'] 

for c in columns_list: 

    print(c, df_target[c].sum()) 

 

# target 264 

# ADHD_evidence 17 

# other_Neuro 4 

# SFARI 41 

 

df_use = df_target[['Gene_name'] + columns_list] 

df_use.to_csv('/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609genes.tsv', sep='\t', index=None) 

 

dc_genes = {} 

dc_genes['target'] = set(df_target['Gene_name']) 

dc_genes['ADHD_evidence'] = genes_ADHD 

dc_genes['other_Neuro'] = genes_other_Neuro 

dc_genes['SFARI'] = genes_SFARI 

 

genes_all = set([i for j in dc_genes.values() for i in j])#1726 

df_genes = pd.DataFrame(index=sorted(genes_all)) 

for c in columns_list: 

    df_genes[c] = df_genes.index.isin(dc_genes[c]).astype(int) 

 

df_genes.index.name = 'gene_symbol' 

df_genes = df_genes.reset_index() 

 

df_genes.to_excel('/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609genes_all.xlsx', index=None) 
