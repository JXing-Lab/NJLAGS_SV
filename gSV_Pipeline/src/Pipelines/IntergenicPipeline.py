import pandas as pd

# define a function to remove blacklisted genes from each row
def remove_blacklist(row, blacklist):
    return [gene for gene in row['brain_eQTL_genes'] if gene not in blacklist]

candidates_nonbenign=pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/nonbenign_candidates.csv")
# Getting summarized family data
All = [x for x in candidates_nonbenign.columns if "family_id" in x]
ASD_only = [x for x in candidates_nonbenign.columns if "family_id_ASD_only" in x]
ASD_LI = [x for x in candidates_nonbenign.columns if "family_id_ASD_LI" in x]
ASD_RI = [x for x in candidates_nonbenign.columns if "family_id_ASD_RI" in x]

phenos = ["All", "ASD_only", "ASD_LI", "ASD_RI"]
for ph in phenos:
    ph_list = eval(ph)

    s = candidates_nonbenign[ph_list[0]].fillna('').str.replace(r'\[|\]', ',')
    for col in ph_list[1:]:
        s += candidates_nonbenign[col].fillna('').str.replace(r'\[|\]', ',')

    unique_fams = s.str.split(',').apply(lambda x: ','.join(set(x) - {''}))
    nb_unique_fams = [len(x.split(",")) for x in unique_fams]

    candidates_nonbenign[ph+"_unique_fams"] = unique_fams
df_expression = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/df_expression.csv")
df_primarymotor = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/df_primarymotor.csv")
df_cerebral = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/df_cerebral.csv")
impc_final = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/impc_final.csv")

## Tourettes list NDDS
ndds = pd.read_excel('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211201_NDDs/20190319TouretteRelatedGenes.all.xlsx', engine='openpyxl')
# Formatting messed up gene names
ndds.at[678,'Gene'] = 'SEPTIN10'
ndds.at[1061,'Gene'] = 'SEPTIN5'
ndds.at[1319,'Gene'] = 'MARCH11'
## SFARI
sfari = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211201_NDDs/SFARI.csv')
## BIOTYPES
# gene_biotypes = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211228_Biotypes/gene_biotypes.csv')
gene_biotypes = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/gencode_biotypes.csv')


#filt for intergenic sites
intergenics = candidates_nonbenign[pd.isna(candidates_nonbenign['Gene_name'])==True]

# #filtering for eqtl contianign variants
# intergenics_eqtl = intergenics[intergenics['isCausal_for_brain_tissue']==True]
### I prefiltered eQTLs to be brain eQTLs, so just need to check for possession of eQTLs. 
intergenics_eqtl = intergenics[intergenics['tissue_count_ME_isCausal']>0]

# #Get brain annotations and other annotations associated 
# all_genes = intergenics_eqtl['eQTL_genes'].tolist()
# all_brain_genes = list()
# for gene_list in all_genes:
#     raw_list = gene_list.split('| ')
#     brain_genes = list()
#     for gene in raw_list:
#         if 'Brain' in gene:
#         # appending gene name to list
#             brain_genes.append(gene.split(' ')[0])
#     all_brain_genes.append(brain_genes)
    
# brain_gene_list = list()
# for i in all_brain_genes:
#     for j in i:
#         brain_gene_list.append(j)

# brain_gene_df = pd.DataFrame(brain_gene_list, columns=['Gene_name'])
# brain eqtl associations and annotations
all_genes = intergenics_eqtl['eQTL_genes'].tolist()
all_ann_ids = intergenics_eqtl['AnnotSV_ID'].tolist()
all_brain_genes = list()
for gene_list, annot in zip(all_genes,all_ann_ids): # Change
# for gene_list in all_genes:
    raw_list = gene_list.split('| ')
    brain_genes = list()
    for gene in raw_list:
        if 'Brain' in gene:
        # appending gene name to list
            brain_genes.append(gene.split(' ')[0])
    # brain_genes.append(brain_genes)
    all_brain_genes.append((brain_genes, annot)) # Change

brain_gene_list = list()
corresp_annot_ids = list() # Change
for i in all_brain_genes:
    # for j in i:
    for j in i[0]: # Change
        brain_gene_list.append(j)
        corresp_annot_ids.append(i[1]) # Change
brain_gene_list

brain_gene_df = pd.DataFrame({'Gene_name': brain_gene_list, 'AnnotSV_ID': corresp_annot_ids}) # Change
non_dups = list(brain_gene_df['Gene_name'].drop_duplicates().index) # Change
brain_gene_df.iloc[non_dups].reset_index(drop=True) # Change

#tpm1
brain_gene_tpm1 = brain_gene_df.merge(df_expression, how = 'left', left_on = "Gene_name", right_on ='GENE')
del brain_gene_tpm1['GENE']

#tpm2
brain_gene_tpm2 = brain_gene_tpm1.merge(df_primarymotor, how = 'left', left_on = "Gene_name", right_on ='geneSymbol') 
del brain_gene_tpm2['geneSymbol']

#tpm3
brain_gene_tpm3 = brain_gene_tpm2.merge(df_cerebral, how = 'left', left_on = "Gene_name", right_on ='Gene_Name') 
del brain_gene_tpm3['Gene_Name']

#tpm filt
tpm_filt = ((brain_gene_tpm3['GTEx_brain_max'] > 5) | (brain_gene_tpm3['brainspan_maxTPM'] > 5)) | (brain_gene_tpm3['brainPrenatal_maxTPM'] > 5)
brain_gene_tpm = brain_gene_tpm3[tpm_filt]

#impc
brain_gene_impc = brain_gene_tpm.merge(impc_final, how='left', left_on='Gene_name', right_on='Symbol')
del brain_gene_impc['Symbol']

#adding a column for the actual brain genes corresponding to each variant
# intergenics_eqtl['brain_eQTL_genes'] = all_brain_genes
intergenics_eqtl['brain_eQTL_genes'] = [genes for genes,annot in all_brain_genes] # Change

#removing sites whose brain eQTL genes don't pass TPM filt (looking at introns_eqtl['brain_eQTL_genes'] and intron_brain_gene_impc to do this manually)
intergenics_eqtl_tpm = intergenics_eqtl.copy()
blacklist = brain_gene_tpm3[tpm_filt==False]["Gene_name"].tolist()
intergenics_eqtl_tpm['brain_eQTL_genes'] = intergenics_eqtl_tpm.apply(remove_blacklist, args=(blacklist,), axis=1)
intergenics_eqtl_tpm = intergenics_eqtl_tpm[intergenics_eqtl_tpm["brain_eQTL_genes"].apply(lambda x: len(x) > 0)]

#look at these candidates and annotations for brain eqtl genes (these are kept seperate since it will be more confusing if they're merged)
# brain_gene_impc.rename(columns={'Gene_name': 'eQTL_Gene_name', 'GTEx_brain_max': 'eQTL_GTEx_brain_max', 'brainspan_maxTPM':'eQTL_brainspan_maxTPM', 'brainPrenatal_maxTPM':'eQTL_brainPrenatal_maxTPM', 'IMPC':'eQTL_IMPC', 'top_level_mp_term_name':'eQTL_top_level_mp_term_name', 'mp_term_name':'eQTL_mp_term_name'}, inplace=True)
brain_gene_impc_final=pd.merge(brain_gene_impc, ndds, how='left', left_on='Gene_name', right_on='Gene')
#candidates_ndd
brain_gene_impc_final = pd.merge(brain_gene_impc_final, sfari, how='left', left_on='Gene_name', right_on='gene-symbol')
#candidates_sfari
brain_gene_impc_final.loc[pd.isna(brain_gene_impc_final['Gene'])==False, 'ndd_tourettes'] = True
brain_gene_impc_final.loc[pd.isna(brain_gene_impc_final['gene-symbol'])==False, 'ndd_sfari'] = True
# brain_gene_impc_final = pd.merge(brain_gene_impc_final, gene_biotypes, how='left', left_on='Gene_name', right_on='hgnc_symbol')
brain_gene_impc_final = pd.merge(brain_gene_impc_final, gene_biotypes, how='left', left_on='Gene_name', right_on='Gene Name')

# #Saving results for intergenic eqtl-focused variants
intergenics_eqtl_tpm[['AnnotSV_ID', 'SV_length', 'Gene_name', 'Location', 'DDD_disease', 'OMIM_phenotype','ExAC_misZ', 'Inheritance_Patterns','GTEx_brain_max','brainspan_maxTPM', 'brainPrenatal_maxTPM', 'IMPC','top_level_mp_term_name','mp_term_name', 'AF', 'tissue_count_ME_isCausal','causal_tissues','eQTL_genes', 'brain_eQTL_genes', 'eQTL_origin', 'ACMG_class', 'psv']].to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/intergenics_eqtl-focused.csv', index=None)
intergenics_eqtl_tpm.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_intergenics_eqtl-focused.csv', index=None)
brain_gene_impc_final[['Gene_name', 'AnnotSV_ID', 'GTEx_brain_max', 'brainspan_maxTPM', 'brainPrenatal_maxTPM', 'IMPC', 'top_level_mp_term_name', 'mp_term_name', 'ndd_tourettes', 'ndd_sfari', 'Gene Type']].to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/intergenics_eqtl-focused_annotations.csv', index=None)



# Making table of counts
# Functions
def get_unique_fams_for_pheno(df, pheno):
    lst = df[pheno+"_unique_fams"].tolist()
    unique_fams = set()
    for item in lst:
        fams = item.split(',')
        for fam in fams:
            fam = fam.strip()
            if fam.isdigit():
                unique_fams.add(int(fam))
    return list(unique_fams)

def get_unique_stats(df, pheno):
    nb_uq_sites = len(df['AnnotSV_ID'].unique())
    uq_genes = df['Gene_name'].unique().tolist()
    nb_uq_genes = len(df['Gene_name'].unique())
    uq_fams = get_unique_fams_for_pheno(df, pheno)
    nb_uq_fams = len(uq_fams)
    brain_eQTL_genes = ["Not_counted"]
    if 'brain_eQTL_genes' in df:
        ex = df['brain_eQTL_genes']
        brain_eQTL_genes = list(set([item for sublist in ex.tolist() for item in sublist]))
    nb_brain_eQTL_genes = len(brain_eQTL_genes)
    
    return uq_genes, nb_uq_genes, nb_uq_sites, uq_fams, nb_uq_fams, brain_eQTL_genes, nb_brain_eQTL_genes

dataframes = ["intergenics", "intergenics_eqtl", "intergenics_eqtl_tpm"] # This line and line below are only lines that need to change across pipelines
pipeline = ["IntergenicEQTL"]*len(dataframes)*4   # 4 because there is a row for each phenotype, plus one more row for all combined
phenotypes = ["All", "ASD_only", "ASD_LI", "ASD_RI"]
final = pd.DataFrame({
    "Pipeline": pipeline, 
    "Dataframe": [item for item in dataframes for _ in range(4)],
    "Phenotypes": phenotypes*len(dataframes),
})
genes = []
nb_genes = []
nb_sites = []
fams = []
nb_fams = []
nb_eQTL_genes = []
eQTL_genes = []

for typ in dataframes:
    df = eval(typ)
    for pheno in phenotypes:
        if pheno=="All":
            subdf = df
        else:
            subdf = df[df["Inheritance_Patterns"].str.contains(pheno)]
        ug, nug, nus, uf, nuf, eg, neg = get_unique_stats(subdf, pheno)
        genes.append(ug)
        nb_genes.append(nug)
        nb_sites.append(nus)
        fams.append(uf)
        nb_fams.append(nuf)
        nb_eQTL_genes.append(neg)
        eQTL_genes.append(eg)

final["Nb_Unique_Genes"] = nb_genes
final["Nb_Unique_Sites"] = nb_sites
final["Nb_Unique_Fams"] = nb_fams
final["Unique_Genes"] = genes
final["Unique_Fams"] = fams
final["Nb_Unique_eQTL_Genes"] = nb_eQTL_genes
final["Unique_eQTL_Genes"] = eQTL_genes
final.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/intergenics_eqtl-focused_counts.csv', index=None)

