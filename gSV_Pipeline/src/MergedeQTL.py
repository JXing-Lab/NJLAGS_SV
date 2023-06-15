import pandas as pd
from biomart import BiomartServer
import io


# # Reading csv generated from getSVeQTLs
# eqtl_df = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/xiaolong_eqtl_table_s2/Supplemental_Dataframe_S2.csv")
# # Applying above filters
# SV_filt=eqtl_df["lead_variant_type"]=="SV"
# Tissue_filt=eqtl_df["tissue"].str.contains("Brain_")
# Caviar_filt=eqtl_df["lead_sv_caviar_prob"] > 0.1
# Caviar_rank_filt=eqtl_df["lead_sv_cariar_rank"] < 6
# eqtl_df = eqtl_df[SV_filt & Tissue_filt & Caviar_filt & Caviar_rank_filt]
# # Getting gene lists
# gene_list = eqtl_df['gene'].unique().tolist()
# no_version_gene_list = [x.split(".")[0] for x in gene_list]

# ## Biomart package, for exchanging ensembleIDs into regular. Not able to handle versions (the decimal in the gene field) so will remove them
# # Connect to the Biomart server - GRCh37 build.
# server = BiomartServer("http://grch37.ensembl.org/biomart")
# # Set the dataset and attributes
# dataset = server.datasets["hsapiens_gene_ensembl"]
# attributes = ["ensembl_gene_id", "external_gene_name"]
# # Define the list of Ensembl gene IDs you want to convert
# ensembl_ids = no_version_gene_list
# # Query the Biomart server to get the gene names
# response = dataset.search({
#     'filters': {"ensembl_gene_id": ensembl_ids},
#     'attributes': attributes
# })
# # Convert the bytes response to a string
# content_str = response.content.decode()
# # Convert the string to a DataFrame
# df = pd.read_csv(io.StringIO(content_str), sep='\t', names=attributes)
# # Print the DataFrame
# df = df.rename(columns={"ensembl_gene_id": "no_version_gene"})

# ### Merging with SV eqtls
# eqtl_df["no_version_gene"] = [x.split(".")[0] for x in eqtl_df["gene"]]
# eqtl_df = eqtl_df.merge(df, on="no_version_gene")
# # Splitting lead_sv_coord into components
# raw_split = [x.split(":") for x in eqtl_df["lead_sv_coord"].tolist()]
# chroms = [x[0] for x in raw_split]
# positions = [x[1].split("-") for x in raw_split]
# starts = [x[0] for x in positions]
# ends = [x[1] for x in positions]
# eqtl_df["chrom"]=chroms
# eqtl_df["start"]=starts
# eqtl_df["end"]=ends
# short_eqtl_df = eqtl_df[["chrom", "start", "end", "external_gene_name", "tissue"]]

# ### Formatting sv_eqtls before merging with mei_eqtls
# ## Combine rows with same chrom and start and end. Gene and tissue union set
# # define aggregation function
# def merge_sets(x):
#     return pd.Series({'external_gene_name': set(x['external_gene_name']),
#                       'tissue': set(x['tissue'])})
# # group the dataframe by chrom, start, and end columns and apply the aggregation function
# merged_eqtls = short_eqtl_df.groupby(['chrom', 'start', 'end']).apply(merge_sets).reset_index()
# # rename columns
# merged_eqtls = merged_eqtls.rename(columns={'external_gene_name': 'gene', 'tissue': 'tissues'})
# # more formatting for merging
# merged_eqtls["isCausal_for_brain_tissue"]=True
# formatted_merged_eqtls = merged_eqtls.rename(columns={"start":"snp_pos", "tissues": "Tissues", "gene":"gene_name", "chrom":"snp_chr"})
# formatted_merged_eqtls["tissue_count_ME_isCausal"] = [len(x) for x in formatted_merged_eqtls["Tissues"].tolist()]
# # changing gene column from sets to strings containing "(Brain) | " after each gene name
# genes = formatted_merged_eqtls['gene_name'].apply(lambda x: str(x)[2:-2]).tolist()
# genes_annot = [f"{x} (Brain) |" if "," not in x else (" (Brain) | ".join(x.split("', '")) + " (Brain)") for x in genes ]
# formatted_merged_eqtls["gene_name"] = genes_annot
# sv_eqtls = formatted_merged_eqtls[["snp_pos", "tissue_count_ME_isCausal","Tissues", "gene_name", "snp_chr"]]
# # changing tissues column from sets to lists
# sv_eqtls["Tissues"] = [list(x) for x in sv_eqtls["Tissues"]]
# sv_eqtls['eqtl_origin'] = 'SV'

# # Save csv of SV eqtls
# sv_eqtl_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/eqtl_sv.csv"
# sv_eqtls.to_csv(sv_eqtl_path, index=False)




# ''' 
# GTEx MEI eQTLs
# '''



# #importing caviar data
# mei_eqtls = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211113_eQTL/20190805GTExv7_q.1.normalizedExpression.CAVIAR.CausalPair.added.csv', sep='\t')


# # Adding tissue data
# tissueColumn = list()
# for i in range(0, 4190):
#     not_nan_values = mei_eqtls.iloc[i, :48].notna().tolist()
#     not_nan_values
#     columns_wout_nan = mei_eqtls.iloc[:, :48].columns[not_nan_values].tolist()
#     tissueColumn.append(columns_wout_nan)

# for tissueList in tissueColumn:
#     for tissue in range(len(tissueList)):
#         tissueList[tissue] = tissueList[tissue].split('|')[0]

# mei_eqtls['Tissues'] = tissueColumn


# #if a variant has 'brain' in Tissues, put 'Brain' next to the gene name in parenthesis
# brain_filter = mei_eqtls['Tissues'].astype(str).str.contains('Brain')
# geneAffectsBrain = brain_filter.tolist()
# for i in range(0, 4190):
#     if (geneAffectsBrain[i]==True):
#         mei_eqtls['gene_name'][i] = mei_eqtls['gene_name'][i] + ' (Brain) '

# mei_eqtls['gene_name'] = mei_eqtls['gene_name'].astype(str) + ' | '
# mei_eqtls = mei_eqtls[['tissue_count_ME_isCausal', 'Tissues', 'snps', 'snp_pos', 'snp_chr', 'gene_name']]
# mei_eqtls = mei_eqtls.groupby('snp_pos').aggregate({'tissue_count_ME_isCausal': 'sum', 'Tissues': 'sum', 'gene_name': 'sum','snps': 'first', 'snp_chr': 'first'}, inplace=True).reset_index()
# mei_eqtls['eqtl_origin'] = 'MEI'

# ### Adding in the sv eqtls (concatting sv eqtls with mei_eqtls)
# # sv_eqtls = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/eqtl_sv.csv")
# # combined_eqtls_df = pd.concat([mei_eqtls, sv_eqtls])
# # combined_eqtls_df["snp_pos"] = combined_eqtls_df["snp_pos"].astype('int64')
# # mei_eqtls = combined_eqtls_df # Concatenation complete.


''' 


GTEx MEI eQTLs
273 unique pos:chrom tuples MEI


'''

#importing caviar data
caviar_df = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211113_eQTL/20190805GTExv7_q.1.normalizedExpression.CAVIAR.CausalPair.added.csv', sep='\t')


# Adding tissue data
tissueColumn = list()
for i in range(0, 4190):
    not_nan_values = caviar_df.iloc[i, :48].notna().tolist()
    not_nan_values
    columns_wout_nan = caviar_df.iloc[:, :48].columns[not_nan_values].tolist()
    tissueColumn.append(columns_wout_nan)

for tissueList in tissueColumn:
    for tissue in range(len(tissueList)):
        tissueList[tissue] = tissueList[tissue].split('|')[0]

caviar_df['Tissues'] = tissueColumn


#if a variant has 'brain' in Tissues, put 'Brain' next to the gene name in parenthesis
brain_filter = caviar_df['Tissues'].astype(str).str.contains('Brain')
geneAffectsBrain = brain_filter.tolist()
for i in range(0, 4190):
    if (geneAffectsBrain[i]==True):
        caviar_df['gene_name'][i] = caviar_df['gene_name'][i] + ' (Brain) '

caviar_df['gene_name'] = caviar_df['gene_name'].astype(str) + ' | '
caviar_df_short = caviar_df[['tissue_count_ME_isCausal', 'Tissues', 'snp_pos', 'snp_chr', 'gene_name']]
caviar_df_short = caviar_df_short.groupby('snp_pos').aggregate({'tissue_count_ME_isCausal': 'sum', 'Tissues': 'sum', 'gene_name': 'sum', 'snp_chr': 'first'}, inplace=True).reset_index()
caviar_df_short['eqtl_origin'] = 'MEI'

BrainTissue = [True if ('Brain' in str(x)) else False for x in caviar_df_short["Tissues"].tolist()]
mei_eqtls = caviar_df_short[BrainTissue]
mei_eqtls



''' 


GTEx SV eQTLs
377 unique pos:chrom tuples MEI

SV selection criteria:
-'lead_variant_type' == 'SV'
-'tissue' == 'Brain_*'
-'lead_sv_caviar_prob' > 0.1
-'lead_sv_cariar_rank' < 6


'''
# Reading csv generated from getSVeQTLs
eqtl_df = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/xiaolong_eqtl_table_s2/Supplemental_Dataframe_S2.csv")
# Applying above filters
SV_filt=eqtl_df["lead_variant_type"]=="SV"
Tissue_filt=eqtl_df["tissue"].str.contains("Brain_")
Caviar_filt=eqtl_df["lead_sv_caviar_prob"] > 0.1
Caviar_rank_filt=eqtl_df["lead_sv_cariar_rank"] < 6
eqtl_df = eqtl_df[SV_filt & Tissue_filt & Caviar_filt & Caviar_rank_filt]
# Getting gene lists
gene_list = eqtl_df['gene'].unique().tolist()
no_version_gene_list = [x.split(".")[0] for x in gene_list]

## Biomart package, for exchanging ensembleIDs into regular. Not able to handle versions (the decimal in the gene field) so will remove them
# Connect to the Biomart server - GRCh37 build.
server = BiomartServer("http://grch37.ensembl.org/biomart")
# Set the dataset and attributes
dataset = server.datasets["hsapiens_gene_ensembl"]
attributes = ["ensembl_gene_id", "external_gene_name"]
# Define the list of Ensembl gene IDs you want to convert
ensembl_ids = no_version_gene_list
# Query the Biomart server to get the gene names
response = dataset.search({
    'filters': {"ensembl_gene_id": ensembl_ids},
    'attributes': attributes
})
# Convert the bytes response to a string
content_str = response.content.decode()
# Convert the string to a DataFrame
df = pd.read_csv(io.StringIO(content_str), sep='\t', names=attributes)
# Print the DataFrame
df = df.rename(columns={"ensembl_gene_id": "no_version_gene"})

### Merging with SV eqtls
eqtl_df["no_version_gene"] = [x.split(".")[0] for x in eqtl_df["gene"]]
eqtl_df = eqtl_df.merge(df, on="no_version_gene")
# Splitting lead_sv_coord into components
raw_split = [x.split(":") for x in eqtl_df["lead_sv_coord"].tolist()]
chroms = [x[0] for x in raw_split]
positions = [x[1].split("-") for x in raw_split]
starts = [x[0] for x in positions]
ends = [x[1] for x in positions]
eqtl_df["chrom"]=chroms
eqtl_df["start"]=starts
eqtl_df["end"]=ends
short_eqtl_df = eqtl_df[["chrom", "start", "end", "external_gene_name", "tissue"]]

### Formatting sv_eqtls before merging with mei_eqtls
## Combine rows with same chrom and start and end. Gene and tissue union set
# define aggregation function
def merge_sets(x):
    return pd.Series({'external_gene_name': set(x['external_gene_name']),
                      'tissue': set(x['tissue'])})
# group the dataframe by chrom, start, and end columns and apply the aggregation function
merged_eqtls = short_eqtl_df.groupby(['chrom', 'start', 'end']).apply(merge_sets).reset_index()
# rename columns
merged_eqtls = merged_eqtls.rename(columns={'external_gene_name': 'gene', 'tissue': 'tissues'})
# more formatting for merging
merged_eqtls["isCausal_for_brain_tissue"]=True
formatted_merged_eqtls = merged_eqtls.rename(columns={"start":"snp_pos", "tissues": "Tissues", "gene":"gene_name", "chrom":"snp_chr"})
formatted_merged_eqtls["tissue_count_ME_isCausal"] = [len(x) for x in formatted_merged_eqtls["Tissues"].tolist()]
# changing gene column from sets to strings containing "(Brain) | " after each gene name
genes = formatted_merged_eqtls['gene_name'].apply(lambda x: str(x)[2:-2]).tolist()
genes_annot = [f"{x} (Brain) |" if "," not in x else (" (Brain) | ".join(x.split("', '")) + " (Brain)") for x in genes ]
formatted_merged_eqtls["gene_name"] = genes_annot
sv_eqtls = formatted_merged_eqtls[["snp_pos", "tissue_count_ME_isCausal","Tissues", "gene_name", "snp_chr"]]
# changing tissues column from sets to lists
sv_eqtls["Tissues"] = [list(x) for x in sv_eqtls["Tissues"]]
sv_eqtls=sv_eqtls.groupby(['snp_pos','snp_chr']).aggregate({'tissue_count_ME_isCausal': 'sum', 'Tissues': 'sum', 'gene_name': 'sum'}, inplace=True).reset_index()
sv_eqtls['eqtl_origin'] = 'SV'
sv_eqtls



''' 


GTEx NEI+SV eQTLs
647 unique pos:chrom tuples MEI

For now, I will not update the gene_name, Tissues, and tissue_count_ME_isCausal columns because I'm not sure how to do it quickly


'''



combined_eqtls_df = pd.concat([mei_eqtls, sv_eqtls])
combined_eqtls_df['snp_pos'] = combined_eqtls_df['snp_pos'].astype(int)
combined_eqtls_df['snp_chr'] = combined_eqtls_df['snp_chr'].astype(str)
combined_eqtls_df=combined_eqtls_df.groupby(['snp_pos','snp_chr']).aggregate({'tissue_count_ME_isCausal': 'sum', 'Tissues': 'sum', 'gene_name': 'sum', 'eqtl_origin': 'sum'}, inplace=True).reset_index()

combined_eqtls_df
combined_eqtls_df.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/eqtl_combined.csv", index=None)
