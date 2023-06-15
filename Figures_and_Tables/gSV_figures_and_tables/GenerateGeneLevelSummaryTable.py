from Bio import Entrez
import mygene
import pandas as pd
import numpy as np

# Defining functions for database querying
def get_gene_info(gene_symbol):
    handle = Entrez.esearch(db='gene', term=f"{gene_symbol}[Gene Name]")
    record = Entrez.read(handle)
    gene_id = record['IdList'][0]

    handle = Entrez.esummary(db='gene', id=gene_id)
    record = Entrez.read(handle)

    # Extract desired information from the record
    name = record['DocumentSummarySet']['DocumentSummary'][0]['Description']
    synonyms = record['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases']

    return name, synonyms

def get_ensembl_gene_id(gene_symbol):
    mg = mygene.MyGeneInfo()
    results = mg.query(gene_symbol, species='human', fields='ensembl.gene', assembly='GRCh37')

    if 'hits' in results and len(results['hits']) > 0:
        gene_info = results['hits'][0]
        if 'ensembl' in gene_info and 'gene' in gene_info['ensembl']:
            ensembl_gene_id = gene_info['ensembl']['gene']
            return ensembl_gene_id
    return None

def get_mgi_id(gene_symbol):
    mouse_symbol = gene_symbol.lower().capitalize()
    return df_homolog[(df_homolog['Symbol']==mouse_symbol) & (df_homolog['Common Organism Name']=='mouse, laboratory')]['Mouse MGI ID'].item()

# Start of script
if __name__=="__main__":
    # Loading previous gene-level summaries
    gene_level_summ_MEISV = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/gene_level_summary_table_MEISV.csv")
    gene_level_summ_CNV = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/gene_level_summary_table_CNV.csv")

    # Creating combined summary table
    gene_level_summ = pd.concat([gene_level_summ_MEISV,gene_level_summ_CNV]).drop_duplicates()

    # Handling unions with an aggregation dictionary
    aggregations = {
        'Name': 'first',
        'Ensembl': 'first',
        'MGD': 'first',
        'Gene Type': 'first',
        'Synonym': 'first',
        'GnomAD_pLI': 'mean',
        'IMPC': 'first',
        'GTEx_brain_max': 'first',
        'brainPrenatal_maxTPM': 'first',
        'brainspan_maxTPM': 'first',
        'ndd_tourettes': 'first',
        'mp_term_name': 'first',
        'ndd_sfari': 'first',
        'DDD_disease': 'first', # Works since first takes the first non-null value
        'OMIM_phenotype': 'first', # Change
        'ExAC_misZ': 'first',
        'ASD_only': 'any',
        'ASD_LI': 'any',
        'ASD_RI': 'any',
        'Pipeline': lambda x: ";".join(x)
    }
    final_gene_level = gene_level_summ.groupby("Gene_name").agg(aggregations).reset_index()

    # Set your email address (required by NCBI)
    Entrez.email = 'your_email@example.com'

    # Query relevant servers/databases
    df_homolog = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20210808_IMPC/HOM_AllOrganism.rpt', sep='\t')
    Name = []
    Ensembl = []
    MGD = []
    Synonym = []
    gene_symbols = final_gene_level['Gene_name'].tolist()
    for gene_symbol in gene_symbols:
        name, synonyms = get_gene_info(gene_symbol)
        ensembl_gene_id = get_ensembl_gene_id(gene_symbol)
        try:
            mgi_gene_id = get_mgi_id(gene_symbol)
        except:
            mgi_gene_id = None

        Name.append(name)
        Ensembl.append(ensembl_gene_id)
        MGD.append(mgi_gene_id)
        Synonym.append(synonyms)

    # Adding the desired columns
    final_gene_level['Name'] = Name
    final_gene_level['Ensembl'] = Ensembl
    final_gene_level['MGD'] = MGD
    final_gene_level['Synonym'] = Synonym

    # PPI network flags taken from: /lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609genes_all.xlsx, and converted into a .csv file
    xiaolong_gene_flags = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/rohan_data/20230609genes_all.csv')
    xiaolong_gene_flags['Known_gene'] = xiaolong_gene_flags.apply(lambda row: ';'.join(row[['ADHD_evidence', 'other_Neuro', 'SFARI']][row[['ADHD_evidence', 'other_Neuro', 'SFARI']] == 1].index), axis=1)
    xiaolong_gene_flags['Known_gene'] = xiaolong_gene_flags['Known_gene'].replace('', np.nan)

    # # Display the updated dataframe
    # print(df[['Gene_name', 'Known_gene']])

    final_final_gene_level = final_gene_level.merge(xiaolong_gene_flags, left_on='Gene_name', right_on='gene_symbol', how='left')
    final_final_gene_level=final_final_gene_level.drop(columns=['target', 'ADHD_evidence', 'other_Neuro', 'SFARI', 'ndd_sfari', 'ndd_tourettes'])
    final_final_gene_level.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/gene_level_summary_table.csv")

