import pandas as pd

# Filtering out these rare variants that have denovo inheritance patterns for multiple people (likley artifcats)
def get_inheritance(sv_id, data):
    denovo_col_famid = [col for col in data if (col.endswith('denovo') & col.startswith('family_id'))]
    columns = data[denovo_col_famid].columns.tolist()
    for column in columns:
        in_multiple_families=False
        # if its a list (so its non nan)...
        if isinstance((data[data['AnnotSV_ID']==sv_id][column].values[0]), str):
            # and multiple have families have denovo...
            if (len(eval(data[data['AnnotSV_ID']==sv_id][column].values[0])) > 1):
                in_multiple_families=True
                break
    return in_multiple_families  

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


# Filtering for introns
# tron_filt2 = (candidates_nonbenign['0'].str.contains('intron')) & (candidates_nonbenign['1'].str.contains('intron'))
tron_filt2 = (pd.isna(candidates_nonbenign['Gene_name'])==False) & (candidates_nonbenign['Location'].str.contains('intron')) & (candidates_nonbenign['0'] == candidates_nonbenign['1'])
introns = candidates_nonbenign[tron_filt2]

# tpm filt
tpm_filt = ((introns['GTEx_brain_max'] > 5) | (introns['brainspan_maxTPM'] > 5)) | (introns['brainPrenatal_maxTPM'] > 5)
introns_tpm = introns[tpm_filt]

# filt by sampAF
introns_sampAF = introns_tpm[introns_tpm['AF']<0.05]
        
introns_sampAF_2=introns_sampAF.copy()
del introns_sampAF_2['0']
del introns_sampAF_2['1']
introns_sampAF_2

# filt by denovo
sv_list = introns_sampAF_2['AnnotSV_ID'].tolist()
sv_list
x = list()
for sv in sv_list:
    in_multiple_families = get_inheritance(sv, introns_sampAF_2)
    x.append(in_multiple_families)
introns_sampAF_2['denovo_in_multiple_families'] = x
introns_denovo = introns_sampAF_2[introns_sampAF_2['denovo_in_multiple_families']==False]

# Filtering by NDD genes
introns_denovo_ndd = introns_denovo[(introns_denovo['ndd_tourettes']==True) | (introns_denovo['ndd_sfari']==True)]

# Filtering by popAF
introns_denovo_ndd['popAF'] = introns_denovo_ndd['popAF'].fillna(0)
introns_denovo_popAF = introns_denovo_ndd[introns_denovo_ndd['popAF']<0.05]
introns_denovo_popAF = introns_denovo_popAF.drop_duplicates(subset=['AnnotSV_ID','SV_length','Gene_name','Location','Inheritance_Patterns'])

# Saving results for Intronic AF-focused candidates
introns_denovo_popAF[['AnnotSV_ID', 'SV_length', 'Gene_name', 'Location', 'DDD_disease', 'OMIM_phenotype','ExAC_misZ', 'GnomAD_pLI', 'Inheritance_Patterns','GTEx_brain_max','brainspan_maxTPM', 'brainPrenatal_maxTPM', 'IMPC','top_level_mp_term_name','mp_term_name', 'AF', 'tissue_count_ME_isCausal','causal_tissues','eQTL_genes', 'ACMG_class', 'ndd_tourettes', 'ndd_sfari', 'Gene Type', 'psv']].to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/introns_af-focused.csv', index=None)
introns_denovo_popAF.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_introns_af-focused.csv', index=None)


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
    
    # uq_families = len(df[])
    # stats = f'Unique Sites: {uq_sites}\nUnique Genes: {uq_genes}\n'
    return uq_genes, nb_uq_genes, nb_uq_sites, uq_fams, nb_uq_fams

dataframes = ["introns", "introns_tpm", "introns_sampAF", "introns_denovo", "introns_denovo_ndd", "introns_denovo_popAF"] # This line and line below are only lines that need to change across pipelines
pipeline = ["IntronicAF"]*len(dataframes)*4   # 4 because there is a row for each phenotype, plus one more row for all combined
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
for typ in dataframes:
    df = eval(typ)
    for pheno in phenotypes:
        if pheno=="All":
            subdf = df
        else:
            subdf = df[df["Inheritance_Patterns"].str.contains(pheno)]
        ug, nug, nus, uf, nuf = get_unique_stats(subdf, pheno)
        genes.append(ug)
        nb_genes.append(nug)
        nb_sites.append(nus)
        fams.append(uf)
        nb_fams.append(nuf)

final["Nb_Unique_Genes"] = nb_genes
final["Nb_Unique_Sites"] = nb_sites
final["Nb_Unique_Fams"] = nb_fams
final["Unique_Genes"] = genes
final["Unique_Fams"] = fams
final.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/introns_af-focused_counts.csv', index=None)