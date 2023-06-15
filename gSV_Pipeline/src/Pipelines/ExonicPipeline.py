import pandas as pd

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

# Filtering for exons
# exon_filt = (pd.isna(candidates_nonbenign['Location']) != True) & ((candidates_nonbenign['Location'].str.contains('exon')) | (candidates_nonbenign['Location'].str.contains('tx')) | (candidates_nonbenign['0'] != candidates_nonbenign['1']))
exon_filt = (pd.isna(candidates_nonbenign['Gene_name'])==False) & (pd.isna(candidates_nonbenign['Location']) != True) & ((candidates_nonbenign['Location'].str.contains('exon')) | (candidates_nonbenign['Location'].str.contains('tx')) | (candidates_nonbenign['0'] != candidates_nonbenign['1']))
exons = candidates_nonbenign[exon_filt]

# Filtering for max TPM > 5
tpm_filt = ((exons['GTEx_brain_max'] > 5) | (exons['brainspan_maxTPM'] > 5)) | (exons['brainPrenatal_maxTPM'] > 5)
exons_tpm = exons[tpm_filt]

# Filtering on sampAF
exons_sampAF = exons_tpm[exons_tpm['AF']<0.05]

# Filtering on popAF
exons_popAF = exons_sampAF[exons_sampAF['popAF']<0.05]

# Saving condensed results
exons_popAF[['AnnotSV_ID', 'SV_length', 'Gene_name', 'Location', 'DDD_disease', 'OMIM_phenotype','ExAC_misZ', 'GnomAD_pLI', 'Inheritance_Patterns','GTEx_brain_max','brainspan_maxTPM', 'brainPrenatal_maxTPM', 'IMPC','top_level_mp_term_name','mp_term_name', 'AF', 'tissue_count_ME_isCausal','causal_tissues','eQTL_genes', 'ACMG_class', 'ndd_tourettes', 'ndd_sfari', 'Gene Type', 'psv']].to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/exons_af-focused.csv', index=None)
# Saving full results
exons_popAF.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/results/full_exons_af-focused.csv', index=None)


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
    # uq_gene_types = df['Gene Type']
    
    # uq_families = len(df[])
    # stats = f'Unique Sites: {uq_sites}\nUnique Genes: {uq_genes}\n'
    return uq_genes, nb_uq_genes, nb_uq_sites, uq_fams, nb_uq_fams

dataframes = ["exons", "exons_tpm", "exons_sampAF", "exons_popAF"] # Only line that needs to change across pipelines
pipeline = ["ExonicAF"]*len(dataframes)*4   # 4 because there is a row for each phenotype, plus one more row for all combined
phenotypes = ["All", "ASD_only", "ASD_LI", "ASD_RI"]
final = pd.DataFrame({
    "Pipeline": pipeline, 
    "Dataframe": [item for item in dataframes for _ in range(4)],
    "Phenotypes": phenotypes*4,
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
final.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/exons_af-focused_counts.csv', index=None)