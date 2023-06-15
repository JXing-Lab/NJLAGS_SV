import pandas as pd

ea = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/exons_af-focused_counts.csv")
ia = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/introns_af-focused_counts.csv")
ie = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/introns_eqtl-focused_counts.csv")
im = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/counts/intergenics_eqtl-focused_counts.csv")

total_counts_df = pd.concat([ea,ia,ie,im])
final_pipeline_dataframes = ["exons_popAF", "introns_denovo_popAF", "introns_eqtl_tpm", "intergenics_eqtl_tpm"]
shortened_counts_df = total_counts_df[total_counts_df["Dataframe"].isin(final_pipeline_dataframes)]

nb_unique_genes_col = []
nb_unique_fams_col = []
unique_genes_col = []
unique_fams_col = []

phenotypes = ["All", "ASD_only", "ASD_LI", "ASD_RI"]
for pheno in phenotypes:
    df = shortened_counts_df[shortened_counts_df["Phenotypes"]==pheno].reset_index(drop=True)

    # li = df["Unique_Genes"].tolist()
    li = pd.concat([df["Unique_Genes"][0:2],df["Unique_eQTL_Genes"][2:4]]).tolist()
    unique_genes = list(set([y for x in li for y in eval(x)]))
    nb_unique_genes = len(unique_genes)

    li = df["Unique_Fams"].tolist()
    unique_fams = list(set([y for x in li for y in eval(x)]))
    nb_unique_fams = len(unique_fams)

    nb_unique_genes_col.append(nb_unique_genes)
    nb_unique_fams_col.append(nb_unique_fams)
    unique_genes_col.append(unique_genes)
    unique_fams_col.append(unique_fams)

counts_by_pheno = pd.DataFrame({
    "Phenotype": phenotypes, 
    "Nb_Unique_Genes": nb_unique_genes_col, 
    "Nb_Unique_Fams":nb_unique_fams_col,
    "Unique_Genes":unique_genes_col,
    "Unique_Fams":unique_fams_col,
})

counts_by_pheno.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/counts_by_phenotype.csv", index=None)
