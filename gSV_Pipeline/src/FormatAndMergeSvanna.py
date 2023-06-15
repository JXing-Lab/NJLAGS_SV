import pandas as pd

svanna_dir="/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/SvAnna"
PSV_df=pd.read_csv(f"{svanna_dir}/Output/SvAnna_ready.SVANNA.csv")
PSV_df = PSV_df.astype({"contig":str, "start":int, "end":int, "vtype":str})
PSV_df.loc[PSV_df['vtype'] == 'SYMBOLIC', 'vtype'] = 'IDP'
# PSV_df

concat_df=pd.read_csv(f"{svanna_dir}/lifted_concat.bed", sep="\t")
concat_df.columns=["contig", "start", "end", "vtype", "concatted_info"]
concat_df = concat_df.astype({"contig":str, "start":int, "end":int, "vtype":str})
concat_df["contig"]=concat_df["contig"].str.replace("chr", "")

merged_PSV_concat_df = pd.merge(PSV_df, concat_df, on=["contig", "start", "end", "vtype"])
merged_PSV_concat_df.drop(columns=["contig", "start", "end", "id", "vtype", "failed_filters"], inplace=True)
# Split the "values" column into separate columns
merged_PSV_concat_df[['SV_chrom', 'SV_start', 'SV_end', 'SV_type']] = merged_PSV_concat_df['concatted_info'].str.split(';', expand=True)
merged_PSV_concat_df.drop(columns=["concatted_info"], inplace=True)
merged_PSV_concat_df["SV_chrom"]=merged_PSV_concat_df["SV_chrom"].str.replace("chr", "")
merged_PSV_concat_df = merged_PSV_concat_df.astype({"SV_chrom":str, "SV_start":int, "SV_end":int, "SV_type":str})

# Read csv
candidates_biotypes = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/biotypes_candidates.csv")

# Add psv column to candidates
candidates_biotypes = candidates_biotypes.astype({"SV_chrom":str, "SV_start":int, "SV_end":int, "SV_type":str})
candidates_svanna = pd.merge(candidates_biotypes, merged_PSV_concat_df, on=['SV_chrom', 'SV_start', 'SV_end', 'SV_type'], how='outer')
candidates_svanna

# Save csv
candidates_svanna_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/svanna_candidates.csv"
candidates_svanna.to_csv(candidates_svanna_path, index=False)