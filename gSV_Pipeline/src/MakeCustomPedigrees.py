import pandas as pd


# Define file paths
pedigree_file = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/entire_pedigree.ped"
pedigree_loc = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI"

# Read pedigree data
pedigree = pd.read_csv(pedigree_file, sep="\t")

# Create new column 'ASD_only' with values from 'phenotype'
pedigree['ASD_only'] = pedigree['phenotype']

# Define a function to filter phenotype values
def filt(x, y):
    if x == 0 or y == 0:
        return 0
    elif x == 2 or y == 2:
        return 2
    else:
        return 1

# Apply filter function to 'phenotype' and 'LI_pheno' columns to create new 'ASD_LI' column
pedigree['ASD_LI'] = [filt(x, y) for x, y in zip(pedigree['phenotype'], pedigree['LI_pheno'])]

# Apply filter function to 'phenotype' and 'RI_pheno' columns to create new 'ASD_RI' column
pedigree['ASD_RI'] = [filt(x, y) for x, y in zip(pedigree['phenotype'], pedigree['RI_pheno'])]

# Define list of pedigree types
pedigrees = ["ASD_LI", "ASD_RI", "ASD_only"]

# Loop through pedigree types and save corresponding data to file
for ped in pedigrees:
    ped_file_loc = f"{pedigree_loc}/{ped}.ped"
    df = pedigree[["family_id", "name", "paternal_id", "maternal_id", "sex", ped]]
    df = df.rename(columns={ped: "phenotype"})
    df.to_csv(ped_file_loc, index=None, sep="\t")
