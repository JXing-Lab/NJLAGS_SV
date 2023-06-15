
import pandas as pd

# Read csv
# candidates_biotypes = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/biotypes_candidates.csv")
candidates_svanna = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/svanna_candidates.csv")

# Inh Pattern Filt
candidates_inh = candidates_svanna[candidates_svanna['Inheritance_Patterns'].isna()==False]

# Save csv
candidates_inh_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/inh_candidates.csv"
candidates_inh.to_csv(candidates_inh_path, index=False)