
import pandas as pd

# Read csv
candidates_inh = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/inh_candidates.csv")

# Benign Filt
benign_filt = (pd.isna(candidates_inh['B_gain_source']) == True) & (pd.isna(candidates_inh['B_loss_source']) == True) & (pd.isna(candidates_inh['B_ins_source']) == True) & (pd.isna(candidates_inh['B_inv_source']) == True)
candidates_nonbenign = candidates_inh[benign_filt]

# Save csv
candidates_nonbenign_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/nonbenign_candidates.csv"
candidates_nonbenign.to_csv(candidates_nonbenign_path, index=False)