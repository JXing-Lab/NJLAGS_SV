import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chisquare, mannwhitneyu

# Extracting counts from vcf
people = []
var_counts = np.array([0]*272)
with open('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/modified_SURVIVOR_merged_callset.vcf') as f:
    for line in f.readlines():
        # Skipping useless header
        fields = line.split('\t')

        if line.startswith("##"):
            continue
        # Getting people names
        elif line.startswith("#"):
            people = fields[9:]
        # No longer in the header
        else:
            add_to_var_counts = np.array([1 if '1' in x[:3] else 0 for x in fields[9:]])
            var_counts+=add_to_var_counts

people_var_counts = {k.strip(): v for k,v in zip(people,var_counts)}

# Extracting which individuals are affected and not
asd = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/ASD_only.ped', sep='\t')
asd_aff_names = asd[asd['phenotype']==2]['name'].tolist()
asd_aff_counts = pd.DataFrame([v for k,v in people_var_counts.items() if k in asd_aff_names]).rename(columns={0: 'Count'})
asd_aff_counts_minus_knome = pd.DataFrame([v for k,v in people_var_counts.items() if (k in asd_aff_names) and (not k.startswith('LP'))]).rename(columns={0: 'Count'})
asd_unaff_counts = pd.DataFrame([v for k,v in people_var_counts.items() if k not in asd_aff_names]).rename(columns={0: 'Count'})

'''
GRAPH + STATS EXCLUDING KNOME BATCH (13 affecteds)
'''

# Read in the data
dfp = asd_aff_counts_minus_knome
dfu = asd_unaff_counts

proband_median = np.median(dfp['Count'])
unaff_median = np.median(dfu['Count'])
print(proband_median)
print(unaff_median)

proband_mean = np.mean(dfp['Count'])
unaff_mean = np.mean(dfu['Count'])
print(proband_mean)
print(unaff_mean)

# Prepare data for graph
unaff = pd.DataFrame({'group': 'ASD_only_unaff', 'value': dfu['Count']})
proband = pd.DataFrame({'group': 'ASD_only_proband', 'value': dfp['Count']})
plot_data = pd.concat([unaff, proband], ignore_index=True)

# Plot data
sns.boxplot(data=plot_data, x='group', y='value', palette='viridis')
sns.stripplot(data=plot_data, x='group', y='value', color='black', size=3, alpha=0.9)
plt.savefig('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/figures/no_knome_burden_analysis.png')
plt.savefig('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/figures/no_knome_burden_analysis.pdf',format='pdf')
plt.show()

# Run chi-squared test
m = pd.crosstab(dfp['Count'], dfu['Count'])
chi2, p_value = chisquare(m.values.flatten())
print(f"Chi-squared test statistic: {chi2}")
print(f"Chi-squared p-value: {p_value}")

# Run Mann-Whitney-Wilcoxon test
unaff_counts = dfu['Count'].to_numpy()
proband_counts = dfp['Count'].to_numpy()
_, p_value = mannwhitneyu(unaff_counts, proband_counts, alternative='two-sided')
print(f"Mann-Whitney-Wilcoxon p-value: {p_value}")

# Saving stats to a log file
with open('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/figures/no_knome_burden_analysis_stats.txt', 'w') as f:
    f.write("MEI+SV Burden Analysis\n\n\n")
    f.write(f"Proband median: {proband_median}\n")
    f.write(f"Unaffected median: {unaff_median}\n")
    f.write("------------------------------\n")
    f.write(f"Proband mean: {proband_mean}\n")
    f.write(f"Unaffected mean: {unaff_mean}\n")
    f.write("------------------------------\n")
    f.write(f"Chi-squared test statistic: {chi2}\n")
    f.write(f"Chi-squared p-value: {p_value}\n")
    f.write(f"Mann-Whitney-Wilcoxon p-value: {p_value}")