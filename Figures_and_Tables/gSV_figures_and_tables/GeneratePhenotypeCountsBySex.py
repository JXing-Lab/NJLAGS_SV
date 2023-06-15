import subprocess
import pandas as pd

# Define the commands
commands = [
    # ASD_only
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_only.db | head -n 1 | awk \'{print "People that have ASD: "}{print NF}\'',
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2 and sex==1) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_only.db | head -n 1 | awk \'{print "Males: "}{print NF}\'',
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2 and sex==2) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_only.db | head -n 1 | awk \'{print "Females: "}{print NF}\'',
    # ASD_LI
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_LI.db | head -n 1 | awk \'{print "People that have ASD_LI: "}{print NF}\'',
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2 and sex==1) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_LI.db | head -n 1 | awk \'{print "Males: "}{print NF}\'',
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2 and sex==2) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_LI.db | head -n 1 | awk \'{print "Females: "}{print NF}\'',
    # ASD_RI
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_RI.db | head -n 1 | awk \'{print "People that have ASD_RI: "}{print NF}\'',
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2 and sex==1) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_RI.db | head -n 1 | awk \'{print "Males: "}{print NF}\'',
    '/lab01/Tools/gemini/bin/gemini query --header -q "SELECT (gts).(phenotype==2 and sex==2) FROM variants" /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/modified_ASD_RI.db | head -n 1 | awk \'{print "Females: "}{print NF}\'',
]

# Define variables to store the numeric outputs
people_with_pheno = []
males = []
females = []

# Execute the commands and capture the outputs
for i, command in enumerate(commands):
    output = subprocess.check_output(command, shell=True).decode('utf-8').strip()
    numeric_output = int(output.split('\n')[-1])
    
    if i%3 == 0:
        people_with_pheno.append(numeric_output)
    elif i%3 == 1:
        males.append(numeric_output)
    elif i%3 == 2:
        females.append(numeric_output)

# Print the stored values
print(f"People: {people_with_pheno}")
print(f"Males: {males}")
print(f"Females: {females}")
phenotypes_by_sex = pd.DataFrame({"Phenotypes": ["ASD_only", "ASD_LI", "ASD_RI"] ,"Males": males, "Females": females, "Total": people_with_pheno})

# Save table
save_dir = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/"
phenotypes_by_sex.to_csv(f"{save_dir}/phenotype_counts_by_sex.csv", index=None)