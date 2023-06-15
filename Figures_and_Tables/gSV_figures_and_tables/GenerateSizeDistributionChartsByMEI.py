import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd

file = 'MEI_merged.vcf'

def take_line_get_svlen(line):
    info = line.split("\t")[7]

    # Extract SVLEN
    svlen_match = re.search(r"SVLEN=(-?\d+)", info)
    if svlen_match:
        svlen = int(svlen_match.group(1))
    else:
        svlen = None

    # Extract SVTYPE_ALT
    svtype_alt_match = re.search(r"SVTYPE_ALT=([^;]+)", info)
    if svtype_alt_match:
        svtype_alt = svtype_alt_match.group(1)
    else:
        svtype_alt = None

    return svlen, svtype_alt

def update_counts(svlen, svtype_alt):
    if 'Alu' in svtype_alt:
        dALU.append(svlen)
    elif 'L1' in svtype_alt:
        dLINE1.append(svlen)
    elif 'SVA_' in svtype_alt:
        dSVA.append(svlen)
    elif 'ALU' in svtype_alt:
        ALU.append(svlen)
    elif 'LINE1' in svtype_alt:
        LINE1.append(svlen)
    elif 'SVA' == svtype_alt:
        SVA.append(svlen)
    else:
        print("I missed a case scenario!")

ALU = []
LINE1 = []
SVA = []
dALU = []
dLINE1 = []
dSVA = []
# Getting svlens and putting them in respective list
full_path = f'/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/starting_callsets/{file}'
with open(full_path) as g:
    for line in g.readlines():
        if line.startswith("#"):
            continue
        svlen, svtype_alt = take_line_get_svlen(line)
        update_counts(svlen, svtype_alt)
            
files_alt = ['ALU', 'LINE1', 'SVA', 'dALU', 'dLINE1', 'dSVA']
range_counts_list = []

for file in files_alt:
    integer_list = eval(file)
    range_counts = {
        '0-50': 0,
        '50-100': 0,
        '100-500': 0,
        '500-1000': 0,
        '1000-10000': 0,
    }

    # Counting the occurrences in each range
    for num in integer_list:
        num = int(num)
        if 0 <= num <= 50:
            range_counts['0-50'] += 1
        elif 50 < num <= 100:
            range_counts['50-100'] += 1
        elif 100 < num <= 500:
            range_counts['100-500'] += 1
        elif 500 < num <= 1000:
            range_counts['500-1000'] += 1
        elif 1000 < num <= 10000:
            range_counts['1000-10000'] += 1
        else:
            # print('yo', num)
            continue
            # there are no values greater than 10,000. So next.

    range_counts_list.append(range_counts)

# Creating the stacked bar chart
N = 5
ind = np.arange(N)  # the x locations for the groups
width = 0.8

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1])

bottom = np.zeros(N)  # variable to keep track of the bottom positions for each bar

for i, file in enumerate(files_alt):
    range_counts = range_counts_list[i]
    vals = list(range_counts.values())

    ax.bar(ind, vals, width, bottom=bottom, label=file)
    
    bottom += vals

ax.set_xlabel('Size')
ax.set_ylabel('Count')
ax.set_title('Size Distribution of MEIs')
keys = range_counts.keys()
ax.set_xticks(ind, keys)
ax.set_yscale('log')
ax.set_ylim(1,1000000)
ax.legend(labels=['ALU INS', 'LINE1 INS', 'SVA INS', 'ALU DEL', 'LINE1 DEL', 'SVA DEL'])

# Save figure, then show
fig.savefig('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/figures/size_distribution_stacked.jpg', bbox_inches='tight', dpi=150)
plt.savefig('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/figures/size_distribution_stacked.pdf', format='pdf', bbox_inches='tight', dpi=150)




# Make table for clarity
data = range_counts_list
df = pd.DataFrame(data)
hd = ['ALU', 'LINE1', 'SVA', 'dALU', 'dLINE1', 'dSVA']
df["Type"] = hd
df = df[['Type', '0-50', '50-100', '100-500', '500-1000', '1000-10000']]
df.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/size_distribution_MEIs.csv", index=None)



'''

Below code was used to make 6 individual charts, not 1 stacked chart

'''




# import numpy as np
# import matplotlib.pyplot as plt

# files = '''
# 20200623ALU.final_comp.vcf
# 20200623LINE1.final_comp.vcf
# 20200623SVA.final_comp.vcf
# 20200623ALU.DEL.final_comp.vcf
# 20200623LINE1.DEL.final_comp.vcf
# 20200623SVA.DEL.final_comp.vcf
# '''.split()
# files_alt = ['ALU', 'LINE1', 'SVA', 'dALU', 'dLINE1', 'dSVA']

# def take_line_get_svlen(line):
#     if line.startswith("#"):
#         return -1
#     info = line.split("\t")[7]
#     svlen_index = info.find("SVLEN=")  # Find the starting index of "SVLEN="
#     if svlen_index != -1:
#         svlen_end_index = info.find(";", svlen_index)  # Find the ending index of the SVLEN value
#         if svlen_end_index != -1:
#             svlen_value = info[svlen_index + len("SVLEN="):svlen_end_index]  # Extract the SVLEN value
#             return svlen_value
#         else:
#             print("SVLEN value is not properly terminated.")
#             return -1
#     else:
#         print("SVLEN key not found in the string.")
#         return -1
    
# def update_counts(svlen, file_index):
#     if svlen != -1:
#         eval(files_alt[file_index]).append(svlen)
        
# ALU = []
# LINE1 = []
# SVA = []
# dALU = []
# dLINE1 = []
# dSVA = []
# # Getting svlens and putting them in respective list
# for i, file in enumerate(files):
#     full_path = f'/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20200623_RawData/{file}'
#     with open(full_path) as g:
#         for line in g.readlines():
#             svlen = take_line_get_svlen(line)
#             update_counts(svlen, i)
            

# files_alt = ['ALU', 'LINE1', 'SVA', 'dALU', 'dLINE1', 'dSVA']
# for file in files_alt:
#     integer_list = eval(file)
#     range_counts = {
#         '0-50': 0,
#         '50-100': 0,
#         '100-1000': 0,
#         '1000-10000': 0,
#         '10000+': 0
#     }

#     # Creating bins based off svlen size
#     for num in integer_list:
#         num = int(num)
#         if 0 <= num <= 50:
#             range_counts['0-50'] += 1
#         elif 50 < num <= 100:
#             range_counts['50-100'] += 1
#         elif 100 < num <= 1000:
#             range_counts['100-1000'] += 1
#         elif 1000 < num <= 10000:
#             range_counts['1000-10000'] += 1
#         else:
#             range_counts['10000+'] += 1

#     # Making the chart for the current file
#     N = 5
#     vals = range_counts.values()
#     ind = np.arange(N) # the x locations for the groups
#     width = 0.80
#     fig = plt.figure()
#     ax = fig.add_axes([0,0,1,1])
#     ax.bar(ind, vals, width, color='cornflowerblue')
#     ax.set_xlabel('Size')
#     ax.set_ylabel('Count')
#     ax.set_title(f'Size Distribution of {file} MEIs')
#     keys = range_counts.keys()
#     ax.set_xticks(ind, keys)
#     plt.yscale('log')

#     # Save figure, then show
#     fig.savefig(f'/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/figures/size_distribution_{file}.jpg',bbox_inches='tight', dpi=150)
#     plt.show()