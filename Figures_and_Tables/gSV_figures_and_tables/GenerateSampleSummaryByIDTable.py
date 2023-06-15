# Prior to running this, ran GetMEIpassGenotypes.sh in this dir

import pandas as pd

def return_counts(types):
    counts = [0]*272
    for typ in types:
        with open(f'/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/MEI_pass_genotypes/{typ}.csv') as f:
            lines = f.readlines()
            for line in lines:
                all_cols = line.split(' ')
                for col_ind in range(272):
                    col = all_cols[col_ind]
                    if '1' in col[:3]:
                        counts[col_ind]+=1
    return counts

mei_pass_ins = ['ALU', 'LINE1', 'SVA']
mei_pass_del = ['dALU', 'dLINE1', 'dSVA']
ins = return_counts(mei_pass_ins)
dell = return_counts(mei_pass_del)

with open('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/MEI_pass_genotypes/WGS_ID_order.txt') as f:
    names = f.readline().split(' ')[:-1]
    short_names = [x.split('_')[0] for x in names]

ins_dict = {x:y for x,y in zip(short_names, ins)}
del_dict = {x:y for x,y in zip(short_names, dell)}
pass_df = pd.read_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/2023_05_08_sample_summary_by_id.csv')
pass_df['MEI_pass'] = pass_df['WGS ID'].map(ins_dict)
pass_df['dMEI_pass'] = pass_df['WGS ID'].map(del_dict)
pass_df['MEI_total'] = pass_df['MEI_pass'] + pass_df['dMEI_pass']
pass_df.to_csv('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/2023_05_09_sample_summary_by_id.csv')