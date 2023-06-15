import pandas as pd

sv = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/svanna_candidates.csv")
l = sv[~pd.isna(sv['Inheritance_Patterns'])][['AnnotSV_ID', 'Inheritance_Patterns']].drop_duplicates()['Inheritance_Patterns'].tolist()

df_cols=['ASD_only', 'ASD_LI', 'ASD_RI']
auto_rec=[0]*3
auto_dom=[0]*3
x_rec=[0]*3
x_dom=[0]*3
x_den=[0]*3
den=[0]*3

def count_inh(sublist):
    for p in sublist:

        ind=0
        if 'ASD_only' in p:
            ind=0
        elif 'ASD_LI' in p:
            ind=1
        elif 'ASD_RI' in p:
            ind=2

        if 'auto_rec' in p:
            auto_rec[ind]+=1
        elif 'auto_dom' in p:
            auto_dom[ind]+=1
        elif 'x_rec' in p:
            x_rec[ind]+=1
        elif 'x_dom' in p:
            x_dom[ind]+=1
        elif 'x_den' in p:
            x_den[ind]+=1
        elif 'den' in p:
            den[ind]+=1

[count_inh(x.split(";")) for x in l]

tot = [auto_rec, auto_dom, x_rec, x_dom, x_den, den]
df = pd.DataFrame(tot).rename(
    columns={0:'ASD_only', 1: 'ASD_LI', 2: 'ASD_RI'}
)
df['Inheritance_Pattern'] = ['auto_rec', 'auto_dom', 'x_rec', 'x_dom', 'x_den', 'den']
df = df[['Inheritance_Pattern','ASD_only','ASD_LI','ASD_RI']]
# df = df.set_index(df.columns[0], drop=True)
df.to_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/tables/inheritance_counts_by_phenotype.csv", index=None)