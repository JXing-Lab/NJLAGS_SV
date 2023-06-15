import pandas as pd

if __name__ == "__main__":
    mod_merged_callset = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/modified_SURVIVOR_merged_callset.vcf"
    candidates_impc = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/impc_candidates.csv")

    site = 0
    hetero = [0] * 100000
    homo = [0] * 100000
    homoref = [0] * 100000
    nocall = [0] * 100000

    ls = []
    columns = [] #use to store the order for the dataframe
    for line in open(mod_merged_callset):
        if not line.startswith('#'):
            es = line.split()
            dc = {}
            dc['SV_chrom'] = es[0].replace("chr", "")
            dc['SV_start'] = es[1]
            column = 10
            while column < len(es):
                if (es[column][0] == '0') and (es[column][2] == '0'):
                    homoref[site] += 1
                elif (es[column][0] == '.') or (es[column][2] == '.'):
                    nocall[site] += 1
                elif es[column][0] == es[column][2]:
                    homo[site] += 1
                elif es[column][0] != es[column][2]:
                    hetero[site] += 1
                column += 1
            dc['HET'] = hetero[site]
            dc['HOM'] = homo[site]
            dc['HOMREF'] = homoref[site]
            dc['NOCALL'] = nocall[site]
            # 272 stands for the sample number, according to this formula: 
            # Alternative Allele Frequency = (2 "1/1" + 1 "1/0") / (2 sample_number - 2 './.')
            # Which equals = ("1/1" + "1/0"/2) / (sample_number - "./.")
            if dc['NOCALL']==272:
                dc['AF']=0
            else:
                dc['AF'] = (dc['HOM'] + dc['HET']/2) / (272-dc['NOCALL'])
            dc['ALTCOUNT'] = 2*dc['HOM'] + dc['HET']
            ls.append(dc)
        site+=1
    dfGT = pd.DataFrame(ls)
    dfGT.drop(columns=['HET', 'HOM', 'HOMREF', 'NOCALL', 'ALTCOUNT'], inplace=True)
    dfGT = dfGT.astype({'SV_start': int, 'SV_chrom': str})
    dfGT.drop_duplicates(subset=['SV_chrom', 'SV_start'], keep='first', inplace=True)
    candidates_impc = candidates_impc.astype({'SV_start': int, 'SV_chrom': str})
    candidates_af = pd.merge(candidates_impc, dfGT, on=['SV_start', 'SV_chrom'])

    candidates_af_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/af_candidates.csv"
    candidates_af.to_csv(candidates_af_path, index=False)