def get_shorthand(base_name):
    replacements = [("x_linked", "x"),("autosomal", "auto"), ("dominant", "dom"), ("recessive", "rec"), ("de_novo", "denovo")]
    shorthand = base_name
    for old,rep in replacements:
        shorthand=shorthand.replace(old, rep)
    return shorthand

    
def formatting(full_path, shorthand):
    dfs = pd.read_csv(full_path, sep="\t")

    dfs['SV_start'] = dfs['start']
    dfs['SV_end'] = dfs['end']

    # Handling insertions with duplicate start and end fields from the VCF.
    # AnnotSV automatically makes the end field +1 for these types of variants.
    # Therefore, we have to add 1 to the end field before merging with AnnotSV output.
    duplicate_rows = dfs[dfs['SV_start'] == dfs['SV_end']]
    for index, row in duplicate_rows.iterrows():
        row['SV_end'] += 1
        dfs.loc[index, 'SV_end'] = row['SV_end']
        print(row)

    dfs = dfs.astype({'SV_start': str, 'SV_end': str})
    dfs['AnnotSV_ID'] = dfs['chrom'].str.slice(start=3) + '_' + dfs['SV_start'] + '_' + dfs['SV_end'] + '_' + dfs['sub_type'] + '_1'
    dfs.drop(columns=['chrom', 'start', 'end', 'sub_type', 'variant_id', 'gene'], inplace=True)
    
    aggregations = {
        'family_id': lambda x: list(x),
        'family_members': lambda x: list(x),
        'family_genotypes': lambda x: list(x),
        'samples': lambda x: list(x),
        'family_count': 'first',
    }
    dfs = dfs.groupby('AnnotSV_ID').agg(aggregations)
    dfs[f'Inheritance_Type_{shorthand}'] = f'{shorthand};'
    dfs.rename(columns={'family_id': f'family_id_{shorthand}', 'family_members': f'family_members_{shorthand}', 'family_genotypes': f'family_genotypes_{shorthand}', 'samples': f'samples_{shorthand}', 'family_count': f'family_count_{shorthand}'}, inplace=True)
    return dfs

def merge_inh_patterns(candidates):
    inheritance_filt = candidates.filter(regex='^Inheritance_Type_')
    candidates["Inheritance_Patterns"] = inheritance_filt.fillna("").sum(axis=1)
    candidates.drop(columns=inheritance_filt, inplace=True)
    return candidates
    
if __name__=="__main__":
    import pandas as pd
    import os

    data_directory = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/"
    annotations_path = os.path.join(data_directory, "annotations/annotated.tsv")
    annotations_df = pd.read_csv(annotations_path, sep='\t')

    inh_directory = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/GEMINI/Inheritance/"
    for filename in sorted(os.listdir(inh_directory)):
        # ADHD, LI, RI, SRS, then ASD
        # auto_dom, auto_rec, denovo, x_denovo, x_dom, then x_rec
        base_name = filename.split(".")[0]
        if not base_name.startswith("ASD_"):
            continue
        full_path = os.path.join(inh_directory,filename)
        shorthand = get_shorthand(base_name)
        try:
            pd.read_csv(full_path, sep="\t")
            dfs = formatting(full_path, shorthand)
            annotations_df = annotations_df.join(dfs, on="AnnotSV_ID")
        except:
            print(f"{filename} is empty!")
            
    seg_candidates = merge_inh_patterns(annotations_df)
    seg_candidates_path = os.path.join(data_directory, "seg_candidates.csv")
    print(seg_candidates_path)
    seg_candidates.to_csv(seg_candidates_path, index=False)
    
