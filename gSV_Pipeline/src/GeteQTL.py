#Function to read in VCFs as dataframes
def readVCF2df(filename): 
    ''' 
    filename is a vcf file, return a dataframe with dtype str and a string of the annotation part of the vcf file 

    ''' 
    import gzip 
    if filename.endswith('.gz'): 
        fo = gzip.open(filename,'rt') 
    else: 
        fo = open(filename) 

    header = [] 
    for line in fo: 
        if line.startswith('##'): 
            header.append(line) 
        else: 
            break 
    fo.close() 

    columns = line.strip().split('\t') 
    df = pd.read_csv(filename, sep='\t', comment='#', header=None,na_filter = False, dtype=str) 
    df.columns = columns 

    return df, header 


# Function to find eQTL data
#If a match is found, return the tissue_count_ME_isCausal value
def findeQTL(ME, error = 100): 
    ''' 
    ME is like 3_50879175_INS_ME_ALU 
    return most similar ID in df_1KGP, and allele frequency 
    site difference allowed is error 
    ''' 
    ME_chr, ME_pos, other = ME.split('_', maxsplit=2) 
    
    #casting
    #caviar_df_short has, by default, a snp_chr column with mixed datatypes.
    combined_eqtls_df['snp_chr'] = combined_eqtls_df['snp_chr'].astype('str')
    ME_pos = int(ME_pos) 
    
    currentChrom = combined_eqtls_df[combined_eqtls_df['snp_chr'] == ME_chr] 
    positions = set(currentChrom['snp_pos']) 
    
    if ME_pos in positions: 
        ME_position = ME_pos 
    else: 
        ME_position = None 
        for i in range(error): 
            
            if ME_pos + 1 + i in positions: 
                ME_position = ME_pos + 1 + i 
                break 

            if ME_pos -1 - i in positions: 
                ME_position = ME_pos -1 -i 
                break 

    if ME_position is None: 
        return 0, None, None, None
    
    location = combined_eqtls_df[combined_eqtls_df['snp_pos'] == ME_position]
    tissueCount = location['tissue_count_ME_isCausal'].values[0]
    tissueNames = location['Tissues'].values[0]
    tissueGenes = location['gene_name'].values[0]
    eqtlOrigin = location['eqtl_origin'].values[0]
    return tissueCount, tissueNames, tissueGenes, eqtlOrigin

# Function to find population AF (from same database)
# If a match is found, return the tissue_count_ME_isCausal value
def findAF(ME, error = 100): 
    ME_chr, ME_pos, other = ME.split('_', maxsplit=2) 
    
    ME_pos = int(ME_pos) 
    
    if ME_chr=='Y':
        return None
    currentChrom = dc_1KGB[ME_chr]
    positions = set(currentChrom['POS']) 
    
    if ME_pos in positions: 
        ME_position = ME_pos 
    else: 
        ME_position = None 
        for i in range(error): 
            
            if ME_pos + 1 + i in positions: 
                ME_position = ME_pos + 1 + i 
                break 

            if ME_pos -1 - i in positions: 
                ME_position = ME_pos -1 -i 
                break 

    if ME_position is None: 
        return None
    
    row = currentChrom[currentChrom['POS'] == ME_position]
    alleles = list(row.iloc[0,9:])
    dc = {'0|0':0,'0|1':1,'1|0':1,'1|1':2} 
    popAF = sum([dc[i] if i in dc else 0 for i in alleles]) / len([i for i in alleles if i in dc]) /2 
    
    return popAF


if __name__ == "__main__":
    import pandas as pd

    ### Import the merged eQTL file
    combined_eqtls_df=pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/eqtl_combined.csv")

    df_1KGP, header_1KGP = readVCF2df('/lab01/Projects/Sammy_Projects/NJLAGS_MEIs/data/20211113_eQTL/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf') 
    df_1KGP['POS'] = df_1KGP['POS'].astype(int) 
    dc_1KGB = {k:v for k,v in df_1KGP.groupby('#CHROM')}

    # Merging eQTL data with rest of annotaitons
    candidates_af = pd.read_csv("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/af_candidates.csv")
    candidates_eqtl = candidates_af.copy()


    # Creating lists for adding columns
    candList = candidates_eqtl['AnnotSV_ID'].tolist()
    eQTLcounts = list()
    eQTLnames = list()
    eQTLgenes = list()
    eQTLorigin = list()

    for ME in candList:
        count, names, genes, origin = findeQTL(ME)
        eQTLcounts.append(count)
        eQTLnames.append(names)
        eQTLgenes.append(genes)
        eQTLorigin.append(origin)
    popAFs = list()

    for ME in candList:
        AF = findAF(ME)
        popAFs.append(AF)

    candidates_eqtl['tissue_count_ME_isCausal'] = eQTLcounts
    candidates_eqtl['causal_tissues'] = eQTLnames
    candidates_eqtl['eQTL_genes'] = eQTLgenes
    candidates_eqtl['eQTL_origin'] = eQTLorigin
    candidates_eqtl['popAF'] = popAFs


    # Marking candidates whose eQTLs are causal in some kind of brain tissue
    # BrainTissue = list()
    # for i in (candidates_eqtl['causal_tissues'].str.join(sep=',').astype(str)):
    #     #print(type(i))
    #     if 'Brain' in i:
    #         BrainTissue.append('True')
    #     else:
    #         BrainTissue.append('False')
    # Using this new version because I think old version is messed up
    
    candidates_eqtl["causal_tissues"] = candidates_eqtl["causal_tissues"].fillna('')

    ### Don't need this since I'm prefiltering eqtls for brain tissue
    # BrainTissue = ["True" if ('Brain' in str(x)) else "False" for x in candidates_eqtl["causal_tissues"].tolist()]
    # candidates_eqtl['isCausal_for_brain_tissue'] = BrainTissue


    # Adding a few needed columns
    columns = candidates_eqtl['Location'].str.split('-', expand=True)
    candidates_eqtl = pd.concat([candidates_eqtl, columns], axis=1)


    # Save csv
    candidates_eqtl_path = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/eqtl_candidates.csv"
    candidates_eqtl.to_csv(candidates_eqtl_path, index=False)