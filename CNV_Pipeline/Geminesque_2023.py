"""
Geminesque.py
> Takes in a .vcf file annotated by AnnotSV
> Finds variants that fit a given inheritance pattern
> Uses logic based off of GEMINI
"""

import pandas as pd
from collections import defaultdict

def Geminesque_main(phenotype, AnnotSV_input):
    
    print("calling main Geminesque method")
    
    # import file
    file_genotype = AnnotSV_input
    
    df_genotype = pd.read_csv(file_genotype, sep='\t')#8670
    
    # extract individuals from reference file
    header_list = df_genotype.columns.values.tolist()
    indv_list = header_list[14:497] # nuber of individuals is HARDCODED
    
    df_genotype['sample_name'] = df_genotype['Samples_ID']
    
    
    # add column "affectCDS"
    df_genotype['affectCDS'] = df_genotype['Location2'].apply(lambda x:x in ["5'UTR-3'UTR","5'UTR-CDS","CDS-3'UTR","CDS"])#4016, total 5288
    df_genotype = df_genotype[df_genotype['Annotation_mode'] == 'split']#6024
    
    
    # load in pedigree
    fped = "../Data/" + phenotype + "/" + phenotype + ".ped"
    df_ped = pd.read_csv(fped,  sep='\t', dtype=str)#304, 47 families
    columns_ped = ['family_id','individual_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
    df_ped.columns = columns_ped
    df_ped = df_ped.iloc[:-1,:]
    # convert to integers
    #for column in df_ped:
        #df_ped[column] = pd.to_numeric(df_ped[column])
        
    df_ped['sex'] = pd.to_numeric(df_ped['sex'])
    df_ped['phenotype'] = pd.to_numeric(df_ped['phenotype'])
    
    # individuals
    #individuals = '1002001 1002002 1002003 1002004 1002005 1002006 1003001 1003002 1003003 1003004 1003005 1004001 1004002 1004003 1004005 1005001 1005002 1005003 1005004 1005005 1005006 1006001 1006002 1006003 1007001 1007002 1007003 1007004 1007005 1007006 1008001 1008002 1008003 1009002 1009003 1009005 1010001 1010003 1010004 1010005 1011001 1011002 1011003 1011004 1012001 1012002 1012003 1012004 1012005 1012007 1012008 1012009 1012010 1012011 1012012 1012013 1014001 1014002 1014003 1014004 1014005 1019001 1019002 1019003 1019004 1019005 1019007 1019008 1019009 1020001 1020002 1020003 1020004 1020008 1020009 1020010 1020011 1020012 1021002 1021003 1021004 1022001 1022002 1022003 1022004 1022005 1025001 1025002 1025003 1025004 1026001 1026002 1026003 1026004 1026005 1027002 1027003 1027004 1027005 1027006 1030001 1030002 1030003 1030004 1033001 1033002 1033003 1034001 1034002 1034003 1035001 1035002 1035003 1035004 1035005 1036001 1036003 1036004 1036005 1037001 1037002 1037003 1037004 1037005 1037006 1038001 1038002 1038003 1038005 1039002 1039003 1039004 1042001 1042002 1042003 1043001 1043002 1043003 1043004 1043005 1044002 1044003 1044005 1045001 1045002 1045003 1045005 1046001 1046002 1046003 1046004 1046005 1048001 1048003 1048004 1048005 1049001 1049002 1049003 1051001 1051002 1051003 1051004 1051005 1051006 1051007 1051008 1051009 1051010 1052001 1052002 1052003 1052004 1054002 1054003 1054004 1058001 1058002 1058003 1058004 1058005 1060001 1060002 1060003 1060004 1060005 1061001 1061002 1061003 1061004 1061005 1062001 1062003 1063001 1063002 1063003 1063005 1063006 1065001 1065002 1065003 1065004 1065005 1065006 1066001 1066002 1066003 1069001 1069002 1069003 1070001 1070002 1070003 1070004 1071001 1071002 1071003 1071005 1071006 1072001 1072002 1072003 1072004 1072005 1075003 1075004 1075005 1077001 1077002 1077003 1077004 1078001 1078002 1078003 1078004 1079002 1079003 1080001 1080002 1080003 1080004 2000001 2000002 2000003 2002001 2002003 2002004 2002005 2002006 2006001 2006002 2006003 2006004 2006005 2008001 2008002 2008003 2008004 2008005 2008006 2009001 2009002 2009003 2009004 2010003 2012001 2012002 2012003 2012004 2012005 2012006 2013001 2013002 2013003 2013004 2013006 2013007 2013008 2013009 2015001 2015002 2015003 2015004 2015005 2016001 2016002 2016003 2016004 2016005 2016006 2019001 2019002 2019003 2019004 2019005 2023001 2023002 2023003 2026001 2026002 2026003 2026004 2026005 2026006 2028001 2028002 2028003 2028005 2035001 2035002 2035003 2036001 2036002 2036003 2039001 2039002 2039003 2039004 2039005 2039007 2040001 2040002 2040003 2040004 2040005 2040006 2043001 2043002 2043003 2043004 2043005 2045001 2045002 2045003 2045005 2045007 2045008 2045009 2045010 2045011 2045012 2045013 2046001 2046002 2046003 2046004 2047001 2047002 2047003 2047004 2048001 2048002 2048003 2048004 2048005 2054002 2054003 2054005 2055001 2055002 2055003 2055004 2055005 2056002 2056004 2056005 2057003 2057004 2057005 2058002 2058003 2058004 2058005 2058007 2058008 2058009 2061002 2061003 2061004 2062002 2062003 2062004 2066001 2066002 2066003 2066004 2066005 2069001 2069002 2069003 2070001 2070002 2070003 2070004 2070005 2070006 2071001 2071002 2071003 2071004 2071005 2071006 2071007 2073001 2073002 2073003 2073004 2073005 2073006 2073007 2074001 2074002 2074003 2074004 2079001 2079002 2079003 2079004 2079005 2080001 2080002 2080003 2080004 2084001 2084002 2084003 2084004 2085001 2085002 2085004 2085005 2086001 2086002 2086003 2086004 2086005 2087001 2087002 2087003 2087004 2088001 2088002 2088003 2088004 2088005 2091002 2091003 2091004 2091005 2092001 2092002 2092003 2092004 2098001 2098002 2098003 2098004 2099001 2099002 2099003 2099004 2102001 2102002 2102003 2102004 2105001 2105002 2105003 2105004 2105005 2105006 2105007 2105008 2105009 2105010 2106001 2106002 2106003 2106004 2106005 2106006 2106007 2106008 2107001 2107002 2107003 2107004 2112001 2112002 2112003 2112004 2112005 2113001 2113002 2113003 2113004 2113005 2113006 2113007 2113008 2113009 2115001 2115002 2115003 2115005 2117002 2117004 2117005 2118001 2118002 2118003 2118004 2118005 2118006 2118007 2118009 2118010 2118011 2118012 2118013 2118014 2118015 2118016 2118017'.split()
    individuals = indv_list
    
    columns_nochange = [e for e in df_genotype.columns if e not in individuals]
    columns_change = [e for e in df_genotype.columns if e in individuals]
    
    # melt them together
    df_genotype = pd.melt(df_genotype, id_vars=columns_nochange, value_vars=individuals,var_name='sample_name1',value_name='SV_genotype')
    df_genotype = df_genotype[df_genotype['SV_genotype'] != './.:0:0:0:0:0:null:0'].copy()#2040
    df_genotype = df_genotype[df_genotype['SV_genotype'].apply(lambda x:x.split(':')[0] != "./.")]
                                                   
    columns_ped = ['family_id','individual_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
    
    # only add family_id to df_genotype
    dc_individual2family = {}
    for c in ['individual_id','paternal_id', 'maternal_id']:
        for tdf in [df_ped]:
            tdf = tdf[['family_id',c]]
            tdf = tdf[tdf[c] != '0']
            dc_individual2family.update(dict(zip(tdf[c], tdf['family_id'])))
    
    # EVERYTHING IS FUNCTIONAL UP TO THIS POINT
    
    # add family_id
    # df_genotype['family_id'] = df_genotype['sample_name'].map(dc_individual2family)
    
    # of course there's a built-in function for this
    df_genotype = df_genotype.loc[:,~df_genotype.columns.duplicated()]
    df_genotype = df_genotype.rename(columns={'sample_name':'sample_name_ori', 'sample_name1':'sample_name'})
    
    # convert the column sample_name into a series
    family_id_ser = df_genotype['sample_name'].squeeze()
    
    # map the dc_individual2family dictionary onto the series
    family_id_ser = family_id_ser.map(dc_individual2family)
    
    # assign the series family_id_ser to a new column family_id in the dataframe df_genotype
    
    df_genotype = df_genotype.assign(family_id = family_id_ser)
    df_genotype = df_genotype[df_genotype['family_id'].notnull()]
    
    df_genotype['family_id'] = df_genotype['sample_name'].map(dc_individual2family)
    
    # add sex
    dc_individual2sex = {}
    for tdf in [df_ped]:
        dc_individual2sex.update(dict(zip(tdf['individual_id'], tdf['sex'])))
    
    # convert the column sample_name into a series
    sex_ser = df_genotype['sample_name'].squeeze()
        
    # df_genotype['sex'] = df_genotype['sample_name'].map(dc_individual2sex)
    
    # map the dc_individual2family dictionary onto the series
    sex_ser = sex_ser.map(dc_individual2sex)
    
    # assign the series family_id_ser to a new column family_id in the new dataframe df_geno
    df_genotype = df_genotype.assign(sex = sex_ser)
    
    # simplify df_genotype to keep only several columns
    df_geno = df_genotype[['AnnotSV_ID', 'SV_chrom', 'SV_start', 'SV_end', 'SV_length', 'SV_type','Gene_name', 'sample_name','SV_genotype','ACMG_class','family_id','sex']].copy()
    
    
    # check how many genes were affected by multiple SVs
    l = [e for e in df_geno.groupby(['Gene_name','sample_name']) if e[1].shape[0] >1]#45 ('Gene_name','sample_name') were affected by multiple SVs
    df_multipleSV_gene = pd.DataFrame([e[0] for e in l])
    df_multipleSV_gene.columns = ['Gene_name', 'sample_name']
    df_multipleSV_gene['AnnotSV_IDs'] = [','.join(e[1]['AnnotSV_ID']) for e in l]
    df_multipleSV_gene['AnnotSV_ID_N'] = [len(set(e[1]['AnnotSV_ID'])) for e in l]
    
    # print out preliminary results for genes affected by multiple SVs
    #print(df_multipleSV_gene['Gene_name'].value_counts().to_dict())
    #print(df_multipleSV_gene['sample_name'].value_counts().to_dict())
    
# =============================================================================
#     # output SV counts for genes
#     gene_value_counts = df_multipleSV_gene['Gene_name'].value_counts().to_dict()
#     file_out = open("gene_value_counts.tsv", "w")
#     for key, value in gene_value_counts.items():
#         file_out.write(key + "\t" + str(value) + "\n")
#     file_out.close()
#     
#     # output SV counts for individuals
#     sample_value_counts = df_multipleSV_gene['sample_name'].value_counts().to_dict()
#     file_out = open("sample_value_counts.tsv", "w")
#     for key, value in sample_value_counts.items():
#         file_out.write(key + "\t" + str(value) + "\n")
#     file_out.close()
# =============================================================================
    
    # create a dataframe to save results
    df_genes = pd.DataFrame(index=sorted(set(df_geno['Gene_name'])))# 361 genes
    
    # add column chr, CDSaffected 
    tdc = dict(zip(df_geno['Gene_name'], df_geno['SV_chrom'])) 
    df_genes['chrom'] = [tdc[i] for i in df_genes.index] 
    genes_CDSaffected = set(df_genotype[df_genotype['affectCDS']]['Gene_name']) 
    df_genes['CDSaffected'] = [i in genes_CDSaffected for i in df_genes.index] 
    tdc = df_geno.groupby('Gene_name')['AnnotSV_ID'].apply(lambda x:','.join(sorted(set([str(e) for e in x])))).to_dict() 
    df_genes['SVs'] = [tdc[i] if i in tdc else '' for i in df_genes.index] 
    df_genes['SVsN'] = df_genes['SVs'].apply(lambda x:0 if x == '' else x.count(',')+1) # 137 genes found in more than one SVs 
    
    # add Inds, IndN, Fams, FamN 
    tdc = df_geno.groupby('Gene_name')['family_id'].apply(lambda x:','.join(sorted(set([str(e) for e in x])))).to_dict() 
    df_genes['Fams'] = [tdc[i] if i in tdc else '' for i in df_genes.index] 
    df_genes['FamN'] = df_genes['Fams'].apply(lambda x:0 if x == '' else x.count(',')+1) # 193 genes found in more than one family 
    tdc = df_geno.groupby('Gene_name')['sample_name'].apply(lambda x:','.join(sorted(set([str(e) for e in x])))).to_dict() 
    df_genes['Inds'] = [tdc[i] if i in tdc else '' for i in df_genes.index] 
    df_genes['IndN'] = df_genes['Inds'].apply(lambda x:0 if x == '' else x.count(',')+1) # 213 genes found in more than one individual 
    
    # add ADHD affectedInds,affectedFams, affectedIndN,affectedFamN, column
    key = 'ASD'
    tdf_ped = df_ped
    
    tdf_ped = tdf_ped[tdf_ped['phenotype'] == '2']#68
    tF = set(tdf_ped['family_id'])#47
    tI = set(tdf_ped['individual_id'])#68
    tdf_geno = df_geno[df_geno['sample_name'].isin(tI)]#1440
    tdc = tdf_geno.groupby('Gene_name')['sample_name'].apply(lambda
    x:','.join(sorted(set([str(e) for e in x])))).to_dict()
    df_genes[key + '|affectedInds'] = [tdc[i] if i in tdc else '' for i
    in df_genes.index]
    df_genes[key + '|affectedIndN'] = df_genes[key +
    '|affectedInds'].apply(lambda x:0 if x == '' else x.count(',')+1)
    tdc = tdf_geno.groupby('Gene_name')['family_id'].apply(lambda
    x:','.join(sorted(set([str(e) for e in x])))).to_dict()
    df_genes[key + '|affectedFams'] = [tdc[i] if i in tdc else '' for i
    in df_genes.index]
    df_genes[key + '|affectedFamN'] = df_genes[key +
    '|affectedFams'].apply(lambda x:0 if x == '' else x.count(',')+1)
    
    #
    df_geno['SV_genotype'] = df_geno['SV_genotype'].apply(lambda x:x.split(':')[0])
    
    
    def checkGeno(individual, tdc_genotype): 
        '''tdc_genotype is a dict like 
        {'10C117899': '0/1', '10C117871': '0/1', '10C117873': '0/1'} 
        individual is an individual, which can be a key of tdc_genotype. 
        if the individual is not in tdc_genotype, it's Genotype is homref. 
        if it is in tdc_genotype, '0/0','0|0': homref 
            '0/1','0|1': het 
            '1/1', '1|1': homalt 
            './.', '.|.': missing 
        ''' 
        if individual in ['0',0]: 
            return 'missing' 
        # individuals with no annotation, consider as homref by default 
        if individual not in tdc_genotype: 
            return 'homref' 
        geno = tdc_genotype[individual] 
        # print(individual)
        if geno in ['./.','.|.']: 
            return 'missing' 
        if geno in ['0/0','0|0']: 
            return 'homref' 
        if geno in ['0/1','0|1']: 
            return 'het' 
        if geno in ['1/1', '1|1']: 
            return 'homalt' 
        print(geno,'undefined case') 
    
        return 'missing' 
    
    def getPhenoFromCode(i):
        if i in [1,'1']:
            return 'unaffected'
        if i in [2,'2']:
            return 'affected'
        return 'missing'
    
    def ped2str(tped): 
        '''return a string description of a ped dataframe 
    
         family_id individual_id paternal_id maternal_id  sex  phenotype  paternal_pheno  maternal_pheno paternal_geno maternal_geno individual_geno 
    
    263       2057     10C117874           0           0    1        1.0             0.0             0.0        homref        homref          homref 
    
    264       2057     10C117872           0           0    2        1.0             0.0             0.0        homref        homref          homref 
    
    265       2057     10C117873   10C117874   10C117872    1        1.0             1.0             1.0        homref        homref             het 
    
    266       2057     10C117899   10C117874   10C117872    2        2.0             1.0             1.0        homref        homref             het 
    
    267       2057     10C117871   10C117874   10C117872    1        1.0             1.0             1.0        homref        homref             het 
    
        ''' 
        family_ids = list(tped['family_id']) 
        tls = [] 
        for family_id in family_ids: 
            ttped = tped[tped['family_id'] == family_id] 
            for _, r in ttped.iterrows(): 
                tls.append('variant_pos:{individual_SVs},family:{family_id}[child:{individual_id}({individual_geno}|{phenotype}),mom:{maternal_id}({maternal_geno}|{maternal_pheno}),dad:{paternal_id}({paternal_geno}|{paternal_pheno})]'.format(individual_SVs=r['individual_SVs'], family_id=family_id, individual_id=r['individual_id'], individual_geno=r['individual_geno'], phenotype=getPhenoFromCode(r['phenotype']), maternal_id=r['maternal_id'], maternal_geno=r['maternal_geno'], maternal_pheno=getPhenoFromCode(r['maternal_pheno']), paternal_id=r['paternal_id'], paternal_geno=r['paternal_geno'], paternal_pheno=getPhenoFromCode(r['paternal_pheno'])))
        return ';'.join(tls)
    
    def get_ped_withAllGenoPheno(gene, tdf_ped, df_geno): 
        '''filter by gene. 
        for ped file, add pheno and geno for mum and dad 
        ''' 
        tdf_ped = tdf_ped.copy() 
        # add paternal_pheno and maternal_pheno column 
        tdc_ind2pheno = dict(zip(tdf_ped['individual_id'], tdf_ped['phenotype'])) 
        tdf_ped['paternal_pheno'] = tdf_ped['paternal_id'].apply(lambda x:tdc_ind2pheno[x] if x in tdc_ind2pheno else 0) 
        tdf_ped['maternal_pheno'] = tdf_ped['maternal_id'].apply(lambda x:tdc_ind2pheno[x] if x in tdc_ind2pheno else 0) 
        tdf_geno =df_geno[df_geno['Gene_name'] == gene] 
        tdc_geno = dict(zip(tdf_geno['sample_name'], tdf_geno['SV_genotype'])) 
        tFams = set(tdf_geno['family_id'])# only keep families with selected gene affected.  
        tdf_ped = tdf_ped[tdf_ped['family_id'].isin(tFams)].copy() 
        tdf_ped['paternal_geno'] = tdf_ped['paternal_id'].apply(lambda x:checkGeno(x, tdc_geno)) 
        tdf_ped['maternal_geno'] = tdf_ped['maternal_id'].apply(lambda x:checkGeno(x, tdc_geno)) 
        tdf_ped['individual_geno'] = tdf_ped['individual_id'].apply(lambda x:checkGeno(x, tdc_geno)) 
        # add columns 'paternal_SVs','maternal_SVs' and 'individual_SVs' 
        tdc_sample2SV = tdf_geno.groupby('sample_name')['AnnotSV_ID'].apply(lambda x:','.join(x)).to_dict() 
        tdf_ped['individual_SVs'] = tdf_ped['individual_id'].apply(lambda x:tdc_sample2SV[x] if x in tdc_sample2SV else None) 
        tdf_ped['maternal_SVs'] = tdf_ped['maternal_id'].apply(lambda x:tdc_sample2SV[x] if x in tdc_sample2SV else None) 
        tdf_ped['paternal_SVs'] = tdf_ped['paternal_id'].apply(lambda x:tdc_sample2SV[x] if x in tdc_sample2SV else None) 
        tdf_ped['individual_multipleSV'] = tdf_ped['individual_SVs'].apply(lambda x:',' in str(x)) 
        tdf_ped['gene'] = gene 
        tdf_ped['gene_chrom'] = tdf_geno.iloc[0]['SV_chrom'] 
        return tdf_ped, tdf_geno 
    
    def check_de_novo(gene, df_ped, df_geno, strict=False, allow_unaffected=False, single_affected=False): 
        '''return '' if does not follow de novo pattern. 
            return families following the pattern 
            all affecteds must be het 
            [affected] all unaffected must be homref or homalt 
            at least 1 affected kid must have unaffected parents 
            [strict] if an affected has affected parents, it’s not de_novo 
            [strict] all affected kids must have unaffected (or no) parents 
            [strict] warning if none of the affected samples have parents. 
        ''' 
    
        tdf_ped, tdf_geno = get_ped_withAllGenoPheno(gene, df_ped, df_geno) 
        tped_affected = tdf_ped[tdf_ped['phenotype'] == 2] 
        # return '' if no affected 
        if tped_affected.shape[0] == 0: 
            return '' 
        affected_genos = set(tped_affected['individual_geno']) 
        if single_affected is False: 
            # all affected must be het 
            if not affected_genos == {'het'}: 
                return '' 
        else: 
            # at leaste one affected is required 
            if not 'het' in affected_genos: 
                return ''
        # all unaffected must be homref or homalt. use if allow_unaffected is False 
        unaffected_genos = set(tdf_ped[tdf_ped['phenotype'] == 1]['individual_geno']) 
        if allow_unaffected is False: 
            # unaffected cannot be 'het' 
            if not (len(unaffected_genos - {'homref'}) == 0 or len(unaffected_genos - {'homalt'}) == 0):
                return '' 
        # individuals following inheritance model 
        tped_good = tped_affected[tped_affected['individual_geno'] == 'het'] 
        if tped_good.shape[0] == 0:# if no individual with good inheritance model 
            return '' 
        #at least 1 affected kid must have unaffected parents 
        tped_denovo = tped_good[(tped_good['paternal_pheno'] == 1) & (tped_good['maternal_pheno'] == 1)] 
        if tped_denovo.shape[0] == 0: 
            return '' 
        #[strict] all affected kids must have unaffected (or no) parents. The parents cannot pheno cannot be 2 
        if strict: 
            tdf_check = tped_good[(tped_good['paternal_pheno'] != 2) & (tped_good['maternal_pheno'] != 2)] 
            if tdf_check.shape[0] != tped_good.shape[0]: 
                return '' 
            else: 
                return ped2str(tdf_check) 
        
        return ped2str(tped_good)+tdf_check["individual_SVs"] 
    
    def check_x_linked_de_novo(gene, df_ped, df_geno, single_affected=False, allow_unaffected=False): 
        '''return '' if does not follow de novo pattern. 
        return families following the pattern 
            affected female child must be het 
            affected male child must be hom_alt (or het) 
            parents should be unaffected and hom_ref 
        ''' 
        tdf_ped, tdf_geno = get_ped_withAllGenoPheno(gene, df_ped, df_geno) 
        chrom = tdf_geno.iloc[0]['SV_chrom'] 
        if chrom not in ['X', 'chrX', 'ChrX']: 
            return '' 
        tped_affected = tdf_ped[(tdf_ped['phenotype'] == 2) ] 
        if tped_affected.shape[0] == 0:# return '' if no affected 
            return '' 
        # check affected female 
        tped_affected_female = tped_affected[(tped_affected['sex'] == 2)]  
        if tped_affected_female.shape[0] != 0:# check if there is at least one affected female 
            affected_genos = set(tped_affected_female['individual_geno']) 
            if single_affected is False: 
                #affected female child must be het 
                if len(affected_genos - {'het'}) != 0: 
                    return '' 
            else: 
                # at least one affected female child is het 
                if 'het' not in affected_genos: 
                    return '' 
        # check affected male 
        tped_affected_male = tped_affected[(tped_affected['sex'] == 1)] 
        if tped_affected_male.shape[0] != 0:#check if there is at least one affected male 
            affected_genos = set(tped_affected_male['individual_geno']) 
            if single_affected is False: 
                # affected male child must be hom_alt (or het) 
                if len(affected_genos - {'het', 'homalt'}) != 0: 
                    return '' 
            else: 
                # at least one affected male is hom_alt or het 
                if ('het' not in affected_genos) and ('homalt' not in affected_genos): 
                    return '' 
        # [added by XC] unaffected cannot be het or homalt, activate if allow_unaffected is False 
        unaffected_genos = set(tdf_ped[tdf_ped['phenotype'] == 1]['individual_geno']) 
        if allow_unaffected is False: 
            if 'homalt' in unaffected_genos or 'het' in unaffected_genos: 
                return '' 
        # individuals following inheritance model 
        tped_good_female = tped_affected_female[tped_affected_female['individual_geno'] == 'het'] 
        tped_good_male = tped_affected_male[tped_affected_male['individual_geno'].isin(['het','homalt'])] 
        tped_good = pd.concat([tped_good_female,tped_good_male]) 
        # parents should be unaffected and hom_ref 
        tped_x_novo = tped_good[(tped_good['maternal_pheno'] != 2) & (tped_good['paternal_pheno'] != 2) & (tped_good['paternal_geno'] == 'homref') & (tped_good['maternal_geno'] == 'homref')] 
        if tped_x_novo.shape[0] == 0: 
            return '' 
        return ped2str(tped_x_novo)
    
    def check_autosomal_recessive(gene, df_ped, df_geno, strict=False, allow_unaffected=False, single_affected=False): 
        '''return '' if does not follow de novo pattern. 
        return families following the pattern 
            all affecteds must be hom_alt 
            [affected] no unaffected can be hom_alt (can be unknown) 
            [strict] if parents exist they must be unaffected and het for all affected kids 
            [strict] if there are no affecteds that have a parent, a warning is issued. 
        ''' 
        tdf_ped, tdf_geno = get_ped_withAllGenoPheno(gene, df_ped, df_geno) 
        tped_affected = tdf_ped[tdf_ped['phenotype'] == 2] 
        affected_genos = set(tped_affected['individual_geno']) 
        if single_affected is False: 
            # all affected must be hom_alt 
            if not affected_genos == {'homalt'}: 
                return '' 
        else: 
            # at least one affected is homalt 
            if not 'homalt' in affected_genos: 
                return '' 
    
        # [affected] no unaffected can be hom_alt (can be unknown). use if allow_unaffected is False 
        unaffected_genos = set(tdf_ped[tdf_ped['phenotype'] != 2]['individual_geno']) 
        if allow_unaffected is False: 
            if len(unaffected_genos & {'homalt'}) != 0: 
                return '' 
        
        # individuals following inheritance model 
        tped_good = tped_affected[tped_affected['individual_geno'] == 'homalt'] 
        if not strict: 
            return ped2str(tped_good) 
        #[strict] if parents exist they must be unaffected and het for all affected kids 
        To_check = tped_good.apply(lambda x:x['paternal_pheno']!=2 and x['maternal_pheno'] !=2 and (x['paternal_geno'] in ['het', 'missing']) and (x['maternal_geno'] in ['het', 'missing']), axis=1) 
        tdf_check = tped_good[To_check] 
        if tdf_check.shape[0] == 0: 
            return '' 
    
        return ped2str(tdf_check) 
    
    def check_autosomal_dominant(gene, df_ped, df_geno, strict=False, log=False, allow_unaffected=False, single_affected=False): 
        '''return '' if does not follow de novo pattern. 
        return families following the pattern 
            All affecteds must be het 
            [affected] No unaffected can be het or homalt (can be unknown) 
            de_novo mutations are not auto_dom (at least not in the first generation) 
            At least 1 affected must have 1 affected parent (or have no parents). 
            If no affected has a parent, a warning is issued. 
            [strict] All affecteds must have parents with known phenotype. 
            [strict] All affected kids must have at least 1 affected parent 
        ''' 
        tdf_ped, tdf_geno = get_ped_withAllGenoPheno(gene, df_ped, df_geno) 
        tped_affected = tdf_ped[tdf_ped['phenotype'] == 2] 
        affected_genos = set(tped_affected['individual_geno']) 
        if single_affected is False: 
            # all affected must be het 
            if not affected_genos == {'het'}: 
                if log:print(gene, 'affected_genos not het', affected_genos,'\n', tped_affected) 
                return '' 
        else: 
            # at leaste one affected is required 
            if not 'het' in affected_genos: 
                if log:print(gene, 'affected_genos not het', affected_genos,'\n', tped_affected) 
                return '' 
        # [affected] No unaffected can be het or homalt (can be unknown). only use if allow_unaffected is False 
        unaffected_genos = set(tdf_ped[tdf_ped['phenotype'] != 2]['individual_geno']) 
        if allow_unaffected is False: 
            if len(unaffected_genos & {'homalt','het'}) != 0: 
                if log:print(gene, 'unaffected_genos not good', unaffected_genos,'\n', tped_affected) 
                return '' 
        # de_novo mutations are not auto_dom (at least not in the first generation) 
        a = check_de_novo(gene, tdf_ped, df_geno, strict=strict, allow_unaffected=allow_unaffected, single_affected=single_affected) 
        if a != '': 
            if log:print(gene, 'de_novo mutation', '\n', tped_affected) 
            return '' 
        # individuals following inheritance model 
        tped_good = tped_affected[tped_affected['individual_geno'] == 'het'] 
        # At least 1 affected must have 1 affected parent (parent need to be het or homalt) (or have no parents). 
        To_check = tped_good.apply(lambda x:(x['paternal_pheno'] in [2] and x['paternal_geno'] in ['het','homalt']) or (x['maternal_pheno']  in [2] and x['maternal_geno'] in ['het','homalt']), axis=1) 
        tdf_check = tped_good[To_check] 
        if tdf_check.shape[0] == 0: 
            return '' 
        if not strict: 
            return ped2str(tdf_check) 
        #[strict] All affecteds must have parents with known phenotype. 
        #[strict] All affected kids must have at least 1 affected parent 
        if strict: 
            To_check = tped_good.apply(lambda x:(x['paternal_pheno'] == 2 and x['paternal_geno'] == 'het') or (x['maternal_pheno'] == 2 and x['maternal_geno'] == 'het'), axis=1) 
            tdf_check = tped_good[To_check] 
            if tdf_check.shape[0] == tped_good.shape[0]: 
                return ''
            return ped2str(tdf_check) 
        
        
    def check_x_linked_recessive(gene, df_ped, df_geno, single_affected=False): 
        tdf_ped, tdf_geno = get_ped_withAllGenoPheno(gene, df_ped, df_geno) 
        chrom = tdf_geno.iloc[0]['SV_chrom'] 
        if chrom not in ['X', 'chrX', 'ChrX']: 
            return '' 
        tped_affected = tdf_ped[(tdf_ped['phenotype'] == 2) ] 
        if tped_affected.shape[0] == 0: 
            return '' 
    
        # check affected female 
        tped_affected_female = tped_affected[(tped_affected['sex'] == 2)] 
        if tped_affected_female.shape[0] != 0: 
            affected_genos = set(tped_affected_female['individual_geno']) 
            if single_affected is False: 
                #Affected females must be HOM_ALT 
                if len(affected_genos - {'homalt'}) != 0: 
                    return '' 
            else: 
                # at least one affected female is HOM_ALT 
                if 'homalt' not in affected_genos: 
                    return ''
    
        # check affected male 
        tped_affected_male = tped_affected[(tped_affected['sex'] == 1)]  
        if tped_affected_male.shape[0] != 0: 
            affected_genos = set(tped_affected_male['individual_geno']) 
            if single_affected is False: 
                #Affected males are not HOM_REF 
                if 'homref' in affected_genos: 
                    return '' 
            else: 
                # at least one male is not HOM_REF 
                if len(affected_genos -{'homref'}) == 0: 
                    return '' 
                
        # Unaffected females are HET or HOM_REF 
        tped_unaffected = tdf_ped[(tdf_ped['phenotype'] == 1) ] 
        tped_unaffected_female = tped_unaffected[tped_unaffected['sex'] == 2] 
        tped_unaffected_male = tped_unaffected[tped_unaffected['sex'] == 1] 
        unaffected_genos = set(tped_unaffected_female['individual_geno']) 
        if 'homalt' in unaffected_genos: 
            return '' 
        #Unaffected males are HOM_REF 
        unaffected_genos = set(tped_unaffected_male['individual_geno']) 
        if len(unaffected_genos - {'homref', 'missing'}) != 0: 
            return '' 
    
        tped_good_female = tped_affected_female[(tped_affected_female['individual_geno'] == 'homalt')] 
        tped_good_male = tped_affected_male[tped_affected_male['individual_geno'].isin('het','homalt')] 
        tped_good = pd.concat([tped_good_female, tped_good_male]) 
        return ped2str(tped_good) 
    
    def check_x_linked_dominant(gene, df_ped, df_geno, single_affected=False): 
        tdf_ped, tdf_geno = get_ped_withAllGenoPheno(gene, df_ped, df_geno) 
        chrom = tdf_geno.iloc[0]['SV_chrom'] 
        if chrom not in ['X', 'chrX', 'ChrX']: 
            return '' 
        tped_affected = tdf_ped[(tdf_ped['phenotype'] == 2) ] 
        if tped_affected.shape[0] == 0: 
            return '' 
        
        # check affected female 
        tped_affected_female = tped_affected[(tped_affected['sex'] == 2)]  
        if tped_affected_female.shape[0] != 0: 
            affected_genos = set(tped_affected_female['individual_geno']) 
            if single_affected is False: 
                #Affected females must be HET 
                if not len(affected_genos - {'het'}) == 0: 
                    return '' 
            else: 
                # at least one affected female is HET 
                if 'het' not in affected_genos: 
                    return ''
                
        # check affected male 
        tped_affected_male = tped_affected[(tped_affected['sex'] == 1)]  
        if tped_affected_male.shape[0] != 0: 
            affected_genos = set(tped_affected_male['individual_geno']) 
            if single_affected is False: 
                #Affected males are HET or HOM_ALT 
                if len(affected_genos - {'het', 'homalt'}) != 0: 
                    return '' 
            else: 
                if ('het' not in affected_genos) and ('homalt' not in affected_genos): 
                    return '' 
                
        # Unaffecteds must be HOM_REF 
        tped_unaffected = tdf_ped[(tdf_ped['phenotype'] == 1) ] 
        unaffected_genos = set(tped_unaffected['individual_geno']) 
        if len(unaffected_genos - {'homref', 'missing'} != 0): 
            return '' 
    
        #girls of affected dad must be affected 
        affected_dads = set(tdf_ped[tdf_ped['paternal_pheno'] == 2]['paternal_id']) - {'0',0,'-9',-9} 
        girls_of_affected_dads = tdf_ped[tdf_ped['paternal_id'].isin(affected_dads) & tdf_ped['sex'].eq(2)] 
        girls_of_affected_dads_phenos = set(girls_of_affected_dads['phenotype']) 
        if 1 in girls_of_affected_dads_phenos: 
            return '' 
        
        #boys of affected dad must be unaffected 
        boys_of_affected_dads = tdf_ped[tdf_ped['paternal_id'].isin(affected_dads) & tdf_ped['sex'].eq(1)] 
        boys_of_affected_dads_phenos = set(boys_of_affected_dads['phenotype']) 
        if 2 in boys_of_affected_dads_phenos: 
            return '' 
        
        tped_good_female = tped_affected_female[tped_affected_female['individual_geno'] == 'het'] 
        tped_good_male = tped_affected_male[tped_affected_male['individual_geno'].isin(['het','homalt'])] 
        tped_good = pd.concat([tped_good_female, tped_good_male]) 
        
        # mothers of affected males must be het (and affected) [added in 0.19.1] 
        mothers_of_affected_males_genos = set(tped_good_male['maternal_geno']) 
        if len(mothers_of_affected_males_genos - {'het','missing'}) != 0: 
            return '' 
        
        # at least 1 parent of affected females must be het (and affected). [added in 0.19.1] 
        parent_of_affected_females_genos = set(tped_good_female['paternal_geno']) | set(tped_good_female['maternal_geno']) 
        if 'het' not in parent_of_affected_females_genos: 
            return '' 
        
        return ped2str(tped_good)
    
    # COMMANDS
    df_genes['de_novo__single_affected'] = [check_de_novo(gene, df_ped = df_ped, df_geno=df_geno, strict=True, allow_unaffected=False, single_affected=True) for gene in df_genes.index]
    df_genes['autosomal_recessive__single_affected'] = [check_autosomal_recessive(gene, df_ped = df_ped, df_geno=df_geno, strict=True, allow_unaffected=False, single_affected=True) for gene in df_genes.index]#None 
    df_genes['autosomal_dominant__single_affected'] = [check_autosomal_dominant(gene, df_ped = df_ped, df_geno=df_geno, strict=True, allow_unaffected=False, single_affected=True) for gene in df_genes.index]#None
    
    # output the column names
    df_genes_columns_out = open("df_genes_column.txt", "w")
    df_genes_columns = df_genes.columns.values.tolist()
    df_genes_columns_string = "\t".join(df_genes_columns)
    df_genes_columns_out.write(df_genes_columns_string)
    df_genes_columns_out.close()
     
    
    # create a separate dataframe to store the actually output I want with the columns that I want
    df_output = df_genes[["chrom","CDSaffected", "SVs","SVsN","de_novo__single_affected","autosomal_recessive__single_affected","autosomal_dominant__single_affected"]]
    
    print("Main Geminesque method  completed. Outputting results now.")
    
    output_filepath = "../Data/" + phenotype + "/GMN_out_" + phenotype + ".tsv"
    df_output.to_csv(output_filepath, sep='\t')
    
    print("Results outputted to " + output_filepath)


# add in OVerlapped_CDS_percent column from AnnotSV
# define function to hold process
def add_Overlapped(phenotype, AnnotSV_input):
    print("Running add_Overlapped() method")
    
    # import pandas
    import pandas as pd
    
    # open relevant files
    GMN_in = open("../Data/" + phenotype + "/GMN_out_" + phenotype + ".tsv", "r")
    AnnotSV_in = open(AnnotSV_input, "r")
    GMN_out = open("../Data/" + phenotype + "/GMN_out_" + phenotype + "_final.tsv", "w")

    # read out geneIDs from the GMN_out file
    GMN_list = GMN_in.read().split("\n")
    anno_list = AnnotSV_in.read().split("\n")
    
    # insert new information to header
    header_list = GMN_list[0].split("\t")
    header_list.insert(3,"Overlapped_CDS_percent")
    header_string = "\t".join(header_list)
    
    # write header out to file    
    GMN_out.write(header_string + "\n")
    
    # get location of anno_OCDS
    dfa = pd.read_csv(AnnotSV_input, sep = "\t")
    gene_name_index = dfa.columns.get_loc("Gene_name")
    OCDS_index = dfa.columns.get_loc("Overlapped_CDS_percent")
    
    # iterate through GMN_list and add AnnotSV info
    for i in range(1, len(GMN_list)):
        line_list = GMN_list[i].split("\t")
        var_gene_ID = line_list[0]        
    
        # use geneIDs to search through AnnotSV 
        for j in range(1, len(anno_list)):
            # only search through lines in the "split" annotation mode
            if("split" not in anno_list[j]):
                continue
            else:
                row_list = anno_list[j].split("\t")                
                anno_gene_ID = row_list[gene_name_index] # geneID
                if(var_gene_ID == anno_gene_ID):
                    anno_OCDS = row_list[OCDS_index]
                    break
        
        line_list.insert(3, anno_OCDS)
        output_string = "\t".join(line_list)
        GMN_out.write(output_string + "\n")
        
        
        
    
    # close files
    GMN_in.close()
    AnnotSV_in.close()
    GMN_out.close()

