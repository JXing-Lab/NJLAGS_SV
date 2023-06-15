"""
CNV_Prioritizer.py
> calls Prioritizer scripts to run on CNV Pipeline

Created: 2023.03.03
Updated: 2023.03.10

"""

"""
╒◖═════════════════════◗╕
        GMNP Step
╘◖═════════════════════◗╛
"""

"""
GMN_Prioritizer.py
> Takes in an annotated CNV tab-delimited file
> Calculates internal frequency of each variant and output it to a separate file
> Takes in Geminesque output
> Uses positions from Geminesque to search through the AnnotSV file
> Removes benign variants
> Checks against frequency file and removes common genotypes
"""


# VARIANT BASED PROCESSOR METHOD 
def var_process(gmn_input, annotsv_input, ped_input, output_filename):
    """
    var_process_2 ()
    > imports an annotated .tsv and Geminesque output
    > sorts through each variant and outputs relevant columns
    > incorporates data from Geminesque output
    > filter on basis of phenotype
    > modified to also accept autosomal dominant and autosomal recessive calls
    """
    print("Running GMN_Prioritizer.py var_process() method...")
    
    import pandas as pd

    # import files
    gmn_in = open(gmn_input, "r")
    annotsv_in = open(annotsv_input, "r")
    var_out = open(output_filename, "w")
    
    # read annotsv file into list and separate a header line 
    annotsv_list = annotsv_in.read().split("\n")
    header_list = annotsv_list[0].split("\t")
    
    # read Geminesque file into list
    gmn_list = gmn_in.read().split("\n")
    
    # write out output header line
    output_header_list = ["chrom", "start", "stop", "sv_type", "overlap_CDS_percent", "transcript_id", "genes", "de_novo_anno", "auto_rec_anno", "auto_dom_anno", "patho_gain_list[phen|hpo|source|coord]", "patho_loss_list[phen|hpo|source|coord]", "benign_gain_list[source|coord]", "benign_loss_list[source|coord]", "cohort_freq", "individuals"]
    output_header_string = "\t".join(output_header_list)
    var_out.write(output_header_string + "\n")
    
    # get indices
    dfa = pd.read_csv(annotsv_input, sep = "\t")
    gene_name_index = dfa.columns.get_loc("Gene_name")
    anno_mode_index = dfa.columns.get_loc("Annotation_mode")
    OCDS_index = dfa.columns.get_loc("Overlapped_CDS_percent")
    # individual indices
    indv_start_index = dfa.columns.get_loc("FORMAT") + 1
    indv_stop_index = dfa.columns.get_loc("Annotation_mode") - 1
    # pathogenic gain/loss/ins
    pgain_start_index = dfa.columns.get_loc("P_gain_phen")
    pgain_stop_index = dfa.columns.get_loc("P_gain_coord")
    ploss_start_index = dfa.columns.get_loc("P_loss_phen")
    ploss_stop_index = dfa.columns.get_loc("P_loss_coord")
    pins_start_index = dfa.columns.get_loc("P_ins_phen")
    pins_stop_index = dfa.columns.get_loc("P_ins_coord")
    # benign gain/loss/ins
    bgain_start_index = dfa.columns.get_loc("B_gain_source")
    bgain_stop_index = dfa.columns.get_loc("B_gain_coord")
    bloss_start_index = dfa.columns.get_loc("B_loss_source")
    bloss_stop_index = dfa.columns.get_loc("B_loss_coord")
    bins_start_index = dfa.columns.get_loc("B_ins_source")
    bins_stop_index = dfa.columns.get_loc("B_ins_coord")
    # transcript id
    transcript_index = dfa.columns.get_loc("Tx")
    
    
    # link individual IDs to column number
    indv_dict = {}
    for z in range(indv_start_index, indv_stop_index + 1):
        indv_dict[str(z)] = header_list[z]
    
    # calculate count of individuals
    indv_count = len(indv_dict)
    print("# of individuals: " + str(indv_count))

    # iterate through AnnotSV file
    for a in range(1, len(annotsv_list)-1):
        # retain relevant columns
        var_list = annotsv_list[a].split("\t")
        var_svid = var_list[0]
        sv_chrom = int(var_list[1])
        sv_start = int(var_list[2])
        sv_stop = int(var_list[3])
        sv_length = int(var_list[4])
        sv_type = var_list[5]
        gene_id = var_list[gene_name_index]
        de_novo_anno = ""
        anno_mode = var_list[anno_mode_index]
        transcript_id = var_list[transcript_index]
        # account for events wherein the overlap % column is empty:
        if(var_list[OCDS_index] == ""):
            overlap_CDS_percent = "N/A"
        else:
            overlap_CDS_percent = int(var_list[OCDS_index])
            
        pathogenic = "False"
        benign = "False"
        
        # GATHER INFORMATION FOR PATHOGENIC AND BENIGN ANNOTATION        
        # declare list for pathogenic gain columns
        p_gain_list = var_list[pgain_start_index:pgain_stop_index + 1] # the +1 is to ensure inclusivity
        for x in range(0, len(p_gain_list)):
            if(p_gain_list[x] == ""):
                p_gain_list[x] = "."
        p_gain_string = "|".join(p_gain_list)
                
        # declare list for pathogenic loss columns
        p_loss_list = var_list[ploss_start_index:ploss_stop_index + 1]
        for x in range(0, len(p_loss_list)):
            if(p_loss_list[x] == ""):
                p_loss_list[x] = "."
        p_loss_string = "|".join(p_loss_list)
        
        # declare list for pathogenic ins columns
        p_ins_list = var_list[pins_start_index:pins_stop_index + 1]
        for x in range(0, len(p_ins_list)):
            if(p_ins_list[x] == ""):
                p_ins_list[x] = "."
        p_ins_string = "|".join(p_ins_list)
                
        
        # declare list for benign gain columns
        b_gain_list = var_list[bgain_start_index:bgain_stop_index + 1]
        for x in range(0, len(b_gain_list)):
            if(b_gain_list[x] == ""):
                b_gain_list[x] = "."
        b_gain_string = "|".join(b_gain_list)
                
        # declare list for benign loss columns
        b_loss_list = var_list[bloss_start_index:bloss_stop_index + 1]
        for x in range(0, len(b_loss_list)):
            if(b_loss_list[x] == ""):
                b_loss_list[x] = "."
        b_loss_string = "|".join(b_loss_list)
        
        # declare list for benign loss columns
        b_ins_list = var_list[bins_start_index:bins_stop_index + 1]
        for x in range(0, len(b_ins_list)):
            if(b_ins_list[x] == ""):
                b_ins_list[x] = "."
        b_ins_string = "|".join(b_ins_list)
        
        # compile list of individuals the variant is present in 
        indvs_present = []
        
        # calculate cohort frequency
        indv_list = []
        fam_list = []
        for c in range(indv_start_index, indv_stop_index + 1):
            if var_list[c] != "0/0:0:0:0:0:0:null:0":
                indv_list.append(var_list[c])
                # isolate individual ID to determine family relation
                indv_ID = indv_dict[str(c)]
                indvs_present.append(indv_ID) # keep track of list of individuals the variant is present in
                fam_ID = indv_ID[0:4]
                if fam_ID not in fam_list:
                    fam_list.append(fam_ID)
                else:
                    continue
                
        # change addition depending on whether it's present in more than one family
        if(len(fam_list) <= 1):
            cohort_freq = str(len(indv_list)) + "/" + str(indv_count)
        else:
            cohort_freq = str(len(indv_list)) + "/" + str(indv_count)
            
        # add family info from the Geminesque file using SV_id as index
        for b in range(1, len(gmn_list)-2):
            gene_list = gmn_list[b].split("\t")
            gene_svid = gene_list[4]
            
            # compile list of individuals into a single column
            indvs_string = "|".join(indvs_present)
            
            
        
            """
            dfp = pd.read_csv(ped_input, sep = "\t")
            for z in range(0, len(indvs_present)):
                indv_index = dfp.index[dfp["PARTID"] == int(indvs_present[z])].tolist()[0]
                indv_pheno = int(dfp.at[indv_index, "Phenotype"])
                if indv_pheno == 2:
                    in_affected = True
                    break
                else:
                    continue       
            """
            # variables to test for segregation analysis
            de_novo_field = gene_list[6]
            auto_rec_field = gene_list[7]
            auto_dom_field = gene_list[8]
            
            # retain if SVID matches AND there is de novo annotation AND annotation mode is full AND is in an affected individual
            hasSegAnno = False
            if de_novo_field != "" or auto_rec_field != "" or auto_dom_field != "":
                hasSegAnno = True
            if(var_svid in gene_svid and hasSegAnno == True and anno_mode == "full" and indvs_string != ""): # also eliminate variants with no individuals
                # assign Geminesque data
                if de_novo_field != "":
                    de_novo_anno = de_novo_field
                else:
                    de_novo_anno = ""
                if auto_rec_field != "":
                    auto_rec_anno = auto_rec_field
                else:
                    auto_rec_anno = ""
                if auto_dom_field != "":
                    auto_dom_anno = auto_dom_field
                else:
                    auto_dom_anno = ""
                
                # test for affected individual
                in_affected = False
                dfp = pd.read_csv(ped_input, sep = "\t")
                for z in range(0, len(indvs_present)):
                    indv_index = dfp.index[dfp["PARTID"] == int(indvs_present[z])].tolist()[0]
                    indv_pheno = int(dfp.at[indv_index, "Phenotype"])
                    if indv_pheno == 2:
                        in_affected = True
                        break
                    else:
                        continue   
                
                # calculate number of de novo calls
                #de_novo_count = len(de_novo_field.split(";"))
                
                # test for cohort rarity in de novo
                rare_enough = False
                if auto_dom_anno != "" or auto_rec_anno != "": # de novo not a criteria here
                    rare_enough = True
                elif len(indvs_present) < 3:
                    rare_enough = True
                else:
                    rare_enough = False
                
                
                if in_affected and rare_enough: # and de_novo_count < 3
                    # write out whole list
                    output_list = [str(sv_chrom), str(sv_start), str(sv_stop), sv_type, str(overlap_CDS_percent), transcript_id, gene_id, de_novo_anno, auto_rec_anno, auto_dom_anno, p_gain_string, p_loss_string, b_gain_string, b_loss_string, cohort_freq, indvs_string]
                    output_string = "\t".join(output_list)
                    var_out.write(output_string + "\n")
                else:
                    continue
                
                    # old troubleshooting script, no longer necessary, I don't think
# =============================================================================
#                 if(sv_chrom == 9 and sv_start == 134064728 and sv_stop == 139351737):
#                     print("CHECK")
#                     print("ID_bool: " + str(var_svid in gene_svid))
#                     print("de_novo_bool: " + str(de_novo_field != ""))
#                     print("anno_bool: " + str(anno_mode == "full"))
# =============================================================================
            
                break
            else:
                continue 


# CHECK OVERLAP SUBMETHOD
def check_overlap(sv_chrom, sv_start, sv_stop, test_chrom, test_start, test_stop):
    final_outer_start = 0
    final_outer_stop = 0
    test_length = test_stop - test_start
    sv_length = sv_stop - sv_start
    
    # must first be from the same chromosome
    if(sv_chrom == test_chrom):
        # determine overlap (inner) and outer positions
        if(test_start >= sv_start): # OVERLAP START
            overlap_start = test_start # greater start is rightmost or inner start
            outer_start = sv_start # lesser start is leftmost or outer start
        else:
            overlap_start = sv_start # greater start is rightmost or inner start
            outer_start = test_start # lesser start is leftmost or outer start
        if(test_stop <= sv_stop): # OVERLAP STOP
            overlap_stop = test_stop # lesser stop is the leftmost or inner stop
            outer_stop = sv_stop # greater stop is rightmost or outer stop
        else:
            overlap_stop = sv_stop # lesser stop is the leftmost or inner stop
            outer_stop = test_stop # greater stop is the rightmost or outer stop
        
        # adjust final_outer variables
        if(final_outer_start == 0 or final_outer_start >= sv_start):
            final_start = sv_start
        else:
            final_start = final_start
        if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
            final_outer_stop = outer_stop
        else:
            final_outer_stop = final_outer_stop
        
        # calculate overlap length
        overlap_length = overlap_stop - overlap_start
        outer_length = outer_stop - outer_start
    
        # perform overlap 70% calculation
        test_70 = round(test_length * .7)
        sv_70 = round(sv_length * .7)
        if(overlap_length >= test_70):
            test_overlap = True
        else:
            test_overlap = False
    else:
        test_overlap = False

    return(test_overlap)

# CHROMOSOME VALUE CONVERTER SUBMETHOD
def chrom_value(input_list):
    """
    chrom_value
    > read in a chrom from the first column of a .vcf line
    > translate it into a numerical value for the sake of sorting
    """    
    
    if(input_list == ""):
        return(0)
    
    input_list = input_list.split("\t")
    
    chrom_val = input_list[1]
    if(chrom_val == "X"):
        chrom_val = 24
    elif(chrom_val == "Y"):
        chrom_val = 25
    else:
        chrom_val = int(chrom_val)
    return(chrom_val)
    
# VARIANT SORTER METHOD
def var_sort(var_input, var_output):
    """
    var_sort()
    > imports GMN_variants output and then sorts them
    > simple bubble sort
    """
    print("Running GMN_Prioritizer.py var_sort() method...")

    # use pandas to acquire positions
    import pandas as pd
    dfv = pd.read_csv(var_input, sep = "\t")
    chrom_index = dfv.columns.get_loc("chrom")
    start_index = dfv.columns.get_loc("start")
    stop_index = dfv.columns.get_loc("stop")

    # import files
    var_file = open(var_input, "r")
    
    # declare export file
    var_out = open(var_output, "w")
    
    # PERFORM THE SORTING (bubble sort)
    unsorted_list = var_file.read().split("\n")
    
    # set holding list
    sorted_holding_list = []    
    
    # set sorting boolean
    swapped = True
    while swapped == True:
        swapped = False
        for i in range(1, len(unsorted_list) - 1):
            # skip empty lines
            if(unsorted_list[i] == "\n"):
                continue
            
            # break up each line in unsorted_list into columns
            line1_list = unsorted_list[i].split("\t")
            line2_list = unsorted_list[i + 1].split("\t")
            
            # sort by chromosome
            if(chrom_value(unsorted_list[i]) > chrom_value(unsorted_list[i + 1])):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by start position
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[start_index]) > int(line2_list[start_index]))):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by stop position
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[start_index]) == int(line2_list[start_index])) and (int(line1_list[stop_index]) > int(line2_list[stop_index]))):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
                
    # append to sorted holding list
    for i in range(len(unsorted_list)):
         # skip any empty lines
        if(len(unsorted_list[i]) <= 1):
            continue
        
        # write to file
        sorted_holding_list.append(unsorted_list[i])
    
    # write to file
    sorted_string = "\n".join(sorted_holding_list)
    var_out.write(sorted_string)
        
    # close files
    var_file.close()
    var_out.close()

    print("var_sort() method has completed running.")

# VARIANT PRIORITIZER METHOD
def var_prioritize(input_file, output_file, disp_file, ndd_file, leniency):
    """
    var_prioritize()
    > imports sorted GMN variants and prioritizes them
    > sorts through each variant, filters/ranks
    > outputs passing variants
    """
    
    print("Running GMN_Prioritizer.py var_prioritize() method...")
    
    # import pandas dataframe for column indices
    import pandas as pd
    dfv = pd.read_csv(input_file, sep = "\t")
    
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_file + "_" + leniency, "w")
    disp_in = open(disp_file, "r")
    ndd_in = open(ndd_file, "r")
    
    # build list of variants
    var_list = file_in.read().split("\n")
    
    # build list of dispensable genes
    disp_genes = disp_in.read().split("\n")
    
    # build list of neurodevelopmental disease genes
    ndd_list = ndd_in.read().split("\n")
    ndd_genes = []
    for a in range(1, len(ndd_list)):
        ndd_line_list = ndd_list[a].split("\t")
        ndd_genes.append(ndd_line_list[0])
    
    # output header
    header_string = var_list[0] + "\t" + "prior_annotation"
    file_out.write(header_string + "\n")
    
    # iterate through list of variants and remove ones with matching benign annotation and no corresponding pathogenic
    for variant in range(1, len(var_list)):
        line_list = var_list[variant].split("\t")
        
        # define columns
        transcript_id = line_list[dfv.columns.get_loc("transcript_id")]
        chrom = line_list[dfv.columns.get_loc("chrom")]
        start = line_list[dfv.columns.get_loc("start")]
        stop = line_list[dfv.columns.get_loc("stop")]
        sv_type = line_list[dfv.columns.get_loc("sv_type")]
        overlap_CDS = line_list[dfv.columns.get_loc("overlap_CDS_percent")]
        de_novo = line_list[dfv.columns.get_loc("de_novo_anno")]
        auto_rec = line_list[dfv.columns.get_loc("auto_rec_anno")]
        auto_dom = line_list[dfv.columns.get_loc("auto_dom_anno")]
        p_gain = line_list[dfv.columns.get_loc("patho_gain_list[phen|hpo|source|coord]")]
        p_loss = line_list[dfv.columns.get_loc("patho_loss_list[phen|hpo|source|coord]")]
        b_gain = line_list[dfv.columns.get_loc("benign_gain_list[source|coord]")]
        b_loss = line_list[dfv.columns.get_loc("benign_loss_list[source|coord]")]
        cohort_freq = line_list[dfv.columns.get_loc("cohort_freq")]
        genes = line_list[dfv.columns.get_loc("genes")]
        individuals = line_list[dfv.columns.get_loc("individuals")]
        
        # FILTER ON PRESENCE OF PRIOR ANNOTATION
        gene_list = genes.split(";") 
        
        disp_tag = ""
        for gene in range(0, len(gene_list)):
            if(gene_list[gene] in disp_genes):
                disp_tag = "dispensable"
            elif(gene_list[gene] in ndd_genes):
                disp_tag = "known NDD gene"
            else:
                disp_tag = "not dispensable"
                break
        line_list.append(disp_tag)
        
        var_list[variant] = "\t".join(line_list)
        
        # FILTER ON PATHOGENIC/BENIGN ANNOTATION
        if leniency == "strict":
            # filter on presence of benign annotation
            if b_gain == ".|." and b_loss == ".|." and disp_tag != "dispensable":
                file_out.write(var_list[variant] + "\n")
            else:
                continue
        elif leniency == "stricter":     
            # filter on presence of benign annotation AND absence of pathogenic annotation
            if b_gain != ".|." and b_loss != ".|." and p_gain == ".|.|.|." and p_loss == ".|.|.|.":
                continue
            else:
                file_out.write(var_list[variant] + "\n")
        else:
            file_out.write(var_list[variant] + "\n")
                
        
    
    # close files
    file_in.close()
    
    print("var_prioritize() method has completed running.")

# ADD OVERLAP_CDS METHOD
def add_overlap_CDS(var_input, CNV_input, var_output):
    """
    add_overlap_cds()
    > imports sorted and prioritized variants
    > replaces blank overlap_cds column with the list of corresponding values from the split columns
    > exports modified variant table
    > also add transcript ID

    """
    
    print("Running GMN_Prioritizer.py add_overlap_cds() method...")
    
    import pandas as pd
    
    # open files
    file_in = open(var_input, "r")
    file_out = open(var_output, "w")
    CNV_in = open(CNV_input, "r")
    
    # build list of variants
    var_list = file_in.read().split("\n")
    # build list of variants from annotation file
    CNV_list = CNV_in.read().split("\n")
    
    # output header
    header_string = var_list[0]
    file_out.write(header_string + "\n")
    
    # acquire correct indices
    dfa = pd.read_csv(CNV_input, sep = "\t")
    dfv = pd.read_csv(var_input, sep = "\t")
    gene_name_index = dfa.columns.get_loc("Gene_name")
    anno_mode_index = dfa.columns.get_loc("Annotation_mode")
    OCDS_index = dfa.columns.get_loc("Overlapped_CDS_percent")
    transcript_index = dfa.columns.get_loc("Tx")
    
    
    # iterate through list of variants and identify full lines as keys, and append split lines to them as values
    pos_list = []
    search_start = 1
    for variant in range(1, len(var_list)-1):
        line_list = var_list[variant].split("\t")
        
        # identify full variants by position
        chrom = line_list[dfv.columns.get_loc("chrom")]
        start = line_list[dfv.columns.get_loc("start")]
        stop = line_list[dfv.columns.get_loc("stop")]
        var_pos = str(chrom) + "_" + str(start) + "_" + str(stop)
        
        # append to pos list
        pos_list.append(var_pos)
        
        # declare holding list for overlap_CDS values and transcript ids
        OCDS_list = []
        tx_list = []
        gene_list = []
        
        # use var_pos to identify the matching split lines in the CNV2_anno file
        merging_splits = False
        for a in range(search_start, len(CNV_list)-1):
            anno_list = CNV_list[a].split("\t")
                
            # identify position currently being examined
            anno_pos = "_".join((anno_list[0].split("_")[0:3]))
            
            # check annotation mode
            anno_mode = anno_list[anno_mode_index]
            
            # compare var_pos to anno_pos
            if(var_pos == anno_pos and anno_mode == "full"):
                merging_splits = True         
                continue
            elif(var_pos == anno_pos and anno_mode == "split"):
                OCDS_list.append(anno_list[OCDS_index])
                tx_list.append(anno_list[transcript_index])
                gene_list.append(anno_list[gene_name_index])
            elif(merging_splits == True):
                # export OCDS list for this variant
                line_list[dfv.columns.get_loc("stop") + 2] = "|".join(OCDS_list)
                line_list[dfv.columns.get_loc("transcript_id")] = "|".join(tx_list)
                line_list[dfv.columns.get_loc("genes")] = "|".join(gene_list)
                merging_splits = False
                break
        
        # export variant
        file_out.write("\t".join(line_list) + "\n")
                
       
        
    
    print("add_overlap_CDS() method has completed running.")

# ADD GTEx EXPRESSION DATA METHOD
def add_gtex(var_input, var_output, gtex_input):
    """
    add_gtex()
    > imports sorted and prioritized variants
    > searches all the genes listed from each variant in the gtex file and returns tissue expression data
    > exports modified variant table

    """
    
    print("Running GMN_Prioritizer.py add_gtex() method...")
    
    # use pandas to get correct indices
    import pandas as pd
    dfv = pd.read_csv(var_input, sep = "\t")
    genes_index = dfv.columns.get_loc("genes")
    
    # open files
    file_in = open(var_input, "r")
    file_out = open(var_output, "w")
    gtex_in = open(gtex_input, "r")
    
    # declare lists
    file_list = file_in.read().split("\n")
    gtex_list = gtex_in.read().split("\n")
    
    # declare header
    header_list = file_list[0].split("\t")
    header_list.append("GTEx [Amygdala | Anterior cingulate cortex (BA24) | Caudate (basal ganglia) | Cerebellar Hemisphere | Cerebellum | Cortex | Frontal Cortex (BA9) | Hippocampus | Hypothalamus | Nucleus accumbens (basal ganglia) | Putamen (basal ganglia) | Spinal cord (cervical c-1) | Substantia nigra]")
    header_string = "\t".join(header_list)
    file_out.write(header_string + "\n")
    
    # sort through each of the variants and extract all genes
    for variant in range(1, len(file_list)-1):
        var_list = file_list[variant].split("\t")
        genes = var_list[genes_index].split(";")
        
        # declare GTEx output list
        gtex_out_list = []
        
        # sort through all the genes and use them as keys to search through the GTEx data
        for gene in range(0, len(genes)):     
            
            
            # search through GTEx file
            for a in range(1, len(gtex_list)):
                
                # isolate columns for brain expression
                brain_cols = gtex_list[a].split("\t")[9:22]
                brain_cols = [round(float(num), 2) for num in brain_cols]
                brain_cols = [str(num) for num in brain_cols]
                
                # if the gene is found in the gtex_list, print out its expression data for all tissues
                if(genes[gene] in gtex_list[a]):
                    gtex_string = genes[gene] + " : " + ";".join(brain_cols)
                    gtex_out_list.append(gtex_string)
            
        # append gtex out list to variant output
        gtex_out_string = "|".join(gtex_out_list)
        var_list.append(gtex_out_string)
        var_out_string = "\t".join(var_list)
        file_out.write(var_out_string + "\n")
                
    # close files
    file_in.close()
    file_out.close()
    gtex_in.close()
    
    print("add_gtex() method has completed running.")

# ADD EXPRESSION DATA FROM ALL THREE SOURCES METHOD
def add_expression(var_input, var_output, tpm1, tpm2, tpm3):
    """
    add_expression()
    > takes in variant-based .tsv input
    > references three tissue expression databases
    > requires every variant to have at least one gene with tpm > 3 for a given brain tissue
    > exports passing variants
    """
    print("Running add_expression() method...")

    import pandas as pd
    import pickle
    
    # acquire correct indices
    dfv = pd.read_csv(var_input, sep = "\t")
    genes_index = dfv.columns.get_loc("genes")
    
    # PROCESS TPM1
    df_tpm1 = pd.read_csv(tpm1, sep = "\t")
    columns_GTEx = list(df_tpm1.columns[1:-1]) # select GTEx columns from dataframe
    columns_GTEx_brain = [e for e in columns_GTEx if "Brain" in e] # retain only the ones from brain 
    df_tpm1['tpm1_max'] = df_tpm1[columns_GTEx_brain].max(axis=1) # declare new column that collects maximum tpm value from all brain columns
    # reduce df_tpm1 to only the relevant columns
    df_tpm1.rename(columns={'Description': 'GENE'}, inplace=True)
    df_tpm1 = df_tpm1[['GENE', 'tpm1_max']]
    
    
    # PROCESS TPM2-
    df_tpm2 = pd.read_csv(tpm2, sep = "\t")
    # treat all columns as necessary for now
    columns_brainspan = list(df_tpm2.columns[2:]) 
    df_tpm2['tpm2_max'] = df_tpm2[columns_brainspan].max(axis=1)

    
    # PROCESS TPM3
    df_tpm3 = pickle.load(open(tpm3,'rb'))#165 tissues 
    df_tpm3['tpm3_max'] = df_tpm3[df_tpm3.columns[2:]].max(axis=1) 
    df_tpm3 = df_tpm3[['Gene_Name', 'tpm3_max'] + [e for e in df_tpm3.columns if "cerebral_cortex" in e]]
    # The second column is what we care about: TPM3
    df_tpm3.rename(columns={'Gene_Name': 'GENE'}, inplace = True)
    df_tpm3 = df_tpm3[['GENE', 'tpm3_max']]
    

    # declare separate dataframe to hold all values from the various tpm files
    df_tpm_1and2 = df_tpm1.merge(df_tpm2, how = 'left', left_on = "GENE", right_on ='geneSymbol')
    df_tpm_all = df_tpm_1and2.merge(df_tpm3, how = 'left', left_on = "GENE", right_on ='GENE') 
    
    
    
    # COMPARE ALL VARIANTS TO DF_TPM_ALL
    # open variant file and export file
    var_in = open(var_input, "r")
    var_list = var_in.read().split("\n")
    var_out = open(var_output, "w")
    
    # write out header list
    header_list = var_list[0].split("\t")
    header_list.append("Expression")
    var_out.write("\t".join(header_list) + "\n")
    
    # cast df_tpm_all to dict
    tpm_dict = df_tpm_all.set_index("GENE").T.to_dict("list")
    
    # iterate through variant file line by line
    for a in range(1, len(var_list)-1):
        line_list = var_list[a].split("\t")
        gene_list = line_list[genes_index].split(";")
        
        # boolean to determine if any tpm exceeds notability threshold (>3 or >5)
        is_expr = False
        
        # iterate through every gene in variant, and if 0 genes have enough expression, drop it
        for b in range(0, len(gene_list)):
            gene_id = gene_list[b]
            
            if gene_id not in tpm_dict:
                continue
            else:
                gene_tpm = tpm_dict[gene_id]
                # iterate through all the tpm information for a given gene
                for c in range(3, len(gene_tpm)):
                    tpm_value = round(float(gene_tpm[c]), 2)
                    if tpm_value >= 5:
                        is_expr = True
                        break
                    else:
                        continue
        
    
        # after having read through all the genes in the variant, decide to write out or not:
        """
        if is_expr == True:
            var_out.write(var_list[a] + "\n")
        else:
            continue
        """
        var_out.write(var_list[a] + "\t" + str(is_expr) + "\n")
          
    print("add_expression() method has completed.")



"""
╒◖═════════════════════◗╕
     StrVCTVRE Step
╘◖═════════════════════◗╛
"""
# BUILD .BED FILE METHOD
def build_bed(var_input, bed_output):
    """
    Build_StrVCTVRE.py
    > takes an annotated .tsv file and extracts positions
    > composes positions into .bed file
    > inputs .bed file into StrVCTVRE 
    """
    
    print("Recomposing annotated .tsv file into .bed file...")
    
    # open files
    var_in = open(var_input, "r")
    bed_out = open(bed_output, "w")
    
    # use pandas to acquire positions
    import pandas as pd
    dfv = pd.read_csv(var_input, sep = "\t")
    chrom_index = dfv.columns.get_loc("chrom")
    start_index = dfv.columns.get_loc("start")
    stop_index = dfv.columns.get_loc("stop")
    
    var_list = var_in.read().split("\n")
    for a in range(1, len(var_list)-1):
        line_list = var_list[a].split("\t")
        
        chrom = line_list[chrom_index]
        start = line_list[start_index]
        stop = line_list[stop_index]
        
        # convert SVTYPE to DEL or DUP
        svtype = line_list[3]
        if(svtype == "<CN3>" or svtype == "<CN4>"):
            svtype = "DUP"
        elif(svtype == "<CN0>" or svtype == "<CN1>"):
            svtype = "DEL"
        
        output_list = [chrom, start, stop, svtype]
        output_string = "\t".join(output_list)
        
        bed_out.write(output_string + "\n")
    
    # close files
    var_in.close()
    bed_out.close()
    
    print(".bed file for StrVCTVRE has been composed.")

def call_StrVCTVRE(phenotype, projects_folder):
    print("Calling StrVCTVRE from prioritizer script...")
    
    # import subprocess package
    import subprocess 

    # call StrVCTVRE
    subprocess.run(projects_folder + '/Scripts/StrVCTVRE_call_' + phenotype + '.sh')
    
    print("Prioritizer script has finished calling StrVCTVRE.") 

def add_STR(STR_input, var_input, var_output):
    print("Adding StrVCTVRE output to GMNP output...")
    import pandas as pd
    
    # open files
    STR_in = open(STR_input, "r")
    var_in = open(var_input, "r")
    var_out = open(var_output, "w")
    
    # read STR_in and var_in as pandas dataframes
    df_str = pd.read_csv(STR_in, sep = "\t", header = None)
    df_str.columns = ["chrom", "start", "stop", "svtype", "StrVCTVRE"]
    df_var = pd.read_csv(var_in, sep = "\t")
    df_var["StrVCTVRE"] = df_str["StrVCTVRE"].copy()
    
    # output
    df_var.to_csv(var_output, sep = "\t", index = False)
    
    # close files
    STR_in.close()
    var_in.close()
    var_out.close()
    
    print("StrVCTVRE output has been added to GMNP output.")



"""
╒◖═════════════════════◗╕
       SvAnna Step
╘◖═════════════════════◗╛
"""

def convert_CNV_to_SvAnna_vcf(vcf_input, vcf_output):
    print("Converting SV .vcf to SvAnna-compatible .vcf format...")
    
    # import datetime to keep track of when the vcf was updated
    from datetime import date
    date_updated = date.today().strftime("%Y%m%d")
    
    # open files
    vcf_in = open(vcf_input, "r")
    vcf_out = open(vcf_output, "w")
    
    # Change header to match SvAnna example
    header_string = """##fileformat=VCFv4.3
##fileDate=""" + date_updated + """
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##FORMAT=<ID=PI,Number=1,Type=Integer,Description='PARTID'>
##FORMAT=<ID=CA,Number=1,Type=String, Description='Caller'>
##FORMAT=<ID=CH,Number=1,Type=Integer,Description='Chromosome'>
##FORMAT=<ID=SA,Number=1,Type=Integer,Description='Start position'>
##FORMAT=<ID=SO,Number=1,Type=Integer,Description='Stop position'>
##FORMAT=<ID=QL,Number=1,Type=String,Description='Quality reported by CNV filtering'>
##FORMAT=<ID=PR,Number=1,Type=Integer,Description='Number of probes'>
##contig=<ID=1,assembly=b37,length=249250621>
##contig=<ID=2,assembly=b37,length=243199373>
##contig=<ID=3,assembly=b37,length=198022430>
##contig=<ID=4,assembly=b37,length=191154276>
##contig=<ID=5,assembly=b37,length=180915260>
##contig=<ID=6,assembly=b37,length=171115067>
##contig=<ID=7,assembly=b37,length=159138663>
##contig=<ID=8,assembly=b37,length=146364022>
##contig=<ID=9,assembly=b37,length=141213431>
##contig=<ID=10,assembly=b37,length=135534747>
##contig=<ID=11,assembly=b37,length=135006516>
##contig=<ID=12,assembly=b37,length=133851895>
##contig=<ID=13,assembly=b37,length=115169878>
##contig=<ID=14,assembly=b37,length=107349540>
##contig=<ID=15,assembly=b37,length=102531392>
##contig=<ID=16,assembly=b37,length=90354753>
##contig=<ID=17,assembly=b37,length=81195210>
##contig=<ID=18,assembly=b37,length=78077248>
##contig=<ID=19,assembly=b37,length=59128983>
##contig=<ID=20,assembly=b37,length=63025520>
##contig=<ID=21,assembly=b37,length=48129895>
##contig=<ID=22,assembly=b37,length=51304566>
##contig=<ID=MT,assembly=b37,length=16569>
##contig=<ID=X,assembly=b37,length=155270560>
##contig=<ID=Y,assembly=b37,length=59373566>
##contig=<ID=hs37d5,assembly=b37,length=35477943>\n"""

    # write out header
    vcf_out.write(header_string)
    column_list = ["#CHROM", "POS", "ID" , "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample\n"]
    column_string = "\t".join(column_list)
    vcf_out.write(column_string)
    
    # iterate through input vcf
    vcf_list = vcf_in.read().split("\n")
    for variant in range(0, len(vcf_list)-1):
        if vcf_list[variant][0] == "#": # skip past old header lines
            continue
        var_list = vcf_list[variant].split("\t")
        info_list = var_list[7].split(";")
        
        # assign variable indices
        chrom = var_list[0]
        start = var_list[1]
        stop = info_list[1].split("=")[1]
        old_alt = var_list[4]
        
        # set all quality scores to passing, since this is actually microarray data
        var_list[5] = "1000"
        
        # Convert <CNV#> tags to DUP or DEL
        if old_alt == "<CN0>" or old_alt == "<CN1>":
            new_alt = "<DEL>"
        else:
            new_alt = "<DUP>"
        var_list[4] = new_alt
        
        # Manually calculate SVLEN
        svlen = int(stop) - int(start)
        if new_alt == "<DEL>":
            svlen = -abs(svlen)
        info_list.append("SVLEN=" + str(svlen))
        var_list[7] = ";".join(info_list)
        
        # output list
        output_string = "\t".join(var_list[0:10]) + "\n"
        vcf_out.write(output_string)
    
    # close all open files
    vcf_in.close()
    
    print(".vcf has been converted to SvAnna-compatible format.")

# convert to SvAnna bed precursor
def VCF_to_BED(input_file, output_file):
    print("Converting VCF to BED file for use in SvAnna")
    
    # open output file
    bed_out = open(output_file, "w")
    
    # open input file
    with open(input_file, "r") as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            else:
                # define current line
                line_list = line.split("\t")
                chrom = line_list[0]
                start = line_list[1]
                stop = line_list[2].split("_")[2]
                # convert SVTYPE to DEL or DUP
                svtype = line_list[4]
                variant_ID = line_list[2]
                if(svtype == "<CN3>" or svtype == "<CN4>"):
                    svtype = "DUP"
                elif(svtype == "<CN0>" or svtype == "<CN1>"):
                    svtype = "DEL"
                # output new line
                output_list = [chrom, start, stop, svtype, variant_ID]
                output_string = "\t".join(output_list) 
                bed_out.write(output_string + "\n")
    
    print("VCF converted to BED")
    
# call liftOver to convert BED from b37 to hg38
def b37_to_hg38(b37_file, hg19_file, hg38_file):
    print("Converting b37 BED to hg38 bed")
    
    # import subprocess package
    import subprocess 
    
    # run b37 to hg19 liftOver
    subprocess.check_call(["/lab01/Tools/liftOver/liftOver", b37_file, "/lab01/Tools/liftOver/b37tohg19.chain",  hg19_file, " unlifted.bed"])
    # run hg19 to hg38 liftOver
    subprocess.check_call(["/lab01/Tools/liftOver/liftOver", hg19_file, "/lab01/Tools/liftOver/hg19ToHg38.over.chain", hg38_file, "unlifted.bed"])
    
    print("b37 BED converted to hg38 BED")

# convert BED back to VCF
def BED_to_VCF(input_file, output_file, ref_file):
    print("Converting BED back to VCF file for use in SvAnna")
    
    # open output file
    write_vcf = open(output_file, "w")
    
    # Get the header for the new vcf
    header_vcf = ref_file
    with open(header_vcf, "r") as f:
        for line in f.readlines():
            if line.startswith("##"):
                write_vcf.write(line)
            elif line.startswith("#"):
                NEW_line = "\t".join(line.split()[0:10])
                write_vcf.write(f"{NEW_line}\n")
            else:
                break
            
    # Get the content for the new vcf
    bed_input = input_file
    with open(bed_input, "r") as f:
        for line in f.readlines():
            fields = line.split()
            CHROM, POS, END, SVTYPE, VARIANT_ID = fields
            #INFO needs svtype, svlen, and end.
            SVLEN=int(END)-int(POS)
            if SVTYPE == "DEL":
                SVLEN = -abs(SVLEN)
            elif SVTYPE == "INS" and SVLEN == 0:
                SVLEN = 1
            
            ID=VARIANT_ID
            REF="N"
            ALT=f"<{SVTYPE}>"
            QUAL="."
            FILTER="PASS"
            INFO=f"SVTYPE={SVTYPE};END={END};SVLEN={SVLEN}"
            FORMAT="GT:DP:AD"
            SAMPLE="0/1:10:5,5" #This is random, since doesn't matter towards SvAnna calculations
            NEW_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE]
            output_string = "\t".join(NEW_line)
            write_vcf.write(output_string + "\n")
    
    
    print("BED converted back to VCF")
    
# call SvAnna
def call_SvAnna(input_file, output_folder, phenotype):
    print("Calling SvAnna...")
    
    import subprocess
    subprocess.run(["java", "-jar", "/lab01/Tools/SvAnna/svanna-cli-1.0.2-distribution/svanna-cli-1.0.2/svanna-cli-1.0.2.jar", "prioritize", "-d", "/lab01/Tools/SvAnna/svanna-cli-1.0.2-distribution/svanna-cli-1.0.2/svanna-data", "--output-format", "csv", "--vcf", input_file, "--out-dir", output_folder])
    subprocess.run(["gunzip", output_folder + "/SvAnna_input_" + phenotype + ".SVANNA.csv.gz"])
    
    print("SvAnna called and results gunzipped.")
    
# adding in SvAnna scores to GMNP output
def SvAnna_assign(SvAnna_input, GMNP_input, results_file):
    print("Exporting SvAnna output into GMNP output...")
    
    # use pandas to acquire positions
    import pandas as pd
    dfp = pd.read_csv(GMNP_input, sep = "\t")
    chrom_index = dfp.columns.get_loc("chrom")
    start_index = dfp.columns.get_loc("start")
    stop_index = dfp.columns.get_loc("stop")
    
    # open output file
    results_out = open(results_file, "w")
    
    # open GMNP to extract positions
    GMNP_in = open(GMNP_input, "r")
    GMNP_list = GMNP_in.read().split("\n")
    # open SvAnna to extract scores
    SvAnna_in = open(SvAnna_input, "r")
    SvAnna_list = SvAnna_in.read().split("\n")
    
    # write out header
    results_out.write(GMNP_list[0] + "\t" + "SvAnna" + "\n")
    
    
    
    # iterate through GMNP_list
    for a in range(1, len(GMNP_list)-1):
        var_list = GMNP_list[a].split("\t")
        gmnp_chrm = var_list[chrom_index]
        gmnp_start = var_list[start_index]
        gmnp_stop = var_list[stop_index]
        gmnp_psv = "null"
        
        # iterate through SvAnna output
        for b in range(1, len(SvAnna_list)-1):
            line_list = SvAnna_list[b].split(",")
            ID_list = line_list[3].split("_")
            svan_chrm = ID_list[0].strip("chr")
            svan_start = ID_list[1]
            svan_stop = ID_list[2]
            svan_psv = line_list[6]
            
            # test for match
            if gmnp_chrm == svan_chrm and gmnp_start == svan_start and gmnp_stop == svan_stop:
                gmnp_psv = svan_psv.strip("\n") # add psv score to corresponding row
                break
            else: 
                continue
        
        # compile output
        var_list.append(gmnp_psv)
        output_list = var_list
        output_string = "\t".join(output_list)
        results_out.write(output_string + "\n")

    print("Combined GMNP, StrVCTVRE, and SvAnna output.")
