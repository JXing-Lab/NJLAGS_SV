"""
CNV_Builder.py
> merges starting call files by caller
> merges files by batch
> reformats file in preparation for annotation
> returns a .bed file to be annotated by the next step in CNV_Master.py

Created: 2023.01.12
Updated: 2023.02.27
"""
print("---CNV_Builder.py has begun running---")

# PHASE 1 - MERGE BY CALLER

"""
CNV_combine
> take in .csv files from cnv calls folder
> combine without merging or sorting
> output a combined .tsv file
"""    
def CNV_combine(input_folder):
    print("CNV_combine script running in folder located at: " + input_folder)
    
    # MERGE AXIOM BATCH
    # open axiom files
    axiom_qsnp_file = open(input_folder + "axiom_run3_combined_20200306_lbf_v2_bpfix_info.csv", "r")
    axiom_penncnv_file = open(input_folder + "axiom_run7_20200302_gc_all_v2_bpfix_info.csv", "r")    
    axiom_unsorted_file = open("../Data/merge_outputs/axiom_unsorted.tsv", "w")
    
    # declare merged list for both axiom callers
    axiom_merge_list = []
    
    # define column positions for axiom_qsnp_file
    for line in axiom_qsnp_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes"  or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from axiom qsnp file
    axiom_qsnp_list = [] # axiom qsnp holding list
    for line in axiom_qsnp_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc]
        qsnp_list = ["qsnp", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        qsnp_string = ":".join(qsnp_list)
            
        # define output list
        if(chrom == "23"):
            chrom = "X"
        if(chrom == "24"):
            chrom = "Y"
        output_list = [chrom, start, stop, partid, qsnp_string]
        output_string = "\t".join(output_list)
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        axiom_qsnp_list.append(output_string)
    
    # turn axiom_list into a string
    axiom_qsnp_string = "\n".join(axiom_qsnp_list)
        
    # write to file
    axiom_unsorted_file.write(axiom_qsnp_string)
    
    # separate qsnp from penncnv
    axiom_unsorted_file.write("\n")
        
    # define column positions for axiom_penncnv_file
    for line in axiom_penncnv_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from axiom penncnv file    
    axiom_penncnv_list = [] # axiom penncnv holding list
    for line in axiom_penncnv_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc].strip()
        penncnv_list = ["penncnv", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        penncnv_string = ":".join(penncnv_list)
            
        # define output list
        if(chrom == "X"):
            continue
        if(chrom == "Y"):
            continue
        output_list = [chrom, start, stop, partid, penncnv_string]
        output_string = "\t".join(output_list)
        
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append to holding list
        axiom_penncnv_list.append(output_string)
    
    # join axiom_penncnv_list to string
    axiom_penncnv_string = "\n".join(axiom_penncnv_list)
        
    # write to file
    axiom_unsorted_file.write(axiom_penncnv_string)
        
    # close intermediate file
    axiom_unsorted_file.close()
    
    # close axiom files
    axiom_qsnp_file.close()
    axiom_penncnv_file.close()
    
    
    
    # MERGE JAN 2015 ILLUMINA BATCH
    # open jan 2015 files
    jan_qsnp_file = open(input_folder + "illm_jan2015_run2_combined_20200306_lbf_v2_bpfix_info.csv", "r")
    jan_penncnv_file = open(input_folder + "psych_run4_20200302_gc_fixed_all_v2_bpfix_info.csv", "r")    
    jan_unsorted_file = open("../Data/merge_outputs/jan_unsorted.tsv", "w")
    
    # declare merged list for both axiom callers
    jan_merge_list = []
    
    # define column positions for 
    for line in jan_qsnp_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "new_stop"):
                stop_loc = i
            elif(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from jan qsnp file  
    jan_qsnp_list = [] # jan qsnp holding list
    for line in jan_qsnp_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc]
        qsnp_list = ["qsnp", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        qsnp_string = ":".join(qsnp_list)
            
        # define output list
        if(chrom == "23"):
            chrom = "X"
        if(chrom == "24"):
            chrom = "Y"
        output_list = [chrom, start, stop, partid, qsnp_string]
        output_string = "\t".join(output_list)
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # jan qsnp list
        jan_qsnp_list.append(output_string)
    
    # join list into outputabble string
    jan_qsnp_string = "\n".join(jan_qsnp_list)
    
    # write to file
    jan_unsorted_file.write(jan_qsnp_string)
        
        
    # define column positions for jan_penncnv_file
    for line in jan_penncnv_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # separate qsnp from penncnv
    jan_unsorted_file.write("\n")
    
    # import variants from axiom penncnv file    
    jan_penncnv_list = [] # jan penncnv holding list
    for line in jan_penncnv_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc].strip()
        penncnv_list = ["penncnv", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        penncnv_string = ":".join(penncnv_list)
            
        # define output list
        if(chrom == "X"):
            continue
        if(chrom == "Y"):
            continue
        output_list = [chrom, start, stop, partid, penncnv_string]
        output_string = "\t".join(output_list)
        
         # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append output string to holding list
        jan_penncnv_list.append(output_string)
        
    # join holding list into a string
    jan_penncnv_string = "\n".join(jan_penncnv_list)
     
    # write to file
    jan_unsorted_file.write(jan_penncnv_string)
        
    
    # close intermediate file
    jan_unsorted_file.close()
    
    
    # close axiom files
    jan_qsnp_file.close()
    jan_penncnv_file.close()
    
    
    
    # MERGE MAY 2017 ILLUMINA BATCH
    # open may 2017 files
    may_qsnp_file = open(input_folder + "illm_may2017_run3_combined_20200306_lbf_v2_bpfix_info.csv", "r")
    may_penncnv_file = open(input_folder + "psych_may2017_run3_20200302_gc_all_v2_bpfix_info.csv", "r")    
    may_unsorted_file = open("../Data/merge_outputs/may_unsorted.tsv", "w")
    
    # declare merged list for both may callers
    may_merge_list = []
    
    # define column positions for may_qsnp_file
    for line in may_qsnp_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from may qsnp file    
    may_qsnp_list = []
    for line in may_qsnp_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc]
        qsnp_list = ["qsnp", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        qsnp_string = ":".join(qsnp_list)
            
        # define output list
        if(chrom == "23"):
            chrom = "X"
        if(chrom == "24"):
            chrom = "Y"
        output_list = [chrom, start, stop, partid, qsnp_string]
        output_string = "\t".join(output_list)
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append output string to holding list
        may_qsnp_list.append(output_string)
    
    # join holding list into output string
    may_qsnp_string = "\n".join(may_qsnp_list)
        
    # write to file
    may_unsorted_file.write(may_qsnp_string)
    
    # separate qsnp from penncnv
    may_unsorted_file.write("\n")
        
    # define column positions for may_penncnv_file
    for line in may_penncnv_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from may penncnv file
    may_penncnv_list = []
    for line in may_penncnv_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc].strip()
        penncnv_list = ["penncnv", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        penncnv_string = ":".join(penncnv_list)
            
        # define output list
        if(chrom == "X"):
            continue
        if(chrom == "Y"):
            continue
        output_list = [chrom, start, stop, partid, penncnv_string]
        output_string = "\t".join(output_list)
        
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append output_string to holding list
        may_penncnv_list.append(output_string)
        
    # join holding list into output string
    may_penncnv_string = "\n".join(may_penncnv_list)
        
    # write to file
    may_unsorted_file.write(may_penncnv_string)
        
    
    # close intermediate file
    may_unsorted_file.close()
    
    
    # close may files
    may_qsnp_file.close()
    may_penncnv_file.close()
    
    print("CNV_combine() has completed running.\n__________________________________")


# CHROMOSOME CONVERTOR SUBMETHOD
"""
chrom_value
> read in a chrom from the first column of a .vcf line
> translate it into a numerical value for the sake of sorting
"""    
# CHROMOSOME VALUE CONVERSION METHOD
def chrom_value(input_list):
    if(input_list == ""):
        return(0)
    
    input_list = input_list.split("\t")
    
    chrom_val = input_list[0]
    if(chrom_val == "X"):
        chrom_val = 24
    elif(chrom_val == "Y"):
        chrom_val = 25
    else:
        chrom_val = int(chrom_val)
    return(chrom_val)


# CNV SORTER METHOD
"""
CNV_sort
> take in a unsorted .tsv file
> use a bubble sort to order calls by chromosome and positions
> output a sorted .tsv file
"""    
def CNV_sort(input_file, projects_folder):
    print("CNV_sort script running on input file: " + input_file)
    
    # parse input filename
    if("axiom" in input_file):
        filehandle = "axiom"
    elif("jan" in input_file):
        filehandle = "jan"
    elif("may" in input_file):
        filehandle = "may"
    elif("all" in input_file):
        filehandle = "all"
    
    # open input and output files
    unsorted_file = open(input_file, "r")
    sorted_file = open(projects_folder + "/Data/merge_outputs/" + filehandle + "_sorted.tsv", "w")
    
    # PERFORM THE SORTING (bubble sort)
    unsorted_list = unsorted_file.read().split("\n")
    
    # set holding list
    sorted_holding_list = []    
    
    # set sorting boolean
    swapped = True
    while swapped == True:
        swapped = False
        for i in range(len(unsorted_list) - 1):
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
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[1]) > int(line2_list[1]))):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by stop position
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[1]) == int(line2_list[1])) and (int(line1_list[2]) > int(line2_list[2]))):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by individual
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[1]) == int(line2_list[1])) and (int(line1_list[2]) == int(line2_list[2])) and (int(line1_list[3]) > int(line2_list[3]))):    
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
    sorted_file.write(sorted_string)
        
    # close files
    unsorted_file.close()
    sorted_file.close()
    
    print("CNV_sort() completed running.\n__________________________________")
    
    
# MERGE BY BATCH METHOD
"""
caller_merge
> take in a sorted .tsv file
> merge adjacent calls together
> output a merged .tsv file
"""    
def caller_merge(input_file, projects_folder):
    print("caller_merge script running on input file: " + input_file)
    
    # parse input filename
    if("axiom" in input_file):
        filehandle = "axiom"
    elif("jan" in input_file):
        filehandle = "jan"
    elif("may" in input_file):
        filehandle = "may"
    elif("test" in input_file):
        filehandle = "test"
    else:
        filehandle = "all"
    
    # open input and output files
    sorted_file = open(projects_folder + "/Data/merge_outputs/" + filehandle + "_sorted.tsv", "r")
    merged_file = open(projects_folder + "/Data/merge_outputs/" + filehandle + "_merged.tsv", "w")
    
    # import information in the sorted file as a list
    sorted_list = sorted_file.read().split("\n")
    info_list = [] # used for storing an increasing amount of info between loops
    qsnp_index = 0
    penncnv_index = 0
    
    # print out header
    header_list = ["chrom", "start", "stop", "partid", "caller_info[caller, chrom, start, stop, copy_number, qc, probes]\n"]
    header_string = "\t".join(header_list)
    merged_file.write(header_string)
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(0, len(sorted_list) - 1):
        
        # declare current line and next line as lists
        curr_line = sorted_list[i].split("\t")
        next_line = sorted_list[i + 1].split("\t")
        
        # skip any empty lines at the end
        if(len(next_line) <= 1 or len(curr_line) <= 1):
            continue
        
        # assign positions and determine length for current line
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_ID = curr_line[3]
        curr_info = curr_line[4]
        curr_length = curr_stop - curr_start
        
        # assign positions and determine length for next line
        next_chrom = next_line[0]
        next_start = int(next_line[1])
        next_stop = int(next_line[2])
        next_ID = next_line[3]
        next_info = next_line[4]
        next_length = next_stop - next_start
        
        # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
        overlap = False        
        
        # must first be from the same individual and chromosome
        if(curr_ID == next_ID and curr_chrom == next_chrom):
            
            # determine overlap (inner) and outer positions
            if(curr_start >= next_start): # OVERLAP START
                overlap_start = curr_start # greater start is rightmost or inner start
                outer_start = next_start # lesser start is leftmost or outer start
            else:
                overlap_start = next_start # greater start is rightmost or inner start
                outer_start = curr_start # lesser start is leftmost or outer start
            if(curr_stop <= next_stop): # OVERLAP STOP
                overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                outer_stop = next_stop # greater stop is rightmost or outer stop
            else:
                overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                outer_stop = curr_stop # greater stop is the rightmost or outer stop
            
            # calculate overlap length
            overlap_length = overlap_stop - overlap_start
            
            # perform overlap 70% calculation
            curr_70 = curr_length * .7
            next_70 = next_length * .7
            if(overlap_length >=  curr_70 and overlap_length >= next_70):
                overlap = True
            else:
                overlap = False
        
        # if chromosome and individual do not match, no overlap is possible
        else:
            overlap = False
            
        # if overlap is true, perform merge
        if(overlap == True):
            
            # PERFORM MERGE
            merge_ongoing = True # initiate merge process
            
            # declare merge list for output
            merge_list = [curr_chrom, str(outer_start), str(outer_stop), curr_ID]
            
            # append caller information
            curr_info_list = curr_info.split(":")
            if("qsnp" in curr_info_list[0]):
                qsnp_index += 1
                curr_info_list[0] = curr_info_list[0] + "_" + str(qsnp_index)
            elif("penncnv" in curr_info_list[0]):
                penncnv_index += 1
                curr_info_list[0] = curr_info_list[0] + "_" + str(penncnv_index)
            curr_info = ":".join(curr_info_list)
            info_list.append(curr_info)
            
        else:
            # if there is an ongoing merge, output it first
            if(merge_ongoing == True):                
                # append caller information for last line in merge
                curr_info_list = curr_info.split(":")
                if("qsnp" in curr_info_list[0]):
                    qsnp_index += 1
                    curr_info_list[0] = curr_info_list[0] + "_" + str(qsnp_index)
                elif("penncnv" in curr_info_list[0]):
                    penncnv_index += 1
                    curr_info_list[0] = curr_info_list[0] + "_" + str(penncnv_index)
                curr_info = ":".join(curr_info_list)
                
                # append caller info
                info_list.append(curr_info)
                info_string = "|".join(info_list)
                
                # build merge string and output
                merge_list.append(info_string)
                merge_string = "\t".join(merge_list)
                merged_file.write(merge_string + "\n")
                
                # reset ongoing merge variables
                merge_ongoing = False
                qsnp_index = 0
                penncnv_index = 0
                info_list = []
                
                
            # if there is no ongoing merge, just output the current line
            else:
                merged_file.write("\t".join(curr_line) + "\n")
        
    # account for last line if unmerged
    if(merge_ongoing == True):
        pass
    else:
        merged_file.write(sorted_list[-1])
    
            
    # close input and output files
    sorted_file.close()
    merged_file.close()
    print("caller_merge() completed running.\n__________________________________")
    

# PHASE 2 - MERGE BY BATCH

import os
import numpy as np

# METHOD TO COMPILE LIST OF INDIVIDUALS
def indv_find(input_folder, indv_output):
    print("indv_find script running...")
    
    # open input files
    axiom_in = open(input_folder + "axiom_merged.tsv", "r")
    jan_in = open(input_folder + "jan_merged.tsv", "r")
    may_in = open(input_folder + "may_merged.tsv", "r")
    
    # open output file
    indv_out = open(indv_output, "w")
    
    # declare list for holding all the individual IDs
    indv_list = []
    
    # iterate through all the axiom individuals
    for line in axiom_in:
        line_list = line.split("\t")
        
        # skip header line
        if(line_list[3] == "partid"):
            continue
        
        indv_ID = int(line_list[3])
        # append only unique individual IDs
        if(indv_ID in indv_list):
            continue
        else:
            indv_list.append(indv_ID)
            
    # iterate through all the jan individuals
    for line in jan_in:
        line_list = line.split("\t")
        
        # skip header line
        if(line_list[3] == "partid"):
            continue
        
        indv_ID = int(line_list[3])
        # append only unique individual IDs
        if(indv_ID in indv_list):
            continue
        else:
            indv_list.append(indv_ID)
    
    # iterate through all the may individuals
    for line in may_in:
        line_list = line.split("\t")
        
        # skip header line
        if(line_list[3] == "partid"):
            continue
        
        indv_ID = int(line_list[3])
        # append only unique individual IDs
        if(indv_ID in indv_list):
            continue
        else:
            indv_list.append(indv_ID)
        
    # sort indv_list
    indv_list.sort()
    
    # write indv_list out to file
    for i in range(0,len(indv_list)):
        indv_list[i] = str(indv_list[i])
    indv_out.write("\n".join(indv_list))

    # close input files
    axiom_in.close()
    jan_in.close()
    may_in.close()
    
    # close output file
    indv_out.close()
    print("indv_find script complete.")

# METHOD TO CHECK IF CALLERS DIFFER IN COPY NUMBER
def check_CNV(input_file, projects_folder):
    file_in = open(input_file, "r")
    
    if("axiom" in input_file):
        filename = "axiom_merged_2"
    if("jan" in input_file):
        filename = "jan_merged_3"
    if("may" in input_file):
        filename = "may_merged_2"
    if("all" in input_file):
        filename = "all_merged_2"
    
    file_out = open(projects_folder + "/Data/merge_outputs/" + filename + ".mismatches.txt", "w")
    
    # iterate through file line-by-line
    for line in file_in:
        
        # declare variable outside inline loop
        CNV = "ERROR"
        
        # separate line into tab-delimited list
        line_list = line.split("\t")
        
        # iterate through line column-by-column
        for i in range(0, len(line_list)):
            if("/" in line_list[i]): # upon finding a genotype filed
                if(CNV == "ERROR"): # replace placeholder with first CNV
                    CNV = line_list[i].split("/")[0]
                elif(CNV == line_list[i].split("/")[0]): # if matching, proceed as normal
                    continue
                else:
                    file_out.write(line) # if mismatched, write out to mismatch file
                    break
    
    file_in.close()
    
# METHOD TO SPLIT
def CNV_split(input_list, CNV_list, caller_mismatch):
    # separate deletion and duplication columns
    del_list = input_list.copy()
    dup_list = input_list.copy()
    
    # in a deletion variant:
    # recalculate outermost positions
    del_start = "0"
    del_stop = "0"
    for j in range(0, len(del_list)):
        # identify genotype column and delete duplications
        if("/" in del_list[j]):
            if(del_list[j].split("/")[0] == "2" or del_list[j].split("/")[0] == "3" or del_list[j].split("/")[0] == "4"):
                del_list[j] = "."
            else:
                # determine new del_start
                if(int(del_start) == 0):
                    del_start = del_list[j].split(":")[4]
                elif(int(del_start) != 0 and int(del_start) >= int(del_list[j].split(":")[4])):
                    del_start = del_list[j].split(":")[4]

                # determine new del_stop
                if(int(del_stop) == 0):
                    del_stop = del_list[j].split(":")[5]
                elif(int(del_stop) != 0 and int(del_stop) <= int(del_list[j].split(":")[5])):
                    del_stop = del_list[j].split(":")[5]
                
                
        # check if callers differ in CNV
        if("|" in del_list[j]):
            del_indv_list = del_list[j].split("|")
            qsnp_CNV = del_indv_list[0][0]
            penncnv_CNV = del_indv_list[1][0]
            if(qsnp_CNV != penncnv_CNV):
                caller_mismatch = True
            
            # assign avg del
            if(caller_mismatch == True):
                for k in range(0, len(del_indv_list)):
                    del_col_list = del_indv_list[k].split("/")
                    avg_del = round(((qsnp_CNV + penncnv_CNV)/2),1)
                    del_col_list[0] = str(avg_del)
                    del_indv_list[k] = "/".join(del_col_list)
            del_list[j] = "|".join(del_indv_list)
            
        # apply new del_start and del_stop
        del_list[1] = del_start
        del_list[2] = del_stop
    
    del_string = "\t".join(del_list)
                                        
    # in a duplication variant:
    # recalculate outermost positions
    dup_start = "0"
    dup_stop = "0"
    
    # calculate average CNV for duplications
    dup_CNVs = CNV_list.copy()
    dup_CNVs[:] = [x for x in dup_CNVs if x != "0"]
    dup_CNVs[:] = [x for x in dup_CNVs if x != "1"]
    dup_CNVs = [ int(x) for x in dup_CNVs ]
    
    
    for j in range(0, len(dup_list)):
        # identify genotype column and delete deletions
        if("/" in dup_list[j]):
            if(dup_list[j].split("/")[0] == "0" or dup_list[j].split("/")[0] == "1"):
                dup_list[j] = "."
            else:
                # determine new dup_start
                if(int(dup_start) == 0):
                    dup_start = dup_list[j].split(":")[4]
                elif(int(dup_start) != 0 and int(dup_start) >= int(dup_list[j].split(":")[4])):
                    dup_start = dup_list[j].split(":")[4]

                # determine new dup_stop
                if(int(dup_stop) == 0):
                    dup_stop = dup_list[j].split(":")[5]
                elif(int(dup_stop) != 0 and int(dup_stop) <= int(dup_list[j].split(":")[5])):
                    dup_stop = dup_list[j].split(":")[5]
            
            # set average CNV for each dup
           
            if("|" in dup_list[j]):
                dup_indv_list = dup_list[j].split("|")
                # check if callers differ in CNV
                qsnp_CNV = int(dup_indv_list[0][0])
                penncnv_CNV = int(dup_indv_list[1][0])
                if(qsnp_CNV != penncnv_CNV):
                    caller_mismatch = True
                
                # assign avg dup
                if(caller_mismatch == True):
                    for k in range(0, len(dup_indv_list)):
                        dup_col_list = dup_indv_list[k].split("/")
                        avg_dup = round(((qsnp_CNV + penncnv_CNV)/2),1)
                        dup_col_list[0] = str(avg_dup)
                        dup_indv_list[k] = "/".join(dup_col_list)
                dup_list[j] = "|".join(dup_indv_list)
                                              
                
        # apply new dup_start and dup_stop
        dup_list[1] = dup_start
        dup_list[2] = dup_stop
    
    # write out strings
    del_list[3] = "SPLIT, DELETIONS"
    del_string = "\t".join(del_list)
    dup_list[3] = "SPLIT, DUPLICATIONS"
    dup_string = "\t".join(dup_list)
    
    output_string_list = [del_string, dup_string]
    
    return(output_string_list)

# METHOD TO MERGE
def batch_merge(input_file, indv_input, output_file, projects_folder):
    print("batch_merge script running...")
    
    # open input files
    file_in = open(input_file, "r")
    indv_in = open(indv_input, "r")
    
    # open output files
    output_filename = input_file.replace("_sorted.tsv","") + output_file + ".tsv"
    file_out = open(output_filename, "w")
    caller_mismatch_out = open(projects_folder + "/Data/merge_outputs/caller_mismatches.tsv","w")
    
    error_out = open("error.txt", "w")
    
    # PRINT OUT HEADER
    # set individual list
    indv_list = indv_in.read().split("\n")
    header_list = ["#Chrom", "Start", "End", "SV_type"]
    header_list.extend(indv_list)
    header_list.append("\n")
    header_string = "\t".join(header_list)
    
    # write out header
    file_out.write(header_string)
    
    # declare ongoing merge variables
    info_list = []
    qsnp_index = 0
    penncnv_index = 0
    final_outer_start = 0
    final_outer_stop = 0
    merge_count = 0 # keeps track of how many merges are undergone
    
    # POPULATE TABLE
    # read in file as newline-separated list
    file_list = file_in.read().split("\n")
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(1, len(file_list) - 1):
        
        # declare current line and next line as lists
        curr_line = file_list[i].split("\t")
        next_line = file_list[i + 1].split("\t")
        
        # skip any empty lines at the end
        if(len(next_line) <= 1 or len(curr_line) <= 1):
            continue
        
        # assign positions and determine length for current line
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_ID = curr_line[3]
        curr_info = curr_ID + ":" + curr_line[4]
        curr_info_list = curr_info.split(":")
        curr_CNV = curr_info_list[5]
        curr_geno = [curr_CNV + "/.", ":".join(curr_info_list[0:5]), ":".join(curr_info_list[6:8])]
        curr_geno = ":".join(curr_geno)
        curr_length = curr_stop - curr_start
        
        # assign positions and determine length for next line
        next_chrom = next_line[0]
        next_start = int(next_line[1])
        next_stop = int(next_line[2])
        next_ID = next_line[3]
        next_info = next_ID + next_line[4]
        next_info_list = next_info.split(":")
        next_CNV = next_info_list[5]
        next_geno = [next_CNV + "/.", ":".join(next_info_list[0:5]), ":".join(next_info_list[6:8])]
        next_geno = ":".join(next_geno)
        next_length = next_stop - next_start
        
        
        # skip sex chromosomes
        if curr_chrom == "X" or curr_chrom == "Y":
            continue
        
        # declare variable to hold information of whether it's a mismatch and what kind
        split_string = "NO MISMATCH"
        
        # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
        overlap = False
        
        
        
        # must first be from the same chromosome
        if(curr_chrom == next_chrom):
            
            # determine overlap (inner) and outer positions
            if(curr_start >= next_start): # OVERLAP START
                overlap_start = curr_start # greater start is rightmost or inner start
                outer_start = next_start # lesser start is leftmost or outer start
            else:
                overlap_start = next_start # greater start is rightmost or inner start
                outer_start = curr_start # lesser start is leftmost or outer start
                
            if(curr_stop <= next_stop): # OVERLAP STOP
                overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                outer_stop = next_stop # greater stop is rightmost or outer stop
            else:
                overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                outer_stop = curr_stop # greater stop is the rightmost or outer stop
            
            # adjust final_outer variables
            if(final_outer_start == 0 or final_outer_start >= outer_start):
                final_outer_start = outer_start
            else:
                final_outer_start = final_outer_start
            if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
                final_outer_stop = outer_stop
            else:
                final_outer_stop = final_outer_stop
            
            # calculate overlap length
            overlap_length = overlap_stop - overlap_start
            
            # perform overlap 70% calculation
            curr_70 = curr_length * .7
            next_70 = next_length * .7
            if(overlap_length >=  curr_70 and overlap_length >= next_70):
                overlap = True
            else:
                overlap = False
        
        # if chromosomes do not match, no overlap is possible
        else:
            overlap = False
            
        # if overlap is true, perform merge
        if(overlap == True):
            
            # PERFORM MERGE
            merge_ongoing = True # initiate merge process
            merge_count += 1
            
            # declare merge list for output
            merge_list = [curr_chrom, str(final_outer_start), str(final_outer_stop), split_string]
            
            # add individual columns
            for k in range(0, len(indv_list)):
                merge_list.append(".")
            
            # append caller information
            curr_geno_list = curr_geno.split(":")
            if("qsnp" in curr_geno_list[2]):
                qsnp_index += 1
                curr_geno_list[2] = curr_geno_list[2] + "_" + str(qsnp_index)
            elif("penncnv" in curr_geno_list[2]):
                penncnv_index += 1
                curr_geno_list[2] = curr_geno_list[2] + "_" + str(penncnv_index)
            curr_geno = ":".join(curr_geno_list)
            info_list.append(curr_geno)            
            
        else:
            # if there is an ongoing merge, output it first
            if(merge_ongoing == True):   
                
                # append caller information for last line in merge
                curr_geno_list = curr_geno.split(":")
                if("qsnp" in curr_geno_list[2]):
                    qsnp_index += 1
                    curr_geno_list[2] = curr_geno_list[2] + "_" + str(qsnp_index)
                elif("penncnv" in curr_geno_list[2]):
                    penncnv_index += 1
                    curr_geno_list[2] = curr_geno_list[2] + "_" + str(penncnv_index)
                curr_geno = ":".join(curr_geno_list)
                info_list.append(curr_geno)
                
                # add info_list from current individual to the corresponding column
                # iterate through each call at this position
                for a in range(0, len(info_list)):
                    # isolate ID
                    info_ID = info_list[a].split(":")[1]
                    # add caller information to the corresponding individual column
                    for b in range(0, len(indv_list)):
                        if(indv_list[b] == info_ID):
                            if(merge_list[b + 4] == "."):
                                merge_list[b + 4] = info_list[a]
                            else:
                                merge_list[b + 4] = merge_list[b + 4] + "|" + info_list[a]
                        else:
                            continue
                
                # build merge string and output
                merged_string = "\t".join(merge_list)
                
                # SPLIT POSITIONS IF CNVS ARE MISMATCHED
                
                # determine if CNVs are mismatched                
                # declare variables outside inline loop
                inmerge_CNV = "ERROR"
                CNV_mismatch = False
                caller_mismatch = False
                
                # separate line into tab-delimited list
                line_list = merged_string.split("\t")
                
                # see how many versions of the inmerge_CNV are there
                CNV_list = []
                
                # iterate through line column-by-column
                for i in range(0, len(line_list)):
                    if("/" in line_list[i]): # upon finding a genotype filed
                        if(inmerge_CNV == "ERROR"): # replace placeholder with first CNV
                            inmerge_CNV = line_list[i].split("/")[0]
                            CNV_list.append(line_list[i].split("/")[0])
                        elif(inmerge_CNV == line_list[i].split("/")[0]): # if matching, proceed as normal
                            CNV_list.append(line_list[i].split("/")[0])
                            continue
                        else:
                            CNV_mismatch = True
                            CNV_list.append(line_list[i].split("/")[0])
                
                
                
                # if CNVs are mismatched, perform split
                if(CNV_mismatch == True):
                    # if it's only duplications, write out as normal. If it's not:
                    if("0" in CNV_list or "1" in CNV_list):        # are there deletions?    
                        if any(int(x) >= 2 for x in CNV_list): # then there are duplications and a split is necessary
                        
                            # call split method to return a given del_string and a given dup_string
                            del_string = (CNV_split(line_list, CNV_list, caller_mismatch))[0]
                            dup_string = (CNV_split(line_list, CNV_list, caller_mismatch))[1]
                            
                            # write out results of split method
                            file_out.write(del_string + "\n")
                            file_out.write(dup_string + "\n")
                            
                            # write out where callers differ on the same individual
                            if(caller_mismatch == True):
                                caller_mismatch_out.write(dup_string + "\n")
                            
                            # reset ongoing merge variables
                            merge_ongoing = False
                            caller_mismatch = False
                            qsnp_index = 0
                            penncnv_index = 0
                            final_outer_start = 0
                            final_outer_stop = 0
                            outer_start = 0
                            outer_stop = 0 
                            info_list = []      
                            
                                
                        # if it's just deletions, output as normal
                        else:
                            
                            del_list = merged_string.split("\t")
                            
                            for j in range(0, len(del_list)):
                                # check if callers differ in CNV               
                                if("|" in del_list[j]):
                                    del_indv_list = del_list[j].split("|")
                                    qsnp_CNV = int(del_indv_list[0][0])
                                    penncnv_CNV = int(del_indv_list[1][0])
                                    if(qsnp_CNV != penncnv_CNV):
                                        caller_mismatch = True
                                        
                                    # assign avg del
                                    if(caller_mismatch == True):
                                        for k in range(0, len(del_indv_list)):
                                            del_col_list = del_indv_list[k].split("/")
                                            avg_del = round(((qsnp_CNV + penncnv_CNV)/2),1)
                                            del_col_list[0] = str(avg_del)
                                            del_indv_list[k] = "/".join(del_col_list)
                                    del_list[j] = "|".join(del_indv_list)
                            
                            del_list[3] = "DELETION ONLY"
                            del_string = "\t".join(del_list)
                            file_out.write(del_string + "\n")
                        
                            # reset ongoing merge variables
                            merge_ongoing = False
                            caller_mismatch = False
                            qsnp_index = 0
                            penncnv_index = 0
                            final_outer_start = 0
                            final_outer_stop = 0
                            outer_start = 0
                            outer_stop = 0 
                            info_list = []    
                            
                            
                    # if it's just duplications, output as normal
                    else:
                        line_list = merged_string.split("\t")

                        # calculate average CNV for duplications
                        dup_CNVs = CNV_list.copy()
                        dup_CNVs[:] = [x for x in dup_CNVs if x != "0"]
                        dup_CNVs[:] = [x for x in dup_CNVs if x != "1"]
                        dup_CNVs = [ int(x) for x in dup_CNVs ]
                
                        dup_list = line_list

                        for j in range(0, len(dup_list)):
                            if("|" in dup_list[j]):
                                dup_indv_list = dup_list[j].split("|")
                                # check if callers differ in CNV
                                qsnp_CNV = int(dup_indv_list[0][0])
                                penncnv_CNV = int(dup_indv_list[1][0])
                                if(qsnp_CNV != penncnv_CNV):
                                    caller_mismatch = True
                                            
                                # write out where callers differ on the same individual
                                if(caller_mismatch == True):
                                    for k in range(0, len(dup_indv_list)):
                                        dup_col_list = dup_indv_list[k].split("/")
                                        avg_dup = round(((qsnp_CNV + penncnv_CNV)/2),1)
                                        dup_col_list[0] = str(avg_dup)
                                        dup_indv_list[k] = "/".join(dup_col_list)
                                dup_list[j] = "|".join(dup_indv_list)
                        
                        dup_list[3] = "DUPLICATION ONLY"
                        dup_string = "\t".join(dup_list)
                        
                        # write out where callers differ on the same individual
                        if(caller_mismatch == True):
                            caller_mismatch_out.write(dup_string + "\n")

                        file_out.write(dup_string + "\n")
                        
                    
                        
                        # reset ongoing merge variables
                        merge_ongoing = False
                        caller_mismatch = False
                        qsnp_index = 0
                        penncnv_index = 0
                        final_outer_start = 0
                        final_outer_stop = 0
                        outer_start = 0
                        outer_stop = 0 
                        info_list = []         
                
                # if there is no CNV mismatch, output as normal
                else:
                    file_out.write(merged_string + "\n")
                    
                    
                    # reset ongoing merge variables
                    merge_ongoing = False
                    caller_mismatch = False
                    qsnp_index = 0
                    penncnv_index = 0
                    final_outer_start = 0
                    final_outer_stop = 0
                    outer_start = 0
                    outer_stop = 0 
                    info_list = []         
                
   
            # if there is no ongoing merge, just output the current line
            else:
                unmerged_list = curr_line[0:3]
                unmerged_list.append(split_string)
                
                # add columns to unmerged_list
                for l in range(0, len(indv_list)):
                    unmerged_list.append(".")
                
                
                # add caller information to the corresponding individual column
                for b in range(0, len(indv_list)):
                    if(indv_list[b] == curr_ID):
                        unmerged_list[b + 4] = curr_geno
                    else:
                        continue
                
                # output as normal
                merged_string = "\t".join(unmerged_list)
                
                file_out.write(merged_string + "\n")
                
                # reset ongoing merge variables
                merge_ongoing = False
                qsnp_index = 0
                penncnv_index = 0
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0

    
    # close input files
    file_in.close()
    indv_in.close()
    error_out.close()
    
    # close output file
    file_out.close()
    caller_mismatch_out.close()
    
    
    print("batch_merge script complete.")
    
    # return count of merges
    return(merge_count)
    
# METHOD TO REMERGE WITHOUT CHANGING FORMAT
def CNV_remerge(input_file, output_label, cycle_count, error_file):
    """
    > cycles in the output from batch_merge()
    > retains format but performs same merge operations
    > continues merging until no more merges can be performed
    """
    # keep track of how many cycles CNV_remerge has undergone
    print("CNV_remerge script is running cycle #" + str(cycle_count) + "...")
    
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_label, "w")
    error_out = open(error_file, "w")
    
    # REPEAT OVERLAP CHECK
    # read in file as newline-separated list
    file_list = file_in.read().split("\n")
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    info_list = []
    final_outer_start = 0
    final_outer_stop = 0
    remerges_needed = 0
    
    # write out header
    file_out.write(file_list[0] + "\n")
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(1, len(file_list) - 1):
        
        # declare current line and next line as lists
        curr_line = file_list[i].split("\t")
        next_line = file_list[i + 1].split("\t")
        
        # skip any empty lines at the end
        if(len(next_line) <= 1 or len(curr_line) <= 1):
            continue
        
        # assign positions and determine length for current line
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_length = curr_stop - curr_start
        curr_type = curr_line[3]
        
        # assign positions and determine length for next line
        next_chrom = next_line[0]
        next_start = int(next_line[1])
        next_stop = int(next_line[2])
        next_length = next_stop - next_start
        next_type = next_line[3]
        
        # skip sex chromosomes
        if curr_chrom == "X" or curr_chrom == "Y":
            continue
        
        # declare variable to hold information of whether it's a mismatch and what kind
        split_string = "NO MISMATCH"
        
        # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
        overlap = False
        
        # must first be from the same chromosome AND be the same type of CNV!
        if(curr_chrom == next_chrom and curr_type == next_type):
            
            # determine overlap (inner) and outer positions
            if(curr_start >= next_start): # OVERLAP START
                overlap_start = curr_start # greater start is rightmost or inner start
                outer_start = next_start # lesser start is leftmost or outer start
            else:
                overlap_start = next_start # greater start is rightmost or inner start
                outer_start = curr_start # lesser start is leftmost or outer start
                
            if(curr_stop <= next_stop): # OVERLAP STOP
                overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                outer_stop = next_stop # greater stop is rightmost or outer stop
            else:
                overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                outer_stop = curr_stop # greater stop is the rightmost or outer stop
            
            # adjust final_outer variables
            if(final_outer_start == 0 or final_outer_start >= outer_start):
                final_outer_start = outer_start
            else:
                final_outer_start = final_outer_start
            if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
                final_outer_stop = outer_stop
            else:
                final_outer_stop = final_outer_stop
            
            # calculate overlap length
            overlap_length = overlap_stop - overlap_start
            
            # perform overlap 70% calculation
            curr_70 = curr_length * .7
            next_70 = next_length * .7
            if(overlap_length >=  curr_70 and overlap_length >= next_70):
                overlap = True
            else:
                overlap = False
        # if chromosomes do not match, no overlap is possible
        else:
            overlap = False
            
            
        # diagnostic step to troubleshoot merge of 152069090
        if(curr_stop >= 151840000 and curr_stop <= 152070000):
            error_out.write(str(curr_chrom) + "\t" + str(curr_start)  + "\t" +  str(curr_stop) + ", the overlap condition is " + str(overlap) + "\n")
            error_out.write("\t".join(next_line) + "\n")
            error_out.write("curr_start = " + str(curr_start) + "\n")
            error_out.write("curr_stop = " + str(curr_stop) + "\n")
            error_out.write("curr_length = " + str(curr_length) + "\n")
            error_out.write("curr_70 = " + str(curr_70) + "\n")
            error_out.write("next_start = " + str(next_start) + "\n")
            error_out.write("next_stop = " + str(next_stop) + "\n")
            error_out.write("next_length = " + str(next_length) + "\n")
            error_out.write("next_70 = " + str(next_70) + "\n")
            error_out.write("overlap_start = " + str(overlap_start) + "\n")
            error_out.write("overlap_stop = " + str(overlap_stop) + "\n")
            error_out.write("overlap_length = " + str(overlap_length) + "\n")
            # calculate overlap percentages
            if(curr_length >= overlap_length):
                overlap_with_curr = overlap_length/curr_length
            else:
                overlap_with_curr = curr_length/overlap_length
            if(next_length >= overlap_length):
                overlap_with_next = overlap_length/next_length
            else:
                overlap_with_next = next_length/overlap_length
            error_out.write("overlap_with_curr = " + str(overlap_with_curr) + "\n")
            error_out.write("overlap_with_next = " + str(overlap_with_next) + "\n\n")

        # IF OVERLAP IS TRUE, PERFORM MERGE
        if(overlap == True):
            remerges_needed += 1
            
            # perform merge
            merge_ongoing = True # initiate merge process
            
            # declare merge list for output
            merge_list = [curr_chrom, str(final_outer_start), str(final_outer_stop), curr_type]
            
            # sort through curr_line and next_line to see where they differ on individual data
            for A in range(4, len(curr_line)):
                if curr_line[A] != next_line[A]:
                    if curr_line[A] == ".": # overwrite null fields with data
                        merge_list.append(next_line[A])
                    elif next_line[A] == ".":
                        merge_list.append(curr_line[A])
                    else:
                        chimaera_string = curr_line[A] + "|" + next_line[A]
                        merge_list.append(chimaera_string)
                else:
                    merge_list.append(curr_line[A]) # copy over the matching field
                
        else: # IF OVERLAP IS FALSE, EITHER OUTPUT ONGOING MERGE OR CURRENT LINE
            
            # if there is an ongoing merge, output it
            if(merge_ongoing == True):   
                # build merge string and output
                merged_string = "\t".join(merge_list) + "\n"
                file_out.write(merged_string)

                # reset ongoing merge variables
                merge_ongoing = False
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0 
                
            else:
                file_out.write(file_list[i] + "\n")
                
                # reset ongoing merge variables
                merge_ongoing = False
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0

            
    print("Remerges needed: " + str(remerges_needed))
            
    
    # close files
    file_in.close()
    file_out.close()
    error_out.close()
    
    print("CNV_remerge has finished running cycle #" + str(cycle_count))

# METHOD TO PERFORM AN ADJACENCY-AGNOSTIC MERGE
def AJA_remerge(input_file, output_label, cycle_count, error_file):
    """
    > sorts though remerge output
    > takes in every position, then searches every OTHER position for a match
    > when a match is found for Variant A with Variant B:
        > writes out merged line in place of Variant A
        > marks Variant B to be skipped when printing non-overlapping lines
    """
    # keep track of how many cycles CNV_remerge has undergone
    print("AJA_remerge script is running cycle #" + str(cycle_count) + "...")
    
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_label, "w")
    error_out = open(error_file, "w")
    
    # REPEAT OVERLAP CHECK
    # read in file as newline-separated list
    file_list = file_in.read().split("\n")
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    info_list = []
    final_outer_start = 0
    final_outer_stop = 0
    remerges_needed = 0
    previously_overlapped = [] # list to store positions that have already been merged and do not need to be listed again
    
    # write out header
    file_out.write(file_list[0] + "\n")
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(1, len(file_list) - 1):
        
        # declare current line and associated variables
        curr_line = file_list[i].split("\t")
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_length = curr_stop - curr_start
        curr_type = curr_line[3]
        
        overlap = False
        
        # compare current line to every other line
        for j in range(1, len(file_list) -1):
            # skip the curr_line
            if file_list[j] == file_list[i]:
                continue
            else:
                # declare variant B's associated variables
                next_line = file_list[j].split("\t")
                next_chrom = next_line[0]
                next_start = int(next_line[1])
                next_stop = int(next_line[2])
                next_length = next_stop - next_start
                next_type = next_line[3]
                
                # declare positions to test for skips
                position_to_test = [str(curr_chrom), str(curr_start), str(curr_stop), curr_type]
                position_to_test_string = "\t".join(position_to_test)
                
                # skip any empty lines at the end
                if(len(next_line) <= 1 or len(curr_line) <= 1):
                    continue
            
                # skip sex chromosomes
                if curr_chrom == "X" or curr_chrom == "Y":
                    continue
                
                # skip positions in previously_overlapped
                if position_to_test_string in previously_overlapped:
                    continue
                
                # declare variable to hold information of whether it's a mismatch and what kind
                split_string = "NO MISMATCH"
                
                # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
                # must first be from the same chromosome AND be the same type of CNV!
                if(curr_chrom == next_chrom and curr_type == next_type):
                    
                    # determine overlap (inner) and outer positions
                    if(curr_start >= next_start): # OVERLAP START
                        overlap_start = curr_start # greater start is rightmost or inner start
                        outer_start = next_start # lesser start is leftmost or outer start
                    else:
                        overlap_start = next_start # greater start is rightmost or inner start
                        outer_start = curr_start # lesser start is leftmost or outer start
                        
                    if(curr_stop <= next_stop): # OVERLAP STOP
                        overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                        outer_stop = next_stop # greater stop is rightmost or outer stop
                    else:
                        overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                        outer_stop = curr_stop # greater stop is the rightmost or outer stop
                    
                    # adjust final_outer variables
                    if(final_outer_start == 0 or final_outer_start >= outer_start):
                        final_outer_start = outer_start
                    else:
                        final_outer_start = final_outer_start
                    if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
                        final_outer_stop = outer_stop
                    else:
                        final_outer_stop = final_outer_stop
                    
                    # calculate overlap length
                    overlap_length = overlap_stop - overlap_start
                    
                    # perform overlap 70% calculation
                    curr_70 = curr_length * .7
                    next_70 = next_length * .7
                    if(overlap_length >=  curr_70 and overlap_length >= next_70):
                        overlap = True
                    else:
                        overlap = False
                # if chromosomes do not match, no overlap is possible
                else:
                    overlap = False
        
                # IF OVERLAP IS TRUE, PERFORM MERGE AND OUTPUT
                if(overlap == True):
                    remerges_needed += 1
                    
                    # declare merge list for output
                    merge_list = [str(curr_chrom), str(final_outer_start), str(final_outer_stop), curr_type]
                    
                    # sort through curr_line and next_line to see where they differ on individual data
                    for A in range(4, len(curr_line)):
                        if curr_line[A] != next_line[A]:
                            if curr_line[A] == ".": # overwrite null fields with data
                                merge_list.append(next_line[A])
                            elif next_line[A] == ".":
                                merge_list.append(curr_line[A])
                            else:
                                chimaera_string = curr_line[A] + "|" + next_line[A]
                                merge_list.append(chimaera_string)
                        else:
                            merge_list.append(curr_line[A]) # copy over the matching field
                    
                    # write out merged line in place of current line
                    file_out.write("\t".join(merge_list) + "\n")
                    position_to_skip = [str(next_chrom), str(next_start), str(next_stop), next_type]
                    position_to_skip_string = "\t".join(position_to_skip)
                    previously_overlapped.append(position_to_skip_string)
            
                    # diagnostic
                    error_out.write(str(curr_chrom) + "\t" + str(curr_start)  + "\t" +  str(curr_stop) + ", the overlap condition is " + str(overlap) + "\n")
                    error_out.write("\t".join(next_line) + "\n")
                    error_out.write("curr_start = " + str(curr_start) + "\n")
                    error_out.write("curr_stop = " + str(curr_stop) + "\n")
                    error_out.write("curr_length = " + str(curr_length) + "\n")
                    error_out.write("curr_70 = " + str(curr_70) + "\n")
                    error_out.write("next_start = " + str(next_start) + "\n")
                    error_out.write("next_stop = " + str(next_stop) + "\n")
                    error_out.write("next_length = " + str(next_length) + "\n")
                    error_out.write("next_70 = " + str(next_70) + "\n")
                    error_out.write("overlap_start = " + str(overlap_start) + "\n")
                    error_out.write("overlap_stop = " + str(overlap_stop) + "\n")
                    error_out.write("overlap_length = " + str(overlap_length) + "\n")
                    # calculate overlap percentages
                    if(curr_length >= overlap_length):
                        overlap_with_curr = overlap_length/curr_length
                    else:
                        overlap_with_curr = curr_length/overlap_length
                    if(next_length >= overlap_length):
                        overlap_with_next = overlap_length/next_length
                    else:
                        overlap_with_next = next_length/overlap_length
                    error_out.write("overlap_with_curr = " + str(overlap_with_curr) + "\n")
                    error_out.write("overlap_with_next = " + str(overlap_with_next) + "\n\n")
            
                    # reset merge variables
                    final_outer_start = 0
                    final_outer_stop = 0
                    outer_start = 0
                    outer_stop = 0
                    
                    # break the j loop
                    break
                else:
                    # reset merge variables
                    final_outer_start = 0
                    final_outer_stop = 0
                    outer_start = 0
                    outer_stop = 0
                    continue
                
            # reset merge variables
            final_outer_start = 0
            final_outer_stop = 0
            outer_start = 0
            outer_stop = 0
        
        # if it has progressed through every other line and overlap is still false, then no known overlaps
        else:
            if position_to_test_string in previously_overlapped:
                continue
            
                # reset merge variables
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0
                
            else:
                file_out.write(file_list[i] + "\n")
                
                # reset merge variables
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0
            
    print("Remerges needed: " + str(remerges_needed))
            
    
    # close files
    file_in.close()
    file_out.close()
    error_out.close()
    
    return(remerges_needed)
    
    print("AJA_remerge has finished running cycle #" + str(cycle_count))



# PHASE 3 - REFORMATTING
"""
"""
def convert_to_vcf(input_file, output_file, fasta_file):
    print("VCFize method running...")
    

    # open input and output files
    file_in = open(input_file, "r")
    file_out = open(output_file, "w")
    
    # declare vcf specification
    file_out.write("##fileformat=VCFv4.3\n")
    file_out.write("##fileDate=20220112\n")
    file_out.write("##INFO=<ID=CNV_TYPE,Number=1,Type=String,Description='Type of CNV'>\n")
    file_out.write("##INFO=<ID=END,Number=1,Type=String,Description='Ending position'>\n")
    
    ALT_list = ["##ALT=<ID=CNV,Description='Copy Number Polymorphism'>",
    "##ALT=<ID=CN9,Description='Copy number allele: 9 copies'>",
    "##ALT=<ID=CN8,Description='Copy number allele: 8 copies'>",
    "##ALT=<ID=CN7,Description='Copy number allele: 7 copies'>",
    "##ALT=<ID=CN6,Description='Copy number allele: 6 copies'>",
    "##ALT=<ID=CN5,Description='Copy number allele: 5 copies'>",
    "##ALT=<ID=CN4,Description='Copy number allele: 4 copies'>",
    "##ALT=<ID=CN3,Description='Copy number allele: 3 copies'>",
    "##ALT=<ID=CN2,Description='Copy number allele: 2 copies'>",
    "##ALT=<ID=CN1,Description='Copy number allele: 1 copies'>",
    "##ALT=<ID=CN0,Description='Copy number allele: 0 copies'>"]
    
    contigs_list = ["##contig=<ID=1,assembly=b37,length=249250621>",
    "##contig=<ID=2,assembly=b37,length=243199373>",
    "##contig=<ID=3,assembly=b37,length=198022430>",
    "##contig=<ID=4,assembly=b37,length=191154276>",
    "##contig=<ID=5,assembly=b37,length=180915260>",
    "##contig=<ID=6,assembly=b37,length=171115067>",
    "##contig=<ID=7,assembly=b37,length=159138663>",
    "##contig=<ID=8,assembly=b37,length=146364022>",
    "##contig=<ID=9,assembly=b37,length=141213431>",
    "##contig=<ID=10,assembly=b37,length=135534747>",
    "##contig=<ID=11,assembly=b37,length=135006516>",
    "##contig=<ID=12,assembly=b37,length=133851895>",
    "##contig=<ID=13,assembly=b37,length=115169878>",
    "##contig=<ID=14,assembly=b37,length=107349540>",
    "##contig=<ID=15,assembly=b37,length=102531392>",
    "##contig=<ID=16,assembly=b37,length=90354753>",
    "##contig=<ID=17,assembly=b37,length=81195210>",
    "##contig=<ID=18,assembly=b37,length=78077248>",
    "##contig=<ID=19,assembly=b37,length=59128983>",
    "##contig=<ID=20,assembly=b37,length=63025520>",
    "##contig=<ID=21,assembly=b37,length=48129895>",
    "##contig=<ID=22,assembly=b37,length=51304566>",
    "##contig=<ID=MT,assembly=b37,length=16569>",
    "##contig=<ID=X,assembly=b37,length=155270560>",
    "##contig=<ID=Y,assembly=b37,length=59373566>",
    "##contig=<ID=hs37d5,assembly=b37,length=35477943>"]
    
    for a in range(0, len(contigs_list)):
        file_out.write(contigs_list[a] + "\n")
    
    file_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype as total number of copy variants'>\n")
    file_out.write("##FORMAT=<ID=PI,Number=1,Type=Integer,Description='PARTID'>\n")
    file_out.write("##FORMAT=<ID=CA,Number=1,Type=String, Description='Caller'>\n")
    file_out.write("##FORMAT=<ID=CH,Number=1,Type=Integer,Description='Chromosome'>\n")
    file_out.write("##FORMAT=<ID=SA,Number=1,Type=Integer,Description='Start position'>\n")
    file_out.write("##FORMAT=<ID=SO,Number=1,Type=Integer,Description='Stop position'>\n")
    file_out.write("##FORMAT=<ID=QL,Number=1,Type=String,Description='Quality reported by CNV filtering'>\n")
    file_out.write("##FORMAT=<ID=PR,Number=1,Type=Integer,Description='Number of probes'>\n")
    
    for a in range(0, len(ALT_list)):
        file_out.write(ALT_list[a] + "\n")
        
    # separate file into lines
    file_list = file_in.read().split("\n")
    
    # change the header
    header_list = file_list[0].split("\t")
    header_list[0] = "#CHROM"
    header_list[1] = "POS"
    header_list[2] = "ID"
    header_list[3] = "REF"
    header_list = pos_insert("ALT",header_list, 5)
    header_list = pos_insert("QUAL",header_list, 6)
    header_list = pos_insert("FILTER",header_list, 7)
    header_list = pos_insert("INFO",header_list, 8)
    header_list = pos_insert("FORMAT",header_list, 9)
    
    # remove trailing tab
    header_list.pop()
    
    header_string = "\t".join(header_list)
    file_out.write(header_string + "\n")
    
    # iterate through input line by line
    for a in range(1, len(file_list)):

        # separate each line into tab-delimited columns
        line_list = file_list[a].split("\t")
        
        # skip empty lines
        if len(line_list) < 2:
            continue
        
        # find CNV genotype to determine del or dup
        for b in range(0, len(line_list)):
            if(":" in line_list[b]):
                CNV_number = line_list[b][0]
                
        # find information for REF allele
        start_pos = int(line_list[1])
        stop_pos = int(line_list[2])
        ref_allele = get_fasta(start_pos, stop_pos, fasta_file)
        
        # find information for ALT allele
        if(line_list[3] == "DELETION ONLY" or line_list[3] == "SPLIT, DELETIONS"):
            alt_allele = "CNV"
        elif(line_list[3] == "DUPLICATION ONLY" or line_list[3] == "SPLIT, DUPLICATIONS"):
            alt_allele = "CNV"
        else:
            if(CNV_number == "0" or CNV_number == "1"):
                alt_allele = "CNV"
            else:
                alt_allele = "CNV"
        
        # fill in empty columns
        line_list[2] = "chr" + line_list[0] + "_" + str(start_pos) + "_" + str(stop_pos)
        line_list[3] = ref_allele # for REF field
        line_list = pos_insert("<CN" + CNV_number + ">", line_list, 5) # for ALT field
        line_list = pos_insert("NULL", line_list, 6) # for QUAL field
        line_list = pos_insert("PASS", line_list, 7) # for FILTER field
        line_list = pos_insert("SVTYPE=" + alt_allele + ";END=" + str(stop_pos), line_list, 8) # for INFO field
        line_list = pos_insert("GT:PI:CA:CH:SA:SO:QL:PR", line_list, 9) # for FORMAT field
        
        # replace alternative genotypes with VCF-correct format
        for c in range(9, len(line_list)):
            if(":" in line_list[c]):
                geno_list = line_list[c].split(":")
                # translate CNV count into genotype
                if(geno_list[0] == "0/."): # homozygous deletion
                    geno_list[0] = "1/1"
                if(geno_list[0] == "1/."): # heterozygous deletion
                    geno_list[0] = "0/1"
                if(geno_list[0] == "3/."): # heterozygous insertion
                    geno_list[0] = "0/1"
                if(geno_list[0] == "4/."): # homozygous insertion
                    geno_list[0] = "1/1"
                line_list[c] = ":".join(geno_list)
                if("|" in line_list[c]):
                    line_list[c] = line_list[c].split("|")[0]
            elif(line_list[c] == "."):
                line_list[c] = "0/0:0:0:0:0:0:null:0"
        
        
        
        # write out line
        line_string = "\t".join(line_list)
        file_out.write(line_string + "\n")
    
    # close input and output files
    file_in.close()
    file_out.close()
    
    print("VCFizing complete.")
    
# INSERT POSITION SUBROUTINE
def pos_insert(x, n_list, pos):
    return n_list[:pos-1]+[x]+n_list[pos-1:]

# FASTA SUBROUTINE
def get_fasta(start_pos, stop_pos, fasta_file):
    # import relevant packages
    import pysam
    
    # query FASTA file    
    genome = pysam.Fastafile(fasta_file)
    sequence = genome.fetch("1", start_pos, start_pos+1)
    return(sequence)
    


# DECLARE GLOBAL VARIABLES
#projects_folder = "/home/rohan/CNV_Project/2023/CNV_Pipeline"
#print("Projects folder located at: " + projects_folder)
"""

# PHASE 1 COMMANDS
# merge starting call files by caller
CNV_combine(projects_folder + "/Data/starting_call_files/")
# sort and merge by batch
CNV_sort(projects_folder + "/Data/merge_outputs/axiom_unsorted.tsv")
caller_merge(projects_folder + "/Data/merge_outputs/axiom_sorted.tsv")
CNV_sort(projects_folder + "/Data/merge_outputs/jan_unsorted.tsv")
caller_merge(projects_folder + "/Data/merge_outputs/jan_sorted.tsv")
CNV_sort(projects_folder + "/Data/merge_outputs/may_unsorted.tsv")
caller_merge(projects_folder + "/Data/merge_outputs/may_sorted.tsv")

# merge to write and form all_sorted
axiom_in = open(projects_folder + "/Data/merge_outputs/axiom_unsorted.tsv", "r")
jan_in = open(projects_folder + "/Data/merge_outputs/jan_unsorted.tsv", "r")
may_in = open(projects_folder + "/Data/merge_outputs/may_unsorted.tsv", "r")
all_out = open(projects_folder + "/Data/merge_outputs/all_unsorted.tsv", "w")
# read in all files as stri
axiom_string = axiom_in.read()
jan_string = jan_in.read()
may_string = may_in.read()
all_string = axiom_string + "\n" + jan_string + "\n" + may_string
# write out combined file, merge, and sort
all_out.write(all_string)
CNV_sort(projects_folder + "/Data/merge_outputs/all_unsorted.tsv")
caller_merge(projects_folder + "/Data/merge_outputs/all_sorted.tsv") # this produces the all_merged.tsv file
# close open files
axiom_in.close()
jan_in.close()
may_in.close()
all_out.close()


# PHASE 2 COMMANDS
batch_merge(projects_folder + "/Data/merge_outputs/axiom_sorted.tsv", "_merge_intermediary_1")
batch_merge(projects_folder + "/Data/merge_outputs/jan_sorted.tsv", "_merge_intermediary_1")
batch_merge(projects_folder + "/Data/merge_outputs/may_sorted.tsv", "_merge_intermediary_1")

# TO DO: Fix merging error
batch_merge(projects_folder + "/Data/merge_outputs/all_sorted.tsv", "_merge_intermediary_1")

# remerge consecutive positions
CNV_remerge(projects_folder + "/Data/merge_outputs/all_merge_intermediary_1.tsv", projects_folder + "/Data/merge_outputs/all_merge_intermediary_2.tsv", 1, projects_folder + "/Data/error_reports/error1.txt") # returns 16 remerges
CNV_remerge(projects_folder + "/Data/merge_outputs/all_merge_intermediary_2.tsv", projects_folder + "/Data/merge_outputs/all_merge_intermediary_3.tsv", 2, projects_folder + "/Data/error_reports/error2.txt") # returns 0 remerges

# adjacency-agnostic murge
remerges_required = 99999
file_counter = 2
cycle_count = 1

# remerger loop
while remerges_required > 0:
    remerge_input = projects_folder + "/Data/merge_outputs/all_merge_intermediary_" + str(file_counter) + ".tsv"
    remerge_output = projects_folder + "/Data/merge_outputs/all_merge_intermediary_" + str(file_counter + 1) + ".tsv"
    error_output = projects_folder + "/Data/error_reports/error_" + str(file_counter + 1) + ".txt"
    remerges_required = AJA_remerge(remerge_input, remerge_output, cycle_count, error_output) # returns 120 remerges
    file_counter += 1
    cycle_count += 1

    
# PHASE 3 COMMANDS
# convert to vcf
convert_to_vcf(projects_folder + "/Data/merge_outputs/all_merge_final.tsv", projects_folder + "/Data/CNV.vcf", projects_folder + "/Data/human_g1k_v37.fasta")
"""