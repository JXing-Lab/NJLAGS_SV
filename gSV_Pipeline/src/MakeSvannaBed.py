import re

original_vcf = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/modified_SURVIVOR_merged_callset.vcf"
bed_file = open("/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/SvAnna/liftOver_ready.bed", "w")

with open(original_vcf, "r") as f:
    for line in f.readlines():
        if line.startswith("#"):
            continue
        
        # Getting easily accessible fields
        fields = line.split()
        CHROM = fields[0]
        POS = fields[1]
        INFO = fields[7]
        
        # Regex to get END and SVTYPE
        END = int(re.search(r"END=(\d+);", INFO).group(1))
        SVTYPE = re.search(r"SVTYPE=(\w+);", INFO).group(1)
        
        # Writing to file
        output_list = [CHROM, str(POS), str(END), SVTYPE]
        output_string = "\t".join(output_list)
        bed_file.write(output_string + "\n")

bed_file.close()