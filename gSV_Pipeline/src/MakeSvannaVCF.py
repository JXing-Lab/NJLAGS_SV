svanna_dir="/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/SvAnna"
new_vcf = f"{svanna_dir}/SvAnna_ready.vcf"
write_vcf = open(new_vcf, "w")

### Get the header for the new vcf
header_vcf = "/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_two/modified_SURVIVOR_merged_callset.vcf"
with open(header_vcf, "r") as f:
    for line in f.readlines():
        if line.startswith("##"):
            write_vcf.write(line)
        elif line.startswith("#"):
            NEW_line = "\t".join(line.split()[0:10])
            write_vcf.write(f"{NEW_line}\n")
        else:
            break
        
### Get the content for the new vcf
lifted_bed = f"{svanna_dir}/lifted.bed"
with open(lifted_bed, "r") as f:
    for line in f.readlines():
        fields = line.split()
        CHROM, POS, END, SVTYPE = fields
        #INFO needs svtype, svlen, and end.
        SVLEN=int(END)-int(POS)
        if SVTYPE == "DEL":
            SVLEN = -abs(SVLEN)
        elif SVTYPE == "INS" and SVLEN == 0:
            SVLEN = 1
        
        ID="."
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
        
write_vcf.close()