import re

def is_match(line1, line2):
    return (line1.split()[:2] == line2.split()[:2])

def get_scores(line):
    string = line.split("\t")[7]
    match = re.search(r"ASSESS=(\d+);.*SR=(\d+)", string)
    if match:
        sr = int(match.group(2))
        assess = int(match.group(1))
        return (sr, assess)
    return (0, 0)

data_dir="/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data"
f2 = open(f"{data_dir}/merge_prep-MEI_callset/MEI_merged_duplicate_cleaned.vcf", "w")
with open(f"{data_dir}/starting_callsets/MEI_merged.vcf", "r") as f1:
    lines = f1.readlines()

    # Find location to start
    for i in range(len(lines)):
        if not lines[i].startswith("#"):
            break
        else:
            f2.write(lines[i])

    # Looking for matches and writing to file
    while i < len(lines)-1:
        # Edge case 1, keep the 55139184 deletion and not the insertion
        if lines[i].split()[1]=='55139184' and lines[i].split()[4]=='<INS:ME:ALU>':
            i+=1
            continue
        # There is another edge case (Remove one row of “chr17   1066028) but the code resolves this automatically.

        if is_match(lines[i], lines[i+1]):
            one = get_scores(lines[i])
            two = get_scores(lines[i+1])
            # print(one[0], two[0])

            if one[0] > two[0]: # Comparing SR values
                f2.write(lines[i])
                i+=1
            # elif one[0] < two[0]:
            #     continue
            elif one[0]==two[0]:
                if one[1] == 5:
                    f2.write(lines[i])
                    i+=1
        else:
            f2.write(lines[i])
        i+=1
    f2.write(lines[i])
