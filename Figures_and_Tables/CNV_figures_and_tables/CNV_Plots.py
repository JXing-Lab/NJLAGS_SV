"""
CNV_Plots.py

"""
# PLOTS

def CN_counts_boxplot(data_file, filter_status, linearity):
    
    # import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon
    from matplotlib import figure
    import pandas as pd
    
    
    # import data
    df = pd.read_csv(data_file, sep = "\t")
    CN0 = df["CN0"].tolist()
    CN1 = df["CN1"].tolist()
    CN3 = df["CN3"].tolist()
    CN4 = df["CN4"].tolist()
    
    # push data to plot
    CN_categories = [CN0, CN1, CN3, CN4]
    
    # define and label axes
    fig, ax = plt.subplots()
    pos = np.arange(len(CN_categories)) + 1
    bp = ax.boxplot(CN_categories, sym='k+', positions=pos,
                    notch=True)
    ax.set_xlabel('CN categories', fontsize = 20)
    ax.set_ylabel('Number of CNVs', fontsize = 20)
    
    # adjust plot to logarithmic scale
    if linearity == "logscale":
        plt.yscale("log")
    
    # set other plot parameters
    plt.title('CNV counts by CN category, ' + filter_status + ", " + linearity)
    plt.setp(bp['whiskers'], color='k', linestyle='-')
    plt.setp(bp['fliers'], markersize=10.0)
    plt.xticks([1, 2, 3, 4], ["CN0", "CN1", "CN3", "CN4"])
    plt.show()
    
    fig.set_size_inches(12, 8)
    fig.savefig("../Results/Plots/CN_counts_" + filter_status + "_" + linearity + ".png")
    fig.savefig("../Results/Plots/CN_counts_" + filter_status + "_" + linearity + ".pdf", format="pdf", bbox_inches="tight") 
    
    
def CN_counter(input_file, output_file):
    print("Running CN_counter() method from CNV_Plots.py...\n")
    
    import pandas as pd
    import numpy as np
    
    # set up dataframe
    dfv = pd.read_csv(input_file, skiprows = 49, sep = "\t")
    partid_list = dfv.columns.tolist()[9:]
    dfi = pd.DataFrame(partid_list, columns = ["PARTID"])
    
    #dfi = dfi.astype({'PARTID':'int'})
    dfi["CN0"] = 0
    dfi["CN1"] = 0
    dfi["CN3"] = 0
    dfi["CN4"] = 0
    dfi.set_index("PARTID")

    # import vcf
    vcf_in = open(input_file, "r")
    vcf_list = vcf_in.read().split("\n")

    # VCF header list
    header_list = vcf_list[49].split("\t")
    
    # iterate through vcf list, incrementing corresponding individual as variant is found
    for b in range(50, len(vcf_list) - 1):
        var_list = vcf_list[b].split("\t")
        CN = var_list[4].strip(">").strip("<") # assign CN number for comparison    
        # iterate through columns until a nonzero genotype is found
        for c in range(9, len(var_list)):
            geno_list = var_list[c].split(":")
            genotype = geno_list[0]
            if genotype != "0/0":
                indv_index = int(dfi[dfi["PARTID"] == header_list[c]].index[0])
                dfi.at[indv_index, CN] += 1

    # calculate CN sum
    CN_columns = ["CN0","CN1", "CN3", "CN4"]
    dfi["CN_sum"] = dfi[CN_columns].sum(axis = 1)

    dfi.to_csv(output_file, sep = "\t", index = False)
    
    
    print("Finished running CN_counter() method from CNV_Plots.py.")
    
def indv_counts_bargraph(data_file):
    
    # import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon
    from matplotlib import figure
    import pandas as pd
    
    # import data
    df = pd.read_csv(data_file, sep = "\t")
    CN0 = sum(df["CN0"].tolist())
    CN1 = sum(df["CN1"].tolist())
    CN3 = sum(df["CN3"].tolist())
    CN4 = sum(df["CN4"].tolist())
    
    # set plot parameters    
    CN_categories = ['CN0', "CN1", "CN3", "CN4"]
    indv_counts = [CN0, CN1, CN3, CN4]
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.bar(CN_categories, indv_counts)
    plt.title('Individual counts by CNV category')
    fig.set_size_inches(12, 8)
    fig.savefig("../Results/Plots/indv_counts.png")
    plt.show()
    
# CALL PLOT COMMANDS

# CN counts boxplot
CN_counter("../Data/CNV3.vcf", "../Results/Plots/CN_counts_postfilter.tsv")
CN_counts_boxplot("../../CNV_Pipeline/Data/Starting_Manifest/CNV_counts_prefilter.tsv", "prefilter", "linear")
CN_counts_boxplot("../../CNV_Pipeline/Data/Starting_Manifest/CNV_counts_prefilter.tsv", "prefilter", "logscale")
CN_counts_boxplot("../Results/Plots/CN_counts_postfilter.tsv", "postfilter", "linear")
CN_counts_boxplot("../Results/Plots/CN_counts_postfilter.tsv", "postfilter", "logscale")

#indv_counts_bargraph("../Data/Starting_Manifest/CNV_counts_prefilter.tsv")