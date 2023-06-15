#!/bin/sh

echo "Running StrVCTVRE.py"
cd /lab01/Tools/StrVCTVRE
source /home/xcao/p/anaconda3/etc/profile.d/conda.sh
conda activate StrVCTVRE_py_2.7
python StrVCTVRE.py -i /lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Data/ASD_only/GMN.bed -o /lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Data/ASD_only/GMN_StrVCTVRE.tsv -f bed -a GRCh37 -l /lab01/Tools/liftOver/liftOver
echo "StrVCTVRE.py has completed running"

