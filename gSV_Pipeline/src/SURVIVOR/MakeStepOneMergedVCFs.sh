cd /lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/merge_step_one/pre-merge_sample_vcfs
mkdir -p ../merged_sample_vcfs

for file in `ls ../sample_sets`; do
    /lab01/Tools/SURVIVOR/Debug/SURVIVOR merge ../sample_sets/$file 100 0 1 1 0 0 ../merged_sample_vcfs/${file}.vcf
done