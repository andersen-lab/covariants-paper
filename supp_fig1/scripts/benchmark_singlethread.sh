set -ex


declare -a read_counts=("1000" "10000" "100000")

for read_count in "${read_counts[@]}"; do
    hyperfine --warmup 3 --runs 5 --export-json ../results/single_thread_performance_${read_count}.json \
        'SAM_Refiner -S ../bygul/jn.1.1_'${read_count}'_merged.trimmed.sam  -r ../NC_045512_Hu-1.fasta --min_count 1 --min_samp_abund 0.0 --ntcover 1 --wgs 1 --AAreport 0 --mp 1' \
        'covar -i ../bygul/jn.1.1_'${read_count}'_merged.trimmed.bam -r ../NC_045512_Hu-1.fasta -a ../NC_045512_Hu-1.gff --min_depth 1 --min_quality 0 --threads 1'
done
