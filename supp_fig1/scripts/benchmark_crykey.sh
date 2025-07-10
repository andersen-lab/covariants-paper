set -ex


declare -a read_counts=("1000" "10000" "100000")

for read_count in "${read_counts[@]}"; do
    hyperfine --warmup 3 --runs 5 --export-json ../results/crykey_performance_${read_count}.json \
        'python crykey_wastewater.py -i ../bygul/jn.1.1_'${read_count}'_merged.trimmed.metadata.tsv -r ../NC_045512_Hu-1.fasta -d ../crykey/crykey_dbs -o ../results/crykey' 
done