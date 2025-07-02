set -ex

hyperfine --warmup 3 --runs 5 --export-json ../results/single_thread_performance.json \
    'freyja covariants ../bygul/jn.1.1.trimmed.bam --annot ../NC_045512_Hu-1.gff --output ../results/freyja_covariants.tsv --min_quality 20 --min_count 0 --threads 1' \
    'SAM_refiner -S ../bygul/jn.1.1.trimmed.sam -r ../NC_045512_Hu-1.fasta --mp 1 --min_count 1 --min_samp_abund 0.0 --min_col_abund 0.0 --ntabund 0.0 --ntcover 0 --max_dist 300 --max_covar 300 --read 1' \
    'covar -i ../bygul/jn.1.1.trimmed.bam -r ../NC_045512_Hu-1.fasta -a ../NC_045512_Hu-1.gff -o ../results/covar.tsv  --min_count 1 --min_quality 20 --threads 1'