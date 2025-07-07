set -ex

hyperfine --warmup 3 --runs 5 --export-json ../results/multi_thread_performance.json \
    'SAM_refiner -S ../bygul/jn.1.1_merged.trimmed.sam  -r ../NC_045512_Hu-1.fasta --min_count 1 --min_samp_abund 0.0 --ntcover 1 --wgs 1 --AAreport 0 --mp 1' \
    'covar -i ../bygul/jn.1.1_merged.trimmed.bam -r ../NC_045512_Hu-1.fasta -a ../NC_045512_Hu-1.gff -o ../results/jn.1.1_coVar.tsv --min_count 0 --min_quality 0 --threads 1' \
    'freyja covariants ../bygul/jn.1.1_merged.trimmed.bam --annot ../NC_045512_Hu-1.gff --output ../results/freyja_covariants.tsv --min_count 0 --min_quality 0 --threads 1'