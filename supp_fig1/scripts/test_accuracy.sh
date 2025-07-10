set -ex


SAM_Refiner -S jn.1.1_10000_merged.trimmed.sam  -r ../NC_045512_Hu-1.fasta --min_count 1 --min_samp_abund 0.01 --ntcover 1 --wgs 1 --AAreport 0 --mp 4
covar -i ../bygul/jn.1.1_10000_merged.trimmed.bam -o ../results/covar/covar_10000.tsv -r ../NC_045512_Hu-1.fasta -a ../NC_045512_Hu-1.gff --min_depth 1 --min_quality 0 --min_frequency 0.01 --threads 4