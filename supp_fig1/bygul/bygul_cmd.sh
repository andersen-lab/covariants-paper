set -ex

bygul simulate-proportions EPI_ISL_18967918.fasta ../ARTIC_V4-1.bed ../NC_045512_Hu-1.fasta  --outdir output --redo --readcnt 5000 --read_length 100 --maxmismatch 0

minimap2 -ax sr ../NC_045512_Hu-1.fasta  output/reads_1.fastq output/reads_2.fastq | samtools view -bS | samtools sort -o jn.1.1_test.bam
ivar trim -i jn.1.1_test.bam -b ARTIC_V4-1.bed | samtools sort -o jn.1.1.trimmed.bam
samtools index jn.1.1.trimmed.bam

samtools view -h jn.1.1.trimmed.bam > jn.1.1.trimmed.sam