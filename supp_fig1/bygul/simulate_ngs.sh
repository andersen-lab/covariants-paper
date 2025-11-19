set -ex

declare -a read_counts=("1000" "10000" "100000") # reads per amplicon

for read_count in "${read_counts[@]}"; do
    bygul simulate-proportions EPI_ISL_18967918.fasta ../ARTIC_V4-1.bed ../NC_045512_Hu-1.fasta  --outdir output --redo --readcnt $read_count --read_length 100 --maxmismatch 0

    vsearch --fastq_mergepairs output/reads_1.fastq --reverse output/reads_2.fastq --fastqout output/merged.fastq
    minimap2 -ax sr ../NC_045512_Hu-1.fasta output/merged.fastq | samtools view -bS | samtools sort -o jn.1.1_${read_count}_merged.bam
    ivar trim -i jn.1.1_${read_count}_merged.bam -b ../ARTIC_V4-1.bed | samtools sort -o jn.1.1_${read_count}_merged.trimmed.bam
    samtools index jn.1.1_${read_count}_merged.trimmed.bam
    samtools view -h jn.1.1_${read_count}_merged.trimmed.bam > jn.1.1_${read_count}_merged.trimmed.sam
    lofreq call --no-default-filter -A -B -a 1 -b 1 --min-bq 0 --min-alt-bq 0 -f ../NC_045512_Hu-1.fasta jn.1.1_${read_count}_merged.trimmed.bam  -o jn.1.1_${read_count}_merged.trimmed.vcf --force-overwrite # relaxed settings for testing

    rm jn.1.1_${read_count}_merged.bam
done
