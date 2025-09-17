set -ex

for file in $(ls ../RSA/fastq); do
    basename=$(basename $file .R1_001.fastq.gz)
    minimap2 -ax sr ../sars2_metadata/NC_045512_Hu-1.fasta ${basename}.R1_001.fastq.gz ${basename}.R2_001.fastq.gz | samtools view -bS | samtools sort -o ../RSA/bam/${basename}.bam
    ivar trim -i ${basename}.bam -b ../sars2_metadata/ARTICv5.3.2.bed | samtools sort -o ../RSA/bam/${basename}.bam
    samtools index ${basename}.bam
done
