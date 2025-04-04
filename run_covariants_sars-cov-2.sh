
for file in $(ls sorted_bam/*.bam); do
    #check if output file exists
    if [ -f covariants/$(basename $file).covariants.tsv ]; then
        continue
    fi
    freyja covariants $file 21563 25384 --annot freyja_metadata/sars-cov-2/NC_045512_Hu-1.gff --output covariants/$(basename $file).covariants.tsv --threads 8 --ref-genome freyja_metadata/sars-cov-2/NC_045512_Hu-1.fasta
done