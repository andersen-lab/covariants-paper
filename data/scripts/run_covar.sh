set -ex

#cargo install --path covar

for file in $(ls ../bam/*.bam); do

    covar -i $file -r ../sars2_metadata/2019-nCoV.fasta -a ../sars2_metadata/NC_045512_Hu-1.gff --output ../covar_out/$(basename $file).covariants.tsv -f 0.0 -d 10 -t 16 -s 21563 -e 25384 || continue
    #covar -i $file -r ../sars2_metadata/NC_045512_Hu-1.fasta -a ../sars2_metadata/NC_045512_Hu-1.gff --output ../covar_out/$(basename $file).covariants.tsv -f 0.0 -d 10 -t 16 -s 21563 -e 25384 || continue
done