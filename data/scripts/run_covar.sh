set -ex

#cargo install --path covar

for file in $(ls ../bam/*.bam); do

    # Check if output file already exists
    if [ -f "../covar_out/$(basename $file).covariants.tsv" ]; then
        echo "Output file ../covar_out/$(basename $file).covariants.tsv already exists. Skipping $file."
        continue
    fi
    #covar -i $file -r ../sars2_metadata/2019-nCoV.fasta -a ../sars2_metadata/NC_045512_Hu-1.gff --output ../covar_out/$(basename $file).covariants.tsv -c 10 -t 16 -s 21563 -e 25384 || continue
    covar -i $file -r ../sars2_metadata/NC_045512_Hu-1.fasta -a ../sars2_metadata/NC_045512_Hu-1.gff --output ../covar_out/$(basename $file).covariants.tsv -c 10 -t 16 -s 21563 -e 25384 || continue
done