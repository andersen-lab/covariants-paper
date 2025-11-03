
for file in bam/*.bam; do
    base=$(basename "$file")
    out="bam_fixed/$base"

    echo "Processing $file -> $out"

    samtools view -h "$file" | \
    awk 'BEGIN {OFS="\t"} 
         # Update @SQ header lines
         /^@SQ/ {sub("SN:2019-nCoV", "SN:NC_045512.2", $0); print; next} 
         # Keep other header lines unchanged
         /^@/ {print; next} 
         # Update RNAME column in alignments
         {if ($3=="2019-nCoV") $3="NC_045512.2"; print}' | \
    samtools view -b -o "$out" -
done
