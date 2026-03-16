# Create header
echo -e "file\tavg_base_quality" > data/filtered_bam_avg_quality.tsv

# Loop through BAM files (excluding .bai index files)
for f in data/filtered_bam/*.bam; do
  [ -f "$f" ] && echo -e "$(basename "$f")\t$(samtools stats "$f" 2>/dev/null | grep '^SN' | grep 'average quality' | awk -F'\t' '{print $4}')"
done >> data/filtered_bam_avg_quality.tsv