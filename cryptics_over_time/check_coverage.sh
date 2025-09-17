set -ex

# Create or overwrite the output TSV file
output_file="../median_coverage_summary.tsv"
echo -e "sample\tMedianCoverage" > "$output_file"

for file in ../bam/*.bam; do
    # Extract coverage values and calculate the median
    median_coverage=$(samtools depth "$file" | awk '{print $3}' | sort -n | awk '{
        count[NR] = $1
    } END {
        if (NR % 2) {
            print count[(NR + 1) / 2]
        } else {
            print (count[NR / 2] + count[(NR / 2) + 1]) / 2
        }
    }')
    echo -e "${file##*/}\t$median_coverage" >> "$output_file"
done