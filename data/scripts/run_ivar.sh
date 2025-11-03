#set -ex

for file in $(ls ../original_bam/*.bam);
    do
    base=$(basename $file .bam)
    # echo base: $base
    # samtools view $file | grep "NM_" || true
    tmp=$(mktemp)
    (samtools mpileup -aa -A -d 0 -B -Q 0 --reference ../sars2_metadata/2019-nCoV.fasta "${file}" | ivar variants -p ../ivar/"${base}" -r ../sars2_metadata/2019-nCoV.fasta) 2> "$tmp"
    if [ -s "$tmp" ]; then
        echo "${base}" >> ../ivar/ivar_stderr_nonempty.txt
    fi
    rm -f "$tmp"
