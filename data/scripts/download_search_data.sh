set -ex

while read -r line; do
  aws s3 cp s3://ucsd-all/${line} ../raw_bam --recursive --exclude "*" --include "*.trimmed.sorted*.bam*"
done < ../search_bucket.txt

cp -r ../raw_bam/*/*.bam* ../bam/
cp ../raw_bam/*/*/*/*/*.bam* ../bam/ 
