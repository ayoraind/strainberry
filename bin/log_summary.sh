#!/bin/bash

FILE=${1}

echo "Filename"$'\t'"contig_number_input"$'\t'"contig_length_input"$'\t'"N50_input"$'\t'"contig_number_output"$'\t'"contig_length_output"$'\t'"N50_output" > ${FILE}.logsummary

contig_number_input=$(cat ${FILE}.log | awk 'NR == 18 {print $2}')
contig_length_input=$(cat ${FILE}.log | awk 'NR == 19 {print $2}')
N50_input=$(cat ${FILE}.log | awk 'NR == 20 {print $2}')
contig_number_output=$(cat ${FILE}.log | awk 'NR == 23 {print $2}')
contig_length_output=$(cat ${FILE}.log | awk 'NR == 24 {print $2}')
N50_output=$(cat ${FILE}.log | awk 'NR == 25 {print $2}')
printf "%s\t%s\t%s\t%s\t%s\t%s\n" $FILE $contig_number_input $contig_length_input $N50_input $contig_number_output $contig_length_output $N50_output >> ${FILE}.logsummary
# remove bam file
#rm -rf /MIGE/01_DATA/02_MAPPING/${FILE}*.{bam,bai} 
