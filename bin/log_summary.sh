#!/bin/bash

FILE=${1}

echo "Filename"$'\t'"contig_number_input"$'\t'"contig_length_input"$'\t'"N50_input"$'\t'"contig_number_iteration1"$'\t'"contig_length_iteration1"$'\t'"N50_iteration1"$'\t'"contig_number_iteration2"$'\t'"contig_length_iteration2"$'\t'"N50iteration2"$'\t'"contig_number_FinalIteration"$'\t'"contig_length_FinalIteration"$'\t'"N50_FinalIteration" > ${FILE}.logsummary

log_file=${FILE}.log

# Check if "Iteration #1" and "Iteration 2" are both present

if grep -q "Iteration #1" ${log_file} && grep -q "Iteration #2" ${log_file}
then

contig_number_input=$(cat ${log_file} | awk 'NR == 32 {print $2}')
contig_length_input=$(cat ${log_file} | awk 'NR == 33 {print $2}')
N50_input=$(cat ${log_file} | awk 'NR == 34 {print $2}')
contig_number_iteration1=$(cat ${log_file} | awk 'NR == 37 {print $2}')
contig_length_iteration1=$(cat ${log_file} | awk 'NR == 38 {print $2}')
N50_iteration1=$(cat ${log_file} | awk 'NR == 39 {print $2}')
contig_number_iteration2=$(cat ${log_file} | awk 'NR == 42 {print $2}')
contig_length_iteration2=$(cat ${log_file} | awk 'NR == 43 {print $2}')
N50_iteration2=$(cat ${log_file} | awk 'NR == 44 {print $2}')
contig_number_FinalIteration=$(cat ${log_file} | awk 'NR == 47 {print $2}')
contig_length_FinalIteration=$(cat ${log_file} | awk 'NR == 48 {print $2}')
N50_FinalIteration=$(cat ${log_file} | awk 'NR == 49 {print $2}')

# Check if only "Iteration #1" is present
elif grep -q "Iteration #1" ${log_file}
then

contig_number_input=$(cat ${log_file} | awk 'NR == 32 {print $2}')
contig_length_input=$(cat ${log_file} | awk 'NR == 33 {print $2}')
N50_input=$(cat ${log_file} | awk 'NR == 34 {print $2}')
contig_number_iteration1=$(cat ${log_file} | awk 'NR == 37 {print $2}')
contig_length_iteration1=$(cat ${log_file} | awk 'NR == 38 {print $2}')
N50_iteration1=$(cat ${log_file} | awk 'NR == 39 {print $2}')
contig_number_iteration2="NA"
contig_length_iteration2="NA"
N50_iteration2="NA"
contig_number_output=$(cat ${log_file} | awk 'NR == 42 {print $2}')
contig_length_output=$(cat ${log_file} | awk 'NR == 43 {print $2}')
N50_output=$(cat ${log_file} | awk 'NR == 44 {print $2}')

# If none of the patterns are present
else

contig_number_input=$(cat ${log_file} | awk 'NR == 18 {print $2}')
contig_length_input=$(cat ${log_file} | awk 'NR == 19 {print $2}')
N50_input=$(cat ${log_file} | awk 'NR == 20 {print $2}')
contig_number_iteration1="NA"
contig_length_iteration1="NA"
N50_iteration1="NA"
contig_number_iteration2="NA"
contig_length_iteration2="NA"
N50_iteration2="NA"
contig_number_FinalIteration=$(cat ${log_file} | awk 'NR == 23 {print $2}')
contig_length_FinalIteration=$(cat ${log_file} | awk 'NR == 24 {print $2}')
N50_FinalIteration=$(cat ${log_file} | awk 'NR == 25 {print $2}')

fi

printf "%s\t%s\t%s\t%s\t%s\t%s\n" $FILE $contig_number_input $contig_length_input $N50_input $contig_number_iteration1 $contig_length_iteration1 $N50_iteration1 $contig_number_iteration2 $contig_length_iteration2 $N50_iteration2 $contig_number_FinalIteration $contig_length_FinalIteration $N50_FinalIteration >> ${FILE}.logsummary
