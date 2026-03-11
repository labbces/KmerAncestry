#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

k=(11 13 15 17 19 21 23 25 27 29 31)

final_table="comp_cost.csv"

mkdir -p kmc_temp

echo "k,sample,reads,unique_kmers,db_size_mb,elapsed_sec,user_sec,sys_sec,max_rss_kb" > "$final_table"

module load kmc 

for i in "${k[@]}"
do
	temp_dir="kmc_temp_${i}"
	mkdir -p "$temp_dir"
	
	for sample in *.fastq.gz
	do
		prefix="${sample%.fastq.gz}"
		rm -rf "${temp_dir}"/*
		log_file="kmc_run.log"
		
		time_line=$( ( /usr/bin/time -f "%e, %U, %S, %M" \
			kmc -k${i} "$sample" "kmers_k${i}_${prefix}" "${temp_dir}" \ 
			> /dev/null ) 2> "$log_file" | tail -n 1 )
		reads=$(grep "Total no. of reads" "$log_file" | awk '{print $NF}')
		unique_kmers=$(grep "No. of unique k-mers" "$log_file" | awk '{print $NF}')
		db_size=$(du -sm "kmers_k${i}_${prefix}"* 2>/dev/null | awk '{sum+=$1}' 'END{print sum}')
		echo "${i},${prefix},${reads},${unique_kmers},${db_size},${time_line}" >> "$final_table"
		
		kmc_tools transform "kmers_k${i}_${prefix}" dump "kmers_k${i}_${prefix}.txt"
	done
done


