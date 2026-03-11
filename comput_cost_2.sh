#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

k=(11 13 15 17 19 21 23 25 27 29 31)

mkdir -p kmc_temp
mkdir -p logs

module load kmc

for i in "${k[@]}"
do
    temp_dir="kmc_temp_${i}"
    mkdir -p "$temp_dir"

    for sample in *.fastq.gz
    do
        prefix="${sample%.fastq.gz}"

        rm -rf "${temp_dir}"/*

        log_file="logs/kmc_k${i}_${prefix}.log"
        time_file="logs/time_k${i}_${prefix}.log"
        size_file="logs/dbsize_k${i}_${prefix}.log"

        echo "k=${i}" > "$log_file"
        echo "sample=${prefix}" >> "$log_file"

        /usr/bin/time -f "elapsed_sec=%e user_sec=%U sys_sec=%S max_rss_kb=%M" \
            -o "$time_file" \
            kmc -k${i} "$sample" "kmers_k${i}_${prefix}" "${temp_dir}" \
            >> "$log_file" 2>&1

        du -sm kmers_k${i}_${prefix}* > "$size_file"

        kmc_tools transform "kmers_k${i}_${prefix}" dump "kmers_k${i}_${prefix}.txt"
    done
done
