#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

set -euo pipefail
shopt -s nullglob

cd /Storage/data1/yuri.dantas/data || exit 1

ks=(4 7 11 13 15 17 19 21 23 25 27 29 31)
final_table="comp_cost.csv"

echo "k,sample,elapsed_sec,user_sec,sys_sec,max_rss_kb,output_size_bytes" > "$final_table"

module load kmc

for k in "${ks[@]}"; do
    temp_dir="kmc_temp_${k}"
    mkdir -p "$temp_dir"

    for sample in *.fastq.gz; do
        prefix="${sample%.fastq.gz}"
        out_prefix="kmers_k${k}_${prefix}"
        log_file="${out_prefix}.time.log"

        rm -rf -- "${temp_dir:?}/"*

        /usr/bin/time -f "%e,%U,%S,%M" -o "$log_file" \
            kmc -k"$k" "$sample" "$out_prefix" "$temp_dir" > /dev/null

        IFS=',' read -r elapsed_sec user sys mem < "$log_file"
        size_pre=$(stat -c%s "${out_prefix}.kmc_pre")
        size_suf=$(stat -c%s "${out_prefix}.kmc_suf")
        size=$((size_pre + size_suf))

        echo "${k},${prefix},${elapsed_sec},${user},${sys},${mem},${size}" >> "$final_table"
    done
done
