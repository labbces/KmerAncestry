#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

module load sourmash/4.9.4

cd /Storage/data1/yuri.dantas/data/samples

samples=("SRR32673054" "SRR12935194" "SRR12934169" "SRR12929232")
ks=(4 7 11 13 15 17 19 21 23 25 27 29 31)

outdir="branchwater_results"
logdir="branchwater_logs"
csv=branchwater_test.csv

echo "sample,k,scaled,elapsed,user,sys,mem,size" > "$csv"

mkdir -p "$outdir" "$logdir"
for sample in "${samples[@]}"
do
	for k in "${ks[@]}"
	do
		infile="manysketch_${sample}.csv"
    		outfile="${outdir}/${sample}_k${k}.zip"
    		logfile="${logdir}/${sample}_k${k}.log"

    		/usr/bin/time -v sourmash scripts manysketch "$infile" \
        	-o "$outfile" \
        	-p k=${k},scaled=1000,abund,dna \
       		2> "$logfile"
    		elapsed=$(grep "Elapsed (wall clock) time" $logfile | awk '{print $8}')
    		user=$(grep "User time (seconds)" $logfile | awk '{print $4}')
    		sys=$(grep "System time (seconds)" $logfile | awk '{print $4}')
    		mem=$(grep "Maximum resident set size" $logfile | awk '{print $6}')

    		elapsed_sec=$(echo $elapsed | awk -F: '{if(NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}')

    		size=$(stat -c%s "$outfile")

    		echo "$sample,$k,1000,$elapsed_sec,$user,$sys,$mem,$size" >> $csv
	done
done
