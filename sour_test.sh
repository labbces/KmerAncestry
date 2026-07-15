#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

module load sourmash/4.9.4

cd /Storage/data1/yuri.dantas/data/

samples=("SRR32673054" "SRR12935194" "SRR12934169" "SRR12929232")
ks=(4 7 11 13 15 17 19 21 23 25 27 29 31)

outdir="sour_test_results"
logdir="sour_test_logs"
csv=sour_test.csv

mkdir -p $outdir $logdir

echo "sample,k,scaled,elapsed_sec,user_sec,sys_sec,max_rss_kb,output_size_bytes" > $csv

for sample in "${samples[@]}"
do 
	for k in "${ks[@]}"
	do 
		outfile="${outdir}/${sample}_k${k}.sig"
		logfile="${logdir}/${sample}_k${k}.log"

		/usr/bin/time -v bash -c "sourmash sketch dna -p k=${k},scaled=1000,abund --name ${sample}_k${k} -o ${outfile} ${sample}_?.fastq.gz" 2> ${logfile}
		
		elapsed=$(grep "Elapsed (wall clock) time" $logfile | awk '{print $8}')
		user=$(grep "User time (seconds)" $logfile | awk '{print $4}')
		sys=$(grep "System time (seconds)" $logfile | awk '{print $4}')
		mem=$(grep "Maximum resident set size" $logfile | awk '{print $6}')
		elapsed_sec=$(echo $elapsed | awk -F: '{if(NF==3) print $1*3600+$2-60+$3; else print $1*60+$2}')
		size=$(stat -c%s "$outfile")
		
		echo "$sample,$k,1000,$elapsed_sec,$user,$sys,$mem,$size">>$csv
	done
done
