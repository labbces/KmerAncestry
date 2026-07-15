#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

module load singularity-ce 

cd /Storage/data1/yuri.dantas/data/

samples=("SRR32673054" "SRR12935194" "SRR12934169" "SRR12929232")
ks=(4 7 11 13 15 17 19 21 23 25 27 29 31)

outdir="frac_test_results"
logdir="frac_test_logs"
csv=frac_test.csv

mkdir -p $outdir $logdir

echo "sample,k,scaled,elapsed_sec,user_sec,sys_sec,max_rss_kb,output_size_bytes" > $csv

for sample in "${samples[@]}"
do 
	for k in "${ks[@]}"
	do 
		outfile="${outdir}/${sample}_k${k}"
		logfile="${logdir}/${sample}_k${k}.log"

		/usr/bin/time -v bash -c "singularity run frackmc.sif ./samples/${sample}_inter.fastq.gz ${outfile} --ksize ${k} --fq --a" 2> ${logfile}
		
		elapsed=$(grep "Elapsed (wall clock) time" $logfile | awk '{print $8}')
		user=$(grep "User time (seconds)" $logfile | awk '{print $4}')
		sys=$(grep "System time (seconds)" $logfile | awk '{print $4}')
		mem=$(grep "Maximum resident set size" $logfile | awk '{print $6}')
		elapsed_sec=$(echo $elapsed | awk -F: '{if(NF==3) print $1*3600+$2-60+$3; else print $1*60+$2}')
		size=$(stat -c%s "$outfile")
		
		echo "$sample,$k,1000,$elapsed_sec,$user,$sys,$mem,$size">>$csv
	done
done
