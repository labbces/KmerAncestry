#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V

mkdir -p kmc_temp


module load kmc

for sample in *.fastq.gz
do
	prefix="${sample%.fastq.gz}"
	/usr/bin/time -v kmc -k17 "$sample"  "kmers_$prefix" kmc_tmp
	kmc_tools transform "kmers_$prefix" dump "kmers_$prefix.txt"
done

