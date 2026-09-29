#!/bin/bash

#$ -q all.q
#$ -cwd 
#$ -V

module load sourmash

# K values

ks=(7 19 31)

# Species list
sp=(officinarum robustum barberi spontaneum)

# Saccharum officinarum
officinarum=(
SRR32744622
SRR32744623
SRR32744568
SRR32744569
SRR32744502
SRR32744503
SRR32744402
SRR32744403
SRR32743801
SRR32743801
SRR32740712
SRR32740713
SRR32740695
SRR32740696
SRR32673054
SRR32636969
SRR32634002
SRR32633743
SRR32632985
SRR32632781
SRR32632778
SRR32632685
SRR13165725
SRR13165728
)

# Saccharum spontaneum
spontaneum=(
SRR32744254
SRR32744255
SRR32744256
SRR32744241
SRR32744242
SRR32743804
SRR32743805
SRR12935194
)

# Saccharum barberi
barberi=(
SRR12929232
SRR32637094
SRR32745130
SRR32745131
SRR32745132
SRR32745133
SRR32745134
SRR32745135
SRR32744643
SRR32744644
SRR32744645
SRR32744570
SRR32744571
SRR32744572
SRR32744526
SRR32744527
SRR32744528
SRR32744504
SRR32744505
SRR32744506
SRR32744468
SRR32744469
SRR32744470
SRR32673125
SRR32637548
SRR32637126
)

# Saccharum robustum
robustum=(
SRR12934169
SRR32763979
SRR32763980
SRR32744547
SRR32744548
SRR32744463
SRR32744464
SRR32744347
SRR32744348
SRR32744339
SRR32744340
SRR32744250
SRR32744251
SRR32744013
SRR32744014
SRR32744016
SRR32744015
SRR32743592
SRR32743593
SRR32743551
SRR32743552
)


cd /Storage/data1/yuri.dantas/data/samples
mkdir -p catalogs

set -euo pipefail
shopt -s nullglob 

for species in "${sp[@]}" 
do
	tmpdir="${species}_tmp" 
	mkdir -p "$tmpdir"

	declare -n samples=$species

	for data in "${samples[@]}" 
	do
		if ! wget -q "https://labbces.cena.usp.br/sra_datasets/${data}_1.fastq.gz" --user "" --password "" -O "$tmpdir/${data}_1.fastq.gz"; then
			echo "Erro ao baixar ${data}_1"
			exit 1
		fi 
		if ! wget -q "https://labbces.cena.usp.br/sra_datasets/${data}_2.fastq.gz" --user "" --password "" -O "$tmpdir/${data}_2.fastq.gz"; then
			echo "Erro ao baixar ${data}_2"
			exit 1
		fi
	done
	echo "name,read1,read2" > "$tmpdir/${species}.csv"
	cd "$tmpdir" || exit 1

	for i in *_1.fastq.gz
	do
		base="$(basename  $i _1.fastq.gz)"
    		r2="${base}_2.fastq.gz"

		echo "${base},${i},${r2}" >> "${species}.csv"
	done 

	for k in "${ks[@]}"
	do
		sourmash scripts manysketch "${species}.csv" -p k="$k",scaled=1000,abund,dna -o "/Storage/data1/yuri.dantas/data/samples/catalogs/${species}_${k}.zip"
	done
	
	unset -n samples
cd /Storage/data1/yuri.dantas/data/samples || exit 1
rm -rf "$tmpdir"
done


