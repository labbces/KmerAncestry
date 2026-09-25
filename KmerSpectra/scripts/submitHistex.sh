#!/usr/bin/env bash
set -euo pipefail

PATH=/home/diriano/FASTK/:$PATH
export PATH

for ktab in FastK_*/*.ktab; do
	base=${ktab/.ktab}
	HIST=${ktab/ktab/hist.txt}
	GSFK=${ktab/ktab/hist.gsfk}
	if [[ ! -s $HIST ]]; then
	  echo "Generating Histogram file $base"
 	  Histex -A -h1:1000 $base > $HIST
	  Histex -G $base > $GSFK
	fi
	if [[ -s "${base}.hist.gsfk" ]]; then
		echo "Generating Kmer Spectra $base"
		eval "$(mamba shell hook --shell bash)"
		mamba activate genomescope2
		Rscript plotKmerSpectra.R $GSFK "${base}_kmerSpectra.png" "${base}_kmerSpectra" > "${base}_kmerSpectra.log" 2> "${base}_kmerSpectra.log"
		mamba deactivate
	fi
done
