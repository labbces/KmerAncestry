# KmerAncestry
Sliding window-based ancestral origin detection using k-mer frequency signatures.

KmerAncestry identifies the species-level ancestral origin of genomic windows through k-mer 
spectra comparison. The tool requires a reference genome sequence and computes k-mer profiles 
across sliding windows, comparing each against species-specific k-mer signatures derived from 
Illumina short-read data. Distance-based assignment determines the most probable ancestral 
contributor for each window.
