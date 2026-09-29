# FastK Kmer Spectra

We are using [FastK](https://github.com/thegenemyers/FASTK) tro compute Kmer spectra, adn to find species-specific sets of kmers, which could be called anecestry informative kmers.

We will compute kmer spectra over the following data (Illumina paired end data, except SP80-3280 whihc is PacBio HiFi data):

- [*Saccharum barberi*](data/barberi.srrs.txt)
- [*Saccharum spontaneum*](data/spontaneum.srrs.txt)
- [*Saccharum robustum*](data/robustum.srrs.txt)
- [*Saccharum officinarum*](data/officinarum.srrs.txt)
- [*Saccharum hubrid cultivar SP80-3280*](data/SP803280.srrs.txt)

For each species I used FastK v1.2 to compute kmer catalogs. This was done for each accession and then all accession-specific catalogues were merged in to a single species catalog ([script](scripts/submitFastK.sh)). Then histograms were computed with Histex ([script](scripts/submitHistex.sh)) and plotted in R ([script](scripts/plotKmerSpectra.R)).  The plot script tries to identify the first valley (or sustained shoulder), to define the limit between error and genomic kmers. In the case of S. robustum, it did not work, and that limit was defined by visual inspection.

I will use the kmer spectra plot to define the minum depth to consider kmers as genomic, this is in order to exclude error kmers upto that depth.

| Species | Kmer spectra (k=19) | min_Depth (Error kmers) | First Peak |
| --- | --- | --- | --- |
| *S. barberi* | <img src="Figs/barberi_kmerSpectra_weighted.png" alt="Weighted k-mer spectrum for S. barberi" width="250"> | 9 | 34 |
| *S. spontaneum* | <img src="Figs/spontaneum_kmerSpectra_weighted.png" alt="Weighted k-mer spectrum for S. spontaneum" width="250"> | 10 | 35 |
| *S. robustum* | <img src="Figs/robustum_kmerSpectra_weighted.png" alt="Weighted k-mer spectrum for S. robustum" width="250"> | 6 | 15 |
| *S. officinarum* | <img src="Figs/officinarum_kmerSpectra_weighted.png" alt="Weighted k-mer spectrum for S. officinarum" width="250"> | 12 | 18 |
| *SP80-3280* | <img src="Figs/SP803280.19_kmerSpectra_weighted.png" alt="Weighted k-mer spectrum for SP80-3280" width="250"> | 13 | 30 |

# Pairwise comparisons

I will use KatComp from [Merqury.FK](https://github.com/thegenemyers/MERQURY.FK), to make pairwise comparisons between the kmer catalogs ([script](scripts/pairwise_comp.sh)). Using k=19. The goal here is to check whether there are species-specific kmers that could perhaps bring ancestry information.


|  | *S. barberi* | *S. spontaneum* | *S. robustum* | *S. officinarum* | *SP80-3280* |
| --- | --- | --- | --- | --- | --- |
| *S. barberi* | | <img src="Figs/KatComp__barberi__vs__spontaneum.fi.png" width="250"> | <img src="Figs/KatComp__barberi__vs__robustum.fi.png" width="250"> | <img src="Figs/KatComp__barberi__vs__officinarum.fi.png" width="250"> | <img src="Figs/KatComp__barberi__vs__SP803280.fi.png" width="250">  |
| *S. spontaneum* | <img src="Figs/KatComp__spontaneum__vs__barberi.fi.png" width="250"> | | <img src="Figs/KatComp__spontaneum__vs__robustum.fi.png" width="250"> | <img src="Figs/KatComp__spontaneum__vs__officinarum.fi.png" width="250"> | <img src="Figs/KatComp__spontaneum__vs__SP803280.fi.png" width="250">  |
| *S. robustum* | <img src="Figs/KatComp__robustum__vs__barberi.fi.png" width="250"> | <img src="Figs/KatComp__robustum__vs__spontaneum.fi.png" width="250"> |  | <img src="Figs/KatComp__robustum__vs__officinarum.fi.png" width="250"> | <img src="Figs/KatComp__robustum__vs__SP803280.fi.png" width="250">  |
| *S. officinarum* | <img src="Figs/KatComp__officinarum__vs__barberi.fi.png" width="250"> | <img src="Figs/KatComp__officinarum__vs__spontaneum.fi.png" width="250"> | <img src="Figs/KatComp__officinarum__vs__robustum.fi.png" width="250"> |  | <img src="Figs/KatComp__officinarum__vs__SP803280.fi.png" width="250">  | 
| *SP80-3280* | <img src="Figs/KatComp__SP803280__vs__barberi.fi.png" width="250"> | <img src="Figs/KatComp__SP803280__vs__spontaneum.fi.png" width="250"> | <img src="Figs/KatComp__SP803280__vs__robustum.fi.png" width="250"> | <img src="Figs/KatComp__SP803280__vs__officinarum.fi.png" width="250"> |  | 


# Species-specific kmer catalogs

I will use Logex (part of FASTK) to generate species-specific catalogs, for the *Saccharum* ancestral species. In Logex I can perform set operations. So, for the set of species, A, B, C and D, in order to compute specific kmer for species A, I will computhe the union of kmers of species B, C and D, and then 'substract' A from that union set, and so on for each species.
For that I need to define a lower and max for kmer depth/multiplicity. As lower (**L**) bound I will use 30% of the first peak or the limit for error kmers from the table above. For upper (**U**) bound I will use 3 times the first peak.

| Species | First Peak | L | U | Observations |
| --- | --- | --- | --- | --- | 
| *S. barberi* |  34| 9 | 102 |  |
| *S. spontanuem* | 35 | 10 | 105 |  |  |
| *S. robustum* | 15 | 6 | 45 |  |
| *S. officinarum* | 18 | 6 | 54 |  |  |
| SP80-3280 | 30 | 10 | 90 |  |  |


```bash
Logex -T10 \
  'barberi_specific=A[9-102]-(B[10-]|.C[6-]|.D[6-])' \
  'spontaneum_specific=B[10-105]-(A[9-]|.C[6-]|.D[6-])' \
  'robustum_specific=C[6-45]-(A[9-]|.B[10-]|.D[6-])' \
  'officinarum_specific=D[6-54]-(A[9-]|.B[10-]|.C[6-])' \
  barberi spontaneum robustum officinarum

```