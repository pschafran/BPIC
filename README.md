# Bayesian Phylogenetic Information Content

## Overview
Script for batch processing of homologous loci for up to 12 taxa to calculate Bayesian phylogenetic information content in the data. Steps:
0. Input data can be a FASTA file for each taxon containing one sequence for each locus, or a FASTA a file for each locus containing one sequence per taxon.
1. Alignment (via MAFFT or Clustal Omega).
2. Tree building (via Mr. Bayes)
3. Analyze information content (via Galax)


## Input Data Formatting
Input data can be a FASTA file for each taxon containing one sequence for each locus, or a FASTA a file for each locus containing one sequence per taxon. In either case, sequence names must be identical across files for each taxon/locus. Example formatting below:
### File for each taxon (e.g. exported from Geneious annotations)
Taxon1.fasta
```
>Gene 1
ATGCAGCTGCTGATCGATGCTAATGCTGATCGTCAGTCGTAGTAGCTAGCTAGT
>Gene 2
GCGATGCGCGATCGATGCTGGCTAGCTGCTAGTCATGCTGTAGCTAGTCGTGCA
>Gene 3
CGATCGATCGTAGTGCTGATGCTGTAGTCGTGTAGCTAGCTGTAGTCAGTAGCT
```
Taxon2.fasta
```
>Gene 1
ATGCAGCTGCTCGTAGCTCGTACGTCAGATGCTAGTCTCGTAGTAGTAGCTAGT
>Gene 2
GCGATGCGCGATCGATCGATGCATGCTATCGATGCTGTGTAGCGCAGTCGTGCA
>Gene 3
CGATCGATCGTAGTACGTCGTAGCTAGCTAGTCTAGCTAGTGTAGTCAGTAGCT
```
### File for each locus (e.g. output from Orthofinder)
Locus1.fasta
```
>Taxon 1
ATGCAGCTGCTGATCGATGCTAATGCTGATCGTCAGTCGTAGTAGCTAGCTAGT
>Taxon 2
GCGATGCGCGATCGATGCTGGCTAGCTGCTAGTCATGCTGTAGCTAGTCGTGCA
>Taxon 3
CGATCGATCGTAGTGCTGATGCTGTAGTCGTGTAGCTAGCTGTAGTCAGTAGCT
```
Locus2.fasta
```
>Taxon 1
ATGCAGCTGCTGATCGATGCTAATGCTGATCGTCAGTCGTAGTAGCTAGCTAGT
>Taxon 2
GCGATGCGCGATCGATGCTGGCTAGCTGCTAGTCATGCTGTAGCTAGTCGTGCA
>Taxon 3
CGATCGATCGTAGTGCTGATGCTGTAGTCGTGTAGCTAGCTGTAGTCAGTAGCT
```

## References

Paul O. Lewis, Ming-Hui Chen, Lynn Kuo, Louise A. Lewis, Karolina Fučíková, Suman Neupane, Yu-Bo Wang, Daoyuan Shi, Estimating Bayesian Phylogenetic Information Content, Systematic Biology, Volume 65, Issue 6, November 2016, Pages 1009–1023, https://doi.org/10.1093/sysbio/syw042

Suman Neupane, Karolina Fučíková, Louise A Lewis, Lynn Kuo, Ming-Hui Chen, Paul O Lewis, Assessing Combinability of Phylogenomic Data Using Bayes Factors, Systematic Biology, Volume 68, Issue 5, September 2019, Pages 744–754, https://doi.org/10.1093/sysbio/syz007

Galax. https://github.com/plewis/galax
