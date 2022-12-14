# Bayesian Phylogenetic Information Content

#### Dev Notes:
* Test files run start to finish. Produces results.html that can be viewed stand-alone
* Locus-based file arrangement untested

#### To Do
* Add CDS option for coding site partitioning
* Ability to continue partially run analyses
* Add options for processing AA seqs
* Fix overlapping genes in results.html

## Overview
Script for batch processing of homologous loci for up to ~12 taxa to calculate Bayesian phylogenetic information content in the data. Steps:

1. Input data can be a FASTA file for each taxon containing one sequence for each locus, or a FASTA a file for each locus containing one sequence per taxon.
2. Alignment (via MAFFT or Clustal Omega).
3. Tree building (via Mr. Bayes)
4. Analyze information content (via Galax)

## Installation
Clone github repo
```
git clone https://github.com/pschafran/BPIC.git
```

#### Dependencies
* Python 3.6
* BioPython [https://biopython.org/](https://biopython.org/)
* ETE Toolkit [http://etetoolkit.org/download/](http://etetoolkit.org/download/)
* MAFFT [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/) and/or Clustal Omega [http://www.clustal.org/omega/](http://www.clustal.org/omega/)
* MrBayes [https://nbisweden.github.io/MrBayes/download.html](https://nbisweden.github.io/MrBayes/download.html)
* Galax [https://github.com/plewis/galax](https://github.com/plewis/galax)

In my experience, using Conda to install MAFFT and Clustal Omega works, MrBayes works best when installed through Homebrew. Installing in a new environment is recommended to avoid conflicts if other conda packages are already installed:
```
conda create -n BPIC -c bioconda -c conda-forge -c etetoolkit mafft clustalo ete3 ete_toolchain python=3.6
pip install biopython
```

The base `BPIC` directory contains an empty directory where you can link dependencies if they are not read directly from the `PATH` variable. In particular, software installed with Conda may not be autodetected by `BPIC`. If you are able to call the program directly in a Terminal, you can do:
```  
cd BPIC/dependencies
ln -s $(which mafft)
ln -s $(which clustalo)
ln -s $(which galax)
ln -s $(which mb)
```
If executables are not in the Terminal's `PATH`, you can link to the absolute path of the file:
```
cd BPIC/dependencies
ln -s /Users/peter/bin/galax-1.1.0-mac/galax1 galax
etc.
```
Alternatively, paths can be provided on the command line when the program is run:
```
BPIC.py -i input_files -o output -f taxon --mafft-path /Users/peter/opt/miniconda3/bin/mafft --mrbayes-path /opt/homebrew/bin/mb ...

```
## Test Run
Once dependencies are installed and/or linked, you can test everything it working with the files included in the `test` directory. It takes about 5 minutes to run with 12 threads on our server, or about 5 minutes to run with 4 threads on a MacBook Air.
```
cd BPIC/test
../BPIC.py -i . -o output -t 4 -f taxon
```

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

## Running
Call the main script `BPIC.py` from command line. At minimum, you must specify the path to the directory containing all the input FASTA files and the format of the files (taxon or locus based as described above). Other parameters are:
```
Required parameters
-i, --input	Directory containing FASTA files
-f, --format	Input file format (either locus or taxon)

Optional parameters (require a value after the flag)
-a, --aligner	Alignment software (either mafft or clustal; default: mafft)
-l, --log	Log file name. File includes more details than screen output (default: printed to screen/STDOUT) // TODO
-o, --output	Output directory name (default: output)
-t, --threads	Maximum number of threads to use (default: 1)

Optional flags (do not require a value after the flag)
--CDS	Partition MrBayes analysis by coding site
--PROT	Analyze sequences as amino acids // TODO
--force Overwrite existing results
--continue	Resume a previously interrupted run // TODO

MrBayes Parameters -v alues must be recognized by MrBayes (see https://nbisweden.github.io/MrBayes/manual.html)
--mrbayes-nst	Substitution model
--mrbayes-rates	Model for among-site rate variation
--mrbayes-ngen	Number of cycles for MCMC
--mrbayes-burninfrac	Proportion of samples to be discarded for convergence calculation (burn-in)
--mrbayes-samplefreq	How often to sample the Markov chain
--mrbayes-nsteps	Number of steps in the stepping-stone analysis

Help
-c, --cite	Show citation information
-h, --help	Show this help menu
-v, --version	Show version number

Advanced
--mrbayes-path	Path to Mr Bayes executable
--mafft-path	Path to MAFFT executable
--clustal-path	Path to Clustal executable
--galax-path	Path to Galax executable
```

## Runtime
Runtime is highly dependent on individual system resources (# CPUS, CPU speed) and the number of taxa and genes in the analysis. In my experience, datasets of ~70 genes with 5-10 taxa will take 24-36 hrs on a 24 CPU server.

## Output
In the output directory,

## MrBayes Bug Workaround
A bug in MrBayes may cause some analyses of *concatenated genes with unlinked branch lengths* to never finish. These processes can be identified in as "task manager" type process viewer such as `top` or `htop` (recommended) and killed when it is clear they're running much longer than other tree analyses. The program will continue until it ends prematurely with an error. Alternatively, you can kill the main BPIC.py process (Ctrl+C in the main terminal window). Rerun your initial `BPIC` command, adding the `--continue mrbayes` flag. The stuck trees will rerun and usually finish normally; rarely, repeating this process a second time is required. Once all trees are completed, `BPIC` will be able to finish and generate `results.html`.

## References

Paul O. Lewis, Ming-Hui Chen, Lynn Kuo, Louise A. Lewis, Karolina Fu????kov??, Suman Neupane, Yu-Bo Wang, Daoyuan Shi, Estimating Bayesian Phylogenetic Information Content, Systematic Biology, Volume 65, Issue 6, November 2016, Pages 1009???1023, https://doi.org/10.1093/sysbio/syw042

Suman Neupane, Karolina Fu????kov??, Louise A Lewis, Lynn Kuo, Ming-Hui Chen, Paul O Lewis, Assessing Combinability of Phylogenomic Data Using Bayes Factors, Systematic Biology, Volume 68, Issue 5, September 2019, Pages 744???754, https://doi.org/10.1093/sysbio/syz007

Galax. https://github.com/plewis/galax

ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Jaime Huerta-Cepas, Francois Serra and Peer Bork. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046
