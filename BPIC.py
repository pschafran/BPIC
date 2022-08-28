#! /usr/bin/env python

import sys
import os
import glob
import subprocess
import shutil
from BPIC.assets.setup import *
from multiprocessing import Pool
from Bio import SeqIO
from collections import Counter


cite = '''
By: Peter W. Schafran
Last Updated: 2022 August 28
https://github.com/pschafran/BPIC
'''
version = "0.1"

help = '''
Required parameters:
-i, --input	Directory containing FASTA files
-f, --format	Input file format (either locus or taxon)

Optional parameters:
-a, --aligner	Alignment software (either mafft or clustal; default: mafft)
--force Force overwrite existing results
-l, --log	Log file name. File includes more details than screen output (default: printed to screen)
-o, --output	Output directory name (default: output)
-t, --threads	Maximum number of threads to use (default: 1)

Help:
-c, --cite	Show citation information
-h, --help	Show this help menu
-v, --version	Show version number

Advanced:
--mrbayes-path	Path to Mr Bayes executable
--mafft-path	Path to MAFFT executable
--clustal-path	Path to Clustal executable
--galax-path	Path to Galax executable
'''


### MAIN ###

if __name__ == "__main__":
	parameterDict = BPIC.assets.setup.checkArgs(sys.argv)
	globals().update(parameterDict)
	logOutput(parameterDict)
	logOutput("Parameter checking passed.\n------------------------")
	checkDependencies(aligner, mrBayesPath, mafftPath, clustalPath, galaxPath)
	checkFastas(inputDir,fileFormat)
