#! /usr/bin/env python

import sys
import os
import glob
import subprocess
import shutil
import BPIC.assets.setup
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
	# Check inputs
	parameterDict = BPIC.assets.setup.checkArgs(sys.argv)
	globals().update(parameterDict)
	logOutput(parameterDict)
	logOutput("Parameter checking passed.\n------------------------")
	checkDependencies(aligner, mrBayesPath, mafftPath, clustalPath, galaxPath)
	fileList = [ file for file in glob.glob("%s/*" % inputDir) if os.path.isfile(file) ]
	locusList = checkFastas(inputDir,fileFormat)

	# Set up files for alignment
	os.mkdir("%s/sequence_files" %(outputDir))
	if fileFormat == "locus":
		for file in fileList:
			shutil.copy(file,"%s/sequence_files/" %(outputDir),follow_symlinks=True)
	elif fileFormat == "taxon":
		for file in filelist:
			taxon = file.split(".")[-1:].join(".")
			record_dict = SeqIO.to_dict(SeqIO.parse("%s/%s" %(inputDir,file), "fasta"))
			for record in record_dict:
				with open("%s/%s.fasta" % (outputDir, record), "a+") as outSeqFile:
					sequence = record_dict[record]
					outSeqFile.write(">%s\n%s\n" % (taxon,sequence))

	# Align files
	
