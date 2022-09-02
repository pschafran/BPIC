#! /usr/bin/env python

import sys
import time
import os
import glob
import subprocess
import shutil
from assets.setup import checkArgs
from assets.setup import logOutput
from assets.setup import checkDependencies
from assets.setup import checkFastas
from assets.alignment import mafftAlign
from assets.alignment import clustalAlign
from assets.alignment import convertFastaToNexus
from assets.alignment import concatenateAlignments
from assets.tree import mrBayes
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
	bpicDir = "/".join(os.path.abspath(__file__).split("/")[0:-1])
	dependencyDir = "%s/dependencies" % bpicDir
	# Check inputs
	#print("Parameter checking...")
	parameterDict = checkArgs(sys.argv)
	globals().update(parameterDict)
	logOutput(log, logFile, parameterDict)
	logOutput(log, logFile, "Parameter checking passed.\n---------------------------")
	mrBayesPath,galaxPath,alignerPath = checkDependencies(aligner, mrBayesPath, mafftPath, clustalPath, galaxPath, dependencyDir, log, logFile)
	logOutput(log, logFile, "Dependency checking passed.\n---------------------------")
	fileList = [ file for file in glob.glob("%s/*" % inputDir) if os.path.isfile(file) ]
	locusList = checkFastas(inputDir,fileFormat, log, logFile)
	logOutput(log, logFile, "File checking passed.\n---------------------------")

	# Set up files for alignment
	os.mkdir("%s" %(outputDir))
	os.mkdir("%s/sequence_files" %(outputDir))
	locusList = []
	if fileFormat == "locus":
		for file in fileList:
			filename = file.split("/")[-1]
			locus = ".".join(filename.split(".")[0:-1])
			if locus not in locusList:
				locusList.append(locus)
			shutil.copy(file,"%s/sequence_files/" %(outputDir),follow_symlinks=True)
	elif fileFormat == "taxon":
		for file in fileList:
			filename = file.split("/")[-1]
			taxon = ".".join(filename.split(".")[0:-1])
			record_dict = SeqIO.to_dict(SeqIO.parse("%s/%s" %(inputDir,file), "fasta"))
			for record in record_dict:
				if record not in locusList:
					locusList.append(record)
				with open("%s/sequence_files/%s.fasta" % (outputDir, record), "a+") as outSeqFile:
					sequence = record_dict[record].seq
					outSeqFile.write(">%s\n%s\n" % (taxon,sequence))

	# Align files
	os.mkdir("%s/alignments" %(outputDir))
	fileList = [ file for file in glob.glob("%s/sequence_files/*" % outputDir) ]
	totalFiles = len(fileList)
	fileCounter = 1
	logOutput(log, logFile, "Aligning files...")
	if aligner == "mafft":
		for file in fileList:
			mafftAlign(alignerPath, file, outputDir, threads)
			#Progress bar
			completionPerc = int(float(100*fileCounter)/float(totalFiles))
			sys.stdout.write('\r')
			sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
			sys.stdout.flush()
			fileCounter+=1
		sys.stdout.write('\n')
	elif aligner == "clustal":
		for file in fileList:
			mafftAlign(alignerPath, file, outputDir, threads)
			#Progress bar
			completionPerc = int(float(100*fileCounter)/float(totalFiles))
			sys.stdout.write('\r')
			sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
			sys.stdout.flush()
			fileCounter+=1
		sys.stdout.write('\n')
	logOutput(log, logFile, "Alignments finished.\n---------------------------")

	#Convert alignments to Nexus
	os.mkdir("%s/alignments/nexus" %(outputDir))
	fileList = [ file for file in glob.glob("%s/alignments/*" % outputDir) if os.path.isfile(file) ]
	logOutput(log, logFile, "Converting alignments...")
	for file in fileList:
		convertFastaToNexus(file, outputDir)
	logOutput(log, logFile, "Conversion finished.\n---------------------------")

	# Concatenate alignments
	logOutput(log, logFile, "Concatenating alignments...")
	for locus1 in locusList:
		for locus2 in locusList:
			if locus1 != locus2:
				concatenateAlignments("%s/alignments/nexus/%s.nex" %(outputDir,locus1), "%s/alignments/nexus/%s.nex" %(outputDir,locus2), outputDir)
	logOutput(log, logFile, "Concatenation finished.\n---------------------------")

	# Generate Trees
	logOutput(log, logFile, "Running MrBayes...")
	fileList = [ file for file in glob.glob("%s/alignments/nexus/*" % outputDir) if os.path.isfile(file) ]
	totalFiles = len(fileList)
	fileCounter = 1
	for file in fileList:
		mrBayes(mrBayesPath, file, outputDir, log, logFile, threads)
		completionPerc = int(float(100*fileCounter)/float(totalFiles))
		sys.stdout.write('\r')
		sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
		sys.stdout.flush()
		fileCounter+=1
	sys.stdout.write('\n')
	logOutput(log, logFile, "MrBayes finished.")
