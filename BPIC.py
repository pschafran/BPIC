#! /usr/bin/env python

import sys
import os
import glob
import subprocess
import shutil
from multiprocessing import Pool


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

def checkArgs(commandline):
	acceptedParameters = ["-i","--input","-f","--format","-a","--aligner","-t","--threads","-o","--output","--force","-h","--help","-v","--version","-c","--cite","--mrbayes-path","--mafft-path","--clustal-path","--galax-path"]

	# Set default parameters
	aligner = "mafft"
	outpurDir = "output"
	threads = 1
	forceOverwrite = FALSE
	logFile = "null"
	mrBayesPath = ""
	mafftPath = ""
	clustalPath = ""
	galaxPath = ""

	if "-h" or "--help" in commandline:
		print(help)
		exit(0)
	if "-v" or "--version" in commandline:
		print(version)
		exit(0)
	if "-c" or "--cite" in commandline:
		print(cite)
		exit(0)
	if "-i" or "--input" not in commandline:
		print("ERROR: No input directory specified.")
		exit(1)
	if "-f" or "--format" not in commandline:
		print("ERROR: Input file format not specified.")
		exit(1)
	if "--force" in commandline:
		forceOverwrite = TRUE
	i = 0
	for parameter in commandline[i:]:
		if parameter in ["-i","--input"]:
			inputDir = commandline[i+1]
			if not os.path.isdir(inputDir):
				print("ERROR: Input directory %s not found." % inputDir)
				exit(1)
		if parameter in ["-f","-format"]:
			fileFormat = commandline[i+1]
			if fileFormat not in ["locus","taxon"]:
				print("ERROR: File format must be 'locus' or 'taxon'.")
				exit(1)
		if parameter in ["-a","--aligner"]:
			aligner = commandline[i+1]
			if aligner not in ["mafft","clustal"]:
				print("ERROR: Aligner must be 'mafft' or 'clustal'.")
				exit(1)
		if parameter in ["-l","--log"]:
			log = TRUE
			logFile = commandline[i+1]
		if parameter in ["-o","--output"]:
			outputDir = commandline[i+1]
			if os.path.isdir(outputDir) and forceOverwrite == FALSE:
				print("ERROR: Output directory already exists. Delete directory or set --force to overwrite.")
				exit(1)
		if parameter in ["-t","--threads"]:
			threads = int(commandline[i+1])
		if parameter == "--mrbayes-path":
			mrBayesPath = commandline[i+1]
		if parameter == "--mafft-path":
			mafftPath = commandline[i+1]
		if parameter == "--clustal-path":
			clustalPath = commandline[i+1]
		if parameter == "--galax-path":
			galaxPath = commandline[i+1]
		i += 1
	parameterDict = {"inputDir" : inputDir, "fileFormat" : fileFormat, "aligner" : aligner, "forceOverwrite" : forceOverwrite, "log" : log, "logFile" : logFile, "outputDir" : outputDir, "threads" : threads, "mrBayesPath" : mrBayesPath, "mafftPath" : mafftPath, "clustalPath" : clustalPath, "galaxPath" : galaxPath}
	return parameterDict

def checkDependencies(aligner,mrBayes,mafft,clustal,galax):
	try:
		mrBayesPath = shutil.which(mrbayes)
	except:
		mrBayesPath = shutil.which("mrbayes")
	else:
		print("ERROR: MrBayes could not be executed.")
	logOutput(mrBayesPath)
	try:
		galaxPath = shutil.which(galax)
	except:
		galaxPath = shutil.which("galax")
	else:
		print("ERROR: Galax could not be executed.")
	logOutput(galaxPath)
	if aligner == "mafft":
		try:
			mafftPath = shutil.which(mafft)
		except:
			mafftPath = shutil.which("mafft")
		else:
			print("ERROR: MAFFT could not be executed.")
		logOutput(mafftPath)
	elif aligner == "clustal":
		try:
			clustalPath = shutil.which(clustal)
		except:
			clustalPath = shutil.which("clustalo")
		else:
			print("ERROR: Clustal Omega could not be executed.")
		logOutput(clustalPath)

def logOutput(logOutput):
	if log == FALSE:
		print(logOutput)
	elif log == TRUE:
		with open(logFile, "a+") as openLogFile:
			openLogFile.write("%s\n" % logOutput)





### MAIN ###

if __name__ == "__main__":
	parameterDict = checkArgs(sys.argv)
	globals().update(parameterDict)
	logOutput(parameterDict)
	logOutput("Parameter checking passed.")
	logOutput("")
	checkDependencies(aligner, mrBayesPath, mafftPath, clustalPath, galaxPath)
