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
from assets.tree import extractMargLik
import multiprocessing
from Bio import SeqIO
from collections import Counter

cite = '''
By: Peter W. Schafran
Last Updated: 2022 August 28
https://github.com/pschafran/BPIC
'''
version = "0.1"

help = '''
Required parameters
-i, --input	Directory containing FASTA files
-f, --format	Input file format (either locus or taxon)

Optional parameters (require a value after the flag)
-a, --aligner	Alignment software (either mafft or clustal; default: mafft)
-l, --log	Log file name. File includes more details than screen output (default: printed to screen)
-o, --output	Output directory name (default: output)
-t, --threads	Maximum number of threads to use (default: 1)

Optional flags (do not require a value after the flag)
--CDS	Partition MrBayes analysis by coding site
--force Overwrite existing results

MrBayes Parameters -v alues must be recognized by MrBayes (see https://nbisweden.github.io/MrBayes/manual.html)
--mrbayes-nst	Substitution model
--mrbayes-rates	Model for among-site rate variation
--mrbayes-ngen	Number of cycles for MCMC
--mrbayes-burninfrac	Proportion of samples to be discarded for convergence calculation (burn-in)
--mrbayes-samplefreq	How often to sample the Markov chain
--mrbayes-nsteps	Number of steps in the stepping-stone analysis
--mrbayes-CDS	Partitions analysis based on codon site

Help
-c, --cite	Show citation information
-h, --help	Show this help menu
-v, --version	Show version number

Advanced
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
	#logOutput(log, logFile, parameterDict)
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
	os.mkdir("%s/alignments/nexus" % outputDir)
	os.mkdir("%s/alignments/nexus/information_content" % outputDir)
	fileList = [ file for file in glob.glob("%s/alignments/*" % outputDir) if os.path.isfile(file) ]
	logOutput(log, logFile, "Converting alignments...")
	for file in fileList:
		convertFastaToNexus(file, outputDir, CDS, mrBayesNST, mrBayesRates, mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesNsteps, threads)
	logOutput(log, logFile, "Conversion finished.\n---------------------------")

	# Concatenate alignments
	logOutput(log, logFile, "Concatenating alignments...")
	concatLocusList = []
	# Iterate over all possible combinations of loci
	for locus1 in locusList:
		for locus2 in locusList:
			if "%s-%s" %(locus2,locus1) not in concatLocusList and locus1 != locus2: # check if the opposite combination of loci has already been done
				concatenateAlignments("%s/alignments/nexus/%s.nex" %(outputDir,locus1), "%s/alignments/nexus/%s.nex" %(outputDir,locus2), outputDir, CDS, mrBayesNST, mrBayesRates, mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesNsteps, threads)
				concatLocusList.append("%s-%s" %(locus1,locus2))
	logOutput(log, logFile, "Concatenation finished.\n---------------------------")

	# Generate Trees
	logOutput(log, logFile, "Generating marginal likelihood trees...")
	fileList = [ file for file in glob.glob("%s/alignments/nexus/*" % outputDir) if os.path.isfile(file) ]
	totalFiles = len(fileList)
	fileCounter = 1
	mbCmdList = []
	for file in sorted(fileList):
		mbCmdList.append("%s %s" %(mrBayesPath, file))
	with multiprocessing.Pool(threads) as pool:
		for count,value in enumerate(pool.imap(mrBayes,mbCmdList)):
			#Progress bar
			completionPerc = int(float(100*count)/float(len(mbCmdList)))
			sys.stdout.write('\r')
			sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
			sys.stdout.flush()
		sys.stdout.write('\n')
	logOutput(log, logFile, "Marginal likelihood trees finished.")
	logOutput(log, logFile, "Generating information content trees...")
	fileList = [ file for file in glob.glob("%s/alignments/nexus/information_content/*" % outputDir) if os.path.isfile(file) ]
	totalFiles = len(fileList)
	fileCounter = 1
	mbCmdList = []
	for file in sorted(fileList):
		mbCmdList.append("%s %s" %(mrBayesPath, file))
	with multiprocessing.Pool(threads) as pool:
		for count,value in enumerate(pool.imap(mrBayes,mbCmdList)):
			#Progress bar
			completionPerc = int(float(100*count)/float(len(mbCmdList)))
			sys.stdout.write('\r')
			sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
			sys.stdout.flush()
		sys.stdout.write('\n')
	logOutput(log, logFile, "Information content trees finished.\n---------------------------")

	# Extract marginal likelihoods from MrBayes logs
	logOutput(log, logFile, "Extracting marginal likelihoods...")
	os.mkdir("%s/tree_info" %(outputDir))
	margLikelihoodDict = {}
	with open("%s/tree_info/marginal_likelihoods.txt" % outputDir, "w") as outfile:
		for locus1 in locusList:
			for locus2 in locusList:
				if "%s-%s" %(locus1,locus2) in concatLocusList:
					 concat_brlenUnlinked_ln = extractMargLik("%s/alignments/nexus/%s-%s_brlen-unlinked.nex.log" %(outputDir,locus1,locus2))
					 concat_brlenLinked_ln = extractMargLik("%s/alignments/nexus/%s-%s_brlen-linked.nex.log" %(outputDir,locus1,locus2))
					 gene1_ln = extractMargLik("%s/alignments/nexus/%s.nex.log" %(outputDir,locus1))
					 gene2_ln = extractMargLik("%s/alignments/nexus/%s.nex.log" %(outputDir,locus2))
					 margLikelihoodDict.update({"%s-%s" : {"gene1": locus1, "gene2": locus2, "concat_brlenUnlinked_ln": concat_brlenUnlinked_ln, "concat_brlenLinked_ln": concat_brlenLinked_ln, "gene1_ln": gene1_ln, "gene2_ln": gene2_ln}})
					 outfile.write("{gene1: %s, gene2: %s, concat_brlenUnlinked_ln: %s, concat_brlenLinked_ln: %s, gene1_ln: %s, gene2_ln: %s},\n" %(locus1,locus2,concat_brlenUnlinked_ln, concat_brlenLinked_ln, gene1_ln, gene2_ln))
	logOutput(log, logFile, "Marginal likelihoods written to %s/tree_info/marginal_likelihoods.txt." % outputDir)

	# Run Galax
	logOutput(log, logFile, "Analyzing trees with Galax...")
	treefileList = [ file for file in glob.glob("%s/alignments/nexus/information_content/*.t" % outputDir) if os.path.isfile(file) ]
	with open("%s/tree_info/galax_treelist.txt" % outputDir, "w") as galaxInputFile:
		for file in treefileList:
			galaxInputFile.write("%s\n" % file)
	galaxCmd = "%s --listfile %s/tree_info/galax_treelist.txt --skip 1000 --outfile %s/tree_info/galax_output" %(galaxPath, outputDir, outputDir)
	process = subprocess.Popen(galaxCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	(out, err) = process.communicate()

	logOutput(log, logFile, "Galax finished.\n---------------------------")





#			#
