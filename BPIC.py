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
from assets.tree import extractIpct
from assets.tree import rfDistance
from assets.html import writeHTML
import multiprocessing
from Bio import SeqIO
from collections import Counter

cite = '''
By: Peter W. Schafran
Last Updated: 2022 September 12
https://github.com/pschafran/BPIC
'''
version = "0.1"

help = '''
Required parameters
-i, --input	Directory containing FASTA files
-f, --format	Input file format (either locus, taxon, or alignment)

Optional parameters (require a value after the flag)
-a, --aligner	Alignment software (Options: mafft or clustal ; Default: mafft)
-l, --log	Log file name. File includes more details than screen output (default: printed to screen/STDOUT)
-o, --output	Output directory name (default: output)
-t, --threads	Maximum number of threads to use (default: 1)
--continue	Resume a previous incomplete run from either alignment or MrBayes stage (Options: align or mrbayes ; Default: off)

Optional flags (do not require a value after the flag)
--CDS	Partition MrBayes analysis by coding site
--force Overwrite existing results

MrBayes Parameters -v alues must be recognized by MrBayes (see https://nbisweden.github.io/MrBayes/manual.html)
--mrbayes-nst	Substitution model (default: 6)
--mrbayes-rates	Model for among-site rate variation (default: gamma)
--mrbayes-ngen	Number of cycles for MCMC (default: 10000000)
--mrbayes-samplefreq	How often to sample the Markov chain (default: 1000)
--mrbayes-nsteps	Number of steps in the stepping-stone analysis (default: 30)

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
	if continueRun == False:
		os.mkdir("%s" %(outputDir))
		os.mkdir("%s/sequence_files" %(outputDir))
		os.mkdir("%s/alignments" %(outputDir))
	#locusList = []
	# If files already organized by locus, just copy to working directory
	if continueRun == False:
		if fileFormat == "locus" or fileFormat == "alignment":
			for file in fileList:
				#filename = file.split("/")[-1]
				#locus = ".".join(filename.split(".")[0:-1])
				#if locus not in locusList:
					#locusList.append(locus)
				shutil.copy(file,"%s/sequence_files/" %(outputDir),follow_symlinks=True)
				if fileFormat == "alignment":
					shutil.copy(file,"%s/alignments/" %(outputDir),follow_symlinks=True)

		# If files organized by taxon, reorganize them into locus format
		elif fileFormat == "taxon":
			for file in fileList:
				filename = file.split("/")[-1]
				taxon = ".".join(filename.split(".")[0:-1])
				record_dict = SeqIO.to_dict(SeqIO.parse("%s/%s" %(inputDir,file), "fasta"))
				for record in record_dict:
					#if record not in locusList:
						#locusList.append(record)
					with open("%s/sequence_files/%s.fasta" % (outputDir, record), "a+") as outSeqFile:
						sequence = record_dict[record].seq
						outSeqFile.write(">%s\n%s\n" % (taxon,sequence))

	# setup lists to track which loci have progressed through which step
	if continueRun == True:
		if continuePoint == "align":
			incompleteAlignments = locusList
		if continuePoint == "mrbayes":
			incompleteMargLikelihood = []
			for locus in locusList:
				incompleteMargLikelihood.append(locus)
			for locus1 in locusList:
				for locus2 in locusList:
					if "%s-%s" %(locus2,locus1) not in incompleteMargLikelihood and locus1 != locus2:
						incompleteMargLikelihood.append("%s-%s" %(locus1,locus2))
			incompleteInformationContent = []
			for locus in locusList:
				incompleteInformationContent.append(locus)

	# Align files
	if fileFormat == "locus" or fileFormat == "taxon":
		alignFileList = [ file for file in glob.glob("%s/sequence_files/*.fasta" % outputDir) ]
		totalFiles = len(alignFileList)
		fileCounter = 1
		if continueRun == False:
			#os.mkdir("%s/alignments" %(outputDir))
			logOutput(log, logFile, "Aligning files...")
			if aligner == "mafft":
				mafftAlign(alignFileList, alignerPath, outputDir, threads, totalFiles, fileCounter)
			elif aligner == "clustal":
				clustalAlign(alignFileList, alignerPath, outputDir, threads, totalFiles, fileCounter)
			logOutput(log, logFile, "Alignments finished.\n---------------------------")
		elif continueRun == True and continuePoint == "align": # If continuing, compare names of alignment files to sequence file, remove
			logOutput(log, logFile, "Resuming aligning files...")
			completeAlignmentFiles = [ file for file in glob.glob("%s/alignments/*.fasta" % outputDir) if os.path.getsize(file) > 0 ]
			for file in completeAlignmentFiles:
				if os.path.getsize(file) == 0: # Mafft creates the empty file before finishing the alignment, need to remove to redo
					os.remove(file)
					completeAlignmentFiles.remove(file)
				else:
					filename = filename = file.split("/")[-1]
					locus = ".".join(filename.split(".")[0:-1])
					if "%s/sequence_files/%s.fasta" %(outputDir, locus) in alignFileList:
						alignFileList.remove("%s/sequence_files/%s.fasta" %(outputDir, locus))
			fileCounter = len(completeAlignmentFiles)
			if aligner == "mafft":
				mafftAlign(alignFileList, alignerPath, outputDir, threads, totalFiles, fileCounter)
			elif aligner == "clustal":
				clustalAlign(alignFileList, alignerPath, outputDir, threads, totalFiles, fileCounter)
			logOutput(log, logFile, "Alignments finished.\n---------------------------")
		elif continueRun == True and continuePoint == "mrbayes":
			logOutput(log, logFile, "Reusing previous alignments.\n---------------------------")
	elif fileFormat == "alignment":
		logOutput(log, logFile, "Input files already aligned.\n---------------------------")

	#Convert alignments to Nexus
	if continueRun == False or (continueRun == True and continuePoint == "align"):
		os.mkdir("%s/alignments/nexus" % outputDir)
		os.mkdir("%s/alignments/nexus/information_content" % outputDir)
	#elif continueRun == True and continuePoint == "mrbayes":
	#	if not os.path.isdir("%s/alignments/nexus" % outputDir):
	#		os.mkdir("%s/alignments/nexus" % outputDir)
	#	if not os.path.isdir("%s/alignments/nexus/information_content" % outputDir):
	#		os.mkdir("%s/alignments/nexus/information_content" % outputDir)
	fileList = [ file for file in glob.glob("%s/alignments/*.fasta" % outputDir) if os.path.isfile(file) ]
	logOutput(log, logFile, "Converting alignments...")
	locusLengthDict = {}
	for file in fileList:
		locus,nchar = convertFastaToNexus(file, outputDir, CDS, mrBayesNST, mrBayesRates, mrBayesNgen, mrBayesSampleFreq, mrBayesNsteps, threads)
		locusLengthDict.update({locus : nchar})
	logOutput(log, logFile, "Conversion finished.\n---------------------------")

	# Concatenate alignments
	logOutput(log, logFile, "Concatenating alignments...")
	concatLocusList = []
	# Iterate over all possible combinations of loci
	for locus1 in locusList:
		for locus2 in locusList:
			if "%s-%s" %(locus2,locus1) not in concatLocusList and locus1 != locus2: # check if the opposite combination of loci has already been done
				concatenateAlignments("%s/alignments/nexus/%s.nex" %(outputDir,locus1), "%s/alignments/nexus/%s.nex" %(outputDir,locus2), outputDir, CDS, mrBayesNST, mrBayesRates, mrBayesNgen, mrBayesSampleFreq, mrBayesNsteps, threads)
				concatLocusList.append("%s-%s" %(locus1,locus2))
	logOutput(log, logFile, "Concatenation finished.\n---------------------------")

	# Generate Trees
	fileList = [ file for file in glob.glob("%s/alignments/nexus/*.nex" % outputDir) if os.path.isfile(file) ]
	totalFiles = len(fileList)
	#fileCounter = 1
	mbCmdList = []
	if continueRun == False or (continueRun == True and continuePoint == "align"):
		logOutput(log, logFile, "Generating marginal likelihood trees...")
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
	# If continuing, check if MrBayes run completed (ran sump), if not delete all partial files except the nexus alignment
	elif continueRun == True and continuePoint == "mrbayes":
		logOutput(log, logFile, "Resuming making marginal likelihood trees...")
		fileCounter = 1
		for file in sorted(fileList):
			if os.path.isfile("%s.pstat" % file):
				fileList.remove(file)
				fileCounter += 1
			else:
				for deleteFile in glob.glob("%s.*" %file):
					if deleteFile.split(".")[-1] != "nex":
						logOutput(log, logFile, "Deleting partial file: %s" % deleteFile)
						os.remove(deleteFile)
		for file in sorted(fileList):
			mbCmdList.append("%s %s" %(mrBayesPath, file))
		filecounter = totalFiles - len(mbCmdList)
		logOutput(log, logFile, "Resuming making marginal likelihood trees...")
		with multiprocessing.Pool(threads) as pool:
			for count,value in enumerate(pool.imap(mrBayes,mbCmdList)):
				#Progress bar
				completionPerc = int(float(100*(filecounter+count))/float(totalFiles))
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
	if continueRun == False or (continueRun == True and continuePoint == "align"):
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
	elif continueRun == True and continuePoint == "mrbayes":
		for file in sorted(fileList):
			if os.path.isfile("%s.pstat" % file):
				fileList.remove(file)
		for file in sorted(fileList):
			mbCmdList.append("%s %s" %(mrBayesPath, file))
		filecounter = totalFiles - len(mbCmdList)
		with multiprocessing.Pool(threads) as pool:
			for count,value in enumerate(pool.imap(mrBayes,mbCmdList)):
				#Progress bar
				completionPerc = int(float(100*(filecounter+count))/float(totalFiles))
				sys.stdout.write('\r')
				sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
				sys.stdout.flush()
			sys.stdout.write('\n')
	logOutput(log, logFile, "Information content trees finished.\n---------------------------")

	# Extract marginal likelihoods from MrBayes logs
	logOutput(log, logFile, "Extracting marginal likelihoods...")
	os.mkdir("%s/tree_info" %(outputDir))
	margLikelihoodDict = {}
	# { locus1-locus2 : {"gene1": locus1, "gene2": locus2, "concat_brlenUnlinked_ln": concat_brlenUnlinked_ln, "concat_brlenLinked_ln": concat_brlenLinked_ln, "gene1_ln": gene1_ln, "gene2_ln": gene2_ln},\
	#   locus1-locus3 : {...},\
	#   locus1-locus4 : {...} }
	with open("%s/tree_info/marginal_likelihoods.csv" % outputDir, "w") as outfile:
		outfile.write("Locus 1,Locus 2,Concatenated BrLen Unlinked (lnL),Concatenated BrLen Linked (lnL),Sum of Separate Trees (lnL),RF Distance,Winner\n")
		for locus1 in locusList:
			for locus2 in locusList:
				if "%s-%s" %(locus1,locus2) in concatLocusList:
					winner = None
					concat_brlenUnlinked_ln = extractMargLik("%s/alignments/nexus/%s-%s_brlen-unlinked.nex.log" %(outputDir,locus1,locus2))
					concat_brlenLinked_ln = extractMargLik("%s/alignments/nexus/%s-%s_brlen-linked.nex.log" %(outputDir,locus1,locus2))
					gene1_ln = extractMargLik("%s/alignments/nexus/%s.nex.log" %(outputDir,locus1))
					gene2_ln = extractMargLik("%s/alignments/nexus/%s.nex.log" %(outputDir,locus2))
					rf = rfDistance("%s/alignments/nexus/%s.nex.con.tre" %(outputDir,locus1), "%s/alignments/nexus/%s.nex.con.tre" %(outputDir,locus2))
					margLikelihoodDict.update({"%s-%s" %(locus1,locus2) : {"gene1": locus1, "gene2": locus2, "concat_brlenUnlinked_ln": concat_brlenUnlinked_ln, "concat_brlenLinked_ln": concat_brlenLinked_ln, "gene1_ln": gene1_ln, "gene2_ln": gene2_ln, "rf_distance" : rf}})
					combinedln = float(gene1_ln) + float(gene2_ln)
					if concat_brlenUnlinked_ln > combinedln or concat_brlenLinked_ln > combinedln:
						winner = "Concatenated"
					elif concat_brlenUnlinked_ln < combinedln and concat_brlenLinked_ln < combinedln:
						winner = "Separate"
					outfile.write("%s,%s,%s,%s,%s,%s,%s\n" %(locus1, locus2, concat_brlenUnlinked_ln, concat_brlenLinked_ln, combinedln, rf, winner))
					#outfile.write("gene1: %s, gene2: %s, concat_brlenUnlinked_ln: %s, concat_brlenLinked_ln: %s, gene1_ln: %s, gene2_ln: %s},\n" %(locus1,locus2,concat_brlenUnlinked_ln, concat_brlenLinked_ln, gene1_ln, gene2_ln))
	logOutput(log, logFile, "Marginal likelihoods written to %s/tree_info/marginal_likelihoods.csv" % outputDir)

	# Run Galax
	logOutput(log, logFile, "Analyzing trees with Galax...")
	treefileList = [ file for file in glob.glob("%s/alignments/nexus/information_content/*.t" % outputDir) if os.path.isfile(file) ]
	with open("%s/tree_info/galax_treelist.txt" % outputDir, "w") as galaxInputFile:
		for file in treefileList:
			galaxInputFile.write("%s\n" % file)
	galaxCmd = "%s --listfile %s/tree_info/galax_treelist.txt --skip 1000 --outfile %s/tree_info/galax_output" %(galaxPath, outputDir, outputDir)
	process = subprocess.run(galaxCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
	#(out, err) = process.communicate()
	# Get Ipct numbers from Galax output file
	ipctDict = extractIpct("%s/tree_info/galax_output.txt" % outputDir)
	logOutput(log, logFile, "Galax finished.\n---------------------------")

	# Format JS data file ### No longer needed
	#with open("%s/tree_info/genedata.js" % outputDir, "w") as outfile:
	#	outfile.write("var pairs = [\n")
	#	for locuspair in margLikelihoodDict:
	#		outfile.write("{gene1:%s, gene2:%s, symtopo:%s, symtree:%s, apotree1:%s, apotree2:%s},\n" % (margLikelihoodDict[locuspair]["gene1"], margLikelihoodDict[locuspair]["gene2"], margLikelihoodDict[locuspair]["concat_brlenUnlinked_ln"], margLikelihoodDict[locuspair]["concat_brlenLinked_ln"], margLikelihoodDict[locuspair]["gene1_ln"], margLikelihoodDict[locuspair]["gene2_ln"]))
	#	outfile.write("]\n\n")
	#	outfile.write("var genes = [\n")
	#	for gene in ipctDict:
	#		outfile.write("{name:%s, chunk:1, ipct:%s, seqlen:100},\n" %(gene, ipctDict[gene]))
	#	outfile.write("]\n\n")
	#	outfile.write("var num_genes = %d\n" %(len(ipctDict.keys())))

	logOutput(log, logFile, "All data written to %s/tree_info/genedata.js" % outputDir)

	# Write HTML
	writeHTML("%s/results.html" % outputDir, margLikelihoodDict, ipctDict, locusLengthDict)
	logOutput(log, logFile, "Final results written to %s/results.html\n---------------------------" % outputDir)
	logOutput(log, logFile, "BPIC finished successfully.\n---------------------------")
#			#
