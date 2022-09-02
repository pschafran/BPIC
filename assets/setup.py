# Setup functions
import os
import shutil
import glob
from Bio import SeqIO


def checkArgs(commandline):
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

	acceptedParameters = ["-i","--input","-f","--format","-a","--aligner","-t","--threads","-o","--output","--force","-h","--help","-v","--version","-c","--cite","--mrbayes-path","--mafft-path","--clustal-path","--galax-path", "--mrbayes-nst", "--mrbayes-rates", "--mrbayes-ngen", "--mrbayes-burninfrac", "--mrbayes-samplefreq", "--mrbayes-printfreq"]

	# Set default parameters
	aligner = "mafft"
	outputDir = "output"
	threads = 1
	forceOverwrite = False
	log = False
	logFile = "null"
	mrBayesPath = ""
	mafftPath = ""
	clustalPath = ""
	galaxPath = ""
	mrBayesNST = 6
	mrBayesRates = "invgamma"
	mrBayesNgen = 10000000
	mrBayesBurninFrac = 0.25
	mrBayesSampleFreq = 1000
	mrBayesPrintFreq = 10000

	if "-h" in commandline or "--help" in commandline:
		print(help)
		exit(0)
	if "-v" in commandline or "--version" in commandline:
		print(version)
		exit(0)
	if "-c" in commandline or "--cite" in commandline:
		print(cite)
		exit(0)
	if "-i" not in commandline and "--input" not in commandline:
		print("ERROR: No input directory specified.")
		exit(1)
	if "-f" not in commandline and "--format" not in commandline:
		print("ERROR: Input file format not specified.")
		exit(1)
	if "--force" in commandline:
		forceOverwrite = True
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
			log = True
			logFile = commandline[i+1]
		if parameter in ["-o","--output"]:
			outputDir = commandline[i+1]
			if os.path.isdir(outputDir) and forceOverwrite == False:
				print("ERROR: Output directory already exists. Delete directory or set --force to overwrite.")
				exit(1)
			elif os.path.isdir(outputDir) and forceOverwrite == True:
				shutil.rmtree(outputDir)
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
		if parameter == "--mrbayes-nst":
			mrBayesNST = commandline[i+1]
		if parameter == "--mrbayes-rates":
			mrBayesRates = commandline[i+1]
		if parameter == "--mrbayes-ngen":
			mrBayesNgen = commandline[i+1]
		if parameter == "--mrbayes-burninfrac":
			mrBayesBurninFrac = commandline[i+1]
		if parameter == "--mrbayes-samplefreq":
			mrBayesSampleFreq = commandline[i+1]
		if parameter == "--mrbayes-printfreq":
			mrBayesPrintFreq = commandline[i+1]

		# Replace base logFile name with outputDir + logFile after whole command line is read
		if log == True:
			logFile = "%s/%s" %(outputDir,logFile)
		i += 1
	parameterDict = {"inputDir" : inputDir, "fileFormat" : fileFormat, "aligner" : aligner, "forceOverwrite" : forceOverwrite, "log" : log, "logFile" : logFile, "outputDir" : outputDir, "threads" : threads, "mrBayesPath" : mrBayesPath, "mafftPath" : mafftPath, "clustalPath" : clustalPath, "galaxPath" : galaxPath, "mrBayesNST" : mrBayesNST, "mrBayesRates" : mrBayesRates, "mrBayesNgen" : mrBayesNgen, "mrBayesBurninFrac" : mrBayesBurninFrac, "mrBayesSampleFreq" : mrBayesSampleFreq, "mrBayesPrintFreq" : mrBayesPrintFreq }
	return parameterDict

def checkDependencies(aligner, mrBayes, mafft, clustal, galax, dependencyDir, log = False, logFile = "null"):
	mrBayesPath = shutil.which(mrBayes)
	if mrBayesPath == None:
		mrBayesPath = shutil.which("mb")
	if mrBayesPath == None:
		mrBayesPath = shutil.which("%s/mb" % dependencyDir)
	if mrBayesPath == None:
		logOutput(log, logFile, "ERROR: MrBayes could not be found.")
		exit(1)
	logOutput(log, logFile, mrBayesPath)

	galaxPath = shutil.which(galax)
	if galaxPath == None:
		galaxPath = shutil.which("galax")
	if galaxPath == None:
		galaxPath = shutil.which("%s/galax" % dependencyDir)
	if galaxPath == None:
		logOutput(log, logFile, "ERROR: Galax could not be executed.")
		exit(1)
	logOutput(log, logFile, galaxPath)
	if aligner == "mafft":
		mafftPath = shutil.which(mafft)
		if mafftPath == None:
			mafftPath = shutil.which("mafft")
		if mafftPath == None:
			mafftPath = shutil.which("%s/mafft" % dependencyDir)
		if mafftPath == None:
			logOutput(log, logFile, "ERROR: MAFFT could not be executed.")
			exit(1)
		logOutput(log, logFile, mafftPath)
		return mrBayesPath,galaxPath,mafftPath
	elif aligner == "clustal":
		clustalPath = shutil.which(clustal)
		if clustalPath == None:
			clustalPath = shutil.which("clustalo")
		if clustalPath == None:
			clustalPath = shutil.which("%s/clustalo" % dependencyDir)
		if clustalPath == None:
			logOutput(log, logFine, "ERROR: Clustal Omega could not be executed.")
			exit(1)
		logOutput(log, logFile, clustalPath)
		return mrBayesPath,galaxPath,clustalPath

def checkFastas(inputDir, fileFormat, log = False, logFile = "null"):
	fileList = [ file for file in glob.glob("%s/*" % inputDir) if os.path.isfile(file) ]
	sequenceDict = {}
	agreementList = []
	disagreementList = []
	duplicateNameDict = {}
	duplicateNames = False
	if fileFormat == "taxon":
		#Parse all files for sequence names, check for duplicates in each file
		for file in fileList:
			locusList = []
			taxon = ".".join(file.split(".")[0:-1])
			iterator = SeqIO.parse("%s/%s" %(inputDir,file), "fasta")
			seqIDs = sorted([record.id for record in iterator])
			for locus in seqIDs:
				if locus in locusList:
					duplicateNames = True
					try:
						duplicateNameDict[file].append(locus)
					except:
						duplicateNameDict.update({file : [locus]})
				else:
					locusList.append(locus)
			sequenceDict.update({taxon : seqIDs})
		# Report any duplicate sequence names found
		if duplicateNames == True:
			logOutput(log, logfile, "ERROR: Duplicate sequence names found.")
			for name in duplicateNameDict:
				logOutput(log, logFile, "%s:\t%s" %(name, duplicateNameDict[name]))
			exit(1)
		# Check that each file has the same sequence names
		for taxon in sequenceDict:
			for taxon2 in sequenceDict:
				if sequenceDict[taxon] != sequenceDict[taxon2]:
					disagreementList.append(taxon)
				else:
					agreementList.append(taxon)
		# Report any files and sequence names that don't match others
		if len(disagreementList) > 0:
			logOutput(log, logFile, "ERROR: At least one locus not present for every taxon. Offending names are:")
			for taxon in sequenceDict:
				for taxon2 in sequenceDict:
					difference = set(sequenceDict[taxon]) - set(sequenceDict[taxon2])
					if len(difference) > 0:
						logOutput(log, logFile, "%s sequences not in %s:\t%s" %(taxon, taxon2, difference))
			exit(1)
		else:
			locusList = sequenceDict[taxon]
			logOutput(log, logFile, "Input file checking passed.\n------------------------")
	elif fileFormat == "locus":
		#Parse all files for sequence names, check for duplicates in each file
		for file in fileList:
			taxonList = []
			locus = ".".join(file.split(".")[0:-1])
			iterator = SeqIO.parse("%s/%s" %(inputDir,file), "fasta")
			seqIDs = sorted([record.id for record in iterator])
			sequenceDict.update({locus : seqIDs})
			for taxon in seqIDs:
				if taxon in taxonList:
					duplicateNames = True
					try:
						duplicateNameDict[file].append(taxon)
					except:
						duplicateNameDict.update({file : [taxon]})
				else:
					taxonList.append(taxon)
			sequenceDict.update({locus : seqIDs})
		# Report any duplicate sequence names found
		if duplicateNames == True:
			logOutput(log, logfile, "ERROR: Duplicate sequence names found.")
			for name in duplicateNameDict:
				logOutput(log, logFile, "%s:\t%s" %(name, duplicateNameDict[name]))
			exit(1)
		#Check that each file has the same sequence names
		for locus in sequenceDict:
			for locus2 in sequenceDict:
				if sequenceDict[locus] != sequenceDict[locus2]:
					disagreementList.append(locus)
				else:
					agreementList.append(locus)
		if len(disagreementList) > 0:
			logOutput(log, logFile, "ERROR: At least one taxon not present for every locus. Offending names are:")
			for locus in sequenceDict:
				for locus2 in sequenceDict:
					difference = set(sequenceDict[locus]) - set(sequenceDict[locus2])
					if len(difference) > 0:
						logOutput(log, logFile, "%s sequences not in %s:\t%s" %(locus, locus2, difference))
			exit(1)
		else:
			locusList = sequenceDict.keys()
			logOutput(log, logFile, "Input file checking passed.\n------------------------")
	return locusList

def logOutput(log, logFile,logOutput):
	if log == False:
		print(logOutput)
	elif log == True:
		with (logFile, "a+") as openLogFile:
			openLogFile.write("%s\n" % logOutput)
