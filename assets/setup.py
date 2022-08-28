# Setup functions

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

def checkFastas(inputDir, fileFormat):
	fileList = [ file for file in glob.glob("%s/*" % inputDir) if os.path.isfile(file) ]
	sequenceDict = {}
	agreementList = []
	disagreementList = []
	if fileFormat == "taxon":
		for file in fileList:
			taxon = file.split(".")[-1:].join(".")
			iterator = SeqIO.parse("%s/%s" %(inputDir,file), "fasta")
			seqIDs = sorted([record.id for record in iterator])
			sequenceDict.update({taxon : seqIDs})
		for taxon in sequenceDict:
			for taxon2 in sequenceDict:
				if sequenceDict[taxon] != sequenceDict[taxon2]:
					disagreementList.append(taxon)
				else:
					agreementList.append(taxon)
		if len(disagreementList) > 0:
			logOutput("ERROR: At least one locus not present for every taxon. Offending names are:")
			for taxon in sequenceDict:
				for taxon2 in sequenceDict:
					difference = sequenceDict[taxon].difference(sequenceDict[taxon2])
					if len(difference) > 0:
						logOutput("%s sequences not in %s:\t%s" %(taxon, taxon2, difference))
			exit(1)
		else:
			logOutput("Input file checking passed.\n------------------------")
	elif fileFormat == "locus":
		for file in fileList:
			taxon = file.split(".")[-1:].join(".")
			iterator = SeqIO.parse("%s/%s" %(inputDir,file), "fasta")
			seqIDs = sorted([record.id for record in iterator])
			sequenceDict.update({taxon : seqIDs})
		for taxon in sequenceDict:
			for taxon2 in sequenceDict:
				if sequenceDict[taxon] != sequenceDict[taxon2]:
					disagreementList.append(taxon)
				else:
					agreementList.append(taxon)
		if len(disagreementList) > 0:
			logOutput("ERROR: At least one taxon not present for every locus. Offending names are:")
			for taxon in sequenceDict:
				for taxon2 in sequenceDict:
					difference = sequenceDict[taxon].difference(sequenceDict[taxon2])
					if len(difference) > 0:
						logOutput("%s taxa not in %s:\t%s" %(taxon, taxon2, difference))
			exit(1)
		else:
			logOutput("Input file checking passed.\n------------------------")

def logOutput(logOutput):
	if log == FALSE:
		print(logOutput)
	elif log == TRUE:
		with open(logFile, "a+") as openLogFile:
			openLogFile.write("%s\n" % logOutput)
