# Tree building functions
import subprocess

def mrBayes(cmd):
	file = cmd.split(" ")[-1]
	logFile = open("%s.log" % file, "w")
	process = subprocess.Popen(cmd, stdout=logFile, stderr=subprocess.PIPE, shell=True, text=True)
	(out, err) = process.communicate() #the stdout and stderr
	logFile.close()

def extractMargLik(mbLogFile):
	linecounter = 1
	captureline = 0
	for line in open(mbLogFile, "r"):
		splitline = line.strip("\n").split(" ")
		if "".join(splitline) == "RunMarginallikelihood(ln)":
			captureline = linecounter + 2
		if linecounter == captureline:
			splitline = line.strip("\n").split(" ")
			cleanline = list(filter(None,splitline))
			ln = cleanline[1]
		linecounter += 1
	return ln

def extractIpct(galaxOutFile):
	linecounter = 1
	startCapture = False
	stopCapture = False
	captureDict = {}
	for line in open(galaxOutFile, "r"):
		splitline = line.strip("\n").split(" ")
		if len(list(filter(None,splitline))) >= 1:
			if list(filter(None,splitline))[0] == "average":
				stopCapture = True
			if startCapture == True and stopCapture == False:
				captureList = list(filter(None,splitline))
				filename = captureList[0].split("/")[-1]
				locus = ".".join(filename.split(".")[0:-2])
				ipct = captureList[6]
				captureDict.update({locus : ipct})
			if "".join(list(filter(None,splitline))) == "treefileuniquecoverageHH*IIpctDDpct":
				startCapture = True
		linecounter += 1
	return captureDict
