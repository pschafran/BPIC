# Tree building functions
import subprocess
from ete3 import Tree
import re
import time
import signal

def mrBayes(cmd, timeout=5):
	file = cmd.split(" ")[-1]
	logFile = open("%s.log" % file, "w")
	try:
		process = subprocess.run(cmd, stdout=logFile, stderr=subprocess.PIPE, timeout = timeout*60, shell=True)
		failed = False
	#(out, err) = process.communicate() #the stdout and stderr
	except subprocess.TimeoutExpired:
		failed = True
	logFile.close()
	return cmd, failed

def extractMargLik(mbLogFile):
	linecounter = 1
	captureline = 0
	ln = None
	for line in open(mbLogFile, "r"):
		splitline = line.strip("\n").split(" ")
		if "".join(splitline) == "RunMarginallikelihood(ln)":
			captureline = linecounter + 2
		if linecounter == captureline:
			splitline = line.strip("\n").split(" ")
			cleanline = list(filter(None,splitline))
			ln = cleanline[1]
		linecounter += 1
	return float(ln)

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

def extractNexusTree(treefile):
	with open(treefile, "r") as infile:
		linecounter = 1
		capture = False
		translateDict = {}
		for line in infile:
			splitline = line.strip("\n").split(" ")
			cleanline = list(filter(None,splitline))
			if len(cleanline) >= 1:
				if cleanline[0] == ";":
					capture = False
				if cleanline[0] == "translate":
					capture = True
				if capture == True:
					translateDict.update({cleanline[0] : cleanline[1].strip(",")})
				if cleanline[0] == "tree":
					tree = cleanline[-1]
					for key in translateDict:
						re.sub("%s:" % key, "%s:" % translateDict[key], tree)
	return tree


def rfDistance(treefile1, treefile2):
	t1 = Tree(extractNexusTree(treefile1))
	t2 = Tree(extractNexusTree(treefile2))
	out = t1.robinson_foulds(t2, unrooted_trees=True)
	rf = out[0]
	rfMax = out[1]
	return rf
