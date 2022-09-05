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
		print(splitline)
