# Tree building functions
import subprocess

def mrBayes(cmd):
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	(out, err) = process.communicate() #the stdout and stderr
