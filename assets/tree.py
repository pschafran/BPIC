# Tree building functions
import subprocess

def mrBayes(mrBayesPath, infile, outputDir, log = False, logfile= None, threads = 1):
	mrBayes_cline = "%s %s" %(mrBayesPath, infile)
	process = subprocess.Popen(mrBayes_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	process.wait()
	(out, err) = process.communicate() #the stdout and stderr
