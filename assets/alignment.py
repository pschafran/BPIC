# Setup functions
import subprocess
import sys
import os
import re
from Bio import SeqIO
from Bio import AlignIO

def clustalAlign(fileList, clustalPath, outputDir, threads, totalFiles, fileCounter):
	for file in fileList:
		filename = file.split("/")[-1]
		locus = ".".join(filename.split(".")[0:-1])
		clustal_cline = "%s --auto --threads %s --in %s --out %s/alignments/%s.fasta" %(clustalPath, threads, file, outputDir, locus)
		process = subprocess.run(clustal_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
		#Progress bar
		completionPerc = int(float(100*fileCounter)/float(totalFiles))
		sys.stdout.write('\r')
		sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
		sys.stdout.flush()
		fileCounter+=1
	sys.stdout.write('\n')
	#process.wait()
	#(out, err) = process.communicate() #the stdout and stderr

def mafftAlign(fileList, mafftPath, outputDir, threads, totalFiles, fileCounter):
	for file in fileList:
		filename = file.split("/")[-1]
		locus = ".".join(filename.split(".")[0:-1])
		mafft_cline = "%s --auto --thread %s %s > %s/alignments/%s.fasta" %(mafftPath, threads, file, outputDir, locus)
		process = subprocess.run(mafft_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
		#Progress bar
		completionPerc = int(float(100*fileCounter)/float(totalFiles))
		sys.stdout.write('\r')
		sys.stdout.write("[%-100s] %d%%" % ('='*completionPerc, completionPerc))
		sys.stdout.flush()
		fileCounter+=1
	sys.stdout.write('\n')
	#process.wait()
	#(out, err) = process.communicate() #the stdout and stderr

def convertFastaToNexus(inputFile, outputDir, CDS, mrBayesNST = 6, mrBayesRates = "gamma", mrBayesNgen = 10000000, mrBayesSampleFreq = 1000, mrBayesNsteps = 30, threads = 1 ):
	filename = inputFile.split("/")[-1]
	locus = ".".join(filename.split(".")[0:-1])
	AlignIO.convert(inputFile, "fasta", "%s/alignments/nexus/%s.nex" %(outputDir,locus), "nexus", molecule_type = "DNA")
	AlignIO.convert(inputFile, "fasta", "%s/alignments/nexus/information_content/%s.nex" %(outputDir,locus), "nexus", molecule_type = "DNA")
	# Some files being converted with apostrophes added around taxon names that break MrBayes. Need to remove them.
	with open("%s/alignments/nexus/%s.nex" %(outputDir,locus), "r") as infile:
		lines = infile.readlines()
	with open("%s/alignments/nexus/%s.nex" %(outputDir,locus), "w") as outfile:
		for line in lines:
			outfile.write(re.sub("\'","", line))
	# End apostrophe replacement block
	with open("%s/alignments/nexus/%s.nex" %(outputDir,locus), "r") as infile:
		linecounter = 1
		for line in infile:
			line = line.strip("\n")
			if linecounter==3:
				splitline = line.strip(';').split(' ')
				nchar = int(splitline[2].split("=")[1])
			linecounter += 1
	with open("%s/alignments/nexus/%s.nex" %(outputDir,locus), "a") as outfile:
		outfile.write("\n")
		outfile.write("begin mrbayes;\n")
		outfile.write("set autoclose=yes nowarn=yes;\n")
		outfile.write("log start filename=%s/alignments/nexus/%s.output.log;\n" %(outputDir, locus))
		outfile.write("set usebeagle = no;\n")
		#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
		#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
		#outfile.write("mcmcp ngen=%s relburnin=yes burninfrac=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesNgen))
		if CDS == True:
			outfile.write("charset first  = 1-%d\\3;\n" % nchar)
			outfile.write("charset second = 2-%d\\3;\n" % nchar)
			outfile.write("charset third  = 3-%d\\3;\n" % nchar)
			outfile.write("partition cds = 3: first, second, third;\n")
			outfile.write("set partition = cds;\n")
			outfile.write("unlink shape=(all) statefreq=(all) revmat=(all);\n")
			outfile.write("prset ratepr=dirichlet(10.0,10.0,10.0);\n")
		else:
			outfile.write("prset ratepr=dirichlet(10.0);\n")
		outfile.write("set seed=5207 swapseed=5207;\n")
		#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
		outfile.write("Lset nst=6 rates=gamma ngammacat=4;\n")
		outfile.write("Prset brlenspr=unconstrained:GammaDir(1.0, 0.1,1.0, 1.0) shapepr=exp(1.0) Revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) statefreqpr=dirichlet(10.0,10.0,10.0,10.0);\n")
		outfile.write("mcmc ngen=10000000 samplefreq=100 printfreq=10000 nchains=1 nruns=1 filename=%s/alignments/nexus/%s.mcmc;\n" %(outputDir,locus))
		#outfile.write("mcmcp ngen=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesSampleFreq, mrBayesNgen))
		outfile.write("sumt filename=%s/alignments/nexus/%s.mcmc;\n" %(outputDir, locus))
		outfile.write("mcmcp filename=%s/alignments/nexus/%s.ss;\n" %(outputDir, locus))
		#outfile.write("ss alpha=0.3 nsteps=%s burninss=-2;\n" % mrBayesNsteps)
		outfile.write("ss alpha=0.3 nsteps=30 burninss=-2;\n")
		outfile.write("sump;\n")
		outfile.write("sumt conformat=simple;\n")
		outfile.write("log stop;\n")
		outfile.write("end;\n")
		outfile.close()
	with open("%s/alignments/nexus/information_content/%s.nex" %(outputDir,locus), "a") as outfile:
		outfile.write("\n")
		outfile.write("begin mrbayes;\n")
		outfile.write("set autoclose=yes nowarn=yes;\n")
		outfile.write("set usebeagle = no;\n")
		outfile.write("log start filename=%s/alignments/nexus/information_content/%s.output.log;\n" %(outputDir, locus))
		#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
		#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
		if CDS == True:
			outfile.write("charset first  = 1-%d\\3;\n" % nchar)
			outfile.write("charset second = 2-%d\\3;\n" % nchar)
			outfile.write("charset third  = 3-%d\\3;\n" % nchar)
			outfile.write("partition cds = 3: first, second, third;\n")
			outfile.write("set partition = cds;\n")
			outfile.write("unlink shape=(all) statefreq=(all) revmat=(all);\n")
			ooutfile.write("prset ratepr=dirichlet(10.0,10.0,10.0);\n")
		else:
			outfile.write("prset ratepr=dirichlet(10.0);\n")
		outfile.write("set seed=5207 swapseed=5207;\n")
		#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
		#outfile.write("set seed=9223 swapseed=9223;")
		outfile.write("Lset nst=6 rates=gamma ngammacat=4;\n")
		outfile.write("Prset brlenspr=unconstrained:GammaDir(1.0, 0.1,1.0, 1.0) shapepr=exp(1.0) Revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) statefreqpr=dirichlet(10.0,10.0,10.0,10.0);\n")
		outfile.write("mcmc ngen=10000000 samplefreq=100 printfreq=10000 nchains=1 nruns=1 filename=%s/alignments/nexus/information_content/%s.mcmc;\n" %(outputDir, locus))
		outfile.write("sumt filename=%s/alignments/nexus/information_content/%s.mcmc;\n" %(outputDir, locus))
		outfile.write("mcmcp filename=%s/alignments/nexus/information_content/%s.ss;\n" %(outputDir, locus))
		#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
		#outfile.write("mcmc ngen=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesSampleFreq, mrBayesNgen))
		outfile.write("sumt conformat=simple;\n")
		outfile.write("sump;\n")
		outfile.write("log stop;\n")
		outfile.write("end;\n")
		outfile.close()
	return locus,nchar

def concatenateAlignments(aln1, aln2, outputDir, CDS, mrBayesNST = 6, mrBayesRates = "gamma", mrBayesNgen = 10000000, mrBayesSampleFreq = 1000, mrBayesNsteps = 30, threads = 1 ):
	filename1 = aln1.split("/")[-1]
	locus1 = ".".join(filename1.split(".")[0:-1])
	filename2 = aln2.split("/")[-1]
	locus2 = ".".join(filename2.split(".")[0:-1])
	outfile = open("%s/alignments/nexus/%s-%s_brlen-linked.nex" %(outputDir,locus1,locus2), 'w')
	masterDict = {}
	#filecounter = 1
	for file in [ aln1, aln2 ]:
		#print "Processing: %s" % file
		infile = open(file, 'r')
		linecounter=1
		taxacounter=1
		taxalist=[]
		seqlist={}
		ntax=[]
		nchar=[]
		masterDict.update({file : {}})
		masterDict[file].update({ "sequences" : {} })
		#infiledict[file] = {}
		for line in infile:
			line = line.strip("\n")
			if linecounter==3:
				splitline = line.strip(';').split(' ')
				splittax = splitline[1].split('=')
				ntax = int(splittax[1])
				masterDict[file]['ntax'] = ntax
				nchar = int(splitline[2].split("=")[1])
				masterDict[file]['nchar'] = nchar
				#print "Number of taxa: %s" % ntax
			if line == "begin mrbayes;":
				break
			if linecounter>5 and len(line.split(" ")) >= 2:
				taxon = line.strip("\n").split(" ")[0]
				sequence = line.strip("\n").split(" ")[-1]
				if taxon not in taxalist:
					taxalist.append(taxon)
				try:
					masterDict[file]['sequences'][taxon].append(sequence)
				except:
					masterDict[file]['sequences'].update({taxon : [sequence]})
			linecounter+=1
	longestTaxon = max(taxalist, key=len)
	totalChar = masterDict[aln1]['nchar'] + masterDict[aln2]['nchar']
	outfile.write("#NEXUS\n")
	outfile.write("begin data;\n")
	outfile.write("\tdimensions ntax=%s nchar=%s;\n" % (masterDict[file]['ntax'], totalChar))
	outfile.write("\tformat datatype=dna missing=? gap=-;\n")
	outfile.write("matrix\n")
	for taxon in taxalist:
		padding = (len(longestTaxon) - len(taxon))+1
		outfile.write("%s%s%s%s\n" %(taxon, " "*padding, "".join(masterDict[aln1]["sequences"][taxon]), "".join(masterDict[aln2]["sequences"][taxon])))
	outfile.write(";\n")
	outfile.write("end;\n")
	outfile.write("\n")
	outfile.write("begin mrbayes;\n")
	outfile.write("set autoclose=yes nowarn=yes;\n")
	outfile.write("set usebeagle = no;\n")
	outfile.write("log start filename=%s/alignments/nexus/%s-%s_brlen-linked.output.log;\n" %(outputDir, locus1, locus2))
	#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
	#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
	#outfile.write("mcmcp ngen=%s relburnin=yes burninfrac=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesNgen))
	if CDS == False:
		outfile.write("charset first = 1-%d;\n" %(masterDict[aln1]['nchar']))
		outfile.write("charset second = %d-%d;\n" %(masterDict[aln1]['nchar']+1, totalChar))
		outfile.write("partition combined = 2: first, second;\n")
		outfile.write("set partition = combined;\n")
		outfile.write("prset ratepr=dirichlet(10.0,10.0);\n")
	elif CDS == True:
		outfile.write("charset first  = 1-%d\\3;\n" % masterDict[aln1]['nchar'])
		outfile.write("charset second = 2-%d\\3;\n" % masterDict[aln1]['nchar'])
		outfile.write("charset third  = 3-%d\\3;\n" % masterDict[aln1]['nchar'])
		outfile.write("charset fourth  = %d-%d\\3;\n" %(masterDict[aln1]['nchar']+1, masterDict[aln2]['nchar']))
		outfile.write("charset fifth = %d-%d\\3;\n" %(masterDict[aln1]['nchar']+1, masterDict[aln2]['nchar']))
		outfile.write("charset sixth  = %d-%d\\3;\n" %(masterDict[aln1]['nchar']+1, masterDict[aln2]['nchar']))
		outfile.write("partition cds = 6: first, second, third, fouth, fifth, sixth;\n")
		outfile.write("set partition = cds;\n")
		outfile.write("prset ratepr=dirichlet(10.0,10.0,10.0,10.0,10.0,10.0);\n")
	#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
	outfile.write("Lset applyto=(all) nst=6 rates=gamma ngammacat=4;\n")
	outfile.write("set seed=5207 swapseed=5207;\n")
	#outfile.write("set seed=9223 swapseed=9223;\n")
	#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
	outfile.write("Prset applyto=(all) brlenspr=unconstrained:GammaDir(1.0, 0.1,1.0, 1.0) shapepr=exp(1.0) Revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) statefreqpr=dirichlet(10.0,10.0,10.0,10.0);\n")
#	outfile.write("mcmcp ngen=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesSampleFreq, mrBayesNgen))
	outfile.write("unlink shape=(all) statefreq=(all) revmat=(all);\n")
	outfile.write("mcmc ngen=10000000 samplefreq=100 printfreq=10000 nchains=1 nruns=1 filename=%s/alignments/nexus/%s-%s_brlen-linked.mcmc;\n" %(outputDir, locus1, locus2))
	outfile.write("sumt filename=%s/alignments/nexus/%s-%s_brlen-linked.mcmc;\n" %(outputDir, locus1, locus2))
	outfile.write("mcmcp filename=%s/alignments/nexus/%s-%s_brlen-linked.ss;\n" %(outputDir, locus1, locus2))
	#outfile.write("ss alpha=0.3 nsteps=%s burninss=-2;\n" % mrBayesNsteps)
	outfile.write("ss alpha=0.3 nsteps=30 burninss=-2;\n")
	outfile.write("sumt conformat=simple;\n")
	outfile.write("sump;\n")
	#outfile.write("sumt burnin=2500;\n")
	outfile.write("log stop;\n")
	outfile.write("end;\n")
	outfile.close()

	# Write second concatenated file with branch lengths unlinked b/w partitions
	outfile = open("%s/alignments/nexus/%s-%s_brlen-unlinked.nex" %(outputDir,locus1,locus2), 'w')
	outfile.write("#NEXUS\n")
	outfile.write("begin data;\n")
	outfile.write("\tdimensions ntax=%s nchar=%s;\n" % (masterDict[file]['ntax'], totalChar))
	outfile.write("\tformat datatype=dna missing=? gap=-;\n")
	outfile.write("matrix\n")
	for taxon in taxalist:
		padding = (len(longestTaxon) - len(taxon))+1
		outfile.write("%s%s%s%s\n" %(taxon, " "*padding, "".join(masterDict[aln1]["sequences"][taxon]), "".join(masterDict[aln2]["sequences"][taxon])))
	outfile.write(";\n")
	outfile.write("end;\n")
	outfile.write("\n")
	outfile.write("begin mrbayes;\n")
	outfile.write("set autoclose=yes nowarn=yes;\n")
	outfile.write("set usebeagle = no;\n")
	outfile.write("log start filename=%s/alignments/nexus/%s-%s_brlen-unlinked.log;\n" %(outputDir, locus1, locus2))
	#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
	#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
	#outfile.write("mcmcp ngen=%s relburnin=yes burninfrac=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesNgen))
	if CDS == False:
		outfile.write("charset first = 1-%d;\n" %(masterDict[aln1]['nchar']))
		outfile.write("charset second = %d-%d;\n" %(masterDict[aln1]['nchar']+1, totalChar))
		outfile.write("partition combined = 2: first, second;\n")
		outfile.write("set partition = combined;\n")
		outfile.write("unlink shape=(all) statefreq=(all) revmat=(all) brlen=(all);\n")
		outfile.write("prset ratepr=dirichlet(10.0,10.0);\n")
	elif CDS == True:
		outfile.write("charset first  = 1-%d\\3;\n" % masterDict[aln1]['nchar'])
		outfile.write("charset second = 2-%d\\3;\n" % masterDict[aln1]['nchar'])
		outfile.write("charset third  = 3-%d\\3;\n" % masterDict[aln1]['nchar'])
		outfile.write("charset fourth  = %d-%d\\3;\n" %(masterDict[aln1]['nchar']+1, masterDict[aln2]['nchar']))
		outfile.write("charset fifth = %d-%d\\3;\n" %(masterDict[aln1]['nchar']+1, masterDict[aln2]['nchar']))
		outfile.write("charset sixth  = %d-%d\\3;\n" %(masterDict[aln1]['nchar']+1, masterDict[aln2]['nchar']))
		outfile.write("partition cds = 6: first, second, third, fouth, fifth, sixth;\n")
		outfile.write("set partition = cds;\n")
		outfile.write("unlink shape=(all) statefreq=(all) revmat=(all) brlen=(all);\n")
		outfile.write("link brlens=(1,2,3);\n")
		outfile.write("link brlens=(4,5,6);\n")
		outfile.write("prset ratepr=dirichlet(10.0,10.0,10.0,10.0,10.0,10.0);\n")
	#outfile.write("lset applyto=(all) nst=%s ngammacat=4 rates=%s;\n" %(mrBayesNST, mrBayesRates))
	outfile.write("Lset  applyto=(all) nst=6  rates=gamma ngammacat=4;\n")
	#outfile.write("set seed=9223 swapseed=9223;")
	outfile.write("set seed=5207 swapseed=5207;\n")
	#outfile.write("prset applyto=(all) statefreqpr=dirichlet(10.0,10.0,10.0,10.0) revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) brlenspr=Unconstrained:GammaDir(1.0,0.100,1.0,1.0) shapepr=exponential(1.0);\n")
	outfile.write("Prset applyto=(all) brlenspr=unconstrained:GammaDir(1.0, 0.1,1.0, 1.0) shapepr=exp(1.0) Revmatpr=dirichlet(10.0,20.0,10.0,10.0,20.0,10.0) statefreqpr=dirichlet(10.0,10.0,10.0,10.0);\n")
	#outfile.write("mcmcp ngen=%s samplefreq=%s printfreq=%s nruns=1 starttree=random nchains=1 savebrlens=yes;\n" %(mrBayesNgen, mrBayesSampleFreq, mrBayesNgen))
	outfile.write("mcmc ngen=10000000 samplefreq=100 printfreq=10000 nchains=1 nruns=1 filename=%s/alignments/nexus/%s-%s_brlen-unlinked.mcmc;\n" %(outputDir, locus1, locus2))
	outfile.write("sumt filename=%s/alignments/nexus/%s-%s_brlen-unlinked.mcmc;\n" %(outputDir, locus1, locus2))
	outfile.write("mcmcp filename=%s/alignments/nexus/%s-%s_brlen-unlinked.ss;\n" %(outputDir, locus1, locus2))
	#outfile.write("ss alpha=0.3 nsteps=%s burninss=-2;\n" % mrBayesNsteps)
	outfile.write("ss alpha=0.3 nsteps=30 burninss=-2;\n")
	outfile.write("sumt conformat=simple;\n")
	outfile.write("sump;\n")
	#outfile.write("sumt burnin=2500;\n")
	outfile.write("log stop;\n")
	outfile.write("end;\n")
	outfile.close()
