# Setup functions
import subprocess
from Bio import SeqIO
from Bio import AlignIO

def clustalAlign(clustalPath, seqFile, outputDir, threads):
	filename = seqFile.split("/")[-1]
	locus = ".".join(filename.split(".")[0:-1])
	clustal_cline = "%s --auto --threads %s --in %s --out %s/alignments/%s.fasta" %(clustalPath, threads, seqFile, outputDir, locus)
	process = subprocess.Popen(clustal_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	process.wait()
	(out, err) = process.communicate() #the stdout and stderr

def mafftAlign(mafftPath, seqFile, outputDir, threads):
	filename = seqFile.split("/")[-1]
	locus = ".".join(filename.split(".")[0:-1])
	mafft_cline = "%s --auto --thread %s %s > %s/alignments/%s.fasta" %(mafftPath, threads, seqFile, outputDir, locus)
	process = subprocess.Popen(mafft_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	process.wait()
	(out, err) = process.communicate() #the stdout and stderr

def convertFastaToNexus(inputFile, outputDir, mrBayesNST = 6, mrBayesRates = "invgamma", mrBayesNgen = 10000000, mrBayesBurninFrac = 0.25, mrBayesSampleFreq = 1000, mrBayesPrintFreq = 10000, threads = 1 ):
	filename = inputFile.split("/")[-1]
	locus = ".".join(filename.split(".")[0:-1])
	AlignIO.convert(inputFile, "fasta", "%s/alignments/nexus/%s.nex" %(outputDir,locus), "nexus", molecule_type = "DNA")
	with open("%s/alignments/nexus/%s.nex" %(outputDir,locus), "a") as outfile:
		outfile.write("\n")
		outfile.write("begin mrbayes;\n")
		outfile.write("set autoclose=yes nowarn=yes;\n")
		#outfile.write("set Usebeagle=Yes Beagleresource=%s;\n" % threads)
		outfile.write("lset nst=%s rates=%s;\n" %(mrBayesNST, mrBayesRates))
		outfile.write("mcmc ngen=%s relburnin=yes burninfrac=%s samplefreq=%s printfreq=%s nchains=4 savebrlens=yes;\n" %(mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesPrintFreq))
		outfile.write("sump burnin=2500;\n")
		outfile.write("sumt burnin=2500;\n")
		#outfile.write("charset %s = 1-%s;\n" %(locus1, masterDict[aln1]['nchar']))
		#outfile.write("charset %s = %s-%s;\n" %(locus2, masterDict[aln1]['nchar']+1, totalChar))
		#outfile.write("partition combined = 2: %s,%s;\n" %(locus1, locus2))
		outfile.write("end;\n")
		outfile.close()



def concatenateAlignments(aln1, aln2, outputDir, mrBayesNST = 6, mrBayesRates = "invgamma", mrBayesNgen = 10000000, mrBayesBurninFrac = 0.25, mrBayesSampleFreq = 1000, mrBayesPrintFreq = 10000, threads = 1 ):
	filename1 = aln1.split("/")[-1]
	locus1 = ".".join(filename1.split(".")[0:-1])
	filename2 = aln2.split("/")[-1]
	locus2 = ".".join(filename2.split(".")[0:-1])
	outfile = open("%s/alignments/nexus/%s-%s.nex" %(outputDir,locus1,locus2), 'w')
	masterDict = {}
	filecounter = 1
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
			#if linecounter >=5 and linecounter == (ntax + 9):
			#	print line
			#	splitline2 = line.strip(';\n').split(' ')
			#	splitchar = splitline2[1].split('=')
			#	nchar = int(splitchar[1])
			#	infiledict[file]['nchar'] = nchar
				#print "Number of characters: %s" % nchar
			#f linecounter >=5 and linecounter >= (ntax + 12):
			#	splitseq = line.strip('\n').split('\t')
			#	if len(splitseq) < 3:
			#		break
			#	seqlist[splitseq[1]] = splitseq[2]
			#infiledict[file]['seqlist'] = seqlist
			#infiledict[file]['taxalist'] = taxalist
			#print infiledict
		#print taxalist
		#print seqlist
		#print infiledict
		filecounter += 1
	#print(masterDict)
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
	#for key in infiledict[filelist[0]]['seqlist']:
	#	if key in infiledict[filelist[1]]['seqlist'].keys():
	#		outfile.write("%s %s%s\n" %(key, infiledict[filelist[0]]['seqlist'][key[0:]], infiledict[filelist[1]]['seqlist'][key[0:]]))
	#	else:
	#		print("ERROR: NAMES IN INPUT FILES DO NOT MATCH!")
	#		exit(1)
			#print key
		#print infiledict[filelist[0]]['seqlist'][key[0:]]
		#print infiledict[filelist[1]]['seqlist'][key[0:]]
		#print ' '
	#print infiledict
	outfile.write(";\n")
	outfile.write("end;\n")
	outfile.write("\n")
	outfile.write("begin mrbayes;\n")
	outfile.write("set autoclose=yes nowarn=yes;\n")
	#outfile.write("set Usebeagle=Yes Beagleresource=%s;\n" % threads)
	outfile.write("lset nst=%s rates=%s;\n" %(mrBayesNST, mrBayesRates))
	outfile.write("mcmc ngen=%s relburnin=yes burninfrac=%s samplefreq=%s printfreq=%s nchains=4 savebrlens=yes;\n" %(mrBayesNgen, mrBayesBurninFrac, mrBayesSampleFreq, mrBayesPrintFreq))
	outfile.write("sump burnin=2500;\n")
	outfile.write("sumt burnin=2500;\n")
	outfile.write("charset %s = 1-%s;\n" %(locus1, masterDict[aln1]['nchar']))
	outfile.write("charset %s = %s-%s;\n" %(locus2, masterDict[aln1]['nchar']+1, totalChar))
	outfile.write("partition combined = 2: %s,%s;\n" %(locus1, locus2))
	outfile.write("end;\n")
	outfile.close()
