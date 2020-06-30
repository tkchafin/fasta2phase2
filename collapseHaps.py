#!/usr/bin/python

import re
import sys
import os
import getopt
from collections import OrderedDict

def main():
	params = parseArgs()

	#Now, get the alignment from the FASTA file (as another dict)
	if params.fasta:

		seqs = OrderedDict()
		for f in read_fasta(params.fasta):
			seqs[f[0]] = f[1]
		
		haps = getHapMap(seqs)
		
		writeHapMap(haps, params.out, params.number)

#Function to write fasta-formatted sequences
def write_fasta(f, aln, width=None):
	with open(f, 'w') as fh:
		try:
			for samp in aln.keys():
				if width:
					ol = ">" + str(samp) + "\n"
					chunks=wrap(aln[samp], width=width, break_on_hyphens=False, drop_whitespace=False)
					for chunk in chunks:
						ol=ol + str(chunk) + "\n"
				else:
					ol = ">" + str(samp) + "\n" + str(aln[samp]) + "\n"
				fh.write(ol)
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()

def write_dict(f, d):
	with open(f, 'w') as fh:
		try:
			for samp in d.keys():
				ol = str(samp) + "\t" + str(d[samp]) + "\n"
				fh.write(ol)
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()

def writeHapMap(haps, prefix, number):
	fasta=dict()
	hapmap=dict()
	count=1
	for h in haps:
		if number: 
			name="Haplotype_"+str(count)
		else:
			name=haps[h][0]
		fasta[name] = h
		for sample in haps[h]:
			hapmap[sample] = name
		count+=1
	write_fasta(str(prefix + ".fasta"), fasta)
	write_dict(str(prefix + ".hapmap"), hapmap)

def getHapMap(seqs):
	haps=dict()
	for ind in seqs:
		s=seqs[ind].upper()
		if s not in haps.keys():
			haps[s] = list()
			haps[s].append(ind)
		else:
			haps[s].append(ind)
	return(haps)


#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
def read_fasta(fas):

	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#print(line)
					if line[0] == ">": #Found a header line
						#If we already loaded a contig, yield that contig and
						#start loading a new one
						if contig:
							yield([contig,seq]) #yield
							contig = "" #reset contig and seq
							seq = ""
						split_line = line.split()
						contig = (split_line[0].replace(">",""))
					else:
						seq += line
				#Iyield last sequence, if it has both a header and sequence
				if contig and seq:
					yield([contig,seq])
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'o:hf:n', \
			["help", "out=", "fasta=", "number"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.out="out"
		self.fasta=None
		self.number=False

		#First pass to see if help menu was called
		for o, a in options:
			if o =="-h" or o=="--help":
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt=="help":
				pass
			elif opt == "o" or opt=="out":
				self.out = arg
			elif opt == "f" or opt=="fasta":
				self.fasta = arg
			elif opt=="n" or opt=="number":
				self.number=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasta:
			self.display_help("Error: Missing required input file (-f).")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ncollapseHaps.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("Description: Collapses redundant sequences from a FASTA and outputs haplotype fasta + table of samples matching each haplotype")

		print("""
	Arguments:
		-f,--fasta	: Path to FASTA file. 
		-o,--out	: Prefix for output (default=out)
		-n,--number	: Number haplotypes (default=names after 1st sample observed in)
		-h,--help	: Displays help menu
		""")
		print()
		sys.exit()


#Function to check if a file path is valid
def fileCheck(f):
	return (os.path.isfile(f))


#Call main function
if __name__ == '__main__':
    main()
