#!/usr/bin/python

import sys
import os
import getopt
import pandas as pd

def main():
	params = parseArgs()
	
	#read fasta
	seqs = dict()
	for f in read_fasta(params.infile):
		seqs[f[0]] = f[1]
	
	#read hap table
	freq = pd.read_csv(params.freq, sep="\t", header=0)
	
	#make new dict of sequences
	expanded=dict()
	for index, row in freq.iterrows():
		hap=row[0]
		idx=1
		for pop in row[1:]:
			if int(pop) > 0:
				for c in range(0,pop):
					#print(c+1, " of ", pop)
					key=str(freq.columns[idx]) + "_" + hap + "_" + str(c+1)
					expanded[key] = seqs[hap]
			idx+=1

	#write new fasta file
	write_fasta(params.out, expanded)

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
			options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:f:', \
			["help", "out=", "freq=", "infile="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		self.infile=None
		self.freq=None
		self.out="out.fas"


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt == "help":
				continue
			elif opt=="infile" or opt=="i":
				self.infile=arg
			elif opt=="freq" or opt=="f":
				self.freq=arg
			elif opt=="out" or opt=="o":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.infile or not self.freq:
			self.display_help("Inputs not provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nreverseHapFreq.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Un-collapses haplotypes to redundant sequences given a table of haplotype frequencies")
		print("""
		-i,--infile	: Input FASTA of haplotypes
		-f,--freq	: tsv Input table of haplotype frequencies
			NOTE: Rows should be haplotypes, and columns should be populations
			      Header should contain population names, and values are counts
				  Column 1 should contain (as row names) matching IDs from the FASTA 
		-o,--out	: Output file name (default=out.fas)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
