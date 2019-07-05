#!/usr/bin/python

import re
import sys
import os
import getopt
import Bio
import alignment_tools as aln
from Bio import AlignIO

def main():
	params = parseArgs()

	#Now, get the alignment from the FASTA file (as another dict)
	if params.fasta:
		pop_assign = dict()
		seqs = dict()
		#print('Reading alignment from FASTA...')

		#check if there is a certain way user wants haps to be coded
		coding=dict()
		if params.code_table:
			coding=parsePopmap(params.code_table)

		for f in read_fasta(params.fasta):
			if coding[f[1]]:
				seqs[f[0]] = coding[f[1]]




	else:
		print("No input provided")
		sys.exit(1)



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



#function reads a tab-delimited popmap file and return dictionary of assignments
def parsePopmap(popmap):

	ret = dict()
	if os.path.exists(popmap):
		with open(popmap, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					else:
						stuff = line.split()
						ret[stuff[0]] = stuff[1]
					#print(line)
				return(ret)
			except IOError:
				print("Could not read file ",popmap)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%popmap)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'o:hf:c:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.out=None
		self.fasta=None

		self.code_table=None

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
			if opt == "p":
				self.popmap= arg
			elif opt == "h":
				pass
			elif opt == "o":
				self.out = arg
			elif opt == "f":
				self.fasta = arg
			elif opt == "c":
				self.code_table = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasta:
			self.display_help("Error: Missing required input file (-f).")
		if not self.code_table:
			self.display_help("Error: Missing required input file (-c).")

		if self.out:
			self.out = self.out + ".txt"
		else:
			self.out = "out.txt"



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nmakeHapMap.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("Description: Writes a tab-delimited table of haplotype codes from a fasta and tsv hap codes list")

		print("""
	Arguments:
		-f	: Path to FASTA file. If diplotypes, separate with _A and _B
		-o	: Prefix for output file <default = ./out>
		-c	: Path to a tab-delimited table for coding haplotypes
		-a	: (Optional) table to append hap codes to (as a new column)
		-h	: Displays help menu
		""")
		print()
		sys.exit()


#Function to check if a file path is valid
def fileCheck(f):
	return (os.path.isfile(f))


#Call main function
if __name__ == '__main__':
    main()
