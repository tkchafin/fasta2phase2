#!/usr/bin/python

import re
import sys
import os
import getopt
import Bio
import alignment_tools as aln
from itertools import product
from Bio import AlignIO

def main():
	params = parseArgs()

	#If fasta provided:
	if params.fasta:
		print("Reading FASTA...")
		#Initialize empty AlignIO object
		locus = Bio.Align.MultipleSeqAlignment([])

		#Add each sequence to the AlignIO
		for seq in read_fasta(params.fasta):
			#Add to AlignIO alignment
			locus.add_sequence(seq[0], seq[1])

		#Generate catalog of variable positions
		alignment = aln.consensAlign(locus, threshold=1.0, mask=1.0)

		#Print header information to output file
		out_fh = open(params.out, "w")

		print("Calculating variable columns...")
		#For each variable column, create outputs for positions and types:
		positions = []
		types = []
		for var in alignment.alnVars:
			positions.append(var.position+1) #add 1 because these are 0-based
			this_type = ""
			if var.value.upper() in ['B', 'D', 'H', 'V']:
				types.append("M")
			else:
				types.append("S")

		out_fh.write(str(len(locus)))
		out_fh.write("\n")
		out_fh.write(str(len(positions))) #write length of polymorphic alignment only
		out_fh.write("\n")

		posline = "P " + " ".join(str(x) for x in positions) + "\n"
		out_fh.write(posline)
		out_fh.write(str("".join(types)))
		out_fh.write("\n")

		#If there are non-diallelic SNPs, print warning:
		if "M" in types:
			print("Warning: There are non-diallelic SNPs in your dataset. PHASE requires these are output in a different format.")
			print("Coding multi-allelic SNPs as integer: [A=0, G=1, C=2, T=3,N=-1]")

		print("Expanding sample sequences...")
		#For each sample, create output lines
		for samp in locus:
			line1 = []
			line2 = []
			for idx, pos in enumerate(positions):
				if types[idx] == "M":
					expand = get_iupac_mult(samp.seq[pos-1]) #convert back to 0-based indexing
				else:
					expand = get_iupac_dip(samp.seq[pos-1])
				line1.append(expand[0])
				line2.append(expand[1])
			#Outputs for sample
			out1 = " ".join(line1) + "\n"
			out2 = " ".join(line2) + "\n"
			name = samp.description + "\n"
			out_fh.write(name)
			out_fh.write(out1)
			out_fh.write(out2)
		print("Done! Output can be found in",params.out)

	else:
		sys.exit("No input provided.")


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'f:o:h', \
			["fasta=","out=","help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.fasta=None
		self.out=None

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
			if opt in ('f', 'fasta'):
				self.fasta = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasta:
			self.display_help("Error: Missing required input file (-f,--fasta).")

		if self.out:
			self.out = self.out + ".inp"
		else:
			self.out = "out.inp"


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfasta2phase.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-f </path/to/fasta  <-o out_prefix>\n")
		print ("Description: Formats a FASTA sequence alignment for haplotype reconstruction in PHASE2")

		print("""
	Input options:
		-f,--fasta	: FASTA file with one consensus sequencer per infividual
		-o,--out	: Prefix for output file <default = ./out>
		-h,--help	: Displays help menu""")
		print()
		sys.exit()

#Function to split character to IUPAC codes, assuming diploidy
def get_iupac_dip(char):
	lower = False
	iupac = {
		"A"	: ["A","A"],
		"G"	: ["G","G"],
		"C"	: ["C","C"],
		"T"	: ["T","T"],
		"N"	: ["?","?"],
		"-"	: ["?","?"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: ["?","?"],
		"D"	: ["?","?"],
		"H"	: ["?","?"],
		"V"	: ["?","?"]
	}
	ret = iupac[char]
	return ret

#Function to split character to IUPAC codes, assuming diploidy
def get_iupac_mult(char):
	lower = False
	iupac = {
		"A"	: ["0","0"],
		"G"	: ["1","1"],
		"C"	: ["2","2"],
		"T"	: ["3","3"],
		"N"	: ["-1","-1"],
		"-"	: ["-1","-1"],
		"R"	: ["0","1"],
		"Y"	: ["2","3"],
		"S"	: ["1","2"],
		"W"	: ["0","3"],
		"K"	: ["1","3"],
		"M"	: ["0","2"],
		"B"	: ["-1","-1"],
		"D"	: ["-1","-1"],
		"H"	: ["-1","-1"],
		"V"	: ["-1","-1"]
	}
	ret = iupac[char]
	return ret

#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
def read_fasta(fas):
	if not fileCheck(fas):
		raise FileNotFoundError("Fatal exception, file %s not found."%fas)

	fh = open(fas)
	try:
		with fh as file_object:
			contig = ""
			seq = ""
			for line in file_object:
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
	finally:
		fh.close()

#Function to check if a file path is valid
def fileCheck(f):
	return (os.path.isfile(f))

#Function to expand ambiguous sequences
#Generator function
def expandAmbiquousDNA(sequence):
   for i in product(*[get_iupac_caseless(j) for j in sequence]):
      yield("".join(i))


#Call main function
if __name__ == '__main__':
    main()
