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
	full_seq = False
	cons_sequences = dict()
	positions = []
	
	#If fasta provided, parse it for consensus sequences and variable positions 
	if params.fasta:
		print("Haplotypes will be exported as full sequences")
		full_seq = True
		#If fasta provided:
		if params.fasta:
			
			print("Reading FASTA...")
			#Initialize empty AlignIO object
			locus = Bio.Align.MultipleSeqAlignment([])

			#Add each sequence to the AlignIO
			for seq in read_fasta(params.fasta):
				#Add to AlignIO alignment
				locus.add_sequence(seq[0], seq[1])
				#also population our dictionary 
				cons_sequences[seq[0]] = seq[1]

			#Generate catalog of variable positions
			alignment = aln.consensAlign(locus, threshold=1.0, mask=1.0)
			
			#Get variable positions from alignment
			for var in alignment.alnVars:
				positions.append(var.position) 
		
	else:
		print("Haplotypes will be exported as variable columns only.")
		
	#parse pairs file 
	if params.pairs:
		print("Reading pairs file...")
		#Open pairs file for reading 
		if os.path.exists(params.pairs):
			with open(params.pairs, 'r') as f:
				try:
					ofh = open(params.out, "w")
					this_ind = None
					assigns = []
					for line in f:
						line = line.strip()
						if not line:
							continue 
						#Check if this is an individual header, or hap data
						if line[0] == "I":
							if this_ind:
								#If an ind already assigned, parse it's list
								best = chooseDiplotype(assigns, params.minp)
								if best:
									if full_seq:
										if this_ind in cons_sequences:
											full_dip = getFullDiplotype(cons_sequences[this_ind], best, positions)
											if full_dip:
												out1 = ">" + this_ind + "_A" + "\n" + full_dip[0] + "\n"
												out2 = ">" + this_ind + "_B" + "\n" + full_dip[1] + "\n"
												ofh.write(out1)
												ofh.write(out2)
											else:
												print("Oh no! Something went wrong with getFullDiplotype()! My creator was too lazy to implement better error check though, so I don't know what went wrong!")
										else:
											print("Warning: Individual %s doesn't seem to be in the FASTA file... Skipping it."%this_ind)
									else:
										#print("PRINTING IND",this_ind)
										#print best to file 
										out1 = ">" + this_ind + "_A" + "\n" + best[0] + "\n"
										out2 = ">" + this_ind + "_B" + "\n" + best[1] + "\n"
										ofh.write(out1)
										ofh.write(out2)
								else:
									print("Warning: No diplotypes of probability >= %s for sample %s"%(params.minp, this_ind))
								#Then, empty assigns and this_ind variables 
								words = line.split()
								this_ind = words[1] #Set new individual, reset assigned haps 
								del assigns[:]
							else:
								words = line.split()
								this_ind = words[1]
								del assigns[:]
								continue 
						else:
							#This line must contain data 
							stuff = line.replace(" ","").split(",")
							assigns.append(stuff)
					if this_ind:
						if len(assigns) > 0:
							best = chooseDiplotype(assigns, params.minp)
							if best:
								if full_seq:
									if this_ind in cons_sequences:
										full_dip = getFullDiplotype(cons_sequences[this_ind], best, positions)
										if full_dip:
											out1 = ">" + this_ind + "_A" + "\n" + full_dip[0] + "\n"
											out2 = ">" + this_ind + "_B" + "\n" + full_dip[1] + "\n"
											ofh.write(out1)
											ofh.write(out2)
										else:
											print("Oh no! Something went wrong with getFullDiplotype()! My creator was too lazy to implement better error check though, so I don't know what went wrong!")
									else:
										print("Warning: Individual %s doesn't seem to be in the FASTA file... Skipping it."%this_ind)
								else:
									#print("PRINTING IND",this_ind)
									#print best to file 
									out1 = ">" + this_ind + "_A" + "\n" + best[0] + "\n"
									out2 = ">" + this_ind + "_B" + "\n" + best[1] + "\n"
									ofh.write(out1)
									ofh.write(out2)
					ofh.close()
				except IOError: 
					print("Could not read file ",pairs)
					sys.exit(1)
				finally:
					f.close()

	else:
		print("No input provided")
		sys.exit(1)

#Returns full diplotype sequence from consensus sequence, and variable-only diplotye 
def getFullDiplotype(seq, v, pos):
	hap1 = seq
	hap2 = seq
	
	index = 0
	for i, p in enumerate(pos):
		p = int(p)
		#hap1 = hap1[:p] + "[" + v[0][i] + "]" + hap1[p+1:]
		#hap2 = hap2[:p] + "[" + v[1][i] + "]" + hap2[p+1:]
		hap1 = hap1[:p] + v[0][i] + hap1[p+1:]
		hap2 = hap2[:p] + v[1][i] + hap2[p+1:]
	return([hap1, hap2])


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


#Parses a list of diplotype assignments, and returns the best one 
#Only returns if diplotype has probabiltiy greater than thresh 
def chooseDiplotype(dips, thresh):
	best = None 
	current = 0.0
	if len(dips) == 1:
		return(dips[0])
	for dip in dips:
		if float(dip[2]) >= 0.0:
			current = float(dip[2])
			best = dip 
			
	if current < thresh:
		return(None)
	else:
		#print("chosen:",best)
		return(best)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'p:o:m:hf:', \
			["pairs=","out=","minp=","help","fasta="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.pairs=None
		self.out=None
		self.minp=0.5
		
		self.fasta=None

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
			if opt in ('p', 'pairs'):
				self.pairs = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ('m','minp'):
				self.minp = float(arg)
			elif opt in ("f","fasta"):
				self.fasta = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.pairs:
			self.display_help("Error: Missing required input file (-p,--pairs).")

		if self.out:
			self.out = self.out + "_pairs.fasta"
		else:
			self.out = "out_pairs.fasta"


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nparsePairs.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-p </path/to/.pairs  <-o out_prefix> <-m min_probability\n")
		print ("Description: Extracts phased haplotypes from PHASE output, using the '.pairs' file")

		print("""
	Arguments:
		-p,--pairs	: Path to .pairs output file of PHASE
		-o,--out	: Prefix for output file <default = ./out>
		-m,--minp	: Minimum posterior probabiltiy to retain diplotype <default=0.5>
		-h,--help	: Displays help menu
		
	Haplotype expansion:
	--Provide if you want haplotypes exported with full sequences.
	--Otherwise, they'll be output as ONLY the variable columns.
		-f,--fasta	: FASTA input file used with fasta2phase.py""")
		print()
		sys.exit()


#Function to check if a file path is valid
def fileCheck(f):
	return (os.path.isfile(f))


#Call main function
if __name__ == '__main__':
    main()
