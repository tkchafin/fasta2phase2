# haploTools
Python scripts to create input files for haplotype reconstruction in Phase, parse outputs, and perform some other manipulations with haplotypes/ diplotypes

### Dependencies

These scripts require Python3 and Biopython. The easiest way to install these is through either pip3:

```
$ pip3 install biopython
```
or conda:
```
$ conda install -c conda-forge biopython
```

### Usage 
All scripts are accessible via a command-line interface. You can view options by calling the script with the -h flag, which displays a help menu:
```
$ python3 ./parsePairs.py

Exiting because help menu was called.

parsePairs.py

Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu

Usage:  parsePairs.py -p </path/to/.pairs  <-o out_prefix> <-m min_probability

Description: Extracts phased haplotypes from PHASE output, using the '.pairs' file

	Arguments:
		-p,--pairs	: Path to .pairs output file of PHASE
		-o,--out	: Prefix for output file <default = ./out>
		-m,--minp	: Minimum posterior probabiltiy to retain diplotype <default=0.5>
		-h,--help	: Displays help menu
		
	Haplotype expansion:
	--Provide if you want haplotypes exported with full sequences.
	--Otherwise, they'll be output as ONLY the variable columns.
		-f,--fasta	: FASTA input file used with fasta2phase.py
    
  ```
  
  Details on the functions and usage of each individual script can be found below.
  
  ### fasta2phase.py
  This script converts a FASTA sequence to the .inp input file required by the [PHASE](http://stephenslab.uchicago.edu/phase/download.html) haplotype phasing software.
  
  
