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
  This script converts a FASTA sequence to the .inp input file required by the [PHASE](http://stephenslab.uchicago.edu/phase/download.html) haplotype phasing software. Viewing the menu, you will see the following options:
  ```
  $ python2 ./fasta2phase.py
  
Exiting because help menu was called.

fasta2phase.py

Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu

Usage:  ./fasta2phase.py -f </path/to/fasta  <-o out_prefix>

Description: Formats a FASTA sequence alignment for haplotype reconstruction in PHASE2

	Input options:
		-f,--fasta	: FASTA file with one consensus sequencer per infividual
		-o,--out	: Prefix for output file <default = ./out>
		-h,--help	: Displays help menu
  ```
  
  Using the script is very simple. All you need to do is provide a FASTA-formatted input file, and an optional prefix for the output:
  ```
  $ python3 ./fasta2phase.py -f test.fasta -o phase_input
  ```
  This will produce a PHASE-formatted output called "phase_input.inp":
  ```
  $ cat phase_input.inp
  1437
12
P 29 122 212 254 255 293 336 347 407 468 524 645
SSSSSSSSSSSS
sample1
C C T A G A G G C A T C
C C T A G A G G C A T C
sample2
C C T A G A G G C A C C
C C T A G A G G C A T C
sample3
C C T A A A G G C A C C
C C T A G A G G C A T C
...
...
...
```
This .inp file can then be provided directly to PHASE.

### parsePairs.py
This script parses the ".pairs" output file from PHASE to extract diplotypes passing a user-determined posterior probability critical threshold. 

To give you an idea of what this script does, here is what the .pairs output from PHASE will look like:
```
IND: sample1
CCTAGAGGCATC , CCTAGAGGCATC , 1.000
IND: sample2
CCTAGAGGCATC , CCTAGAGGCACC , 1.000
IND: sample3
CCTAGAGGCATC , CCTAAAGGCACC , 0.048
CCTAGAGGCACC , CCTAAAGGCATC , 0.952
IND: sample4
CCTAAAGGCATC , CCTAAAGGCATC , 1.000
IND: sample5
TCTAGAGGCACC , TCTAGAGGCACC , 1.000
IND: sample6
CCTAGAGGCATC , CCTAAAGGCACC , 0.048
CCTAGAGGCACC , CCTAAAGGCATC , 0.952
...
...
...
```
Here, PHASE has written diplotype calls, as well as alternative possibilities, each provided with an associated posterior probability. To have parsePairs.py search this file, and keep all samples with diplotypes over PP 0.95, you could call it like so (assuming you named your output file like 'phase.output_pairs:

```
$ python3 ./parsePairs.py -p phase.output_pairs -m 0.95 -o phased_pairs.fasta
```
This will produce an output FASTA file, with haplotypes for each diplotype pair named with _A or _B:

### hapFrequency.py 

  
