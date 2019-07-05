# haploTools
Python scripts to facilitate haplotype phasing, parse outputs of PHASE, and calculate haplotype and/or diplotype frequencies for populations or sample groups.

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

This script parses FASTA file of haplotypes (or diplotypes formatted as in the output of parsePairs.py). 
```
$ python3 ./hapFrequency.py -h
Exiting because help menu was called.

hapFrequency.py

Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu

Usage:  ./hapFrequency.py -f </path/to/fasta>  <-o out_prefix> <-p /path/to/popmap>

Description: Extracts haplotype frequencies from a FASTA file

	Arguments:
		-f	: Path to FASTA file. If diplotypes, separate with _A and _B
		-o	: Prefix for output file <default = ./out>
		-p	: Path to tab-delimited population map, if you want to group
			  --NOTE: Output table will include column "TOTAL" for all samples (inc. not in popmap)
			          and a column "ALLPOPS" which will include only samples from the popmap.
			          This is intended to make it easy to get frequencies/ counts for subsets.
		-d	: Output frequencies/ counts of diplotypes
		-c	: Path to a tab-delimited table for coding haplotypes
		-h	: Displays help menu
```
As you can see, there are options provided for you to provide an input fasta file <-f> and output prefix <-o>. The script will produce two tables from the input FASTA: a counts file containing absolute haplotype counts, and a freq file containing haplotype frequencies. You can also optionally calculate counts/ frequencies of diplotypes, provided that A and B designations are given to associate haplotype pairs. 

You can also optionally provide a tab-delimited file coding samples (haplotypes) into populations or groups for which you want to calculate frequencies. This file should be formatted like so:
```
sample1_A	Pop1
sample1_B	Pop1
sample2_A	Pop3
...
...
```

The "counts.tsv" file will look like this:
```
Haplotype       TOTAL   4-4.n   3-3.n   fawn    1-1.n   5-5.5n ...    ALLPOPS
ATGTGTGAA      657     17      15      14      20      6       52     ...  197
ATGGTGAAT       435     3       20      14      8       9       14     ...   95
ATGGGTGAC       426     3       9       8       5       6       18     ...   59
AGGTCAATA      309     4       13      5       9       4       13     ...   65
...
...
...
```
Here, the TOTAL column shows counts for ALL samples, including those which were excluded from the input population map, and an ALLPOPS column which just shows totals for the included samples. This is intended so that you can easily subset samples while also comparing numbers to the population counts. 

The freqs.tsv file follows a similar format:
```
Haplotype       TOTAL   4-4.n   3-3.n   fawn    1-1.n   5-5.5n  2-2.n   ...     ALLPOPS
ATGTGTGAA       0.2292  0.3864  0.2027  0.2000  0.3125  0.1579  0.3662  0.3571  ...  0.3137
ATGGTGAAT       0.1518  0.0682  0.2703  0.2000  0.1250  0.2368  0.0986  0.1429  ...  0.1513
ATGGGTGAC       0.1486  0.0682  0.1216  0.1143  0.0781  0.1579  0.1268  0.1429  ...  0.0939
...
...
...
```
To simlify the output, you can also provide a tab-delimited file specifying names for haplotypes:
```
ATGGGTGAC	HapA
ATGGTGAAT 	HapB
...
...
```

In which case, the output would reflect the specified coding:
```
Haplotype       TOTAL   4-4.n   3-3.n   fawn    1-1.n   5-5.5n  2-2.n   ...     ALLPOPS
HapB       0.2292  0.3864  0.2027  0.2000  0.3125  0.1579  0.3662  0.3571  ...  0.3137
HapC       0.1518  0.0682  0.2703  0.2000  0.1250  0.2368  0.0986  0.1429  ...  0.1513
HapA       0.1486  0.0682  0.1216  0.1143  0.0781  0.1579  0.1268  0.1429  ...  0.0939
...
...
...
```
