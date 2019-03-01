FaQCs: Quality Control of Next Generation Sequencing Data . [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/faqcs/README.html)
===========================================================
![3D QC plot](http://oi61.tinypic.com/n36p9x.jpg)

#### FaQCs versions 2.x was rewritten by Jason Gans (jgans at lanl.gov) wiht C++. It speeds up more than 10x compared to version 1.x in Perl script.

-------------
PREREQUISITES
-------------

1. zlib version >= 1.2.8
   (https://zlib.net)
2. R for ploting                 
   (http://www.r-project.org/)                             
3. Jellyfish for kmer counting   (Optional) 
   (http://www.cbcb.umd.edu/software/jellyfish/) 


----------------------
COMPILE / INSTALLATION
----------------------
*  The FaQCs is written in C++ and zlib library is required to complie from source. To complie, `cd` into the source direcotry and type `make`

```
   $ git clone https://github.com/LANL-Bioinformatics/FaQCs.git
   $ cd FaQCs
   $ make
 
   # A FaQCs binary executable will be ready to use and it can 
   # be moved to user's or system PATH environment.
```

* Alternatively, use conda to install
```   
   $ conda install -c bioconda -c conda-forge faqcs
```
** Trimming only comparison (--trim_only)

![comparison](https://github.com/LANL-Bioinformatics/EDGE/blob/gh-pages/images/FaQCs_performance.png)



-----------
BASIC USAGE
-----------

* Trimming by quality 5 and filtering reads with any ambiguous base or low complexity.

  $ FaQCs -p reads1.fastq reads2.fastq -d out_directory

* Quailty check only on subsamples of input, no trimming and filtering. 

  $ FaQCs -p reads1.fastq reads2.fastq -d out_directory -qc_only 



-----------
FULL USAGE
-----------

```     
Usage: FaQCs [options] [-u unpaired.fastq] -p reads1.fastq reads2.fastq -d out_directory
Version 2.09
Input File(s):
	-u			<File> Unpaired reads
	-1			<File> First paired read file
	-2			<File> Second paired read file
Trim:
	--mode			"HARD" or "BWA" or "BWA_plus" (default BWA_plus)
				BWA trim is NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)
	-q			<INT> Targets # as quality level (default 5) for trimming
	--5end			<INT> Cut # bp from 5 end before quality trimming/filtering
	--3end			<INT> Cut # bp from 3 end before quality trimming/filtering
	--adapter		<bool> Trim reads with illumina adapter/primers (default: no)
	--rate			<FLOAT> Mismatch ratio of adapters' length (default: 0.2, allow 20% mismatches)
	--polyA			<bool>  Trim poly A ( > 15 )
	--artifactFile		<File> additional artifact (adapters/primers/contaminations) reference file in fasta format
Filters:
	--min_L			<INT> Trimmed read should have to be at least this minimum length (default:50)
	--avg_q			<NUM> Average quality cutoff (default:0, no filtering)
	-n			<INT> Trimmed read has greater than or equal to this number of continuous base "N" will be discarded.
				(default: 2, "NN")
	--lc			<FLOAT> Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)
	--phiX			<bool> Filter phiX reads (slow)
Q_Format:
	--ascii			Encoding type: 33 or 64 or autoCheck (default)
				Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)
	--out_ascii		Output encoding. (default: 33)
Output:
	--prefix		<TEXT> Output file prefix. (default: QC)
	--stats			<File> Statistical numbers output file (default: prefix.stats.txt)
	-d			<PATH> Output directory.
Options:
	-t			<INT > # of CPUs to run the script (default:2 )
	--split_size		<INT> Split the input file into several sub files by sequence number (default: 1000000)
	--qc_only		<bool> no Filters, no Trimming, report numbers.
	--kmer_rarefaction	<bool>
				Turn on the kmer calculation. Turn on will slow down ~10 times. (default:Calculation is off.)
				(meaningless if -subset is too small)
	-m			<INT>     kmer for rarefaction curve (range:[2,31], default 31)
	--subset		<INT>   Use this nubmer x split_size for qc_only and kmer_rarefaction
				(default: 10,  10x1000000 SE reads, 20x1000000 PE reads)
	--discard		<bool> Output discarded reads to prefix.discard.fastq (default: 0, not output)
	--substitute		<bool> Replace "N" in the trimmed reads with random base A,T,C ,or G (default: 0, off)
	--trim_only		<bool> No quality report. Output trimmed reads only.
	--replace_to_N_q	<INT> Replace base G to N when below this quality score (default:0, off)
	--5trim_off		<bool> Turn off trimming from 5'end.
	--debug			<bool> Keep intermediate files
	--version		<bool> Print the version and exit
```

------------
OUTPUT FILES
------------
Expected output files
- QC.1.trimmed.fastq
- QC.2.trimmed.fastq
- QC.unpaired.trimmed.fastq
- [QC.stats.txt](https://raw.githubusercontent.com/LANL-Bioinformatics/FaQCs/master/example/output/QC.stats.txt)
- QC_qc_report.pdf

[Example qc_report.pdf file](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4246454/bin/12859_2014_366_MOESM1_ESM.pdf)

--------
CITATION
--------

Chienchi Lo, PatrickS.G. Chain (2014) Rapid evaluation and Quality Control of Next Generation Sequencing Data with FaQCs. [BMC Bioinformatics. 2014 Nov 19;15 ](http://www.ncbi.nlm.nih.gov/pubmed/25408143)

---------
COPYRIGHT
---------
Copyright (2018).  Triad National Security, LLC. All rights reserved.
 
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National 
Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National 
Nuclear Security Administration.
 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National 
Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, 
paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to 
the public, perform publicly and display publicly, and to permit others to do so.

This is open source software; you can redistribute it and/or modify it under the terms of the GPLv3 License. If software 
is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with 
the version available from LANL. Full text of the [GPLv3 License](https://github.com/LANL-Bioinformatics/FaQCs/blob/v2/LICENSE)  can be found in the License file in the main development 
branch of the repository.
