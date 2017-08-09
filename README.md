FaQCs: Quality Control of Next Generation Sequencing Data
===========================================================
![3D QC plot](http://oi61.tinypic.com/n36p9x.jpg)
-------------
PREREQUISITES
-------------

1. The main program is developed in Perl v 5.8.8.
2. Parallel::ForkManager module from CPAN   
   (http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.9/lib/Parallel/ForkManager.pm)
3. String::Approx module from CPAN   
   (http://search.cpan.org/~jhi/String-Approx-3.27/Approx.pm)
4. R for ploting                 
   (http://www.r-project.org/)                             
5. Jellyfish for kmer counting   (Optional) 
   (http://www.cbcb.umd.edu/software/jellyfish/) 

Note: The two Perl modules can be installed by INSTALL.sh script in the lib directory.

    cd lib
    ./INSTALL.sh

-----------
BASIC USAGE
-----------

* Trimming by quality 5 and filtering reads with any ambiguous base or low complexity.

  $ perl FaQCs.pl -p reads1.fastq reads2.fastq -d out_directory

* Quailty check only on subsamples of input, no trimming and filtering. 

  $ perl FaQCs.pl -p reads1.fastq reads2.fastq -d out_directory -qc_only 

-----------
Full USAGE
-----------
     
    Usage: perl FaQCs.pl [options] [-u unpaired.fastq] -p reads1.fastq reads2.fastq -d out_directory
    Version 1.36
    Input File: (can use more than once)
            -u            <Files> Unpaired reads
            
            -p            <Files> Paired reads in two files and separate by space
    Trim:
            -mode         "HARD" or "BWA" or "BWA_plus" (default BWA_plus)
                          BWA trim is NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)

            -q            <INT> Targets # as quality level (default 5) for trimming
    
            -5end         <INT> Cut # bp from 5 end before quality trimming/filtering 
      
            -3end         <INT> Cut # bp from 3 end before quality trimming/filtering 

            -adapter      <bool> Trim reads with illumina adapter/primers (default: no)
                          -rate   <FLOAT> Mismatch ratio of adapters' length (default: 0.2, allow 20% mismatches)
                          -polyA  <bool>  Trim poly A ( > 15 ) 
                          -keepshort  turn on this will keep short portion of reads instead of keep longer portion of reads
            					
            -artifactFile  <File>    additional artifact (adapters/primers/contaminations) reference file in fasta format 
    Filters:
            -min_L        <INT> Trimmed read should have to be at least this minimum length (default:50)

            -avg_q        <NUM> Average quality cutoff (default:0, no filtering)
            
            -n            <INT> Trimmed read has more than this number of continuous base "N" will be discarded. 
                          (default: 2, "NN") 

            -lc           <FLOAT> Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)

            -phiX         <bool> Filter phiX reads (slow)
            
    Q_Format:
            -ascii        Encoding type: 33 or 64 or autoCheck (default)
                          Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)

            -out_ascii    Output encoding. (default: 33)
    Output:
            -prefix       <TEXT> Output file prefix. (default: QC)

            -stats        <File> Statistical numbers output file (default: prefix.stats.txt)

            -d            <PATH> Output directory.
    Options:
            -t            <INT > # of CPUs to run the script (default:2 )

            -split_size   <INT> Split the input file into several sub files by sequence number (default: 1000000) 

            -qc_only      <bool> no Filters, no Trimming, report numbers.

            -kmer_rarefaction     <bool>   
                          Turn on the kmer calculation. Turn on will slow down ~10 times. (default:Calculation is off.)
                          (meaningless if -subset is too small)
                          -m  <INT>     kmer for rarefaction curve (range:[2,31], default 31)

            -subset       <INT>   Use this nubmer x split_size for qc_only and kmer_rarefaction  
                                  (default: 10,  10x1000000 SE reads, 20x1000000 PE reads)

            -discard      <bool> Output discarded reads to prefix.discard.fastq (default: 0, not output)
 
            -substitute   <bool> Replace "N" in the trimmed reads with random base A,T,C ,or G (default: 0, off)
 
            -trim_only    <bool> No quality report. Output trimmed reads only.
 
            -replace_to_N_q  <INT>  For NextSeq data, to replace base G to N when below this quality score (default:0, off)

            -5trim_off    <bool> Turn off trimming from 5'end.

            -debug        <bool> keep intermediate files

------------
Output Files
------------
Expected output files
- QC.1.trimmed.fastq
- QC.2.trimmed.fastq
- QC.unpaired.trimmed.fastq
- [QC.stats.txt](https://raw.githubusercontent.com/LANL-Bioinformatics/FaQCs/master/example/output/QC.stats.txt)
- QC_qc_report.pdf

[Example qc_report.pdf file](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4246454/bin/12859_2014_366_MOESM1_ESM.pdf)

---------------
VERSION HISTORY
---------------
======== Version 1.36
- add option "-replace_to_N_q" (replace base G to N when below this quality score)
- auto detect pacbio CCS reads quality encoding with score > 41
- fix major bugs: trim 5 end and trim 3 end does not work

======== Version 1.35
- add option "-polyA"  Trim polyA 
- update adapter trimming algorithm. If more than one adapter in a read, the whole read is filtered.

======== Version 1.34
- add option "-5trim_off    <bool> Turn off trimming from 5'end."
- add INSTALL.sh script for two requried perl modules installations.

======== Version 1.33
- input paired no need quote for exploit the autocomplete feature
- add trim_only option
- mode with  "HARD" or "BWA" or "BWA_plus" (default BWA_plus)

======== Version 1.32
- add -5end and -3end to cut x number base from 5' end or 3' end before quality trimming/filtering
- fix bug on phiX filtering with reverse complementary strand hit
- fix error when all reads in subsample are filtered/trimmed.

======== Version 1.31
- report raw subsample graphic results side-by-side with qc results for comparison.  

======== Version 1.3
- add -phiX to filter phiX reads
- add -substitute to replace "N" in the trimmed reads with random base A,T,C ,or G
- change -adapter behavior from filtering to trimming
- change -n behavior from # of tolerance to number of continuous base "N" filtering

======== Version 1.2
- add -adapter and -artifactFile for filtering reads with Adapters/Primers and other contaminations
- require String::Approx module from CPAN for above function

======== Version 1.1
New features and changes in illumina_fastq_qc  version 1.1 with respect to version 1.0:
- add -qc_only option for quick quality check without trimming and filtering
- add -discard option to output discarded reads

======== Version 1.0
Stable function release.
Features:
- trim bidirection
- minimium length filtering after trim
- "N" base filter
- low complexity filter
- average read quality filter
- autocheck quality encoding and quality encoding coversion
- multi-threads  (required Parallel::ForkManager)
- input paired end reads aware

-------------
CITATION
-------------

Chienchi Lo, PatrickS.G. Chain (2014) Rapid evaluation and Quality Control of Next Generation Sequencing Data with FaQCs. [BMC Bioinformatics. 2014 Nov 19;15 ](http://www.ncbi.nlm.nih.gov/pubmed/25408143)

