FaQCs: Quality Control of Next Generateion Sequenceing Data
=======

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

-----------
BASIC USAGE
-----------

* Trimming by quality 5 and filtering reads with any ambiguous base or low complexity.

  $ perl FaQCs.pl -p 'reads1.fastq reads2.fastq' -d out_directory

* Quailty check only on subsamples of input, no trimming and filtering. 

  $ perl FaQCs.pl -p 'reads1.fastq reads2.fastq' -d out_directory -qc_only 

---------------
VERSION HISTORY
---------------
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

