#include "FaQCs.h"

#include <iostream>
#include <deque>
#include <algorithm>

#include <getopt.h>
#include <limits.h>
#include <string.h> // for strncmp
#include <zlib.h>

using namespace std;

Options::Mode parse_mode(string m_mode /*copy*/);
unsigned int strtou(const string &m_buffer);
void parse_artifact_file(const string &m_filename, vector< pair<string, string> > &m_adapter);
string complement(string m_seq /* make a copy */);

struct sort_by_seq_length
{
	inline bool operator()(const pair<string, string> &m_a, const pair<string, string> &m_b)
	{
		// Sort in *descending* order
		return m_a.second.size() > m_b.second.size();
	};
};

// Where possible, keep the command line options the same as the perl script (version 1.35):
// Input File: (can use more than once)
//	-u			<Files> Unpaired reads    
//	-p			<Files> Paired reads in two files and separate by space
// Trim:
//	-mode			"HARD" or "BWA" or "BWA_plus" (default BWA_plus)
//				BWA trim is NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)
//	-q			<INT> Targets # as quality level (default 5) for trimming
//	-5end			<INT> Cut # bp from 5 end before quality trimming/filtering 
//	-3end			<INT> Cut # bp from 3 end before quality trimming/filtering 
//	-adapter		<bool> Trim reads with illumina adapter/primers (default: no)
//	-rate			<FLOAT> Mismatch ratio of adapters' length (default: 0.2, allow 20% mismatches)
//	-polyA			<bool>  Trim poly A ( > 15 ) 			
//	-artifactFile		<File> additional artifact (adapters/primers/contaminations) reference file in fasta format 
// Filters:
//	-min_L			<INT> Trimmed read should have to be at least this minimum length (default:50)
//	-avg_q			<NUM> Average quality cutoff (default:0, no filtering)
//	-n			<INT> Trimmed read has greater than or equal to this number of continuous base "N" will be discarded. 
//				(default: 2, "NN") 
//	-lc			<FLOAT> Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)
//	-phiX			<bool> Filter phiX reads (slow)
// Q_Format:
//	-ascii			Encoding type: 33 or 64 or autoCheck (default)
//				Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)
//	-out_ascii		Output encoding. (default: 33)
// Output:
//	-prefix			<TEXT> Output file prefix. (default: QC)
//	-stats			<File> Statistical numbers output file (default: prefix.stats.txt)
//	-d			<PATH> Output directory.
// Options:
//	-t			<INT > # of CPUs to run the script (default:2 )
//	-split_size		<INT> Split the input file into several sub files by sequence number (default: 1000000) 
//	-qc_only		<bool> no Filters, no Trimming, report numbers.
//	-kmer_rarefaction	<bool>   
//				Turn on the kmer calculation. Turn on will slow down ~10 times. (default:Calculation is off.)
//				(meaningless if -subset is too small)
//	-m			<INT>     kmer for rarefaction curve (range:[2,31], default 31)
//	-subset			<INT>   Use this nubmer x split_size for qc_only and kmer_rarefaction  
//				(default: 10,  10x1000000 SE reads, 20x1000000 PE reads)
//	-discard		<bool> Output discarded reads to prefix.discard.fastq (default: 0, not output)
//	-substitute		<bool> Replace "N" in the trimmed reads with random base A,T,C ,or G (default: 0, off)
//	-trim_only		<bool> No quality report. Output trimmed reads only.
//	-5trim_off		<bool> Turn off trimming from 5'end.
//	-debug			<bool> keep intermediate files
Options::Options(int argc, char* argv[])
{

	// Check to see if we need to handle the depricated '-p "file1 file2"' argument (which does not match the
	// getopt syntax). This test is included to maintain backwards compatibility with the origianl FaQCs.pl
	// perl script.
	for(int i = 1;i < argc;++i){
		
		if(strncmp(argv[i], "-p", 2 /*strlen("-p")*/) == 0){
			
			// We have a match!
			if( (i + 2) >= argc ){
				throw __FILE__ ":Options: Unable to extract file names after depricated '-p' flag";
			}
			
			// Make sure that there at least two arguments to the '-p' flag (do *not* test here to make
			// sure that these arguments are valid file names, as opposed to other command line flags).
			// Since Unix files can start with a '-', differentiating a flag from a file name is tricky!
			input_read1_file = argv[i + 1];
			input_read2_file = argv[i + 2];
			
			cerr << "The paired read flag (-p) has been deprecated. " 
				<< "Please specify paired read files with -1 <file> and -2 <file>" << endl;
		}
	}
	
	// Set default boolean parameter values
	print_usage = (argc == 1);
	protect_5 = false;
	replace_N = false;
	kmer_rarefaction = false;
	discard_output = false;
	qc_only = false;
	trim_only = false;
	filter_adapter = false;
	filter_phiX = false;
	debug = false;
	
	// Set the default enumerated parameter
	mode = BWA_plus;
	
	// Set default string parameter values
	prefix = "QC";
	
	// Set default floating point parameter values
	average_quality = 0.0;
	low_complexity_cutoff_ratio = 0.85;
	filterAdapterMismatchRate = 0.2;
	
	// Set default integer parameter values
	input_quality_offset = AUTO_DETECT_QUALITY_OFFSET;
	output_quality_offset = 33;
	num_thread = 0; // Use all threads
	quality = 5;
	min_read_length = 50;
	max_num_poly_N = 2;
	kmer = 31;
	num_subsample = 10;
	trim_5 = 0;
	trim_3 = 0;
	split_size = 1000000;
	replace_to_N_q = 0;
	
	bool version = false;
	bool trim_polyA = false;
	
	const char* options = "d:t:n:1:2:p:q:u:?h";
	int config_opt = 0;
	int long_index = 0;

	// Use (mostly) the same command line options as the
	// original FaQCs perl script to maintain backwards
	// compatibility
	struct option long_opts[] = {
		{"mode", true, &config_opt, 1},
		{"5end", true, &config_opt, 2},
		{"3end", true, &config_opt, 3},
		{"adapter", false, &config_opt, 4},
		{"rate", true, &config_opt, 5},
		{"polyA", false, &config_opt, 6},
		{"artifactFile", true, &config_opt, 7},
		{"min_L", true, &config_opt, 8},
		{"avg_q", true, &config_opt, 9},
		{"lc", true, &config_opt, 10},
		{"phiX", false, &config_opt, 11},
		{"ascii", true, &config_opt, 12},
		{"out_ascii", true, &config_opt, 13},
		{"prefix", true, &config_opt, 14},
		{"stats", true, &config_opt, 15},
		{"split_size", true, &config_opt, 16},
		{"qc_only", false, &config_opt, 17},
		{"kmer_rarefaction", false, &config_opt, 18},
		{"subset", true, &config_opt, 19},
		{"discard", false, &config_opt, 20},
		{"substitute", false, &config_opt, 21},
		{"trim_only", false, &config_opt, 22},
		{"5trim_off", false, &config_opt, 23},
		{"debug", false, &config_opt, 24},
		{"version", false, &config_opt, 25},
		{"R1", true, &config_opt, 26},
		{"R2", true, &config_opt, 27},
		{"Ru", false, &config_opt, 28},
		{"Rd", false, &config_opt, 29},
		{"QRpdf", false, &config_opt, 30},
		{"replace_to_N_q", true, &config_opt, 31},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;
	
	// Use getopt_long_only() instead of getopt_long() to enable long (word, as opposed to 
	// character base) command line options that start with '-' in addition to '--'. This is
	// to match the behavior of the original FaQCs.pl script.
	while( (opt_code = getopt_long_only( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:

				if(config_opt == 1){ // mode
				
					mode = parse_mode(optarg);
					break;
				}
				
				if(config_opt == 2){ // 5end
				
					trim_5 = strtou(optarg);
					break;
				}
				
				if(config_opt == 3){ // 3end
				
					trim_3 = strtou(optarg);
					break;
				}
				
				if(config_opt == 4){ // adapter
				
					filter_adapter = true;
					break;
				}
				
				if(config_opt == 5){ // rate
				
					filterAdapterMismatchRate = atof(optarg);
					break;
				}
				
				if(config_opt == 6){ // polyA
				
					trim_polyA = true;
					break;
				}
				
				if(config_opt == 7){ // artifactFile
				
					artifact_file = optarg;
					
					// Following the FaQCs.pl script, set the filter_adapter flag if
					// a user specifies an artifact file
					filter_adapter = true;
					break;
				}
				
				if(config_opt == 8){ // min_L
				
					min_read_length = strtou(optarg);
					break;
				}
				
				if(config_opt == 9){ // ave_q
				
					average_quality = atof(optarg);
					break;
				}
				
				if(config_opt == 10){ // lc
				
					low_complexity_cutoff_ratio = atof(optarg);
					break;
				}
				
				if(config_opt == 11){ // phiX
				
					filter_phiX = true;
					break;
				}
				
				if(config_opt == 12){ // ascii
				
					const int local = atoi(optarg);
					
					if( (local <= SCHAR_MIN) || (local > SCHAR_MAX) ){
						throw __FILE__ ":Options::Options: ascii out of bounds!";
					}
					
					input_quality_offset = char(local);
					break;
				}
				
				if(config_opt == 13){ // out_ascii
				
					const int local = atoi(optarg);
					
					if( (local <= SCHAR_MIN) || (local > SCHAR_MAX) ){
						throw __FILE__ ":Options::Options: out_ascii out of bounds!";
					}
					
					output_quality_offset = char(local);
					break;
				}
				
				if(config_opt == 14){ // prefix
				
					prefix = optarg;
					break;
				}
				
				if(config_opt == 15){ // stats
				
					stats_file = optarg;
					break;
				}
				
				if(config_opt == 16){ // split_size
				
					split_size = strtou(optarg);
					break;
				}
				
				if(config_opt == 17){ // qc_only
				
					qc_only = true;
					break;
				}
				
				if(config_opt == 18){ // kmer_rarefaction
				
					kmer_rarefaction = true;
					break;
				}
				
				if(config_opt == 19){ // subset
			
					num_subsample = strtou(optarg);
					break;
				}
				
				if(config_opt == 20){ // discard
			
					discard_output = true;
					break;
				}
				
				if(config_opt == 21){ // substitute
			
					replace_N = true;
					break;
				}
				
				if(config_opt == 22){ // trim_only
			
					trim_only = true;
					break;
				}
				
				if(config_opt == 23){ // 5trim_off
			
					protect_5 = true;
					break;
				}
				
				if(config_opt == 24){ // debug
			
					debug = true;
					break;
				}
				
				if(config_opt == 25){ // version
			
					version = true;
					break;
				}
				
				if(config_opt == 26){ // R1
				
					input_read1_file = optarg;
					break;
				}
				
				if(config_opt == 27){ // R2
				
					input_read2_file = optarg;
					break;
				}
				
				if(config_opt == 28){ // Ru
				
					input_unpaired_file = optarg;
					break;
				}
				
				if(config_opt == 29){ // Rd
				
					trimmed_discard_file = optarg;
					break;
				}
				
				if(config_opt == 30){ // QRpdf
				
					plots_file = optarg;
					break;
				}
				
				if(config_opt == 31){ // replace_to_N_q
				
					replace_to_N_q = strtou(optarg);
					break;
				}
				
				cerr << "Unknown flag!" << endl;
				break;
			case '1':
				input_read1_file = optarg;
				break;
			case '2':
				input_read2_file = optarg;
				break;
			case 'u':
				input_unpaired_file = optarg;
				break;
			case 'd':
				output_dir = optarg;
				break;
			case 'm':
				kmer = strtou(optarg);
				break;
			case 'n':
				max_num_poly_N = strtou(optarg);
				break;
			case 'q':
				{
					const int local = atoi(optarg);
					
					if( (local <= SCHAR_MIN) || (local > SCHAR_MAX) ){
						throw __FILE__ ":Options::Options: q out of bounds!";
					}
					
					quality = char(local);
				}
				break;
			case 't':
				num_thread = strtou(optarg);
				break;
			case 'h':
			case '?':
				print_usage = true;
				break;
			case 'p':
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				break;
		};
	}
	
	if(print_usage){
	
		cerr << "FaQCs version " << FaQCs_VERSION << endl;
		cerr << "Input File(s):" << endl;
		cerr << "\t-u\t\t\t<File> Unpaired reads" << endl;
		cerr << "\t-1\t\t\t<File> First paired read file" << endl;
		cerr << "\t-2\t\t\t<File> Second paired read file" << endl;
		cerr << "Trim:" << endl;
		cerr << "\t--mode\t\t\t\"HARD\" or \"BWA\" or \"BWA_plus\" (default BWA_plus)" << endl;
		cerr << "\t\t\t\tBWA trim is NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)" << endl;
		cerr << "\t-q\t\t\t<INT> Targets # as quality level (default 5) for trimming" << endl;
		cerr << "\t--5end\t\t\t<INT> Cut # bp from 5 end before quality trimming/filtering" << endl;
		cerr << "\t--3end\t\t\t<INT> Cut # bp from 3 end before quality trimming/filtering" << endl;
		cerr << "\t--adapter\t\t<bool> Trim reads with illumina adapter/primers (default: no)" << endl;
		cerr << "\t--rate\t\t\t<FLOAT> Mismatch ratio of adapters' length (default: 0.2, allow 20% mismatches)" << endl;
		cerr << "\t--polyA\t\t\t<bool>  Trim poly A ( > 15 )" << endl;
		cerr << "\t--artifactFile\t\t<File> additional artifact (adapters/primers/contaminations) reference file in fasta format" << endl;
		cerr << "Filters:" << endl;
		cerr << "\t--min_L\t\t\t<INT> Trimmed read should have to be at least this minimum length (default:50)" << endl;
		cerr << "\t--avg_q\t\t\t<NUM> Average quality cutoff (default:0, no filtering)" << endl;
		cerr << "\t-n\t\t\t<INT> Trimmed read has greater than or equal to this number of continuous base \"N\" will be discarded." << endl;
		cerr << "\t\t\t\t(default: 2, \"NN\")" << endl;
		cerr << "\t--lc\t\t\t<FLOAT> Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)" << endl;
		cerr << "\t--phiX\t\t\t<bool> Filter phiX reads (slow)" << endl;
		cerr << "Q_Format:" << endl;
		cerr << "\t--ascii\t\t\tEncoding type: 33 or 64 or autoCheck (default)" << endl;
		cerr << "\t\t\t\tType of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)" << endl;
		cerr << "\t--out_ascii\t\tOutput encoding. (default: 33)" << endl;
		cerr << "Output:" << endl;
		cerr << "\t--prefix\t\t<TEXT> Output file prefix. (default: QC)" << endl;
		cerr << "\t--stats\t\t\t<File> Statistical numbers output file (default: prefix.stats.txt)" << endl;
		cerr << "\t-d\t\t\t<PATH> Output directory." << endl;
		cerr << "Options:" << endl;
		cerr << "\t-t\t\t\t<INT > # of CPUs to run the script (default:2 )" << endl;
		cerr << "\t--split_size\t\t<INT> Split the input file into several sub files by sequence number (default: 1000000)" << endl;
		cerr << "\t--qc_only\t\t<bool> no Filters, no Trimming, report numbers." << endl;
		cerr << "\t--kmer_rarefaction\t<bool>" << endl;
		cerr << "\t\t\t\tTurn on the kmer calculation. Turn on will slow down ~10 times. (default:Calculation is off.)" << endl;
		cerr << "\t\t\t\t(meaningless if -subset is too small)" << endl;
		cerr << "\t-m\t\t\t<INT>     kmer for rarefaction curve (range:[2,31], default 31)" << endl;
		cerr << "\t--subset\t\t<INT>   Use this nubmer x split_size for qc_only and kmer_rarefaction" << endl;
		cerr << "\t\t\t\t(default: 10,  10x1000000 SE reads, 20x1000000 PE reads)" << endl;
		cerr << "\t--discard\t\t<bool> Output discarded reads to prefix.discard.fastq (default: 0, not output)" << endl;
		cerr << "\t--substitute\t\t<bool> Replace \"N\" in the trimmed reads with random base A,T,C ,or G (default: 0, off)" << endl;
		cerr << "\t--trim_only\t\t<bool> No quality report. Output trimmed reads only." << endl;
		
		// Not currently checking read platform
		//cerr << "\t--replace_to_N_q\t<INT> For NextSeq data, to replace base G to N when below this quality score (default:0, off)" << endl;
		cerr << "\t--replace_to_N_q\t<INT> Replace base G to N when below this quality score (default:0, off)" << endl;
		
		cerr << "\t--5trim_off\t\t<bool> Turn off trimming from 5'end." << endl;
		cerr << "\t--debug\t\t\t<bool> Keep intermediate files" << endl;
		cerr << "\t--version\t\t<bool> Print the version and exit" << endl;
		
		return;
	}
	
	if(version){
	
		cerr << "Version: " << FaQCs_VERSION << endl;
		
		// Set print_usage to true to that the program exits before
		// attempting to process any files.
		print_usage = true;
		return;
	}
	
	if( input_read1_file.empty() ^ input_read2_file.empty() ){
		
		if( input_read1_file.empty() ){
			cerr << "Please specify a read one file (-1 <file>)" << endl;
		}
		
		if( input_read2_file.empty() ){
			cerr << "Please specify a read one file (-2 <file>)" << endl;
		}
		
		print_usage = true;
		return;
	}
	else{
		// If we are performing kmer rarefaction with paired end inputs, we need to double 
		// num_subsample in orde to plot the correct number of points in the rarefaction curve
		num_subsample *= 2;
	}
	
	if( input_unpaired_file.empty() ){
		
		if( input_read1_file.empty() && input_read2_file.empty() ){
			
			cerr << "Please specify either a pair of fastq files (-1 <file> -2 <file>) "
				<< " or a single fastq file of unpaired reads (-u <file>)" << endl;
			
			print_usage = true;
			return;	
		}
	}
	
	if( (kmer < 2) || (kmer > 31) ){
		
		cerr << "Please specify a kmer value (-m) in the range 2 <= kmer <= 31" << endl;
		print_usage = true;
		return;
	}

	if( (low_complexity_cutoff_ratio > 1.0) || (low_complexity_cutoff_ratio < 0.0) ){
		
		cerr << "Please specify a low complexity cutoff ration (-lc) in the range 0 <= ratio <= 1.0" << endl;
		print_usage = true;
		return;
	}
	
	if( (filterAdapterMismatchRate > 1.0) || (filterAdapterMismatchRate < 0.0) ){
		
		cerr << "Please specify an adapter mismatch rate (-adapter) in the range 0 <= rate <= 1.0" << endl;
		print_usage = true;
		return;
	}
	
	if(split_size == 0){
		
		cerr << "Please specify a split_size (--split_size) value greater than 0" << endl;
		print_usage = true;
		return;
	}

	if(num_subsample == 0){
		
		cerr << "Please specify a subset (--subset) value greater than 0" << endl;
		print_usage = true;
		return;
	}
	
	if(replace_N){
		cerr << "**Warning** \"-substitue\" is not currently implemented" << endl;
	}
	
	if(filter_adapter){
		
		// Load the default adapters and check for any additional adapter sequences
		// Note that only the actual adapter sequences (and not the reverse complement
		// of the adapter sequences) are searched against the read data. This is an 
		// intentional design choice to be consistent with the behaviour of the 
		// original FaQCs.pl script.
		adapter.push_back( make_pair(
			"cre-loxp-forward", "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"
			) );
			
		adapter.push_back( make_pair(
			"cre-loxp-reverse", "AGCATATTGAAGCATATTACATACGATATGCTTCAATAATGC"
			) );
		
		adapter.push_back( make_pair(
			"TruSeq-adapter-1", "GGGGTAGTGTGGATCCTCCTCTAGGCAGTTGGGTTATTCTAGAAGCAGATGTGTTGGCTGTTTCTGAAACTCTGGAAAA"
			) );
		
		adapter.push_back( make_pair(
			"TruSeq-adapter-3", "CAACAGCCGGTCAAAACATCTGGAGGGTAAGCCATAAACACCTCAACAGAAAA"
			) );
		
		adapter.push_back( make_pair(
			"PCR-primer-1", "CGATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"
			) );
		
		adapter.push_back( make_pair(
			"PCR-primer-2", "GCATAACTTCGTATAGCATACATTATACGAAGTTATACGA"
			) );
		
		adapter.push_back( make_pair(
			"Nextera-primer-adapter-1", "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
			) );
		
		adapter.push_back( make_pair(
			"Nextera-primer-adapter-2", "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
			) );
		
		adapter.push_back( make_pair(
			"Nextera-junction-adapter-1", "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"
			) );
	}
	
	if(trim_polyA){
	
		adapter.push_back( make_pair(
			"polyA", "AAAAAAAAAAAAAAAAAAAA"
			) );
	}
	
	if(filter_phiX){
		
		// For simplicity, we treat the phiX sequence (and its complement) as another
		// adapter for the sake of trimming. This will require a little bit of data
		// manipulation when the time comes to report trimming statistics. This sequence
		// us taken from the FaQCs.pl script and is listed as: PhiX174 NC_001422
		const string phi_x_seq = 
			"GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGAC"
			"TGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAG"
			"TGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGAT"
			"GCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCG"
			"AAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCT"
			"TTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCG"
			"TCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAG"
			"TGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATT"
			"TTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCG"
			"TCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACAT"
			"TTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAAGGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGC"
			"CGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGC"
			"TAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTC"
			"TACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAAC"
			"CTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACA"
			"GACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGG"
			"TGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGC"
			"TGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCG"
			"CCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTC"
			"GTGATAAAAGATTGAGTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTT"
			"CAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTA"
			"TATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGAT"
			"GCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTGAATGGTCGCCATGATGGTGGTTATTATACC"
			"GTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGCCGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAAT"
			"CAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAG"
			"GCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACC"
			"CTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTG"
			"GACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGG"
			"TTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTT"
			"CACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGC"
			"AGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATC"
			"AGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTTCTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAG"
			"CTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTC"
			"CAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTA"
			"CTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATG"
			"TTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTA"
			"TAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATT"
			"ACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTG"
			"ATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTTATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTG"
			"CCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTC"
			"CTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGG"
			"ACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATA"
			"TGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAA"
			"TGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGC"
			"TGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA";

		adapter.push_back( make_pair(
			PHI_X, phi_x_seq
			) );
			
		adapter.push_back( make_pair(
			PHI_X_COMPLEMENT, complement(phi_x_seq)
			) );
	}
	
	if(artifact_file != ""){

		// Read the artifact file (in fasta format) and add the results to
		// the adapter sequence.
		parse_artifact_file(artifact_file, adapter);
	}
	
	if(!input_read1_file.empty() && !input_read2_file.empty() ){
		
		if( trimmed_read1_file.empty() ){
			trimmed_read1_file = output_dir + "/" + prefix + ".1.trimmed.fastq";
		}
		
		if( trimmed_read2_file.empty() ){
			trimmed_read2_file = output_dir + "/" + prefix + ".2.trimmed.fastq";
		}
		
		if( trimmed_unpaired_file.empty() ){
			trimmed_unpaired_file = output_dir + "/" + prefix + ".unpaired.trimmed.fastq";
		}
		
		if( trimmed_discard_file.empty() ){
			trimmed_discard_file = output_dir + "/" + prefix + ".discard.trimmed.fastq";
		}
	}
	
	if( !input_unpaired_file.empty() ){
		
		if( trimmed_unpaired_file.empty() ){
			trimmed_unpaired_file = output_dir + "/" + prefix + ".unpaired.trimmed.fastq";
		}
		
		if( trimmed_discard_file.empty() ){
			trimmed_discard_file = output_dir + "/" + prefix + ".discard.trimmed.fastq";
		}
	}
	
	if( plots_file.empty() ){
		plots_file = output_dir + "/" + prefix + "_qc_report.pdf";
	}
	
	if(!discard_output){
		trimmed_discard_file = "";
	}
	else{
		if( trimmed_discard_file.empty() ){
			trimmed_discard_file = output_dir + "/" + prefix + ".discard.fastq";
		}
	}
	
	if( stats_file.empty() ){
		stats_file = output_dir + "/" + prefix + ".stats.txt";
	}

	switch(mode){
		case Options::HARD:
		
			if(!qc_only){
				cerr << "Hard trimming is used." << endl;
			}
			break;
		case Options::BWA:
			
			if(!qc_only){
				cerr << "Bwa trimming is used." << endl;
			}
			break;
		case Options::BWA_plus:
		
			if(!qc_only){
				cerr << "Bwa extension trimming is used." << endl;
			}
			break;
		case Options::UNDEFINED_MODE:
			
			mode = Options::BWA_plus;
			
			if(!qc_only){
				cerr << "Not recognized mode. Bwa extension trimming algorithm is used." << endl;
			}
			
			break;
		default:
			throw __FILE__ ":Options: Undefined mode!";
	};
}

Options::Mode parse_mode(string m_mode /*copy*/)
{
	// Make the string lower case
	for(string::iterator i = m_mode.begin();i != m_mode.end();++i){
		*i = tolower(*i);
	}
	
	if(m_mode == "hard"){
		return Options::HARD;
	}
	
	if(m_mode == "bwa"){
		return Options::BWA;
	}
	
	if(m_mode == "bwa_plus"){
		return Options::BWA_plus;
	}
	
	return Options::UNDEFINED_MODE;
}

unsigned int strtou(const string &m_buffer)
{
	size_t ret = 0;
	int p = 1;
	
	for(string::const_reverse_iterator i = m_buffer.rbegin();i != m_buffer.rend();++i){
		
		if( !isdigit(*i) ){
			throw __FILE__ ":strtou: Invalid character";
		}
		
		ret += (*i - '0')*p;
		p *= 10;
	}
	
	if(ret > UINT_MAX){
		throw __FILE__ ":strtou: Overflow!";
	}
	
	return (unsigned int)(ret);
}

void parse_artifact_file(const string &m_filename, vector< pair<string, string> > &m_adapter)
{
	// Allow for compressed fasta files
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		cerr << "Unable to open " << m_filename << " for loading artifact sequences" << endl;
		throw "I/O error";
	}
	
	const int buffer_len = 4096;
	char buffer[buffer_len];
	
	string defline;
	deque<char> data;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		char* ptr = strchr(buffer, '>');
		
		if(ptr != NULL){
			
			if( !data.empty() ){
				
				m_adapter.push_back( make_pair( string(), string() ) );
				
				m_adapter.back().first = defline;
				m_adapter.back().second.assign( data.begin(), data.end() );
			}
			
			data.clear();
			
			// Skip the fasta header prefix ('>')
			++ptr;
			
			// Remove any end of line symbols
			for(char* p = ptr;*p != '\0';++p){
				if( (*p == '\n') || (*p == '\r') ){
					*p = '\0';
				}
			}

			// In cases where the defline is longer than the buffer, indicate truncation
			// with an ellipsis
			if( strlen(ptr) == (buffer_len - 1) ){
				defline = buffer + string("...");
			}
			else{
				defline = ptr;
			}

		}
		else{
			for(char* p = buffer;*p != '\0';++p){
				
				if( !isspace(*p) ){
					data.push_back(*p);
				}
			}
		}
	}
	
	if( !data.empty() ){
				
		m_adapter.push_back( make_pair( string(), string() ) );

		m_adapter.back().first = defline;
		m_adapter.back().second.assign( data.begin(), data.end() );
	}
	
	gzclose(fin);
}

// Return the reverse complement while preserving the case of the input sequence string
string complement(string m_seq /* make a copy */)
{
	for(string::iterator i = m_seq.begin();i != m_seq.end();i++){
		
		switch(*i){
			case 'A':
				*i = 'T';
				break;
			case 'a':
				*i = 't';
				break;
			case 'T':
				*i = 'A';
				break;
			case 't':
				*i = 'a';
				break;
			case 'G':
				*i = 'C';
				break;
			case 'g':
				*i = 'c';
				break;
			case 'C':
				*i = 'G';
				break;
			case 'c':
				*i = 'g';
				break;
			// Degenerate bases
			case 'M': // A or C -> T or G = K
				*i = 'K';
				break;
			case 'm':
				*i = 'k';
				break;
			case 'R': // G or A -> C or T = Y
				*i = 'Y';
				break;
			case 'r':
				*i = 'y';
				break;
			case 'S': // G or C -> C or G = S
				*i = 'S';
				break;
			case 's':
				*i = 's';
				break;
			case 'V': // G or C or A -> C or G or T = B
				*i = 'B';
				break;
			case 'v':
				*i = 'b';
				break;
			case 'W': // A or T -> T or A = W
				*i = 'W';
				break;
			case 'w':
				*i = 'w';
				break;
			case 'Y': // T or C -> A or G = R
				*i = 'R';
				break;
			case 'y':
				*i = 'r';
				break;
			case 'H': // A or C or T -> T or G or A = D
				*i = 'D';
				break;
			case 'h':
				*i = 'd';
				break;
			case 'K': // G or T -> C or A = M
				*i = 'M';
				break;
			case 'k':
				*i = 'm';
				break;
			case 'D': // G or A or T -> C or T or A = H
				*i = 'H';
				break;
			case 'd':
				*i = 'h';
				break;
			case 'B': // G or T or C -> C or A or G = V
				*i = 'V';
				break;
			case 'b':
				*i = 'v';
				break;
			case 'N': // G or T or C or A -> C or A or G or T = N
				*i = 'N';
				break;
			case 'n':
				*i = 'n';
				break;
		};
	}
	
	reverse( m_seq.begin(), m_seq.end() );
	
	return m_seq;
}
