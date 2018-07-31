// Validate read trimming results by comparing three file:
//	1) The original FASTQ file
//	2) The trimmed FASTQ file
//	3) The "ground truth" alignment of the original reads
//	   to the appropriate reference genome(s)
//
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Fri Apr 28 10:09:16 2017

#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include <zlib.h>
#include <string.h>
#include "fastq.h"

using namespace std;

#define	VTRIM_VERSION	"0.1"

void evaluate_trimming(const string &m_raw_seq, const string &m_trimmed_seq, 
	const size_t &m_align_start, const size_t &m_align_stop,
	size_t &m_over_5, size_t &m_over_3,
	size_t &m_under_5, size_t &m_under_3,
	size_t &m_correct_5, size_t &m_correct_3);

void evaluate_trimming(const string &m_raw_seq, const string &m_trimmed_seq, 
	size_t &m_correct_5, size_t &m_correct_3);

bool next_align(gzFile m_fin, string &m_def, size_t &m_start, size_t &m_stop);
size_t str_to_size_t(const string &m_str);
string truncate_defline(const string &m_str);

int main(int argc, char *argv[])
{
	try{
		
		string read_file;
		string trimmed_file;
		string alignment_file;
		bool print_usage = (argc == 1);
		
		const char* options = "r:t:a:?h";
		//int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		while( (opt_code = getopt_long_only( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:
					cerr << "Unknown flag!" << endl;
					break;
				case 'r':
					read_file = optarg;
					break;
				case 't':
					trimmed_file = optarg;
					break;
				case 'a':
					alignment_file = optarg;
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
			
			cerr << "VTrim version " << VTRIM_VERSION << endl;
			cerr << "Usage for " << argv[0] << endl;
			cerr << "\t-r <raw fastq read file>" << endl;
			cerr << "\t-t <trimmed fastq read file>" << endl;
			cerr << "\t-a <Tabular BLASTn (-m 8) alignment of raw reads to reference genome>" << endl;
			cerr << "\t[-h|-?] (Usage)" << endl;
			
			return EXIT_FAILURE;
		}
		
		if( read_file.empty() ){
		
			cerr << "Please specify a raw read file in fastq format (-r)" << endl;
			return EXIT_FAILURE;
		}
		
		if( trimmed_file.empty() ){
		
			cerr << "Please specify a trimmed read file in fastq format (-t)" << endl;
			return EXIT_FAILURE;
		}
		
		if( alignment_file.empty() ){
		
			cerr << "Please specify an alignment file in tabular BLASTn format (-a)" << endl;
			return EXIT_FAILURE;
		}
		
		gzFile fraw = gzopen( read_file.c_str(), "r");
		
		if(fraw == NULL){

			cerr << "Unable to open " << read_file << " for loading raw read sequences" << endl;
			throw "I/O error";
		}
		
		gzFile ftrim = gzopen( trimmed_file.c_str(), "r");
		
		if(ftrim == NULL){

			cerr << "Unable to open " << trimmed_file << " for loading trimmed read sequences" << endl;
			throw "I/O error";
		}
		
		gzFile falign = gzopen( alignment_file.c_str(), "r");
		
		if(falign == NULL){

			cerr << "Unable to open " << alignment_file << " for loading BLAST alignment" << endl;
			throw "I/O error";
		}
		
		bool update_raw = true;
		string raw_def;
		string raw_seq;
		string raw_qual;
		
		bool update_trimmed = true;
		string trimmed_def;
		string trimmed_seq;
		string trimmed_qual;
		
		bool update_alignment = true;
		string align_def;
		size_t align_start = 0;
		size_t align_stop = 0;
		
		size_t total_correct_5 = 0;
		size_t total_correct_3 = 0;
		size_t total_over_5 = 0;
		size_t total_under_5 = 0;
		size_t total_over_3 = 0;
		size_t total_under_3 = 0;
		
		// For reads that are completely trimmed away, but
		// still have valid alignments
		size_t total_complete_over = 0;
		
		// For reads that are *not* trimmed away and
		// do not have valid alignments
		size_t total_complete_under = 0;
		
		// For reads that are completely trimmed away, and
		// do not have a valid alignment
		size_t total_complete_correct = 0;
		
		size_t num_read = 0;
		size_t num_bp = 0;
		
		size_t num_read_A = 0; // Aligned and trimmed
		size_t num_read_B = 0; // Not aligned but trimmed
		size_t num_read_C = 0; // Aligned but trimmed away
		size_t num_read_D = 0; // Not aligned and trimmed away
		
		while(true){

			if(update_raw){
				
				if(next_read(fraw, raw_def, raw_seq, raw_qual) == false){
				
					// All done!
					break;
				}
				
				// Since the BLASTn alignment does not allow white space, truncate
				// the defline at the first space
				raw_def = truncate_defline(raw_def);
				
				++num_read;
				num_bp += raw_seq.size();
			}
			
			if(update_trimmed){
			
				next_read(ftrim, trimmed_def, trimmed_seq, trimmed_qual);
				
				// Since the BLASTn alignment does not allow white space, truncate
				// the defline at the first space
				trimmed_def = truncate_defline(trimmed_def);
			}

			if(update_alignment){
		
				// Throw away alignments with the same query defline (this are reads
				// that have multiple alignments to the reference genome). We assume that
				// the first alignment is the correct alignment.
				while(true){
					
					string local_def;
					size_t local_start = 0;
					size_t local_stop = 0;
				
					if(next_align(falign, local_def, local_start, local_stop) == false){
						break;
					}
					
					if(local_def != align_def){
						
						align_def = local_def;
						align_start = local_start;
						align_stop = local_stop;
						
						break;
					}
				}
			}
			
			update_raw = false;
			update_trimmed = false;
			update_alignment = false;
			
			//cerr << "raw: " << raw_def << endl;
			//cerr << "trimmed: " << trimmed_def << endl;
			//cerr << "align: " << align_def << endl << endl;
			
			// All deflines must agree in order to make a comparison
			if(raw_def == trimmed_def){
				
				if(raw_def == align_def){
					
					// DEBUG
					//cerr << "A" << endl;
					
					// This raw read appears in the trimmed output and
					// has a valid alignment
					size_t over_5 = 0;
					size_t over_3 = 0;
					size_t under_5 = 0;
					size_t under_3 = 0;
					size_t correct_5 = 0;
					size_t correct_3 = 0;

					// Evaluate the quality of this trimming					
					evaluate_trimming(raw_seq, trimmed_seq, 
						align_start, align_stop,
						over_5, over_3,
						under_5, under_3,
						correct_5, correct_3);
					
					total_over_5 += over_5;
					total_under_5 += under_5;
					
					total_over_3 += over_3;
					total_under_3 += under_3;
					
					total_correct_5 += correct_5;
					total_correct_3 += correct_3;
					
					update_raw = true;
					update_trimmed = true;
					update_alignment = true;
					
					++num_read_A;
				}
				else{ // raw_def != align_def
					
					// DEBUG
					//cerr << "\tB" << endl;
					
					// This does *not* have a valid alignment
					// but still has trimmed sequence (we've under trimmed)					
					size_t correct_5 = 0;
					size_t correct_3 = 0;
					
					evaluate_trimming(raw_seq, trimmed_seq, 
						correct_5, correct_3);
					
					total_correct_5 += correct_5;
					total_correct_3 += correct_3;
					
					total_complete_under += raw_seq.size() - (correct_5 + correct_3);
					
					update_raw = true;
					update_trimmed = true;
					
					++num_read_B;
				}
			}
			else{ // (raw_def != trimmed_def)
				
				// The trimmed file does *not* contain this raw read.
				if(raw_def == align_def){
					
					// DEBUG
					//cerr << "\t\tC" << endl;
					//cout << ">" << raw_def << ' ';
					//cout << "aligned -> [" << align_start 
					//	<< ", " << align_stop << "]" << endl;
					//	
					//cout << raw_seq << endl;
					
					//cout << raw_def << endl;
					//cout << raw_seq << endl;
					//cout << '+' << endl;
					//cout << raw_qual << endl;
					
					// This read has been completely trimmed away
					// (but *has* an alignment to the reference genome).
					total_complete_over += (align_stop - align_start) + 1;
					
					total_correct_5 += align_start - 1;
					total_correct_3 += raw_seq.size() - align_stop;
					
					update_raw = true;
					update_alignment = true;
					
					++num_read_C;
				}
				else{ // raw_def != align_def
					
					// DEBUG
					//cerr << "\t\t\tD" << endl;
					
					// This read has been completely trimmed away
					// and does *not* have a valid alignment -- good job
					// read trimmer!
					
					total_complete_correct += raw_seq.size();
					
					update_raw = true;
					
					++num_read_D;
				}
			}
		}
		
		gzclose(fraw);
		gzclose(ftrim);
		gzclose(falign);
		
		// DEBUG
		//cerr << "total_correct_5 = " << total_correct_5 << endl;
		//cerr << "total_correct_3 = " << total_correct_3 << endl;
		//cerr << "total_over_5 = " << total_over_5 << endl;
		//cerr << "total_under_5 = " << total_under_5 << endl;
		//cerr << "total_over_3 = " << total_over_3 << endl;
		//cerr << "total_under_3 = " << total_under_3 << endl;
		//cerr << "total_complete_over = " << total_complete_over << endl;
		//cerr << "total_complete_under = " << total_complete_under << endl;
		//cerr << "total_complete_correct = " << total_complete_correct << endl;
		
		
		// Write the performance statistics
		size_t tp = total_correct_5 + total_correct_3 + total_complete_correct;
		size_t fp = total_over_5 + total_over_3 + total_complete_over;
		size_t fn = total_under_5 + total_under_3 + total_complete_under;
		
		cout << "Raw read file: " << read_file << endl;
		cout << "Trimmed read file: " << trimmed_file << endl;
		cout << "Tabular BLASTn alignment file: " << alignment_file << endl;
		
		cout << "Processed " << num_bp << " bp in " << num_read << " reads" << endl;
		
		cout << "\tNumber of reads aligned and trimmed = " << num_read_A << " (" 
			<< (100.0*num_read_A)/num_read << "%)" << endl;
		cout << "\tNumber of reads not aligned but trimmed = " << num_read_B << " (" 
			<< (100.0*num_read_B)/num_read << "%)" << endl;
		cout << "\tNumber of reads aligned but trimmed away = " << num_read_C << " (" 
			<< (100.0*num_read_C)/num_read << "%)" << endl;
		cout << "\tNumber of reads not aligned and trimmed away = " << num_read_D << " (" 
			<< (100.0*num_read_D)/num_read << "%)" << endl;
		
		cout << "Whole read performance:" << endl;
		cout << "\tTp = " << tp << endl;
		cout << "\tFp = " << fp << endl;
		cout << "\tFn = " << fn << endl;
		cout << "\tPrecision [tp/(tp + fp)] = " << double(tp)/(tp + fp) << endl;
		cout << "\tRecall [tp/(tp + fn)] = " << double(tp)/(tp + fn) << endl;
		
		tp = total_correct_5;
		fp = total_over_5;
		fn = total_under_5;
		
		cout << "5' performance:" << endl;
		cout << "\tTp = " << tp << endl;
		cout << "\tFp = " << fp << endl;
		cout << "\tFn = " << fn << endl;
		cout << "\tPrecision [tp/(tp + fp)] = " << double(tp)/(tp + fp) << endl;
		cout << "\tRecall [tp/(tp + fn)] = " << double(tp)/(tp + fn) << endl;
		
		tp = total_correct_3;
		fp = total_over_3;
		fn = total_under_3;
		
		cout << "3' performance:" << endl;
		cout << "\tTp = " << tp << endl;
		cout << "\tFp = " << fp << endl;
		cout << "\tFn = " << fn << endl;
		cout << "\tPrecision [tp/(tp + fp)] = " << double(tp)/(tp + fp) << endl;
		cout << "\tRecall [tp/(tp + fn)] = " << double(tp)/(tp + fn) << endl;
		
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}

bool next_align(gzFile m_fin, string &m_def, size_t &m_start, size_t &m_stop)
{
	if( gzeof(m_fin) ){
		return false;
	}

	const int buffer_len = 4096;
	char buffer[buffer_len];
	
	if(gzgets(m_fin, buffer, buffer_len) == NULL){

		if( gzeof(m_fin) ){
			return false;
		}

		throw __FILE__ ":next_align: Unable to read buffer";
	}
	
	if(strpbrk(buffer, "\n\r") == NULL){
		throw __FILE__ ":next_align: Line buffer overflow";
	}
		
	// Split the BLASTn tabular alignment line into columns based on
	// tab delimiters
	size_t col = 0;
	string record;
	const char delim = '\t';
	
	const size_t QUERY_COL = 0;
	const size_t QUERY_START = 6;
	const size_t QUERY_STOP = 7;
	
	for(char* ptr = buffer;(*ptr != '\n') && (*ptr != '\r');++ptr){
		
		if(*ptr == delim){
			
			switch(col){
				case QUERY_COL:
				
					m_def = record;
					break;
				case QUERY_START:
				
					m_start = str_to_size_t(record);
					break;
				case QUERY_STOP:
				
					m_stop = str_to_size_t(record);
					
					// Stop parsing once we have the query stop
					return true;
			};
			
			record.clear();
			++col;
		}
		else{
			record.push_back(*ptr);
		}
	}
	
	return false;
}

size_t str_to_size_t(const string &m_str)
{
	size_t ret = 0;
	size_t power = 1;
	
	for(string::const_reverse_iterator i = m_str.rbegin();i != m_str.rend();++i){
		
		if( !isdigit(*i) ){
			throw __FILE__ ":str_to_size_t: Illegal digit";
		}
		
		ret += power*(*i - '0');
		
		power *= 10;
	}
	
	return ret;
}

void evaluate_trimming(const string &m_raw_seq, const string &m_trimmed_seq, 
	const size_t &m_align_start, const size_t &m_align_stop,
	size_t &m_over_5, size_t &m_over_3,
	size_t &m_under_5, size_t &m_under_3,
	size_t &m_correct_5, size_t &m_correct_3)
{
	const size_t len = m_raw_seq.size();
	
	// Align the trimmed read to the raw read
	size_t start = m_raw_seq.find(m_trimmed_seq);
	
	if(start == string::npos){
		throw __FILE__ ":evaluate_trimming: Unable to align trimed to raw";
	}
	
	const size_t stop = start + m_trimmed_seq.size();
	
	// BLASTn uses 1's based indicies
	++start;
	
	if(start > m_align_start){
		
		// Over trimmed
		m_over_5 += start - m_align_start;
		m_correct_5 += m_align_start - 1;
	}
	else{
		
		if(start < m_align_start){
		
			// Under trimmed
			m_under_5 += m_align_start - start;
			m_correct_5 += start - 1;
		}
		else{
			m_correct_5 += start - 1;
		}
	}
	
	if(stop < m_align_stop){
		
		// Over trimmed
		m_over_3 += m_align_stop - stop;
		m_correct_3 += len - m_align_stop;
	}
	else{
		
		if(stop > m_align_stop){
		
			// Under trimmed
			m_under_3 += stop - m_align_stop;
			m_correct_3 += len - stop;
		}
		else{
			m_correct_3 += len - stop;
		}
	}
}

// For reads that *should* have been completely trimmed away (i.e. no alignment)
// but were only partially trimmed
void evaluate_trimming(const string &m_raw_seq, const string &m_trimmed_seq, 
	size_t &m_correct_5, size_t &m_correct_3)
{
	const size_t len = m_raw_seq.size();
	
	// Align the trimmed read to the raw read
	size_t start = m_raw_seq.find(m_trimmed_seq);
	
	if(start == string::npos){
		throw __FILE__ ":evaluate_trimming: Unable to align trimed to raw";
	}
	
	const size_t stop = start + m_trimmed_seq.size();
	
	// BLASTn uses 1's based indicies
	++start;
	
	m_correct_5 += start - 1;
	m_correct_3 += len - stop;
}

string truncate_defline(const string &m_str)
{
	const size_t pos = m_str.find(' ');
	
	if(pos == string::npos){
		return m_str;
	}
	
	return m_str.substr(0, pos);
}
