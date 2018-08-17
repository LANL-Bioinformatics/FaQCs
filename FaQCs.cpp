// A C++/OpenMP version of Chien-Chi Lo's FaQCs.pl perl scripts
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Wed Mar  1 10:19:26 2017

#include <iostream>
#include <iomanip>
#include <deque>
#include <algorithm>

#include <stdlib.h>

#ifdef _OPENMP
	#include <omp.h>
#endif // _OPENMP

#include "FaQCs.h"
#include "fastq.h"
#include "file_util.h"

using namespace std;

void process_paired(vector<size_t> &m_filter_stats,
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	PlotInfo &m_info, Options &m_opt);
void process_unpaired(vector<size_t> &m_filter_stats,
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	PlotInfo &m_info, Options &m_opt);
void write_stats(const vector<size_t> &m_filter_stats, 
	const MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	const Options &m_opt);
void remove_existing_files(const Options &m_opt);
void remove_file(const string &m_filename);

int main(int argc, char *argv[])
{
	try{
		
		// Parse the command line options
		Options opt(argc, argv);
		
		if(opt.print_usage){
			return EXIT_FAILURE;
		}
		
		#ifdef _OPENMP
		if(opt.num_thread > 0){
			
			// Limit the maximum number of threads that OpenMP is
			// allowed to use
			omp_set_num_threads(opt.num_thread);
		}
		#endif // _OPENMP
		
		// If the output directory does not exist, create it
		if( directory_exits(opt.output_dir) ){
		}
		else{
			if( !create_directory(opt.output_dir) ){
				
				cerr << "Unable to create requested output directory: \"" << opt.output_dir << '"' << endl;
				return EXIT_FAILURE;
			}
		}
		
		vector<size_t> filter_stats(FilterStat::NUM_STAT);
		MAP< string, pair<size_t /*read*/, size_t /*base*/> > adaptor_stats;
		PlotInfo info;
		
		// Follow the example of the FaQCs.pl script, which removes any existing
		// output files before writing new ones.
		remove_existing_files(opt);
		
		if( opt.has_paired() ){
			
			// Process paired reads
			process_paired(filter_stats, adaptor_stats, 
				info, opt);
		}		
		
		if( opt.has_unpaired() ){
			
			// Process unpaired reads
			process_unpaired(filter_stats, adaptor_stats, 
				info, opt);
		}
		
		// The phiX sequences are treated as adapters. For reporting the trimming statistics,
		// we need to remove their read and base counts from the adaptor_stats structure and
		// put those values into the filter_stats vector.
		if(opt.filter_phiX){
			
			// The phiX sequence is one adaptor ...
			MAP< string, pair<size_t /*read*/, size_t /*base*/> >::const_iterator iter = 
				adaptor_stats.find(PHI_X);
			
			if( iter != adaptor_stats.end() ){
				
				filter_stats[FilterStat::READ_PHIX] += iter->second.first;
				filter_stats[FilterStat::BASE_PHIX] += iter->second.second;
				
				adaptor_stats.erase(iter);
			}
			
			// ... and the complement of the phiX sequence is another adaptor.
			iter = adaptor_stats.find(PHI_X_COMPLEMENT);
			
			if( iter != adaptor_stats.end() ){
				
				filter_stats[FilterStat::READ_PHIX] += iter->second.first;
				filter_stats[FilterStat::BASE_PHIX] += iter->second.second;
				
				adaptor_stats.erase(iter);
			}
		}
		
		if(opt.filter_adapter){
			
			// Collect the adapter stats
			for(MAP< string, pair<size_t, size_t> >::const_iterator i = adaptor_stats.begin();
				i != adaptor_stats.end();++i){
				
				filter_stats[FilterStat::READ_ADAPTER] += i->second.first;
				filter_stats[FilterStat::BASE_ADAPTER] += i->second.second;
			}
		}
		
		write_stats(filter_stats, adaptor_stats, opt);
		
		if(!opt.trim_only){
			plot(info, filter_stats, opt);
		}
		
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

void process_paired(vector<size_t> &m_filter_stats, 
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	PlotInfo &m_info, Options &m_opt)
{

	// Open the input files for reading
	gzFile fin1 = gzopen( m_opt.input_read1_file.c_str(), "r");
		
	if(fin1 == NULL){

		cerr << "Unable to open " << m_opt.input_read1_file << " for loading read one sequences" << endl;
		throw "I/O error";
	}
	
	gzFile fin2 = gzopen( m_opt.input_read2_file.c_str(), "r");
		
	if(fin2 == NULL){

		cerr << "Unable to open " << m_opt.input_read2_file << " for loading read two sequences" << endl;
		
		gzclose(fin1);
		
		throw "I/O error";
	}
	
	// Open the output files for writing (as uncompressed files for now)
	gzFile fout1 = NULL;
	gzFile fout2 = NULL;
	gzFile fout_unpaired = NULL;
	gzFile fout_discard = NULL;
	
	if(!m_opt.qc_only){
	
		fout1 = gzopen( m_opt.trimmed_read1_file.c_str(), "wT");
	
		if(fout1 == NULL){

			cerr << "Unable to open " << m_opt.trimmed_read1_file 
				<< " for writing read one sequences" << endl;
			throw "I/O error";
		}
		
		fout2 = gzopen( m_opt.trimmed_read2_file.c_str(), "wT");
	
		if(fout2 == NULL){

			cerr << "Unable to open " << m_opt.trimmed_read2_file 
				<< " for writing read two sequences" << endl;
			throw "I/O error";
		}
		
		fout_unpaired = gzopen( m_opt.trimmed_unpaired_file.c_str(), "wT");
	
		if(fout_unpaired == NULL){

			cerr << "Unable to open " << m_opt.trimmed_unpaired_file 
				<< " for writing unpaired sequences" << endl;
			throw "I/O error";
		}
		
		if( !m_opt.trimmed_discard_file.empty() ){
		
			fout_discard = gzopen( m_opt.trimmed_discard_file.c_str(), "wT");
	
			if(fout_discard == NULL){

				cerr << "Unable to open " << m_opt.trimmed_discard_file 
					<< " for writing discarded sequences" << endl;
				throw "I/O error";
			}
		}
	}
		
	vector<Read> buffer1;
	vector<Read> buffer2;
	
	vector<Read> discard_buffer1;
	vector<Read> discard_buffer2;
	
	const size_t buffer_size = 32768;
	bool check_for_next_seq = true;
	
	MAP<Word, size_t> kmer_table;

	buffer1.reserve(buffer_size);
	buffer2.reserve(buffer_size);
	
	while(true){
		
		// Avoid additional overhead of copying read data by
		// storing the data directly into the buffer.	
		buffer1.push_back( Read() );
		buffer2.push_back( Read() );
	
		Read &r1 = buffer1.back();
		Read &r2 = buffer2.back();
		
		bool ret1 = next_read(fin1, r1.def, r1.seq, r1.qual);
		bool ret2 = next_read(fin2, r2.def, r2.seq, r2.qual);
	
		if(!ret1 && !ret2){

			// We have reached the end of the sequence data, remove the
			// the last buffer element that we added in anticipation of 
			// additional data to store.
			buffer1.pop_back();
			buffer2.pop_back();
 	
			if(m_opt.input_quality_offset == AUTO_DETECT_QUALITY_OFFSET){
				
				m_opt.input_quality_offset = auto_detect_quality_offset(buffer1);
				
				if( m_opt.input_quality_offset != auto_detect_quality_offset(buffer2) ){
					
					cerr << "Inconsistent quality offset detection between reads one and two" << endl;
					throw __FILE__ ":process_paired: I/O Error";
				}
			}
			
			if( (m_opt.quality < DEFAULT_NEXTSEQ_QUALITY_SCORE) && auto_detect_next_seq(buffer1) ){
					
				cerr << "The input looks like NextSeq data and the quality level (-q) is adjusted to "
					<< DEFAULT_NEXTSEQ_QUALITY_SCORE << " for trimming." << endl;
				m_opt.quality = DEFAULT_NEXTSEQ_QUALITY_SCORE;
			}
			
			if(fout_discard != NULL){
				
				// Make a copy of all input reads so we can write
				// discarded reads to disk
				discard_buffer1 = buffer1;
				discard_buffer2 = buffer2;
			}
			
			trim(buffer1, m_filter_stats, m_adapter_stats, 
				kmer_table, m_info, m_opt);
			
			trim(buffer2, m_filter_stats, m_adapter_stats, 
				kmer_table, m_info, m_opt);
			
			const size_t len = buffer1.size();
			
			// Write the trimmed reads to disk
			for(size_t i = 0;i < len;++i){
				
				const Read& curr1 = buffer1[i];
				const Read& curr2 = buffer2[i];
				
				const bool valid1 = curr1.valid();
				const bool valid2 = curr2.valid();
				
				if(valid1 && valid2){
			
					m_filter_stats[FilterStat::PAIRED_READ_NUMBER] += 2;
					m_filter_stats[FilterStat::PAIRED_BASE_LENGTH] += curr1.seq.size() + curr2.seq.size();
				}
				
				if(!m_opt.qc_only){
				
					// For reads that have not been trimmed away ...
					if(valid1 && valid2){

						write_read(fout1, 
							curr1.def, 
							curr1.seq, 
							curr1.qual);

						write_read(fout2, 
							curr2.def, 
							curr2.seq, 
							curr2.qual);
					}
					else{
						if(valid1 && !valid2){
							write_read(fout_unpaired, 
								curr1.def, 
								curr1.seq, 
								curr1.qual);
						}
						else{
							if(valid2 && !valid1){
								write_read(fout_unpaired, 
									curr2.def, 
									curr2.seq, 
									curr2.qual);
							}
						}
						
						if(fout_discard != NULL){
							
							if(!valid1){
							
								write_read(fout_discard, 
									discard_buffer1[i].def,
									discard_buffer1[i].seq,
									discard_buffer1[i].qual);
							}
							
							if(!valid2){
							
								write_read(fout_discard, 
									discard_buffer2[i].def,
									discard_buffer2[i].seq,
									discard_buffer2[i].qual);
							}
						}
					}
				}
			}

			buffer1.clear();
			buffer2.clear();

			// We're done reading the sequence data
			break;
		}

		if(ret1 != ret2){

			if(ret1){
				cerr << "Did not find a match to read one: " << r1.def << endl;
			}
			else{
				cerr << "Did not find a match to read two: " << r2.def << endl;
			}

			throw __FILE__ "I/O error";
		}

		// Make sure that the read ids match
		if(parse_id(r1.def) != parse_id(r2.def)){

			cerr << "Read one id (" << parse_id(r1.def) << ")\n"
				<< "does not match\n"
				<< "read two id (" << parse_id(r2.def) << ")" << endl;
			throw __FILE__ ":trim: I/O error";
		}

		if(buffer1.size() == buffer_size){

			if(m_opt.input_quality_offset == AUTO_DETECT_QUALITY_OFFSET){
			
				m_opt.input_quality_offset = auto_detect_quality_offset(buffer1);
				
				if( m_opt.input_quality_offset != auto_detect_quality_offset(buffer2) ){
					
					cerr << "Inconsistent quality offset detection between reads one and two" << endl;
					throw __FILE__ ":process_paired: I/O Error";
				}
			}
			
			if(check_for_next_seq){
				
				if( (m_opt.quality < DEFAULT_NEXTSEQ_QUALITY_SCORE) && auto_detect_next_seq(buffer1) ){
					
					cerr << "The input looks like NextSeq data and the quality level (-q) is adjusted to "
						<< DEFAULT_NEXTSEQ_QUALITY_SCORE << " for trimming." << endl;
					m_opt.quality = DEFAULT_NEXTSEQ_QUALITY_SCORE;
				}
				
				check_for_next_seq = false;
			}
			
			if(fout_discard != NULL){
				
				// Make a copy of all input reads so we can write
				// discarded reads to disk
				discard_buffer1 = buffer1;
				discard_buffer2 = buffer2;
			}
			
			trim(buffer1, m_filter_stats, m_adapter_stats, 
				kmer_table, m_info, m_opt);
			
			trim(buffer2, m_filter_stats, m_adapter_stats, 
				kmer_table, m_info, m_opt);
			
			// Write the trimmed reads to disk
			for(size_t i = 0;i < buffer_size;++i){
				
				const Read& curr1 = buffer1[i];
				const Read& curr2 = buffer2[i];
				
				const bool valid1 = curr1.valid();
				const bool valid2 = curr2.valid();
				
				if(valid1 && valid2){
			
					m_filter_stats[FilterStat::PAIRED_READ_NUMBER] += 2;
					m_filter_stats[FilterStat::PAIRED_BASE_LENGTH] += curr1.seq.size() + curr2.seq.size();
				}

				if(!m_opt.qc_only){

					// For reads that have not been trimmed away ...
					if(valid1 && valid2){

						write_read(fout1, 
							curr1.def, 
							curr1.seq, 
							curr1.qual);

						write_read(fout2, 
							curr2.def, 
							curr2.seq, 
							curr2.qual);
					}
					else{
						if(valid1 && !valid2){
							write_read(fout_unpaired, 
								curr1.def, 
								curr1.seq, 
								curr1.qual);
						}
						else{
							if(valid2 && !valid1){
								write_read(fout_unpaired, 
									curr2.def, 
									curr2.seq, 
									curr2.qual);
							}
						}
						
						if(fout_discard != NULL){
							
							if(!valid1){
							
								write_read(fout_discard, 
									discard_buffer1[i].def,
									discard_buffer1[i].seq,
									discard_buffer1[i].qual);
							}
							
							if(!valid2){
							
								write_read(fout_discard, 
									discard_buffer2[i].def,
									discard_buffer2[i].seq,
									discard_buffer2[i].qual);
							}
						}
					}
				}
			}

			buffer1.clear();
			buffer2.clear();
		}
	}
	
	gzclose(fin1);
	gzclose(fin2);
		
	if(!m_opt.qc_only){

		gzclose(fout1);
		gzclose(fout2);

		gzclose(fout_unpaired);
		
		if(fout_discard != NULL){
			gzclose(fout_discard);
		}
	}
	
	// Collect the kmer frequency table
	for(MAP<Word, size_t>::const_iterator i = kmer_table.begin();i != kmer_table.end();++i){
		++m_info.kmer_frequency_histogram[i->second];
	}

	if( m_opt.kmer_rarefaction && m_info.kmer_rarefaction.empty() ){

		// Make sure that the rarefaction curve will have at least one point
		Rarefaction local;

		local.num_seq = m_filter_stats[FilterStat::TOTAL_NUMBER];
		local.distinct_kmer = kmer_table.size();
		local.total_kmer = 0;

		for(MAP<Word, size_t>::const_iterator i = kmer_table.begin();i != kmer_table.end();++i){
			local.total_kmer += i->second;
		}

		m_info.kmer_rarefaction.push_back(local);
        }
}

void process_unpaired(vector<size_t> &m_filter_stats, 
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	PlotInfo &m_info, Options &m_opt)
{
	// Open the input files for reading
	gzFile fin = gzopen( m_opt.input_unpaired_file.c_str(), "r");
		
	if(fin == NULL){

		cerr << "Unable to open " << m_opt.input_unpaired_file 
			<< " for loading unpaired read sequences" << endl;
		throw "I/O error";
	}
	
	// Open the output files for writing (as uncompressed files for now)
	gzFile fout = NULL;
	gzFile fout_discard = NULL;
	
	if(!m_opt.qc_only){
		
		fout = gzopen( m_opt.trimmed_unpaired_file.c_str(), "wT");
	
		if(fout == NULL){

			cerr << "Unable to open " << m_opt.trimmed_unpaired_file 
				<< " for writing unpaired read sequences" << endl;
			throw "I/O error";
		}
		
		if( !m_opt.trimmed_discard_file.empty() ){
		
			fout_discard = gzopen( m_opt.trimmed_discard_file.c_str(), "wT");
	
			if(fout_discard == NULL){

				cerr << "Unable to open " << m_opt.trimmed_discard_file 
					<< " for writing discarded sequences" << endl;
				throw "I/O error";
			}
		}
	}
	
	vector<Read> buffer;
	vector<Read> discard_buffer;
	
	const size_t buffer_size = 32768;
	bool check_for_next_seq = true;
	
	MAP<Word, size_t> kmer_table;

	buffer.reserve(buffer_size);
	
	while(true){
	
		// Avoid additional overhead of copying read data by
		// storing the data directly into the buffer.	
		buffer.push_back( Read() );
	
		Read &r = buffer.back();
		
		bool ret = next_read(fin, r.def, r.seq, r.qual);
	
		if(!ret){
		
			// We have reached the end of the sequence data, remove the
			// the last buffer element that we added in anticipation of 
			// additional data to store.
			buffer.pop_back();
 	
			if(m_opt.input_quality_offset == AUTO_DETECT_QUALITY_OFFSET){				
				m_opt.input_quality_offset = auto_detect_quality_offset(buffer);
			}
			
			if( (m_opt.quality < DEFAULT_NEXTSEQ_QUALITY_SCORE) && auto_detect_next_seq(buffer) ){
					
				cerr << "The input looks like NextSeq data and the quality level (-q) is adjusted to " 
					<< DEFAULT_NEXTSEQ_QUALITY_SCORE << " for trimming." << endl;
					
				m_opt.quality = DEFAULT_NEXTSEQ_QUALITY_SCORE;
			}
			
			if(fout_discard != NULL){
				
				// Make a copy of all input reads so we can write
				// discarded reads to disk
				discard_buffer = buffer;
			}
			
			trim(buffer, m_filter_stats, m_adapter_stats, 
				kmer_table, m_info, m_opt);
			
			const size_t len = buffer.size();
			
			// Write the trimmed reads to disk
			for(size_t i = 0;i < len;++i){
				
				const Read& curr = buffer[i];
				
				const bool valid = curr.valid();
				
				// For reads that have not been trimmed away ...
				if(!m_opt.qc_only){
					
					if(valid){

						write_read(fout,
							curr.def,
							curr.seq,
							curr.qual);
					}
					else{
						if(fout_discard != NULL){
							write_read(fout_discard,
								discard_buffer[i].def,
								discard_buffer[i].seq,
								discard_buffer[i].qual);
						}
					}
				}
			}

			buffer.clear();

			// We're done reading the sequence data
			break;
		}

		if(buffer.size() == buffer_size){

			if(m_opt.input_quality_offset == AUTO_DETECT_QUALITY_OFFSET){
				m_opt.input_quality_offset = auto_detect_quality_offset(buffer);				
			}
			
			if(check_for_next_seq){
				
				if( (m_opt.quality < DEFAULT_NEXTSEQ_QUALITY_SCORE) && auto_detect_next_seq(buffer) ){
					
					cerr << "The input looks like NextSeq data and the quality level (-q) is adjusted to "
						<< DEFAULT_NEXTSEQ_QUALITY_SCORE << " for trimming." << endl;
					m_opt.quality = DEFAULT_NEXTSEQ_QUALITY_SCORE;
				}
				
				check_for_next_seq = false;
			}
			
			if(fout_discard != NULL){
				
				// Make a copy of all input reads so we can write
				// discarded reads to disk
				discard_buffer = buffer;
			}
			
			trim(buffer, m_filter_stats, m_adapter_stats, 
				kmer_table, m_info, m_opt);
			
			// Write the trimmed reads to disk
			for(size_t i = 0;i < buffer_size;++i){
				
				const Read& curr = buffer[i];
				const bool valid = curr.valid();
				
				// For reads that have not been trimmed away ...
				if(!m_opt.qc_only){
					
					if(valid){

						write_read(fout,
							curr.def,
							curr.seq,
							curr.qual);
					}
					else{
						if(fout_discard != NULL){
							write_read(fout_discard,
								discard_buffer[i].def,
								discard_buffer[i].seq,
								discard_buffer[i].qual);
						}
					}
				}
			}

			buffer.clear();
		}
	}
	
	gzclose(fin);
	
	if(!m_opt.qc_only){
	
		gzclose(fout);
		
		if(fout_discard != NULL){
			gzclose(fout_discard);
		}
	}

	// Collect the kmer frequency table
	for(MAP<Word, size_t>::const_iterator i = kmer_table.begin();i != kmer_table.end();++i){
		++m_info.kmer_frequency_histogram[i->second];
	}

	if( m_opt.kmer_rarefaction && m_info.kmer_rarefaction.empty() ){

		// Make sure that the rarefaction curve will have at least one point
		Rarefaction local;

		local.num_seq = m_filter_stats[FilterStat::TOTAL_NUMBER];
		local.distinct_kmer = kmer_table.size();
		local.total_kmer = 0;

		for(MAP<Word, size_t>::const_iterator i = kmer_table.begin();i != kmer_table.end();++i){
			local.total_kmer += i->second;
		}

		m_info.kmer_rarefaction.push_back(local);
        }
}

void write_stats(const vector<size_t> &m_filter_stats, 
	const MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	const Options &m_opt)
{
	ofstream fout( m_opt.stats_file.c_str() );
	
	if(!fout){
	
		// Do throw an error if we are unable to write the statistics file, just
		// print an error an skip this step.
		cerr << "Unable to open " << m_opt.stats_file 
			<< " for writing filtering statistics" << endl;
		return;
	}
	
	using namespace FilterStat;
	
	fout << fixed << setprecision(2); // Used fixed precision up to hundreths (unless noted below).
	
	if(m_opt.qc_only){
	
		fout << "\n";
		fout << "Reads #: " << m_filter_stats[TOTAL_COUNT] << "\n";
      		fout << "Total bases: " << m_filter_stats[TOTAL_LENGTH] << "\n";
		fout << "Reads Length: " 
			<< float(m_filter_stats[TOTAL_LENGTH])/m_filter_stats[TOTAL_COUNT] 
			<< "\n";
		fout << "Processed " << m_filter_stats[TOTAL_NUMBER] << " reads for quality check only\n";
      
		fout << "  Reads length < " 
			<< m_opt.min_read_length 
			<< " bp: " << m_filter_stats[READ_LENGTH] 
			<< " (" << (100.0*m_filter_stats[READ_LENGTH])/m_filter_stats[TOTAL_NUMBER]
			<< " %)\n";
		fout << "  Reads have " << m_opt.max_num_poly_N 
			<< " continuous base \"N\": " << m_filter_stats[READ_NN] 
			<< " (" << (100.0*m_filter_stats[READ_NN])/m_filter_stats[TOTAL_NUMBER] << " %)\n";
		fout << "  Low complexity Reads  (>" << m_opt.low_complexity_cutoff_ratio*100.0
			<< "% mono/di-nucleotides): " << m_filter_stats[READ_LOW_COMPLEXITY] 
			<< " (" << (100.0*m_filter_stats[READ_LOW_COMPLEXITY])/m_filter_stats[TOTAL_NUMBER] << " %)\n";
		fout << "  Reads < average quality " << m_opt.average_quality
			<< ": " << m_filter_stats[READ_AVG_Q] << " (" 
			<< (100.0*m_filter_stats[READ_AVG_Q])/m_filter_stats[TOTAL_NUMBER] << " %)\n";
			
		if(m_opt.filter_phiX){
			fout << "  Reads hits to phiX sequence: " 
				<< m_filter_stats[READ_PHIX] << " (" 
				<< (100.0*m_filter_stats[READ_PHIX])/m_filter_stats[TOTAL_NUMBER] << " %)\n";
		}

		if(m_opt.filter_adapter){
		
			fout << "  Reads with Adapters/Primers: " << m_filter_stats[READ_ADAPTER] 
				<< " (" << (100.0*m_filter_stats[READ_ADAPTER])/m_filter_stats[TOTAL_NUMBER] << " %)\n";
			
			
			// Write the adapters in order of the number of reads matches
			deque< pair<size_t, string> > adapter;
			
			for(MAP< string, pair<size_t, size_t> >::const_iterator i = m_adapter_stats.begin();
				i != m_adapter_stats.end();++i){
				
				adapter.push_back( make_pair(i->second.first, i->first) ); // {reads, adapter name}
			}
			
			// Sort in *ascending* order ...
			sort( adapter.begin(), adapter.end() );
			
			// ... and iterate from the back to print in *descending* order
			for(deque< pair<size_t, string> >::const_reverse_iterator i = adapter.rbegin();
				i != adapter.rend();++i){
				
				MAP< string, pair<size_t, size_t> >::const_iterator iter = 
					m_adapter_stats.find(i->second);
					
				if( iter == m_adapter_stats.end() ){
				
					cerr << "Unable to find adapter statistics!" << endl;
					continue;
				}
				
				const size_t &affected_reads = iter->second.first;
				const size_t &affected_bases = iter->second.second;
				
				fout << "    " << iter->first << " " << affected_reads 
					<< " reads (" << (100.0*affected_reads)/m_filter_stats[TOTAL_NUMBER]
					<< " %) " << affected_bases << " bases ("
					<< (100.0*affected_bases)/m_filter_stats[TOTAL_LENGTH]
					<< " %)\n";
			}
		}
	}
	else{
	
		fout << "Before Trimming\n";
		fout << "Reads #: " << m_filter_stats[TOTAL_NUMBER] << "\n";
		fout << "Total bases: " << m_filter_stats[TOTAL_LENGTH]<< "\n";
		fout << "Reads Length: " 
			<< float(m_filter_stats[TOTAL_LENGTH])/m_filter_stats[TOTAL_NUMBER] 
			<< "\n";

		fout << "\nAfter Trimming\n";
		fout << "Reads #: " << m_filter_stats[TOTAL_TRIMMED_NUMBER] 
			<< " (" << (100.0*m_filter_stats[TOTAL_TRIMMED_NUMBER])/m_filter_stats[TOTAL_NUMBER] 
			<< " %)\n";
		fout << "Total bases: " << m_filter_stats[TOTAL_TRIMMED_LENGTH] 
			<< " (" << (100.0*m_filter_stats[TOTAL_TRIMMED_LENGTH])/m_filter_stats[TOTAL_LENGTH] << " %)\n";

		if(m_filter_stats[TOTAL_TRIMMED_NUMBER] > 0){
		
			fout << "Mean Reads Length: " 
				<< float(m_filter_stats[TOTAL_TRIMMED_LENGTH])/m_filter_stats[TOTAL_TRIMMED_NUMBER]
				<< "\n";
		}
		else{
			fout << "Mean Reads Length: 0\n";
		}
  
		if( m_opt.has_paired() ){
			
			fout << "  Paired Reads #: " << m_filter_stats[PAIRED_READ_NUMBER]
				<< " (" << (100.0*m_filter_stats[PAIRED_READ_NUMBER])/m_filter_stats[TOTAL_TRIMMED_NUMBER] 
				<< " %)\n";
			fout << "  Paired total bases: " <<  m_filter_stats[PAIRED_BASE_LENGTH] 
				<< " (" << (100.0*m_filter_stats[PAIRED_BASE_LENGTH])/m_filter_stats[TOTAL_TRIMMED_LENGTH]  << " %)\n";
			fout << "  Unpaired Reads #: " 
				<< m_filter_stats[TOTAL_TRIMMED_NUMBER] - m_filter_stats[PAIRED_READ_NUMBER] 
				<< " (" 
				<< ( 100.0*(m_filter_stats[TOTAL_TRIMMED_NUMBER]- m_filter_stats[PAIRED_READ_NUMBER]) )/m_filter_stats[TOTAL_TRIMMED_NUMBER] 
				<< " %)\n"; 
			fout << "  Unpaired total bases: " 
				<< m_filter_stats[TOTAL_TRIMMED_LENGTH] - m_filter_stats[PAIRED_BASE_LENGTH]
				<< " (" 
				<< ( 100.0*(m_filter_stats[TOTAL_TRIMMED_LENGTH] - m_filter_stats[PAIRED_BASE_LENGTH]) )/m_filter_stats[TOTAL_TRIMMED_LENGTH] 
				<< " %)\n";
		}
    
		fout << "\nDiscarded reads #: " << m_filter_stats[TOTAL_NUMBER] - m_filter_stats[TOTAL_TRIMMED_NUMBER]
			<< " (" 
			<< ( 100.0*(m_filter_stats[TOTAL_NUMBER] - m_filter_stats[TOTAL_TRIMMED_NUMBER]) )/m_filter_stats[TOTAL_NUMBER]
			<< " %)\n";
		fout << "Trimmed bases: " << m_filter_stats[TOTAL_LENGTH] - m_filter_stats[TOTAL_TRIMMED_LENGTH] 
			<< " ("
			<< ( 100.0*(m_filter_stats[TOTAL_LENGTH] - m_filter_stats[TOTAL_TRIMMED_LENGTH]) )/m_filter_stats[TOTAL_LENGTH]
			<< " %)\n";
		fout << "  Reads Filtered by length cutoff (" 
			<< m_opt.min_read_length << " bp): " << m_filter_stats[READ_LENGTH] 
			<< " (" << (100.0*m_filter_stats[READ_LENGTH])/m_filter_stats[TOTAL_NUMBER] 
			<< " %)\n";
      		fout << "  Bases Filtered by length cutoff: " << m_filter_stats[BASE_LENGTH]
			<< " (" 
			<< (100.0*m_filter_stats[BASE_LENGTH])/m_filter_stats[TOTAL_LENGTH] 
			<< " %)\n";
		fout << "  Reads Filtered by continuous base \"N\" (" << m_opt.max_num_poly_N 
			<< "): " << m_filter_stats[READ_NN] << " (" 
			<< (100.0*m_filter_stats[READ_NN])/m_filter_stats[TOTAL_NUMBER]  << " %)\n";
		fout << "  Bases Filtered by continuous base \"N\": " << m_filter_stats[BASE_NN] 
			<< " (" << (100.0*m_filter_stats[BASE_NN])/m_filter_stats[TOTAL_LENGTH] << " %)\n";
		fout << "  Reads Filtered by low complexity ratio (" 
			<< setprecision(1)
			<< m_opt.low_complexity_cutoff_ratio
			<< setprecision(2)
			<< "): " << m_filter_stats[READ_LOW_COMPLEXITY] 
			<< " ("
			<< (100.0*m_filter_stats[READ_LOW_COMPLEXITY])/m_filter_stats[TOTAL_NUMBER]
			<< " %)\n";
		fout << "  Bases Filtered by low complexity ratio: " 
			<< m_filter_stats[BASE_LOW_COMPLEXITY] << " (" 
			<< (100.0*m_filter_stats[BASE_LOW_COMPLEXITY])/m_filter_stats[TOTAL_LENGTH]
			<< " %)\n";

		if(m_opt.average_quality > 0.0){
		
			fout << "  Reads Filtered by avg quality (" << m_opt.average_quality 
				<< "): " << m_filter_stats[READ_AVG_Q] 
				<< " ("
				<< (100.0*m_filter_stats[READ_AVG_Q])/m_filter_stats[TOTAL_NUMBER]
				<< " %)\n";
			fout << "  Bases Filtered by avg quality: " << m_filter_stats[BASE_AVG_Q] 
				<< " (" 
				<< (100.0*m_filter_stats[BASE_AVG_Q])/m_filter_stats[TOTAL_LENGTH]
				<< " %)\n";
		}

		if(m_opt.filter_phiX){

			fout << "  Reads Filtered by phiX sequence: " 
				<< m_filter_stats[READ_PHIX]
				<< " (" 
				<< (100.0*m_filter_stats[READ_PHIX])/m_filter_stats[TOTAL_NUMBER]
				<< " %)\n";
			fout << "  Bases Filtered by phiX sequence: " << m_filter_stats[BASE_PHIX]
				<< " (" 
				<< (100.0*m_filter_stats[BASE_PHIX])/m_filter_stats[TOTAL_LENGTH] 
				<< " %)\n";
		}
		
		fout << "  Reads Trimmed by quality (" 
			<< setprecision(1)
			<< float(m_opt.quality)
			<< setprecision(2)
			<< "): " 
			<< m_filter_stats[READ_QUAL_TRIM]
			<< " (" 
			<< (100.0*m_filter_stats[READ_QUAL_TRIM])/m_filter_stats[TOTAL_NUMBER]
			<< " %)\n";
		fout << "  Bases Trimmed by quality: " << m_filter_stats[BASE_QUAL_TRIM]
			<< " ("
			<< (100.0*m_filter_stats[BASE_QUAL_TRIM])/m_filter_stats[TOTAL_LENGTH]
			<< " %)\n";

		if(m_opt.trim_5 > 0){
			fout << "  Reads Trimmed with " << m_opt.trim_5 << " bp from 5' end\n";
		}
		
		if(m_opt.trim_3 > 0){
			fout << "  Reads Trimmed with " << m_opt.trim_3 << " bp from 3' end\n";
		}
		
		if(m_opt.filter_adapter){
		
			fout << "  Reads Trimmed with Adapters/Primers: " 
				<< m_filter_stats[READ_ADAPTER]
				<< " ("
				<< (100.0*m_filter_stats[READ_ADAPTER])/m_filter_stats[TOTAL_NUMBER]
				<< " %)\n";				
			fout << "  Bases Trimmed with Adapters/Primers: " 
				<< m_filter_stats[BASE_ADAPTER] 
				<< " (" 
				<< (100.0*m_filter_stats[BASE_ADAPTER])/m_filter_stats[TOTAL_LENGTH]
				<< " %)\n";
			
			// Write the adapters in order of the number of reads matches
			deque< pair<size_t, string> > adapter;
			
			for(MAP< string, pair<size_t, size_t> >::const_iterator i = m_adapter_stats.begin();
				i != m_adapter_stats.end();++i){
				
				adapter.push_back( make_pair(i->second.first, i->first) ); // {reads, adapter name}
			}
			
			// Sort in *ascending* order ...
			sort( adapter.begin(), adapter.end() );
			
			// ... and iterate from the back to print in *descending* order
			for(deque< pair<size_t, string> >::const_reverse_iterator i = adapter.rbegin();
				i != adapter.rend();++i){
				
				MAP< string, pair<size_t, size_t> >::const_iterator iter = 
					m_adapter_stats.find(i->second);
					
				if( iter == m_adapter_stats.end() ){
				
					cerr << "Unable to find adapter statistics!" << endl;
					continue;
				}
				
				const size_t &affected_reads = iter->second.first;
				const size_t &affected_bases = iter->second.second;
				
				fout << "    " << iter->first << " " << affected_reads 
					<< " reads (" << (100.0*affected_reads)/m_filter_stats[TOTAL_NUMBER]
					<< " %) " << affected_bases << " bases ("
					<< (100.0*affected_bases)/m_filter_stats[TOTAL_LENGTH]
					<< " %)\n";
			}
		}
		
		if(m_opt.replace_N){
			fout << "\nN base random substitution: A " << m_filter_stats[N_TO_A]
				<< ", T " << m_filter_stats[N_TO_T]
				<< ", C " << m_filter_stats[N_TO_C]
				<<", G " << m_filter_stats[N_TO_G] << "\n";
		}
	}
}

void remove_existing_files(const Options &m_opt)
{
	remove_file(m_opt.plots_file);
	remove_file(m_opt.stats_file);
	remove_file(m_opt.trimmed_read1_file);
	remove_file(m_opt.trimmed_read2_file);
	remove_file(m_opt.trimmed_unpaired_file);
	remove_file(m_opt.trimmed_discard_file);
}

void remove_file(const string &m_filename)
{
	if( file_exits(m_filename) ){
		
		cerr << "The output " << m_filename << " file exists and will be overwritten." << endl;
		unlink( m_filename.c_str() );
	}
}
