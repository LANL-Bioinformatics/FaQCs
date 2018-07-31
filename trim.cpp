#include "FaQCs.h"
#include "fastq.h"
#include "seq_overlap.h"

#include <iostream>
#include <math.h>

using namespace std;

bool trim_read(Read &m_read, 
	std::vector<size_t> &m_filter_stats, 
	MAP<Word, size_t> &m_kmer_table,
	PlotInfo &m_info, const Options &m_opt);

unsigned int hard_trim(std::string &m_seq, std::string &m_qual, 
	const Options &m_opt);

unsigned int BWA_trim(std::string &m_seq, std::string &m_qual, 
	const Options &m_opt);
	
unsigned int BWA_plus_trim(std::string &m_seq, std::string &m_qual, 
	const Options &m_opt);

unsigned int count_poly_n(const std::string &m_seq);
float average_quality(const std::string &m_quality, const char m_quality_offset);

void update_quality_matrix(matrix<size_t> &m_data, const string &m_qual, 
	const unsigned int &m_offset_5, const char &m_quality_offset);
void update_base_statistics(matrix<size_t> &m_matrix_data, 
	vector<NucleotideCount> &m_comp_data, const string &m_seq, 
	const unsigned int &m_offset_5);
void update_length_histogram(vector<size_t> &m_hist, 
	const unsigned int &m_len);
void update_kmer(MAP<Word, size_t> &m_kmer_table, const string &m_seq, 
	const unsigned int &m_k);
unsigned int trim_adapters_and_phiX(Read &m_read);

void trim_adapters_and_phiX(vector<Read> &m_buffer,
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	const Options &m_opt);
pair<unsigned int, unsigned int> find_mask_range(const vector<bool> &m_mask);

void mask_quality_terminal_N(const string &m_seq, string &m_qual, 
	const Options &m_opt);

// Accumulation for vectors
template <class T>
inline vector<T>& operator+=(vector<T> &m_lhs, const vector<T> &m_rhs)
{
	const size_t len = m_rhs.size();
	
	if( len /*src*/ > m_lhs.size() /*dest*/){
		
		// If the source vector is greater than the destination vector, 
		// grow the destination (left hand) vector to be the same size 
		// as the source (right hand) vector.
		m_lhs.resize(len);
	}
	
	for(size_t i = 0;i < len;++i){
		m_lhs[i] += m_rhs[i];
	}
	
	return m_lhs;
}

void trim(vector<Read> &m_buffer, 
	vector<size_t> &m_filter_stats, 
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	MAP<Word, size_t> &m_kmer_table,
	PlotInfo &m_info,
	Options &m_opt)
{
	const unsigned int buffer_size = m_buffer.size();
		
	#pragma omp parallel
	{
		
		// Accumulate filtering statistics in local buffers that will
		// later be merged to accumulate global statistics
		vector<size_t> local_filter_stats(FilterStat::NUM_STAT);
		MAP< string, pair<size_t /*read*/, size_t /*base*/> > local_adapter_stats;
		MAP<Word, size_t> local_kmer_table;
		
		PlotInfo local_info;
		
		if(m_opt.filter_adapter || m_opt.filter_phiX){
			trim_adapters_and_phiX(m_buffer, local_adapter_stats, m_opt);
		}
		
		#pragma omp for
		for(unsigned int i = 0;i < buffer_size;++i){

			// A reference to the i^th read
			Read &r = m_buffer[i];
			
			const bool valid = trim_read(
				r,
				local_filter_stats,
				local_kmer_table,
				local_info, m_opt);
				
			if(!valid){
				r.seq = r.qual = "";
			}
		}
		
		// Update the statistics
		#pragma omp critical
		{
			m_filter_stats += local_filter_stats;
			
			for(MAP< string, pair<size_t, size_t> >::const_iterator i = local_adapter_stats.begin();
				i != local_adapter_stats.end();++i){
				
				pair<size_t, size_t> &ref = m_adapter_stats[i->first];
				
				ref.first += i->second.first;
				ref.second += i->second.second;
			}
		
			for(MAP<Word, size_t>::const_iterator i = local_kmer_table.begin();i != local_kmer_table.end();++i){
				m_kmer_table[i->first] += i->second;
			}
	
			m_info.pre_quality_matrix += local_info.pre_quality_matrix;
			m_info.post_quality_matrix += local_info.post_quality_matrix;
			
			m_info.pre_base_matrix += local_info.pre_base_matrix;
			m_info.post_base_matrix += local_info.post_base_matrix;
				
			m_info.pre_read_quality_histogram += local_info.pre_read_quality_histogram;
			m_info.post_read_quality_histogram += local_info.post_read_quality_histogram;
			
			m_info.pre_base_quality_histogram += local_info.pre_base_quality_histogram;
			m_info.post_base_quality_histogram += local_info.post_base_quality_histogram;
			
			m_info.pre_nuc_composition += local_info.pre_nuc_composition;
			m_info.post_nuc_composition += local_info.post_nuc_composition;
			
			m_info.pre_length_histogram += local_info.pre_length_histogram;
			m_info.post_length_histogram += local_info.post_length_histogram;
		}
	}

	if(m_opt.kmer_rarefaction){
		
		// Compute the kmer rarefaction curve as "num_subsample" number of samples, where each sample
		// contains "split_size" number of reads.

		const size_t index = m_filter_stats[FilterStat::TOTAL_NUMBER]/m_opt.split_size;
		const size_t num_rarefaction = m_info.kmer_rarefaction.size();
		
		if( ( index > num_rarefaction) && (num_rarefaction < m_opt.num_subsample) ){
	
			Rarefaction local;

			local.num_seq = m_filter_stats[FilterStat::TOTAL_NUMBER];
			local.distinct_kmer = m_kmer_table.size();
			local.total_kmer = 0;

			for(MAP<Word, size_t>::const_iterator i = m_kmer_table.begin();i != m_kmer_table.end();++i){
				local.total_kmer += i->second;
			}
	
			m_info.kmer_rarefaction.push_back(local);
		}

		if(num_rarefaction >= m_opt.num_subsample){
		
			// Stop collecting kmers once we have completed the rarefaction curve
			m_opt.kmer_rarefaction = false;
		}
	}
}

string parse_id(const string &m_def)
{
	const size_t len = m_def.size();
	
	// Return the defline up to the first space
	string::size_type loc = m_def.find(' ');
	
	if(loc == string::npos){
		
		// No space was found, use the entire defline
		loc = len;
	}
	
	// From the FaQCs.pl script:
	// Special parsing is required for fastq variants from the NCBI SRA. These
	// fastq files have a defline that indicates the read (i.e. read 1 or read 2),
	// such as: 
	//	XXXXX.1 or XXXXX.2
	// or
	//	XXXXX/1 or XXXXX/2
	
	// Are there at least *two* characters in the
	// defline (i.e. loc > 1) and is the last symbol a digit?
	if( (loc > 1) && isdigit(m_def[loc - 1]) ){
		
		// Is the next to last symbol a '.' or a '/'?
		if( (m_def[loc - 2] == '.') || (m_def[loc - 2] == '/') ){
			
			// We have a match! Remove the last two characters
			loc -= 2;
		}
	}
	
	return m_def.substr(0, loc);
}

// From the qc_process() perl function
bool trim_read(Read &m_read, 
	vector<size_t> &m_filter_stats, 
	MAP<Word, size_t> &m_kmer_table,
	PlotInfo &m_info, const Options &m_opt)
{
	bool ret = true;
	
	unsigned int len = m_read.seq.size();
	
	// To align quality and base statistics, track any trimming 
	// that occurs at the start (i.e. 5' end) of the read.
	unsigned int offset_5 = 0;
	
	++ m_filter_stats[FilterStat::TOTAL_COUNT];
	++ m_filter_stats[FilterStat::TOTAL_NUMBER];
	m_filter_stats[FilterStat::TOTAL_LENGTH] += len;
	
	// It appears that some sequencing platoforms will produce reads that have an 'N'
	// with a corresponding high quality score! Replace the quality for
	// *terminal* 'N' bases with a quality score of 0
	mask_quality_terminal_N(m_read.seq, m_read.qual, m_opt);
	
	update_quality_matrix(m_info.pre_quality_matrix, m_read.qual, offset_5, 
		m_opt.input_quality_offset);
	update_base_statistics(m_info.pre_base_matrix, m_info.pre_nuc_composition, 
		m_read.seq, offset_5);
	update_length_histogram(m_info.pre_length_histogram, len);
	
	// Compute the average quality and truncate to an integer
	int quality_bin = int( average_quality(m_read.qual, m_opt.input_quality_offset) );
	
	// Use this truncated average quality to update the *pre* trimming histograms
	++ m_info.pre_read_quality_histogram[quality_bin];
	m_info.pre_base_quality_histogram[quality_bin] += len;

	if(m_opt.qc_only && m_opt.kmer_rarefaction){
		update_kmer(m_kmer_table, m_read.seq, m_opt.kmer);
	}
	
	if( len != m_read.qual.size() ){
		
		cerr << m_read.def << " sequence length is no equal to quality string length. It will be filtered." << endl;
		return false;
	}
	
	if(m_opt.filter_adapter || m_opt.filter_phiX){
		
		// The adapter and phiX statistics have *already* been updated.
		// We only need to trim the actual read (and recompute the read length).
		offset_5 += trim_adapters_and_phiX(m_read);
		
		len = m_read.seq.size();
	}
	
	if (m_opt.trim_5 && !m_opt.qc_only){
		
		if(m_opt.trim_5 > len){
		
			m_read.seq = "";
            		m_read.qual = "";

	    		len = 0;
			offset_5 += len;
		}
		else{
		
			m_read.seq = m_read.seq.substr(m_opt.trim_5, len - m_opt.trim_5);
            		m_read.qual = m_read.qual.substr(m_opt.trim_5, len - m_opt.trim_5);

	    		len -= m_opt.trim_5;
			offset_5 += m_opt.trim_5;
		}
        }
	
	if (m_opt.trim_3 && !m_opt.qc_only){
		
		if(m_opt.trim_3 > len){
		
			m_read.seq = "";
            		m_read.qual = "";

	    		len = 0;
		}
		else{
			m_read.seq = m_read.seq.substr(0, len - m_opt.trim_3);
			m_read.qual = m_read.qual.substr(0, len - m_opt.trim_3);

            		len -= m_opt.trim_3;
		}
        }
	
	// Apply length filter
        if( (len < m_opt.min_read_length) || (len == 0) ){
	
		m_filter_stats[FilterStat::BASE_LENGTH] += len;
		++ m_filter_stats[FilterStat::READ_LENGTH];
		
		ret = false;
        }
	
	if( !m_opt.qc_only && ret ){
	
		const unsigned int init_len = len;
		
		// Perform the quality trim
		switch(m_opt.mode){
			case Options::HARD:
				offset_5 += hard_trim(m_read.seq, m_read.qual, m_opt);
				break;
			case Options::BWA:
				offset_5 += BWA_trim(m_read.seq, m_read.qual, m_opt);
				break;
			case Options::BWA_plus:
				offset_5 += BWA_plus_trim(m_read.seq, m_read.qual, m_opt);
				break;
			default:
				throw __FILE__ ":trim_read: Undefined trimming mode!";
		};

		len = m_read.qual.size();
		
		if(init_len != len){
			
			m_filter_stats[FilterStat::BASE_QUAL_TRIM] += init_len - len;
			++ m_filter_stats[FilterStat::READ_QUAL_TRIM];
		}
		
		// Re-apply length filter
        	if( (len < m_opt.min_read_length) || (len == 0) ){
		
			m_filter_stats[FilterStat::BASE_LENGTH] += len;
			++ m_filter_stats[FilterStat::READ_LENGTH];
		
			ret = false;
        	}
	}
	
	// Apply the "N" filter
	if(ret && (count_poly_n(m_read.seq) >= m_opt.max_num_poly_N) ){
		
		m_filter_stats[FilterStat::BASE_NN] += len;
		++ m_filter_stats[FilterStat::READ_NN];
			
		if( !m_opt.qc_only ){
			ret = false;
		}
	}
		
	// Apply the average quality filter
	const float ave_Q = average_quality(m_read.qual, m_opt.input_quality_offset);
	
	if( ret && (ave_Q < m_opt.average_quality) ){
		
		m_filter_stats[FilterStat::BASE_AVG_Q] += len;
		++ m_filter_stats[FilterStat::READ_AVG_Q];
		
		ret = false;
	}
	
	// Update the read length
	len = m_read.qual.size();
	
	// Apply the low complexity filter
	if( ret && (len != 0) ){
		
		if(m_opt.replace_to_N_q > 0){
			
			for(unsigned int i = 0;i < len;++i){
				
				// **Unlike** the FaQCs.pl script, we are not currently testing
				// for sequencing platform (i.e. NextSeq).
				if( (m_read.seq[i] == 'G') && (
					quality_score(m_read.qual[i], m_opt.input_quality_offset) 
						< int(m_opt.replace_to_N_q) ) ){

					m_read.seq[i] = 'N';
				}
			}
		}
	
		// Filter on both single nucletide composition *and* dinucleotide
		// composition
		unsigned int num_A = 0;
		unsigned int num_T = 0;
		unsigned int num_G = 0;
		unsigned int num_C = 0;
		
		// When tracking the dinucleotide composition, we only need to
		// track four bases. {A, T, C, G}. Represent each base with a
		// two bit number: A = 0, T = 1, C = 2, G = 3.
		// A pair of nucleotides can be stored with four bits and can
		// take on 2^4 (=16) possible values (in point of fact, since we
		// don't count the same nucleotide twice, the maximum index to the
		// dinucleotide vector would be (G << 2) | C == 13
		#define	NUM_DINUCLEOTIDE	16 // 2^4
		#define	TWO_BIT_A		0
		#define	TWO_BIT_T		1
		#define	TWO_BIT_C		2
		#define	TWO_BIT_G		3
		#define	INVALID_BASE		4
		
		vector<unsigned int> dc(NUM_DINUCLEOTIDE); // Dinucleotide composition
		
		unsigned char last = INVALID_BASE;
		
		for(string::const_iterator i = m_read.seq.begin();i != m_read.seq.end();++i){
	
			switch(*i){
				case 'A': case 'a':
				
					++num_A;
					
					// Don't count runs of the *same* nucleotide
					if( (TWO_BIT_A != last) && (last != INVALID_BASE) ){
						++ dc[(last << 2) | TWO_BIT_A];
					}
					
					last = TWO_BIT_A;
					break;
				case 'T': case 't':
				
					++num_T;
					
					// Don't count runs of the *same* nucleotide
					if(  (TWO_BIT_T != last) && (last != INVALID_BASE) ){
						++ dc[(last << 2) | TWO_BIT_T];
					}
					
					last = TWO_BIT_T;
					break;
				case 'G': case 'g':
				
					++num_G;
					
					// Don't count runs of the *same* nucleotide
					if( (TWO_BIT_G != last) && (last != INVALID_BASE) ){
						++ dc[(last << 2) | TWO_BIT_G];
					}
					
					last = TWO_BIT_G;
					break;
				case 'C': case 'c':
				
					++num_C;
					
					// Don't count runs of the *same* nucleotide
					if( (TWO_BIT_C != last) && (last != INVALID_BASE) ){
						++ dc[(last << 2) | TWO_BIT_C];
					}
					
					last = TWO_BIT_C;
					break;
				default:
					last = INVALID_BASE;
					break;
			};
		}		
		
		float norm = 1.0/len;
		
		if( (num_A*norm > m_opt.low_complexity_cutoff_ratio) ||
		    (num_T*norm > m_opt.low_complexity_cutoff_ratio) ||
		    (num_G*norm > m_opt.low_complexity_cutoff_ratio) ||
		    (num_C*norm > m_opt.low_complexity_cutoff_ratio) ){

			m_filter_stats[FilterStat::BASE_LOW_COMPLEXITY] += len;
			m_filter_stats[FilterStat::READ_LOW_COMPLEXITY] ++;

			ret = false;
		}
		else{

			// Adjust the normalization to account for the maximum number
			// of dinucleotides in a sequence.
			norm *= 2.0;

			for(vector<unsigned int>::const_iterator i = dc.begin();i != dc.end();++i){

				if((*i)*norm > m_opt.low_complexity_cutoff_ratio){

					m_filter_stats[FilterStat::BASE_LOW_COMPLEXITY] += len;
					m_filter_stats[FilterStat::READ_LOW_COMPLEXITY] ++;

					ret = false;
					break;
				}
			}
		}
	}
	
	// Convert quality score format (if needed)
	if( ret && (m_opt.input_quality_offset != m_opt.output_quality_offset) ){
		
		for(string::iterator i = m_read.qual.begin();i != m_read.qual.end();++i){
			*i = quality_score(*i, m_opt.input_quality_offset) + m_opt.output_quality_offset;
			
			if(*i < 0){
				throw __FILE__ ": quality error!";
			}
		}
	}
	
	if(ret){
	
		m_filter_stats[FilterStat::TOTAL_TRIMMED_LENGTH] += len;
		++ m_filter_stats[FilterStat::TOTAL_TRIMMED_NUMBER];
		
		update_quality_matrix(m_info.post_quality_matrix, m_read.qual, offset_5, 
			m_opt.output_quality_offset);
		update_base_statistics(m_info.post_base_matrix, m_info.post_nuc_composition,
			m_read.seq, offset_5);
		update_length_histogram(m_info.post_length_histogram, len);
		
		// Truncate the trimmed average quality from a float to an integer
		quality_bin = int(ave_Q);
	
		// Use this truncated average quality to update the *post* trimming histograms
		++ m_info.post_read_quality_histogram[quality_bin];
		m_info.post_base_quality_histogram[quality_bin] += len;

		if(!m_opt.qc_only && m_opt.kmer_rarefaction){
                	update_kmer(m_kmer_table, m_read.seq, m_opt.kmer);
        	}
	}
	
	return ret;
}

float average_quality(const string &m_quality, const char m_quality_offset)
{
	int total_qual = 0;
	
	for(string::const_iterator i = m_quality.begin();i != m_quality.end();++i){
	//for(string::const_iterator i = m_quality.begin(), end = m_quality.end();i != end;++i){
	
		// As an optimization, don't call the quality_score function here, since
		// it returns the difference between the encoded quality and the m_quality_offset
		// *for each base*. Just compute the average encoded quality and subtract the 
		// m_quality_offset at the end.

		total_qual += *i;
	}
	
	if( !m_quality.empty() ){
		
		// Clamp the return value to be greater than or equal to zero (to handle the edge case
		// of old quality formats, like solexa+64, that allowed negative quality scores.
		return max( 0.0f, float(total_qual)/m_quality.size() - float(m_quality_offset) );
	}
	
	return 0.0;
}

unsigned int count_poly_n(const string &m_seq)
{
	unsigned int max_poly_n = 0;
	unsigned int curr_poly_n = 0;
	
	for(string::const_iterator i = m_seq.begin();i != m_seq.end();++i){
	//for(string::const_iterator i = m_seq.begin(), end = m_seq.end();i != end;++i){
	
		if(*i == 'N'){
			
			++curr_poly_n;
			max_poly_n = max(max_poly_n, curr_poly_n);
		}
		else{
			curr_poly_n = 0;
		}
	}
	
	return max_poly_n;
}

char auto_detect_quality_offset(const vector<Read> &m_buffer)
{
	for(vector<Read>::const_iterator i = m_buffer.begin();i != m_buffer.end();++i){
	
		for(string::const_iterator j = i->qual.begin();j != i->qual.end();++j){
	
			if(*j > 74){
				return 64; // Solexa/Illumina
			}

			if(*j < 59){
				return 33; // Sanger
			}
		}		
	}
	
	throw __FILE__ ":auto_detect_quality_offset: Unknown quality format!";
	return 0;
}

bool auto_detect_next_seq(const vector<Read> &m_buffer)
{
	for(vector<Read>::const_iterator i = m_buffer.begin();i != m_buffer.end();++i){
		return (i->def.find("@NS") == 0);
	}
	
	return false;
}

// Return the number bases trimmed from the 5' end
unsigned int hard_trim(string &m_seq, string &m_qual, 
	const Options &m_opt)
{
	const int len = m_qual.size();
	
	// Trim the 3' end
	int pos_3 = len - 1;
	int final_pos_5 = 0;
	int final_pos_3 = pos_3;

	while(pos_3 > 0){
	
		if( m_opt.quality < quality_score(m_qual[pos_3], m_opt.input_quality_offset) ){
            		
			final_pos_3 = pos_3;
            		break;
        	}
		
        	--pos_3;
    	}
	
    	// Trim the 5' end
	if (!m_opt.protect_5){

		
    		int pos_5 = final_pos_5;
		
    		while(pos_5 < pos_3){
		
        		if( m_opt.quality < quality_score(m_qual[pos_5], m_opt.input_quality_offset) ){
	    			
				final_pos_5 = pos_5;
	    			break;
			}
			
        		++pos_5;
    		}
	}
	
	m_seq = m_seq.substr(final_pos_5, final_pos_3 - final_pos_5 + 1);
	m_qual = m_qual.substr(final_pos_5, final_pos_3 - final_pos_5 + 1);
	
	return (unsigned int)(final_pos_5);
}

// Return the number bases trimmed from the 5' end
unsigned int BWA_trim(string &m_seq, string &m_qual, 
	const Options &m_opt)
{
	const int len = m_qual.size();
	
	// Cast the quality threshold to a int to avoid the overhead of
	// casting for every comparison.
	const int Q = (int)m_opt.quality;
	
	// Trim the 3' end -- please note that the 5' end is *not* trimmed!
	int pos_3 = len - 1;
	int final_pos_3 = pos_3;
	const int final_pos_5 = 0; 
  
        int area = 0;  
        int maxArea = 0;

	while( (pos_3 > 0) && (area >= 0) ){
	
		area += Q - (int)quality_score(m_qual[pos_3], m_opt.input_quality_offset);
		
		if(area > maxArea){
		
			maxArea = area;
			final_pos_3 = pos_3 - 1;
		}
		
		--pos_3;
	}
	
	m_seq = m_seq.substr(final_pos_5, final_pos_3 - final_pos_5 + 1);
	m_qual = m_qual.substr(final_pos_5, final_pos_3 - final_pos_5 + 1);
	
	return (unsigned int)(final_pos_5);
}
	

// Return the number bases trimmed from the 5' end
// From the bwa_trim_plus() perl function
unsigned int BWA_plus_trim(string &m_seq, string &m_qual, 
	const Options &m_opt)
{
	const int len = m_qual.size();
	
	// Cast the quality threshold to a int to avoid the overhead of
	// casting for every comparison.
	const int Q = (int)m_opt.quality;
	
	int at_least_scan = min(5, len);
   	int num_after_neg = min(2, len);
	int pos_3 = len - 1;
	int final_pos_5 = 0;
	int final_pos_3 = pos_3;
	int area = 0;
	int maxArea = 0;
	
	// Trim 3' end
	while(at_least_scan){
		
		--at_least_scan;
		
		if( (pos_3 > num_after_neg) && (area >= 0) ){
			at_least_scan = num_after_neg;
		}
		
	    	area += Q - (int)quality_score(m_qual[pos_3], m_opt.input_quality_offset);
		
		if(area > maxArea){
		
			maxArea = area;
			final_pos_3 = pos_3 - 1;
		}
		
		--pos_3;
	}
	
	// trim 5' end
	if(!m_opt.protect_5){
    	
		int pos_5 = 0; 
    		
		maxArea = 0;
		area = 0;
    		
		at_least_scan = min(5, len);

    		while(at_least_scan){
		
       			--at_least_scan;
			
        		if( (pos_5 < (final_pos_3 - num_after_neg) ) && (area >= 0) ){
				at_least_scan = num_after_neg;
			}
			
			area += Q - (int)quality_score(m_qual[pos_5], m_opt.input_quality_offset);
			
	    		if(area > maxArea){
			
				maxArea = area;
				final_pos_5 = pos_5 + 1;
	    		}
			
	    		++pos_5;
    		}
	}
	    
	if(final_pos_3 <= final_pos_5){

        	// This read has been trimmed away!
		m_seq = "";
		m_qual = "";
	}
	else{
        	m_seq = m_seq.substr(final_pos_5, final_pos_3 - final_pos_5 + 1);
		m_qual = m_qual.substr(final_pos_5, final_pos_3 - final_pos_5 + 1);
	}
	
	return (unsigned int)(final_pos_5);
}

void update_quality_matrix(matrix<size_t> &m_data, const string &m_qual, 
	const unsigned int &m_offset_5, const char &m_quality_offset)
{
	const unsigned int len = m_qual.size();
	const unsigned int full_len = len + m_offset_5;
	
	if(m_data.get_num_row() < full_len){
		m_data.resize(full_len, MAX_QUALITY_SCORE + 1);
	}
	
	for(unsigned int i = 0;i < len;++i){		
		++ m_data( i + m_offset_5, quality_score(m_qual[i], m_quality_offset) );
	}
}

void update_base_statistics(matrix<size_t> &m_matrix_data, 
	vector<NucleotideCount> &m_comp_data, const string &m_seq, 
	const unsigned int &m_offset_5)
{
	const unsigned int len = m_seq.size();
	const unsigned int full_len = len + m_offset_5;
	
	if(m_matrix_data.get_num_row() < full_len){
		m_matrix_data.resize(full_len, NUM_BASE);
	}
	
	unsigned int num_A = 0;
	unsigned int num_T = 0;
	unsigned int num_C = 0;
	unsigned int num_G = 0;
	unsigned int num_N = 0;
	
	unsigned int index = m_offset_5;
	
	for(string::const_iterator i = m_seq.begin();i != m_seq.end();++i, ++index){
	
		switch(*i){
			case 'A': case 'a':
			
				++ num_A;
				++ m_matrix_data(index, BASE_A);
				break;
			case 'T': case 't':
			
				++ num_T;
				++ m_matrix_data(index, BASE_T);
				break;
			case 'C': case 'c':
			
				++ num_C;
				++ m_matrix_data(index, BASE_C);
				break;
			case 'G': case 'g':
			
				++ num_G;
				++ m_matrix_data(index, BASE_G);
				break;
			case 'N': case 'n':
			
				++ num_N;
				++ m_matrix_data(index, BASE_N);
				break;
		};
	}
	
	const float norm = (len > 0) ? float(NUM_COMPOSITION_BIN - 1)/len : 0.0;
	
	++ m_comp_data[norm*num_A].num_A;

	++ m_comp_data[norm*num_T].num_T;
	
	const unsigned int index_C = norm*num_C;
	++ m_comp_data[index_C].num_C;
	
	const unsigned int index_G = norm*num_G;
	++ m_comp_data[index_G].num_G;
	
	++ m_comp_data[norm*num_N].num_N;

	++ m_comp_data[index_G + index_C].num_GC;
}

void update_length_histogram(vector<size_t> &m_hist, const unsigned int &m_len)
{
	
	if(m_hist.size() <= m_len){
		m_hist.resize(m_len + 1);
	}
	
	++ m_hist[m_len];
}

void update_kmer(MAP<Word, size_t> &m_kmer_table, const string &m_seq, const unsigned int &m_k)
{
	// Count the kmers in m_seq
	const Word comp_shift = 2*(m_k - 1);
	const Word mask = ( 1UL << (2*m_k) ) - 1;

	Word w = 0;
	Word comp = 0;
	unsigned int word_len = 0;

	for(string::const_iterator i = m_seq.begin();i != m_seq.end();++i){
	//for(string::const_iterator i = m_seq.begin(), end = m_seq.end();i != end;++i){

		++word_len;

		switch(*i){
			case 'A': case 'a':
				w = (w << 2) | BASE_A;
				comp = (comp >> 2) | (Word(BASE_T) << comp_shift);
				break;
			case 'T': case 't':
				w = (w << 2) | BASE_T;
				comp = (comp >> 2) | (Word(BASE_A) << comp_shift);
				break;
			case 'G': case 'g':
				w = (w << 2) | BASE_G;
				comp = (comp >> 2) | (Word(BASE_C) << comp_shift);
				break;
			case 'C': case 'c':
				w = (w << 2) | BASE_C;
				comp = (comp >> 2) | (Word(BASE_G) << comp_shift);
				break;
			default:
				word_len = 0;
				break;
		};

		if(word_len >= m_k){
			
			// Count cannonical kmers, where a cannonical kmer is defined
			// as the minimun of the kmer and its complement.
			++m_kmer_table[ min(w & mask, comp & mask) ];
		}
	}
}

// Return the number of bases trimmed from the 5' end of the read
unsigned int trim_adapters_and_phiX(Read &m_read)
{
	const unsigned int len = m_read.seq.size();
	
	if(len == m_read.start_length.second){
		
		// Read has not been modified
		return 0;
	}
		
	m_read.seq = m_read.seq.substr(m_read.start_length.first, m_read.start_length.second);
	m_read.qual = m_read.qual.substr(m_read.start_length.first, m_read.start_length.second);
	
	if(m_read.start_length.second == 0){
	
		// The entire read matches adapter -- all bases have been trimmed
		return len;
	}
	
	return m_read.start_length.first;
}

// ** NOTE **
// The trim_adapters_and_phiX() function is intended to be called within an OpenMP parallel section.
// All of the local variables within this function are therefor *thread local*!
// In addition, it is assumes that the m_adapter_stats variable that is passed to this function
// is also *thread local* (otherwise we will have a race condition on our hands).
void trim_adapters_and_phiX(vector<Read> &m_buffer, 
	MAP< string, pair<size_t /*read*/, size_t /*base*/> > &m_adapter_stats,
	const Options &m_opt)
{
	const unsigned int buffer_size = m_buffer.size();
	const unsigned int num_adapter = m_opt.adapter.size();
	
	// Work in terms of ratios of matching bases, not mismatched bases
	const float filter_adapter_match_rate = 1.0 - m_opt.filterAdapterMismatchRate;
	
	SO::SeqOverlap align(SO::SeqOverlap::SmithWaterman);

	vector< pair<unsigned int, vector<bool> > > mask;
	
	mask.reserve(SO_LEN);
	
	unsigned int current_slot = 0;
	
	#pragma omp for
	for(unsigned int i = 0;i < buffer_size;++i){

		// A reference to the i^th read
		Read &r = m_buffer[i];
		
		const size_t read_len = r.seq.size();
		
		// Set the default valid sequence range for read 1
		r.start_length.first = 0;
		r.start_length.second = read_len;
		
		mask.push_back( make_pair(i, vector<bool>(r.start_length.second, true) ) );
		
		align.pack_query(current_slot, r.seq);
		++current_slot;
		
		if(current_slot == SO_LEN){

			// Track the best matching adapter to each read
			vector< pair<SO::SO_Score, unsigned int> > best( current_slot, make_pair(0, 0) );
			
			// Align all adapters to each read
			for(unsigned int j = 0;j < num_adapter;++j){
			
				// The match threshold is number of bases, computed from the 
				// filter_adapter_match_rate and the *smaller* of read length and
				// adaptor sequence.
				const int match_threshold = filter_adapter_match_rate*
					min( read_len, m_opt.adapter[j].second.size() );
				
				align.pack_target(m_opt.adapter[j].second);
				
				align.align();
				
				for(unsigned int slot = 0;slot < current_slot;++slot){
			
					// Note that the scoring system is currently:
					// 1 = match
					// -1 = mismatch
					// So (match length + score)/2 is the number of perfect matches

					const SO::SO_Score score = align.score(slot);
					const pair<int, int> range = align.alignment_range_query(slot);

					const int match_length = range.second - range.first + 1;
					const int num_match = (match_length + score)/2;

					if(num_match >= match_threshold){

						vector<bool> &m = mask[slot].second;
						
						// We have a match! Mask the matching sequence.
						for(int k = range.first;k <= range.second;++k){
							m[k] = false;
						}

						if(score > best[slot].first){

							best[slot].first = score; // New best alignment score
							best[slot].second = j; // Adapter index
						}
					}
				}
			}
			
			// Update the valid range of each read
			for(unsigned int slot = 0;slot < current_slot;++slot){
				
				if(best[slot].first > 0){ // We have a match for this slot
					
					const pair< unsigned int, vector<bool> > &m = mask[slot];
					
					const pair<unsigned int, unsigned int> range = find_mask_range(m.second);
					
					m_buffer[mask[slot].first].start_length = range;
					
					// Indicate that this read has been trimmed by the best matching adapter
					const string &adapter = m_opt.adapter[best[slot].second].first;
					
					// Record which adapter sequence was responsible for trimming the largest number
					// of bases for this read.
					pair<size_t /*read*/, size_t /*base*/> &stat = m_adapter_stats[adapter];
					
					++ stat.first;
					stat.second += m.second.size() - range.second;
				}
			}
			
			mask.clear();
			current_slot = 0;
		}
	}
	
	// Handle any remaining sequences in the alignment buffer
	if(current_slot > 0){
		
		// Track the best matching adapter to each read
		vector< pair<SO::SO_Score, unsigned int> > best( current_slot, make_pair(0, 0) );

		// Align all adapters to each read
		for(unsigned int j = 0;j < num_adapter;++j){

			const int match_threshold = filter_adapter_match_rate*m_opt.adapter[j].second.size();

			align.pack_target(m_opt.adapter[j].second);

			align.align();

			for(unsigned int slot = 0;slot < current_slot;++slot){

				// Note that the scoring system is currently:
				// 1 = match
				// -1 = mismatch
				// So (match length + score)/2 is the number of perfect matches

				const SO::SO_Score score = align.score(slot);
				const pair<int, int> range = align.alignment_range_query(slot);

				const int match_length = range.second - range.first + 1;
				const int num_match = (match_length + score)/2;

				if(num_match >= match_threshold){
										
					vector<bool> &m = mask[slot].second;

					// We have a match! Mask the matching sequence.
					for(int k = range.first;k <= range.second;++k){
						m[k] = false;
					}

					if(score > best[slot].first){

						best[slot].first = score; // New best alignment score
						best[slot].second = j; // Adapter index
					}
				}
			}
		}

		// Update the valid range of each read
		for(unsigned int slot = 0;slot < current_slot;++slot){

			if(best[slot].first > 0){ // We have a match for this slot

				const pair< unsigned int, vector<bool> > &m = mask[slot];

				const pair<unsigned int, unsigned int> range = find_mask_range(m.second);
				
				m_buffer[mask[slot].first].start_length = range;

				// Indicate that this read has been trimmed by the best matching adapter
				const string &adapter = m_opt.adapter[best[slot].second].first;

				// Record which adapter sequence was responsible for trimming the largest number
				// of bases for this read.
				pair<size_t /*read*/, size_t /*base*/> &stat = m_adapter_stats[adapter];

				++ stat.first;
				stat.second += m.second.size() - range.second;
			}
		}
	}
}

pair<unsigned int, unsigned int> find_mask_range(const vector<bool> &m_mask)
{
	// Keep the *longest* contiguous run of valid (i.e. not masked) bases
	const unsigned int len = m_mask.size();

	unsigned int longest_run_start = 0;
	unsigned int longest_run_length = 0;

	unsigned int run_start = 0;
	unsigned int run_length = 0;

	for(unsigned int i = 0;i < len;++i){

		if(m_mask[i] == false){

			if(run_length > longest_run_length){

				longest_run_length = run_length;
				longest_run_start = run_start;

				run_length = 0;
			}
		}
		else{
			if(run_length == 0){
				run_start = i;
			}

			++run_length;
		}
	}

	if(run_length > longest_run_length){

		longest_run_length = run_length;
		longest_run_start = run_start;
	}

	if(longest_run_length == 0){

		// The entire read matches adapter sequence
		return make_pair(0, 0);
	}

	return make_pair(longest_run_start, longest_run_length);
}

void mask_quality_terminal_N(const string &m_seq, string &m_qual, const Options &m_opt)
{
	string::const_iterator s = m_seq.begin();
	string::iterator q = m_qual.begin();
	
	while( ( s != m_seq.end() ) && (*s == 'N') ){
		
		// By setting the quality to be equal to the input_quality_offset, we
		// effectively set the quality to zero
		*q = m_opt.input_quality_offset;
		++s;
		++q;
	}
	
	string::const_reverse_iterator sr = m_seq.rbegin();
	string::reverse_iterator qr = m_qual.rbegin();
	
	while( ( sr != m_seq.rend() ) && (*sr == 'N') ){
		
		// By setting the quality to be equal to the input_quality_offset, we
		// effectively set the quality to zero
		*qr = m_opt.input_quality_offset;
		++sr;
		++qr;
	}
}
