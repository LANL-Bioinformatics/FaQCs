#ifndef __FAQCS
#define __FAQCS

#include <string>
#include <vector>
#include <unordered_map>
#include <limits.h>
#include "matrix.h"

#define FaQCs_VERSION	"2.08"

#define	AUTO_DETECT_QUALITY_OFFSET	SCHAR_MIN

#define	MIN_QUALITY_SCORE	0
#define	MAX_QUALITY_SCORE	41

#define	DEFAULT_NEXTSEQ_QUALITY_SCORE	20

// The number of bins for measuring base composition -> [0.000, 100.0]
// 10001 bins allows for a four digit number to be represented as an
// integer (i.e. xxxx).
#define	NUM_COMPOSITION_BIN	10001

#define	PATH_SEPARATOR	"/"

// Make it easy to switch between map and unordered_map
#define	MAP	std::unordered_map

// In order to treat the phiX sequence, and its complement, as
// adapters (when trimming), we need to give them unique names
#define	PHI_X			"__PhiX174_NC_001422__"
#define	PHI_X_COMPLEMENT	"__PhiX174_NC_001422_complement__"

typedef unsigned short int Dinucleotide;
typedef size_t Word;

enum {
	BASE_A,
	BASE_T,
	BASE_C,
	BASE_G,
	BASE_N,
	NUM_BASE
};

// Enumerate the types of filter statistics that we will collect
// and report to the user
namespace FilterStat {
	enum{
		TOTAL_COUNT,
		TOTAL_NUMBER,
		TOTAL_LENGTH,
		TOTAL_TRIMMED_NUMBER,
		TOTAL_TRIMMED_LENGTH,
		PAIRED_READ_NUMBER,
		PAIRED_BASE_LENGTH,
		READ_LENGTH,
		BASE_LENGTH,
		READ_NN,
		BASE_NN,
		READ_PHIX,
		BASE_PHIX,
		READ_ADAPTER,
		BASE_ADAPTER,
		READ_AVG_Q,
		BASE_AVG_Q,
		READ_QUAL_TRIM,
		BASE_QUAL_TRIM,
		READ_LOW_COMPLEXITY,
		BASE_LOW_COMPLEXITY,
		N_TO_A,
		N_TO_T,
		N_TO_G,
		N_TO_C,
		NUM_STAT
	};
}

struct Options
{
	bool print_usage;
	bool protect_5;
	bool replace_N;
	bool kmer_rarefaction;
	bool discard_output;
	bool qc_only;
	bool trim_only;
	bool filter_adapter;
	bool filter_phiX;
	bool debug;
	
	typedef enum {
		HARD, 
		BWA, 
		BWA_plus, 
		UNDEFINED_MODE
	}  Mode;
	
	Mode mode;
	
	std::string prefix;
	std::string plots_file;
	std::string stats_file;
	
	std::string input_read1_file;
	std::string input_read2_file;
	std::string input_unpaired_file;
	
	std::string trimmed_read1_file;
	std::string trimmed_read2_file;
	std::string trimmed_unpaired_file;
	std::string trimmed_discard_file;
	std::string output_dir;
	std::string artifact_file;
	
	float average_quality;
	float low_complexity_cutoff_ratio;
	float filterAdapterMismatchRate;
	
	char input_quality_offset;
	char output_quality_offset;
	unsigned int num_thread;
	char quality;
	unsigned int min_read_length;
	unsigned int max_num_poly_N;
	unsigned int kmer;
	unsigned int num_subsample;
	unsigned int trim_5;
	unsigned int trim_3;
	unsigned int split_size;
	unsigned int replace_to_N_q;
	
	std::vector< std::pair<std::string, std::string> > adapter;
	
	Options(int argc, char* argv[]);
	
	inline bool has_paired() const
	{
		return ( !input_read1_file.empty() && !input_read1_file.empty() );
	};
	
	inline bool has_unpaired() const
	{
		return ( !input_unpaired_file.empty() );
	};
};


struct Read
{
	std::string seq;
	std::string qual;
	std::string def;
	
	// The {start, length} after adapter and/or phiX trimming
	std::pair<unsigned int, unsigned int> start_length;
	
	inline bool empty() const
	{
		return def.empty();
	};
	
	inline bool valid() const
	{
		return !seq.empty();
	};
};

struct NucleotideCount
{
	size_t num_A;
	size_t num_T;
	size_t num_C;
	size_t num_G;
	size_t num_N;
	size_t num_GC;
	
	NucleotideCount() :
		num_A(0), num_T(0), num_C(0), num_G(0), num_N(0), num_GC(0)
	{
	};
	
	inline NucleotideCount& operator+=(const NucleotideCount &m_rhs)
	{
		num_A += m_rhs.num_A;
		num_T += m_rhs.num_T;
		num_C += m_rhs.num_C;
		num_G += m_rhs.num_G;
		num_N += m_rhs.num_N;
		num_GC += m_rhs.num_GC;
		
		return (*this);
	};
};

struct Rarefaction
{
	size_t num_seq;
	size_t distinct_kmer;
	size_t total_kmer;
};

// Information needed for plotting summary data
struct PlotInfo
{
	// Position x Quality score
	matrix<size_t> pre_quality_matrix;
	matrix<size_t> post_quality_matrix;
	
	// Position x Base score
	matrix<size_t> pre_base_matrix;
	matrix<size_t> post_base_matrix;
	
	// Quality histogram
	std::vector<size_t> pre_read_quality_histogram;
	std::vector<size_t> pre_base_quality_histogram;
	
	std::vector<size_t> post_read_quality_histogram;
	std::vector<size_t> post_base_quality_histogram;
	
	// Nucleotide composition
	std::vector<NucleotideCount> pre_nuc_composition;
	std::vector<NucleotideCount> post_nuc_composition;
	
	// Read length histogram
	std::vector<size_t> pre_length_histogram;
	std::vector<size_t> post_length_histogram;
	
	// Kmer
	MAP<size_t, size_t> kmer_frequency_histogram;
	std::vector<Rarefaction> kmer_rarefaction;
	
	PlotInfo()
	{
		pre_read_quality_histogram.resize(MAX_QUALITY_SCORE + 1);
		pre_base_quality_histogram.resize(MAX_QUALITY_SCORE + 1);
		
		post_read_quality_histogram.resize(MAX_QUALITY_SCORE + 1);
		post_base_quality_histogram.resize(MAX_QUALITY_SCORE + 1);
		
		pre_nuc_composition.resize(NUM_COMPOSITION_BIN);
		post_nuc_composition.resize(NUM_COMPOSITION_BIN);
	};
};

// In trim.cpp
void trim(std::vector<Read> &m_buffer, 
	std::vector<size_t> &m_filter_stats, 
	MAP< std::string, std::pair<size_t /*read*/, size_t /*base*/> > &m_adaptor_stats,
	MAP<Word, size_t> &m_kmer_table, PlotInfo &m_info, Options &m_opt);

bool auto_detect_next_seq(const std::vector<Read> &m_buffer);
char auto_detect_quality_offset(const std::vector<Read> &m_buffer);
std::string parse_id(const std::string &m_def);

// In plot.cpp
void plot(const PlotInfo &m_info, std::vector<size_t> &m_filter_stats,
	const Options &m_opt);
 
#endif // __ FAQCS
