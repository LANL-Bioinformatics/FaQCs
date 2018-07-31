#ifndef __FASTQ
#define __FASTQ

#include <string>
#include <fstream>
#include <zlib.h>

// Make sure that we have been compiled against a version of zlib that
// supports writing uncompressed files
#if (ZLIB_VERNUM < 0x1280)
	#error "Please compile with a zlib version >= 1.2.8"
#endif

inline char quality_score(const char m_quality, const char m_offset)
{
	// Please note that negative scores *are* allowed! The solexa+64 quality scores
	// start at -5 (granted this is an old format, but we don't know what will the
	// future will hold).
	// Don't allow negative quality scores
	//if(m_offset > m_quality){		
	//	throw __FILE__ ":quality_score: Negative quality score!";
	//}
	
	// Clamp the actual quality score to be greater than zero. This will need to change
	// if a future format defines meaningfull negative quality scores!
	return std::max(0, m_quality - m_offset);
};

bool next_read(gzFile m_fin, std::string &m_def, std::string &m_seq, std::string &m_quality);

void write_read(gzFile m_fout, const std::string &m_def, const std::string &m_seq, 
	const std::string &m_quality);

#endif // __FASTQ
