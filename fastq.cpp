#include <string.h>
#include "fastq.h"

#include <iostream>

using namespace std;

bool next_read(gzFile m_fin, string &m_def, string &m_seq, string &m_quality)
{

	// The FASTQ file format is as follows for each record 
	// (from the Wikipedia entry):
	// Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
    	// Line 2 is the raw sequence letters.
    	// Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
    	// Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

	// Read each line of the fastq file in chunks of size "buffer_len" bytes. As of 2017, most platforms (i.e. Illumina) have a line length 
	// that is less than 256 bytes. However, a notable exception is the PacBio platform (which can have line lengths greater
	// than 1 KB). The value of buffer_len may need to be tuned to the sequencing platform, or, implement an adaptive strategy
	// that grows the buffer size to a platform-optimal value.
	const int buffer_len = 4096000;
	char buffer[buffer_len];
	char *ptr = NULL;
	
	/////////////////////////////////////////////////////////////////////////////////
	// Read the header
	/////////////////////////////////////////////////////////////////////////////////
	
	m_def.clear(); // Start with an emply header
	
	while(true){
	
		if(gzgets(m_fin, buffer, buffer_len) == NULL){

			if( gzeof(m_fin) ){
				return false;
			}

			throw __FILE__ ":next_read: Unable to read header";
		}
		
		if( ( ptr = strpbrk(buffer, "\n\r") ) != NULL ){
			
			*ptr = '\0';
			
			m_def += buffer;
			break;
		}
		
		m_def += buffer;
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// Read the sequence
	/////////////////////////////////////////////////////////////////////////////////
	
	m_seq.clear(); // Start with an emply sequence
	
	while(true){
	
		if(gzgets(m_fin, buffer, buffer_len) == NULL){
			
			cerr << "Error in read: " << m_def << endl;
			throw __FILE__ ":next_read: Unable to read sequence";
		}

		if( ( ptr = strpbrk(buffer, "\n\r") ) != NULL ){
			
			*ptr = '\0';
			
			m_seq += buffer;
			break;
		}
		
		m_seq += buffer;
	}

	// Skip the '+'. Note that we *don't* actually test to make sure we read a single '+' on a line
	// by itself.
	if(gzgets(m_fin, buffer, buffer_len) == NULL){
		
		cerr << "Error in read: " << m_def << endl;
		throw __FILE__ ":next_read: Unable to read '+'";
	}
	
	if( strpbrk(buffer, "\n\r") == NULL ){
	
		cerr << "Error in read: " << m_def << endl;
		throw __FILE__ ":next_read: Error reading '+' delimiter";
	}

	/////////////////////////////////////////////////////////////////////////////////
	// Read the quality
	/////////////////////////////////////////////////////////////////////////////////

	m_quality.clear(); // Start with an emply quality string
	
	while(true){
	
		if(gzgets(m_fin, buffer, buffer_len) == NULL){
			
			cerr << "Error in read: " << m_def << endl;
			throw __FILE__ ":next_read: Unable to read quality";
		}

		if( ( ptr = strpbrk(buffer, "\n\r") ) != NULL ){
			
			*ptr = '\0';
			
			m_quality += buffer;
			break;
		}
		
		m_quality += buffer;
	}
	
	if( m_seq.size() != m_quality.size() ){
		
		cerr << "Error in read: " << m_def << endl;
		throw __FILE__ ":next_read: |Sequence| != |Quality|";
	}
		
	return true;
}

void write_read(gzFile m_fout, const string &m_def, const string &m_seq, 
	const string &m_quality)
{
	//m_fout << m_def << '\n' << m_seq << "\n+\n" << m_quality << endl;
	
	gzputs( m_fout, m_def.c_str() );
	gzputc( m_fout, '\n' );
	gzputs( m_fout, m_seq.c_str() );
	gzputs( m_fout, "\n+\n" );
	gzputs( m_fout, m_quality.c_str() );
	gzputc( m_fout, '\n' );
}
