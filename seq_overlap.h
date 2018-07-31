#ifndef __SEQ_OVERLAP
#define __SEQ_OVERLAP

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <deque>

/////////////////////////////////////////////////////////////////////////////////////////////
//
// SSE instruction defines, macros and datatypes. Note that using SSE4.1 requires
// g++ version >= 4.4 and the compiler flag -msse4.1. Currenlty (4/1/17) the Intel AVX2 instruction
// set (which offers 256 bit registers for the operations we need) will be a good future path 
// forward for modern hardware. From the AVX Wikipedia web page 
// (https://en.wikipedia.org/wiki/Advanced_Vector_Extension) chips supporting AVX2:
//
//    Intel
//        Haswell processor, Q2 2013
//        Haswell E processor, Q3 2014
//        Broadwell processor, Q4 2014
//        Broadwell E processor, Q3 2016
//        Skylake processor, Q3 2015
//        Kaby Lake processor, Q3 2016(ULV mobile)/Q1 2017(desktop/mobile)
//        Cannonlake processor, expected in 2017
//    AMD
//        Carrizo processor, Q2 2015
//        Ryzen processor, Q1 2017
//
//
// Note that ATA32 requires the sse4.1 instruction set.
// Note that ATA16 requires the sse2 instruction set.
//
// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
// appear to be necessary for 64-bit OSes.
/////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ATA32
	#include <immintrin.h>
#elif ATA16
	#include <xmmintrin.h>
#endif

namespace SO {

#define	SSE_ALIGNMENT_SIZE	16

#ifdef ATA32
	// Use a 32 bit integer data type for dynamic programming. Since SSE provides
	// a 128 bit integer data type, we can process 4 elements in parallel.
	#define		SO_LEN				4u
	#define		_MM_SET1			_mm_set1_epi32
	#define		_MM_ADD				_mm_add_epi32
	#define		_MM_SUB				_mm_sub_epi32
	#define		_MM_CMPEQ			_mm_cmpeq_epi32
	#define		_MM_CMPLT			_mm_cmplt_epi32
	#define		_MM_CMPGT			_mm_cmpgt_epi32
	#define		_MM_MAX				_mm_max_epi32
	#define		_MM_MULLO			_mm_mullo_epi32
	#define		_MM_MIN				_mm_min_epi32
	
	typedef 	int				SO_Score;
	typedef		__m128i				SIMDVector;
	
#elif ATA16
	// Use a 16 bit integer data type for dynamic programming. Since SSE provides
	// a 128 bit integer data type, we can process 8 elements in parallel.
	#define		SO_LEN				8u
	#define		_MM_SET1			_mm_set1_epi16
	#define		_MM_ADD				_mm_add_epi16
	#define		_MM_SUB				_mm_sub_epi16
	#define		_MM_CMPEQ			_mm_cmpeq_epi16
	#define		_MM_CMPLT			_mm_cmplt_epi16
	#define		_MM_CMPGT			_mm_cmpgt_epi16
	#define		_MM_MAX				_mm_max_epi16
	#define		_MM_MULLO			_mm_mullo_epi16
	#define		_MM_MIN				_mm_min_epi16
	
	typedef 	short int			SO_Score;
	typedef		__m128i				SIMDVector;
	
#else
	#error "Please specificy a data size!"
#endif

union sse_elem {
	
	sse_elem()
	{
		// Do nothing
	};
	
	sse_elem(const SO_Score &m_scalar)
	{
		sse = _MM_SET1(m_scalar);
	};
	
	sse_elem(const SIMDVector &m_vector)
	{
		sse = m_vector;
	};
	
	~sse_elem()
	{
		// Do nothing
	};
	
	SIMDVector sse;
	SO_Score v[SO_LEN];
};

#define		MINIMUM_ALLOWED_SCORE	-9999

namespace NA {
	
	// Allowed nucleic acid bases
	//	A = adenine
	//	C = cytosine
	//	G = guanine
	//	T = thymine
	//	M = A or C
	//	R = G or A
	//	S = G or C
	//	V = G or C or A
	//	W = A or T
	//	Y = T or C
	//	H = A or C or T
	//	K = G or T
	//	D = G or A or T
	//	B = G or T or C
	//	N = A or T or C or G
	typedef enum {
		A = (1 << 0), 
		C = (1 << 1), 
		G = (1 << 2), 
		T = (1 << 3),
		GAP = (1 << 4),
		M = (A | C),
		R = (G | A),
		S = (G | C),
		V = (G | C | A),
		W = (A | T),
		Y = (T | C),
		H = (A | C | T),
		K = (G | T),
		D = (G | A | T),
		B = (G | T | C),
		N = (A | T | C | G)
	} nucleic_acid;
};

// Both base types (nucleic acid and amino acid) must fit
// into a variable of size base_type
typedef unsigned char base_type;

base_type na_to_bits(const char &m_base);

class SeqOverlap{

	public:
	
	typedef enum {SmithWaterman, Overlap} AlignmentMode;
	
	private:
		
		struct SO_Elem{
			
			// Use the notation of "Biological sequence analysis" by Durbin, Eddy, Krogh and Mitchison
			sse_elem M;
			
			#ifdef ALLOW_GAPS
			sse_elem I_query;  // insertion in query
			sse_elem I_target; // insertion in target
			#endif // ALLOW_GAPS
									
			// The actual alignment is not needed, just the score and the
			// coordinates of the begining and end of the alignment. For this purpose, propagate
			// the alignment start during the dynamic programming.
			
			// The starting row and column of the alignment
			sse_elem M_start_i;
			
			#ifdef INCLUDE_TARGET_RANGE
			sse_elem M_start_j;
			#endif // INCLUDE_TARGET_RANGE
			
			#ifdef ALLOW_GAPS
			sse_elem I_query_start_i;
			sse_elem I_target_start_i;
			
			#ifdef INCLUDE_TARGET_RANGE
			sse_elem I_query_start_j;
			sse_elem I_target_start_j;
			#endif // INCLUDE_TARGET_RANGE
			
			#endif // ALLOW_GAPS
		};
		
		AlignmentMode mode;
		
		SO_Elem *last_row;
		SO_Elem *curr_row;
		SO_Elem max_elem;
		size_t dp_num_elem;
		
		// The location of the max element in the *query*
		sse_elem stop_i;
		
		#ifdef INCLUDE_TARGET_RANGE
		// The location of the max element in the *target*
		sse_elem stop_j;
		#endif // INCLUDE_TARGET_RANGE
		
		sse_elem match;
		
		#ifdef TREAT_N_AS_MASK
		sse_elem mask;
		#endif // TREAT_N_AS_MASK
		
		sse_elem mismatch;
		
		#ifdef ALLOW_GAPS
		sse_elem gap_existance;
		sse_elem gap_extension;
		#endif // ALLOW_GAPS
		
		// The input sequences in 5'-3' orientation
		sse_elem *query;
		sse_elem *target;
		
		sse_elem query_len;
		sse_elem target_len;
		
		// The maximum query and target length for the
		// currently packed sequences.
		SO_Score max_query_len;
		SO_Score max_target_len;
		
		// The current number of elements allocated to the query and target arrays.
		// We only reallocate when we need to *grow* these arrays (so these values may
		// be different from max_query_len and max_target_len, which dynamically grow
		// and shrink as needed).
		unsigned int query_buffer_len;
		unsigned int target_buffer_len;
		
		#ifdef TREAT_N_AS_MASK
		// If mask_N_na == true, then treat 'N' as a masking character that gets a score of mask
		bool mask_N_na;
		#endif // TREAT_N_AS_MASK
		
		void reverse_complement(std::vector<base_type> &m_seq) const
		{			
			size_t len = m_seq.size();

			len = len/2 + (len%2);

			std::vector<base_type>::iterator f = m_seq.begin();
			std::vector<base_type>::reverse_iterator r = m_seq.rbegin();

			for(size_t i = 0;i < len;i++, f++, r++){
				
				const base_type tmp = complement(*f);
				*f = complement(*r);
				*r = tmp;
			}
		};

		inline base_type complement(const base_type &m_base) const
		{
			switch(m_base){
				case NA::A:
					return NA::T;
				case NA::C:
					return NA::G;
				case NA::G:
					return NA::C;
				case NA::T:
					return NA::A;
				case NA::M:
					return NA::K;
				case NA::R:
					return NA::Y;
				case NA::S:
					return NA::S;
				case NA::V:
					return NA::B;
				case NA::W:
					return NA::W;
				case NA::Y:
					return NA::R;
				case NA::H:
					return NA::D;
				case NA::K:
					return NA::M;
				case NA::D:
					return NA::H;
				case NA::B:
					return NA::V;
				case NA::N:
					return NA::N;
				case NA::GAP:
					return NA::GAP;
			};

			throw __FILE__ ":complement: Unknown base";
			return NA::GAP; // Keep the compiler happy
		};

		void align_overlap();
		void align_smith_waterman();
		
		SO_Score trace_back();

	public:
		
		SeqOverlap(const AlignmentMode &m_mode);
			
		~SeqOverlap()
		{
			if(query != NULL){
				_mm_free(query);
			}
			
			if(target != NULL){
				_mm_free(target);
			}
			
			if(last_row != NULL){
				_mm_free(last_row);
			}
			
			if(curr_row != NULL){
				_mm_free(curr_row);
			}
		};

		void align()
		{
			switch(mode){
				case Overlap:
					align_overlap();
					break;
				case SmithWaterman:
					align_smith_waterman();
					break;
				default:
					throw __FILE__ ":align: Unknown mode";
			};
		}
		
		inline void set_mode(const AlignmentMode &m_mode)
		{
			mode = m_mode;
		};
		
		#ifdef TREAT_N_AS_MASK
		inline void enable_mask_N_na()
		{
			mask_N_na = true;
		};
		
		inline void disable_mask_N_na()
		{
			mask_N_na = false;
		};
		#endif // TREAT_N_AS_MASK
		
		// Pack the supplied query into slot m_index
		inline void pack_query(const unsigned int &m_index, const std::string &m_query)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":pack_query: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			const unsigned int len = m_query.size();

			if( ( len > query_buffer_len ) || (query == NULL) ){
				
				// We need to allocate a new query buffer
				sse_elem *temp = (sse_elem*)_mm_malloc(len*sizeof(sse_elem), SSE_ALIGNMENT_SIZE);
				
				if(query != NULL){
					
					// Copy the contents of the current query buffer into the new query buffer
					memcpy( temp, query, query_buffer_len*sizeof(sse_elem) );
					_mm_free(query);
				}
				
				query_buffer_len = len;
				
				query = temp;
				
			}
			
			for(unsigned int i = 0;i < len;i++){
				query[i].v[m_index] = na_to_bits(m_query[i]);
			}

			query_len.v[m_index] = len;
			
			// Recompute the maximum query length in case m_index *used* to store the largest sequence, 
			// in which case we will need to shrink max_query_len.
			max_query_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){			
				max_query_len = std::max(max_query_len, query_len.v[i]);
			}
		};
		
		// Pack the supplied query into *all* slots
		inline void pack_query(const std::string &m_query)
		{
			const unsigned int len = m_query.size();

			if( ( len > query_buffer_len ) || (query == NULL) ){
				
				// We need to allocate a new query buffer
				sse_elem *temp = (sse_elem*)_mm_malloc(len*sizeof(sse_elem), SSE_ALIGNMENT_SIZE);
				
				if(query != NULL){
					
					// Copy the contents of the current query buffer into the new query buffer
					memcpy( temp, query, query_buffer_len*sizeof(sse_elem) );
					_mm_free(query);
				}
				
				query_buffer_len = len;
				
				query = temp;
				
			}
			
			// Update the maximum query length
			max_query_len = len;
			
			for(unsigned int i = 0;i < len;i++){
				query[i].sse = _MM_SET1( na_to_bits(m_query[i]) );
			}
			
			query_len.sse = _MM_SET1(len);
		};
		
		// Clear all queries in all slots
		inline void clear_query()
		{
			// Reset the maximum query length
			max_query_len = 0;
			
			query_len.sse = _MM_SET1(0);
		};
		
		// Clear the query in slot m_index
		inline void clear_query(const unsigned int &m_index)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":clear_query: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			query_len.v[m_index] = 0;
			
			// Recompute the maximum query length
			max_query_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){			
				max_query_len = std::max(max_query_len, query_len.v[i]);
			}
		};
		
		// Pack the supplied target into slot m_index
		inline void pack_target(const unsigned int &m_index, const std::string &m_target)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":pack_target: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			const unsigned int len = m_target.size();
			
			if( ( len > target_buffer_len ) || (target == NULL) ){
				
				// We need to allocate a new query buffer
				sse_elem *temp = (sse_elem*)_mm_malloc(len*sizeof(sse_elem), SSE_ALIGNMENT_SIZE);
				
				if(target != NULL){
					
					// Copy the contents of the current query buffer into the new query buffer
					memcpy( temp, target, target_buffer_len*sizeof(sse_elem) );
					_mm_free(target);
				}
				
				target_buffer_len = len;
				
				target = temp;
				
			}
			
			for(unsigned int i = 0;i < len;i++){
				target[i].v[m_index] = na_to_bits(m_target[i]);
			}

			target_len.v[m_index] = len;
			
			// Recompute the maximum target length in case m_index *used* to store the largest sequence, 
			// in which case we will need to shrink max_target_len.
			max_target_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){			
				max_target_len = std::max(max_target_len, target_len.v[i]);
			}
		};
		
		// Pack the supplied target into *all* slots
		inline void pack_target(const std::string &m_target)
		{
			const unsigned int len = m_target.size();

			if( ( len > target_buffer_len ) || (target == NULL) ){
				
				// We need to allocate a new query buffer
				sse_elem *temp = (sse_elem*)_mm_malloc(len*sizeof(sse_elem), SSE_ALIGNMENT_SIZE);
				
				if(target != NULL){
					
					// Copy the contents of the current query buffer into the new query buffer
					memcpy( temp, target, target_buffer_len*sizeof(sse_elem) );
					_mm_free(target);
				}
				
				target_buffer_len = len;
				
				target = temp;
				
			}
			
			// Update the maximum target length
			max_target_len = len;
			
			for(unsigned int i = 0;i < len;i++){
				target[i].sse = _MM_SET1( na_to_bits(m_target[i]) );
			}
			
			target_len.sse = _MM_SET1(len);
		};
		
		// Clear all targets in all slots
		inline void clear_target()
		{
			
			// Reset the maximum target length
			max_target_len = 0;
			
			target_len.sse = _MM_SET1(0);
		};
		
		// Clear the target in slot m_index
		inline void clear_target(const unsigned int &m_index)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":clear_target: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			target_len.v[m_index] = 0;
			
			// Recompute the maximum target length
			max_target_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){
				max_target_len = std::max(max_target_len, target_len.v[i]);
			}
		};
					
		// The coordinates of the first and last aligned base in the query
		inline std::pair<int, int> alignment_range_query(const unsigned int &m_index) const
		{
			
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":alignment_range_query: Index out of bounds";
			}
			#endif // _DEBUG
			
			return std::make_pair(max_elem.M_start_i.v[m_index], stop_i.v[m_index]);
		};
		
		#ifdef INCLUDE_TARGET_RANGE
		// The coordinates of the first and last aligned base in the query
		inline std::pair<int, int> alignment_range_target(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":alignment_range_target: Index out of bounds";
			}
			#endif // _DEBUG
			
			return std::make_pair(max_elem.M_start_j.v[m_index], stop_j.v[m_index]);
		};
		
		inline void alignment_range(std::pair<int, int> &m_query_range,
			std::pair<int, int> &m_target_range, const unsigned int &m_index) const
		{
		
			m_query_range = alignment_range_query(m_index);
			m_target_range = alignment_range_target(m_index);
		};
		#endif // INCLUDE_TARGET_RANGE
		
		inline SO_Score score(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":score: Index out of bounds";
			}
			#endif // _DEBUG
			
			return max_elem.M.v[m_index];
		};		
};

} // namespace::SO

#endif // __SEQ_OVERLAP
