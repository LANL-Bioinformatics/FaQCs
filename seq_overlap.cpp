#include "seq_overlap.h"

using namespace std;
using namespace SO;

SeqOverlap::SeqOverlap(const AlignmentMode &m_mode)
{
	mode = m_mode;
	
	last_row = curr_row = NULL;
	dp_num_elem = 0;
	
	query = target = NULL;
	
	max_query_len = max_target_len = 0;
	query_buffer_len = target_buffer_len = 0;
	
	query_len.sse = _MM_SET1(0);
	target_len.sse = _MM_SET1(0);
	
	match.sse = _MM_SET1(1);
	mismatch.sse = _MM_SET1(-1);
	
	#ifdef ALLOW_GAPS
	gap_existance.sse = _MM_SET1(-1000); // Could use blastn value?
	gap_extension.sse = _MM_SET1(-1000);
	#endif // ALLOW_GAPS
	
	#ifdef TREAT_N_AS_MASK
	mask_N_na = false;
	mask.sse = _MM_SET1(0);
	#endif // TREAT_N_AS_MASK
};

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqOverlap::align_overlap()
{
	throw __FILE__ ":SeqOverlap::align_overlap: Implement me!";
}

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqOverlap::align_smith_waterman()
{
	// Add one to allow for an initial row and column of zeros
	const size_t num_elem = max_target_len + 1;

	const sse_elem zero(0);
	
	#ifdef TREAT_N_AS_MASK
	const sse_elem all_N(NA::N);
	#endif // TREAT_N_AS_MASK
	
	const sse_elem m_minus_mm( _MM_SUB(match.sse, mismatch.sse) );
	
	// Resize the dynamic programming matrix if needed
	if( !last_row || !curr_row || (num_elem > dp_num_elem) ){

		if(last_row){
			_mm_free(last_row);
		}
		
		if(curr_row){
			_mm_free(curr_row);
		}
		
		dp_num_elem = num_elem;
		
		last_row = (SO_Elem*)_mm_malloc(dp_num_elem*sizeof(SO_Elem), SSE_ALIGNMENT_SIZE);
		
		if(!last_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for last_row";
		}
		
		curr_row = (SO_Elem*)_mm_malloc(dp_num_elem*sizeof(SO_Elem), SSE_ALIGNMENT_SIZE);
		
		if(!curr_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for curr_row";
		}
	}

	// Initialize the dynamic programming matrix to have zeros along the first row 
	for(SO_Score j = 0;j <= max_target_len;++j){

		SO_Elem &elem_ref = last_row[j];

		elem_ref.M.sse = zero.sse;
		
		#ifdef ALLOW_GAPS
		elem_ref.I_query.sse = elem_ref.I_target.sse = gap_existance.sse;
		#endif // ALLOW_GAPS
		
		elem_ref.M_start_i.sse = zero.sse;
		
		#ifdef INCLUDE_TARGET_RANGE
		elem_ref.M_start_j.sse = _MM_SET1(j);
		#endif // INCLUDE_TARGET_RANGE
	}
	
	// Reset the maximum element
	max_elem.M.sse = _MM_SET1(0);
	
	for(SO_Score i = 0;i < max_query_len;++i){
		
		// Initialize the dynamic programming matrix to have zeros along the first column
		SO_Elem &elem_ref = curr_row[0];

		elem_ref.M.sse = zero.sse;
		
		#ifdef ALLOW_GAPS
		elem_ref.I_query.sse = elem_ref.I_target.sse = gap_existance.sse;
		#endif // ALLOW_GAPS
		
		elem_ref.M_start_i.sse = _MM_SET1(i + 1);
		
		#ifdef INCLUDE_TARGET_RANGE
		elem_ref.M_start_j.sse = zero.sse;
		#endif // INCLUDE_TARGET_RANGE
		
		// The dp matrix has max_query_len rows and max_target_len columns
		// A B
		// C X <-- dp[i][j]
		SO_Elem *A_ptr = last_row;
		
		#ifdef ALLOW_GAPS
		SO_Elem *B_ptr = A_ptr + 1;
		SO_Elem *C_ptr = curr_row;
		#endif // ALLOW_GAPS
		
		SO_Elem *X_ptr = curr_row + 1;
		
		const sse_elem q = query[i];
		
		const sse_elem all_i(i);
		const sse_elem valid_query( _MM_CMPLT(all_i.sse, query_len.sse) );
		
		#ifdef TREAT_N_AS_MASK
		const sse_elem query_is_N( _MM_CMPEQ(q.sse, all_N.sse) );
		#endif // TREAT_N_AS_MASK
		
		#ifdef ALLOW_GAPS
		for(SO_Score j = 0;j < max_target_len;++j, ++A_ptr, ++B_ptr, ++C_ptr, ++X_ptr){
		#else
		for(SO_Score j = 0;j < max_target_len;++j, ++A_ptr, ++X_ptr){
		#endif // ALLOW_GAPS
		
			sse_elem tmp_a;
			sse_elem tmp_b;
			sse_elem tmp_c;
			
			const sse_elem all_j(j);
			const sse_elem t = target[j];
			
			tmp_a.sse = _MM_ADD(
				_mm_and_si128(
					_MM_CMPGT(_mm_and_si128(q.sse, t.sse), zero.sse),
					m_minus_mm.sse),
				mismatch.sse);
			
			#ifdef TREAT_N_AS_MASK
			if(mask_N_na){

				tmp_b.sse = _mm_or_si128(query_is_N.sse, _MM_CMPEQ(t.sse, all_N.sse) );

				tmp_a.sse = _mm_or_si128( 
						_mm_and_si128(tmp_b.sse, mask.sse), 
						_mm_andnot_si128(tmp_b.sse, tmp_a.sse) );
			}
			#endif // TREAT_N_AS_MASK
			
			#ifdef ALLOW_GAPS
			
			tmp_c.sse = _MM_MAX( _MM_MAX(A_ptr->M.sse, A_ptr->I_query.sse), A_ptr->I_target.sse );
			
			X_ptr->M.sse = _MM_ADD(
				// Unlike the Overlap alignment, clamp the max value at zero
				_MM_MAX(tmp_c.sse, zero.sse),
				tmp_a.sse);
				
			#else
			
			X_ptr->M.sse = _MM_ADD(
				// Unlike the Overlap alignment, clamp the max value at zero
				_MM_MAX(A_ptr->M.sse, zero.sse),
				tmp_a.sse);
				
			#endif // ALLOW_GAPS
			
			#ifdef ALLOW_GAPS
			//////////////////////////////////////////////////////////////////////////////
			// tmp_a == (A->M < A->I_query) || (A->M < A->I_target)
			tmp_a.sse = _mm_or_si128(
					_MM_CMPLT(A_ptr->M.sse, A_ptr->I_query.sse), 
					_MM_CMPLT(A_ptr->M.sse, A_ptr->I_target.sse) );
			
			// Update M_start if !tmp_a, which means update if A->M is greater than both
			// A->I_query and A->I_target
			X_ptr->M_start_i.sse = _mm_andnot_si128(tmp_a.sse, A_ptr->M_start_i.sse);
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->M_start_j.sse = _mm_andnot_si128(tmp_a.sse, A_ptr->M_start_j.sse);
			#endif // INCLUDE_TARGET_RANGE
					
			//////////////////////////////////////////////////////////////////////////////
			// tmp_a = !(A->I_query < A->I_target) && (A->I_query > A->M)
			//       =  (A->I_query >= A->I_target) && (A->I_query > A->M)
			tmp_a.sse = _mm_andnot_si128(
					_MM_CMPLT(A_ptr->I_query.sse, A_ptr->I_target.sse),
					_MM_CMPGT(A_ptr->I_query.sse, A_ptr->M.sse) );
			
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_query_start_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_query_start_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_and_si128(
					_MM_CMPGT(A_ptr->I_target.sse, A_ptr->M.sse),
					_MM_CMPGT(A_ptr->I_target.sse, A_ptr->I_query.sse) );
			
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_target_start_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_target_start_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
			
			//#else
			// Setting X_ptr->M_start_i = A_ptr->M_start_i and X_ptr->M_start_j = A_ptr->M_start_j is
			// now performed below.
			//X_ptr->M_start_i.sse = A_ptr->M_start_i.sse;
			//X_ptr->M_start_j.sse = A_ptr->M_start_j.sse;
			#endif // ALLOW_GAPS
					
			//////////////////////////////////////////////////////////////////////////////
			// Check for values that were clamped to zero
			#ifdef ALLOW_GAPS
			tmp_a.sse = _MM_CMPGT(zero.sse, tmp_c.sse);
			#else
			tmp_a.sse = _MM_CMPGT(zero.sse, A_ptr->M.sse);
			#endif // ALLOW_GAPS
			
			#ifdef ALLOW_GAPS
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, all_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, all_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
			
			#else
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, A_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, all_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, A_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, all_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
			
			#endif // ALLOW_GAPS
			
			#ifdef ALLOW_GAPS
			//////////////////////////////////////////////////////////////////////////////
			// Unlike the Overlap alignment, clamp the C cell values to zero
			tmp_b.sse = _MM_ADD(_MM_MAX(C_ptr->M.sse, zero.sse), gap_existance.sse);
			tmp_c.sse = _MM_ADD(_MM_MAX(C_ptr->I_query.sse, zero.sse), gap_extension.sse);
			
			X_ptr->I_query.sse = _MM_MAX(tmp_b.sse, tmp_c.sse );

			tmp_a.sse = _MM_CMPLT(tmp_b.sse, tmp_c.sse);
			
			X_ptr->I_query_start_i.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, C_ptr->M_start_i.sse),
							_mm_and_si128(
								tmp_a.sse, C_ptr->I_query_start_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->I_query_start_j.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, C_ptr->M_start_j.sse),
							_mm_and_si128(
								tmp_a.sse, C_ptr->I_query_start_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
			
			//////////////////////////////////////////////////////////////////////////////
			// Unlike the Overlap alignment, clamp the B cell values to zero
			tmp_b.sse = _MM_ADD(_MM_MAX(B_ptr->M.sse, zero.sse), gap_existance.sse);
			tmp_c.sse = _MM_ADD(_MM_MAX(B_ptr->I_target.sse, zero.sse), gap_extension.sse);
								
			X_ptr->I_target.sse = _MM_MAX(tmp_b.sse, tmp_c.sse);
			
			tmp_a.sse = _MM_CMPLT(tmp_b.sse, tmp_c.sse);
			
			X_ptr->I_target_start_i.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, B_ptr->M_start_i.sse),
							_mm_and_si128(
								tmp_a.sse, B_ptr->I_target_start_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE
			X_ptr->I_target_start_j.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, B_ptr->M_start_j.sse),
							_mm_and_si128(
								tmp_a.sse, B_ptr->I_target_start_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
			
			#endif // ALLOW_GAPS
								
			//////////////////////////////////////////////////////////////////////////////
			// Unlike the Overlap alignment, find the maximum scoring element *anywhere*
			// in the dynamic programming matrix
			const sse_elem valid_query_and_target( 
				_mm_and_si128(valid_query.sse, _MM_CMPLT(all_j.sse, target_len.sse) ) );
			
			// tmp_a is true if this element is a new maximum M value
			tmp_a.sse = _mm_andnot_si128( _MM_CMPLT(X_ptr->M.sse, max_elem.M.sse), valid_query_and_target.sse );
			
			max_elem.M.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, X_ptr->M.sse),
						_mm_andnot_si128(tmp_a.sse, max_elem.M.sse) );
			
			max_elem.M_start_i.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, X_ptr->M_start_i.sse),
						_mm_andnot_si128(tmp_a.sse, max_elem.M_start_i.sse) );
			
			stop_i.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, all_i.sse),
						_mm_andnot_si128(tmp_a.sse, stop_i.sse) );
			
			#ifdef INCLUDE_TARGET_RANGE	
			max_elem.M_start_j.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, X_ptr->M_start_j.sse),
						_mm_andnot_si128(tmp_a.sse, max_elem.M_start_j.sse) );
			
			stop_j.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, all_j.sse),
						_mm_andnot_si128(tmp_a.sse, stop_j.sse) );
			#endif // INCLUDE_TARGET_RANGE
		}
		
		// Swap the last_row and curr_row
		swap(last_row, curr_row);
	}
}

base_type SO::na_to_bits(const char &m_base)
{
	switch(m_base){
		case 'A': case 'a':
			return NA::A;
		case 'C': case 'c':
			return NA::C;
		case 'G': case 'g':
			return NA::G;
		case 'T': case 't':
			return NA::T;
		case 'M': case 'm':
			return NA::M;
		case 'R': case 'r':
			return NA::R;
		case 'S': case 's':
			return NA::S;
		case 'V': case 'v':
			return NA::V;
		case 'W': case 'w': 
			return NA::W;
		case 'Y': case  'y':
			return NA::Y;
		case 'H': case 'h':
			return NA::H;
		case 'K': case 'k':
			return NA::K;
		case 'D': case 'd':
			return NA::D;
		case 'B': case 'b':
			return NA::B;
		case 'N': case 'n':
			return NA::N;
		case '-':
			return NA::GAP;
	};

	throw __FILE__ ":na_to_bits: Unknown base!";
	return NA::GAP; // Keep the compiler happy
};

