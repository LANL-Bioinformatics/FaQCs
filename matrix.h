#ifndef __MATRIX
#define	__MATRIX

#include <vector>

template <class T>
class matrix : private std::vector<T>
{
	private:
		size_t num_row;
		size_t num_col;
	public:
		
		matrix()
		{
			num_row = 0;
			num_col = 0;
		};
		
		matrix(const size_t &m_num_row, const size_t &m_num_col)
		{
			matrix::resize(m_num_row, m_num_col);
		};
		
		void resize(const size_t &m_num_row, const size_t &m_num_col)
		{
			num_row = m_num_row;
			num_col = m_num_col;
			
			std::vector<T>::resize(num_row*num_col);
		};
		
		typedef typename std::vector<T>::iterator iterator;
		typedef typename std::vector<T>::const_iterator const_iterator;

		inline iterator begin()
		{
			return std::vector<T>::begin();
		};
		
		inline const_iterator begin() const
		{
			return std::vector<T>::begin();
		};
		
		inline iterator end()
		{
			return std::vector<T>::end();
		};

		inline const_iterator end() const
		{
			return std::vector<T>::end();
		};

		inline T& operator()(const size_t &m_r, const size_t &m_c)
		{
			return (*this)[m_r*num_col + m_c];
		};
		
		inline const T& operator()(const size_t &m_r, const size_t &m_c) const
		{
			return (*this)[m_r*num_col + m_c];
		};
		
		inline size_t get_num_row() const
		{
			return num_row;
		};
		
		inline size_t get_num_col() const
		{
			return num_col;
		};
		
		// Return the selected row
		inline std::vector<T> row(const size_t &m_index) const
		{
			if(m_index >= num_row){
				throw __FILE__ ":matrix::row: index out of bounds";
			}
			
			std::vector<T> ret(num_col);
			
			const size_t offset = m_index*num_col;
			
			for(size_t i = 0;i < num_col;++i){
				ret[i] = (*this)[offset + i];
			}
			
			return ret;
		};
		
		inline std::vector<T> column(const size_t &m_index) const
		{
			if(m_index >= num_col){
				throw __FILE__ ":matrix::column: index out of bounds";
			}
			
			std::vector<T> ret(num_row);
			
			size_t offset = m_index;
			
			for(size_t i = 0;i < num_row;++i, offset += num_col){
				ret[i] = (*this)[offset];
			}
			
			return ret;
		};
		
		matrix& operator+=(const matrix &m_rhs)
		{
			// The right hand side can have fewer rows, but it must
			// have the same number of *columns* (unless the current
			// matrix is empty). The reason for this is that the matricies
			// are packed by column. We rely on the vector::resize function
			// to correctly preserve the old matrix contents.
			if( (num_col == 0) || (num_row == 0) ){ // left hand size is empty
				resize(m_rhs.num_row, m_rhs.num_col);
			}
			else{
				
				// The left hand matrix (i.e. the destination) is allowed to
				// dynamically grow to accomodate additional rows).
				if(num_row < m_rhs.num_row){
					
					if(num_col != m_rhs.num_col){
						throw __FILE__ ":matrix::operator+=: Unequal number of columns!";
					}
				
					resize(m_rhs.num_row, num_col);
				}
			}
			
			const size_t len = m_rhs.size();
			
			for(size_t i = 0;i < len;++i){
				(*this)[i] += m_rhs[i];
			}
			
			return *this;
		};
		
		inline void reserve(const size_t &m_r, const size_t &m_c)
		{
			std::vector<T>::reserve(m_r*m_c);
		};
};

#endif // __MATRIX
