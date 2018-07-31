#ifndef __SHAVE
#define __SHAVE

#include <string>
#include <deque>
#include <unordered_map>
#include <unordered_set>

#define 	SHAVE_VERSION 			"1.1"
#define		DEFAULT_MIN_AVE_READ_QUALITY	20.0f
#define		DEFAULT_MIN_BASE_QUALITY	0
#define		DEFAULT_MIN_TERMINAL_QUALITY	5
#define		DEFAULT_MIN_READ_LENGTH		50
#define		DEFAULT_QUALITY_MAP_HEIGHT	500
#define		DEFAULT_QUALITY_MAP_WIDTH	500
#define		DEFAULT_ADAPTOR_THRESHOLD	0.8f

#define		MAP	std::unordered_map
#define		SET	std::unordered_set

struct Coord
{
	unsigned int x;
	unsigned int y;
	unsigned int tile;
	
	Coord() : x(0), y(0), tile(0)
	{
	};
	
	Coord(const unsigned int &m_x, const unsigned int &m_y, const unsigned int &m_tile = 0) : 
		x(m_x), y(m_y), tile(m_tile)
	{
	};
	
	inline bool operator==(const Coord &m_rhs) const
	{
		return (tile == m_rhs.tile) && (x == m_rhs.x) && (y == m_rhs.y);
	};
	
	// Needed for unordered_map -- currently, tile is not used to compute the hash key
	inline size_t hash_value() const
	{
		return ( ( (size_t)x ) << 32 ) | y;
	};
};

namespace std
{
	template<>
	struct hash<Coord>
	{
	    size_t operator()(const Coord& m_xy) const
	    {
       		return m_xy.hash_value();
	    }
	};
}

struct Read
{
	std::string seq;
	std::string qual;
	std::string def;
	
	inline bool empty() const
	{
		return def.empty();
	};
	
	inline void clear()
	{
		seq.clear();
		qual.clear();
		def.clear();
	};
};

struct Options
{
	bool print_usage;
	bool verbose;
	
	std::string read_one_input;
	std::string read_two_input;
	
	std::string read_one_output;
	std::string read_two_output;
	
	bool read_one_protect_5;
	bool read_one_protect_3;
	
	bool read_two_protect_5;
	bool read_two_protect_3;
	
	std::string read_one_min_qual_5;
	std::string read_one_min_qual_3;
	
	std::string read_two_min_qual_5;
	std::string read_two_min_qual_3;
	
	std::deque<std::string> adaptor;
	float adaptor_threshold;
	
	std::string quality_map_output;
	float min_ave_quality;
	int min_base_quality;
	int min_terminal_quality;
	bool pool_tile_quality;
	
	unsigned int quality_map_width;
	unsigned int quality_map_height;
	
	size_t min_read_length;
	
	Options(int argc, char *argv[]);
};

void write_quality_map(const std::string &m_filename, 
	const MAP<Coord, float> &m_q, 
	const SET<Coord> &m_mask, const unsigned int &m_tile, 
	const Options &m_opt);

#endif // __SHAVE
