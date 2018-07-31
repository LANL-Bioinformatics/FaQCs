#ifndef __FILE_UTIL
#define __FILE_UTIL

#include <string>

bool file_exits(const std::string &m_filename);
bool directory_exits(const std::string &m_dir);
bool create_directory(const std::string &m_dir);

#endif // __FILE_UTIL
