#include "file_util.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <ctype.h>

using namespace std;

bool file_exits(const string &m_filename)
{
	struct stat file_info;

	if(stat(m_filename.c_str(), &file_info) == -1){
		return false;
	}
	
	return S_ISREG(file_info.st_mode);
}

bool directory_exits(const string &m_dir)
{
	struct stat dir_info;

	if(stat(m_dir.c_str(), &dir_info) == -1){
		return false;
	}
	
	return S_ISDIR(dir_info.st_mode);
}

bool create_directory(const string &m_dir)
{
	return (mkdir(m_dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) == 0);
}
