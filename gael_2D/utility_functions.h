#pragma once
#include <functional>
#include <fstream>
#include <time.h>       /* time */
#include <iostream>
#include <string>
#include <sys/stat.h> // For mkdir on POSIX systems

#include "structures.h"

#ifdef _WIN32 // MSC_VER
#define WIN32_LEAN_AND_MEAN
// *sigh* no gettimeofday on Win32/Win64
 inline int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
	 //
	 // TO BE CONTINUED ......
	 //
	 //
	return 0;
}
#else
#include <sys/time.h>
#endif // _WIN32

inline struct def_grad combiner(double a, double b, std::function<def_grad(double, double)> func)
{
	return func(a, b);
}

// functions from main
inline bool fexists(const char *fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

template <typename T, typename U, typename V, typename W>
bool allequal(const T &t, const U &u, const V &v, const W &w)
{
	return (t == u) || (t == v) || (t == w);
}

inline double get_wall_time()
{
	struct timeval time;
	if (gettimeofday(&time, NULL))
	{
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

inline double get_cpu_time()
{
	return (double)clock() / CLOCKS_PER_SEC;
}

std::string get_outputPaht(){
	return "output/"; 
}

std::ofstream openFileInOutputDirectory(const std::string& filename) {
    const std::string outputPath = get_outputPaht();

    struct stat info;
	if (stat(outputPath.c_str(), &info) == 0 && S_ISDIR(info.st_mode))
	{}else{
    // Check if the "output/" directory exists; if not, create it
    #ifdef _WIN32
        if (_mkdir(outputPath.c_str()) != 0) {
    #else
        if (mkdir(outputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
    #endif
            // Handle directory creation error if needed
            throw std::runtime_error("Failed to create the 'output/' directory.");
        }
	}
    // Concatenate the output path and filename
    std::string fullpath = outputPath + filename;

    // Open the file
    std::ofstream file(fullpath);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << fullpath << std::endl;
        throw std::runtime_error("Failed to open the file: " + fullpath);
    }

    // File is successfully opened; return the file object
    return file;
}