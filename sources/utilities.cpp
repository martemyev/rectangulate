#include "utilities.hpp"

#if defined(__linux__) || defined(__APPLE__)
  #include <sys/time.h> // for time measurements
#endif

#include <execinfo.h>
#include <cstdlib>
#include <fstream>
#include <climits>
#include <cerrno>
#include <cstring>

//------------------------------------------------------------------------------
//
// expect and require
//
//------------------------------------------------------------------------------
void requirement_fails(const char *file,
                       unsigned int line,
                       std::string message)
{
  std::string exc = "Exception:\nfile = " + std::string(file) +
                    "\nline = " + d2s(line) +
                    "\nmessage = " + message + "\n";

#if defined(__linux__)
  const int backtraceSize = 20;
  void *array[backtraceSize];
  int size = backtrace(array, backtraceSize);
  char **strings = backtrace_symbols(array, size);

  exc += "backtrace:\nsize = " + d2s<int>(size) + "\n";
  for (int i = 0; i < size; ++i)
    exc += std::string(strings[i]) + "\n";

  free(strings);
#endif

  throw std::runtime_error(exc);
}

//------------------------------------------------------------------------------
//
// Time measurement (wall time)
//
//------------------------------------------------------------------------------
double get_wall_time()
{
#if defined(__linux__) || defined(__APPLE__)
  struct timeval time;
  const int ierr = gettimeofday(&time, NULL);
  require(ierr == 0, "gettimeofday returned an error code " + d2s(ierr));

  return (1.0*time.tv_sec + 1.0e-6*time.tv_usec);

#elif(_WIN32)
  require(false, "Not implemented for Windows");
  return 0.;
#else
  require(false, "Not implemented for this unknown OS");
  return 0.;
#endif
}

//------------------------------------------------------------------------------
//
// Name of a file without a path
//
//------------------------------------------------------------------------------
std::string file_name(const std::string &path)
{
  if (path == "") return path;

  size_t pos = 0;
#if defined(__linux__) || defined(__APPLE__)
  pos = path.find_last_of('/');
#elif defined(_WIN32)
  pos = path.find_last_of('\\');
#endif

  if (pos == std::string::npos)
    return path; // there is no '/' in the path, so this is the filename

  return path.substr(pos + 1);
}

//------------------------------------------------------------------------------
//
// Path of a given file
//
//------------------------------------------------------------------------------
std::string file_path(const std::string &path)
{
  if (path == "") return path;

  size_t pos = 0;
#if defined(__linux__) || defined(__APPLE__)
  pos = path.find_last_of('/');
#elif defined(_WIN32)
  pos = path.find_last_of('\\');
#endif

  if (pos == std::string::npos)
    return ""; // there is no '/' in the path, the path is "" then

  return path.substr(0, pos + 1);
}

//------------------------------------------------------------------------------
//
// Stem of a given file (no path, no extension)
//
//------------------------------------------------------------------------------
std::string file_stem(const std::string &path)
{
  if (path == "") return path;

  // get a file name from the path
  const std::string fname = file_name(path);

  // extract a stem and return it
  size_t pos = fname.find_last_of('.');
  if (pos == std::string::npos)
    return fname; // there is no '.', so this is the stem

  return fname.substr(0, pos);
}

//------------------------------------------------------------------------------
//
// Extension of a given file
//
//------------------------------------------------------------------------------
std::string file_extension(const std::string &path)
{
  if (path == "") return path;

  // extract a file name from the path
  const std::string fname = file_name(path);

  size_t pos = fname.find_last_of('.');
  if (pos == std::string::npos)
    return ""; // there is no '.', so there is no extension

  // extract an extension and return it
  return fname.substr(pos);
}

//------------------------------------------------------------------------------
//
// Check if the given file exists
//
//------------------------------------------------------------------------------
bool file_exists(const std::string &path)
{
  if (path == "") return false; // no file - no existance

  // This is not the fastest method, but it should work on all operating
  // systems. Some people also not that this method check 'availibility' of the
  // file, not its 'existance'. But that's what we actually need. If a file
  // exists, but it's not available (even for reading), we believe, that the
  // file doesn't exist.
  bool exists = false;
  std::ifstream in(path.c_str());
  if (in.good())
    exists = true; // file exists and is in a good state
  in.close();

  return exists;
}

//------------------------------------------------------------------------------
//
// Get an absolute path according to the given relative one
//
//------------------------------------------------------------------------------
std::string absolute_path(const std::string &rel_path)
{
#if defined(__linux__) || defined(__APPLE__)
  char abs_path[PATH_MAX];
  char *res = realpath(rel_path.c_str(), abs_path);
  require(res != NULL, "The function realpath() failed with the input "
          "(relative path) = '" + rel_path + "'. errno is set to " + d2s(errno)+
          " which means '" + std::string(strerror(errno)) + "'");
  return std::string(abs_path);
#else
  require(false, "absolute_path() is not implemented for this OS");
#endif
}

//------------------------------------------------------------------------------
//
// Check endianness
//
//------------------------------------------------------------------------------
bool is_big_endian()
{
  union
  {
    int i;
    char c[sizeof(int)];
  } x;
  x.i = 1;
  return x.c[0] == 1;
}

//------------------------------------------------------------------------------
//
// Get the endianness of the machine
//
//------------------------------------------------------------------------------
std::string endianness()
{
  return (is_big_endian() ? "BigEndian" : "LittleEndian");
}

//------------------------------------------------------------------------------
//
// Check if there is a string in the array of strings
//
//------------------------------------------------------------------------------
int argcheck(int argc, char **argv, const char *arg)
{
  for(int i = 1; i < argc; ++i)
  {
    // strcmp returns 0 if the strings are equal
    if(strcmp(argv[i], arg) == 0)
      return(i);
  }

  return 0;
}

//------------------------------------------------------------------------------
//
// Add some space to the given string so it has the requires length
//
//------------------------------------------------------------------------------
std::string add_space(const std::string &str, int length)
{
  // how spaces need to add (if the string longer than the required length,
  // nothing is added)
  const int n_spaces = std::max(length - (int)str.size(), 0);
  return str + std::string(n_spaces, ' ');
}

