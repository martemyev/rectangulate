#include "parameters.hpp"
#include "utilities.hpp"

#include <vector>
#include <algorithm>



Parameters::Parameters()
  : mesh_filename(DEFAULT_FILE_NAME),
    properties_filename(DEFAULT_FILE_NAME),
    n_rect_elements_x(0),
    n_rect_elements_z(0),
    _parameters(),
    _longest_string_key_len(DEFAULT_PRINT_LEN),
    _longest_string_value_len(DEFAULT_PRINT_LEN)
{
  int p = 0;

  add_option("-meshfile", new OneParam<std::string>("name of file with triangular mesh", &mesh_filename, ++p));
  add_option("-propfile", new OneParam<std::string>("name of file with media properties", &properties_filename, ++p));
  add_option("-nnx", new OneParam<int>("number of rectangular elements in x-direction", &n_rect_elements_x, ++p));
  add_option("-nnz", new OneParam<int>("number of rectangular elements in z-direction", &n_rect_elements_z, ++p));
}




void Parameters::read_command_line(int argc, char **argv)
{
  // Check for help
  if (argc == 1                     ||
      argcheck(argc, argv, "-h")    ||
      argcheck(argc, argv, "-help") ||
      argcheck(argc, argv, "--help"))
  {
    show_options();
    throw 1; // exit with some necessary finishing procedures
  }

  // all the command line entries
  std::vector<std::string> arguments(argv, argv+argc);

  // We start from the 1-st (not 0-th, because arguments[0] is the path to the
  // executable file of this program). Then we consider every second parameter,
  // since the command line goes like this: param0 value0 param1 value1 ...
  // and we need to consider only parameters.
  require((argc-1) % 2 == 0, "The number of command line arguments must be even"
          ", because every parameter is accompanied by a value. But there are "
          "only " + d2s(argc-1) + " of the arguments");
  for (size_t ar = 1; ar < arguments.size(); ar += 2)
  {
    ParaMap::const_iterator iter = _parameters.find(arguments[ar]);
    require(iter != _parameters.end(), "Command line argument '" + arguments[ar]
            + "' wasn't found");
    require(ar+1 < arguments.size(), "Command line argument '" + arguments[ar]
            + "' doesn't have any value");
    iter->second->read(arguments[ar+1]);
  }

  print_parameters(); // show what we work with
  check_parameters(); // check correctness of the input data, if possible
}




void Parameters::show_options() const
{
  update_longest_string_key_len();

  std::cout << "\nAvailable options [default values in brackets]\n\n";

  typedef std::vector<std::pair<std::string, ParamBase*> > ParaVec;

  // sort the map of parameters according to their priority values
  ParaVec sorted_parameters(_parameters.begin(), _parameters.end());
  std::sort(sorted_parameters.begin(),
            sorted_parameters.end(),
            compare_by_parameter_priority);

  ParaVec::const_iterator iter = sorted_parameters.begin();

  for (; iter != sorted_parameters.end(); ++iter)
  {
    const ParamBase *par = iter->second;

    std::cout << add_space(iter->first, _longest_string_key_len + SPACE_BETWEEN)
              << par->get_description()
              << " [" << par->str() << "]\n";
  }

  std::cout << "\n";

  std::cout << add_space("-help",  _longest_string_key_len + SPACE_BETWEEN)
            << "print this menu and exit\n";
  std::cout << "\n";
}




void Parameters::print_parameters() const
{
  update_longest_string_key_len();
  update_longest_string_value_len();

  typedef std::vector<std::pair<std::string, ParamBase*> > ParaVec;

  // sort the map of parameters according to their priority values
  ParaVec sorted_parameters(_parameters.begin(), _parameters.end());
  std::sort(sorted_parameters.begin(),
            sorted_parameters.end(),
            compare_by_parameter_priority);

  ParaVec::const_iterator iter = sorted_parameters.begin();

  for (; iter != sorted_parameters.end(); ++iter)
  {
    const ParamBase *par = iter->second;
    std::cout << add_space(iter->first, _longest_string_key_len + SPACE_BETWEEN)
              << add_space(par->str(), _longest_string_value_len + SPACE_BETWEEN)
              << par->get_description() << "\n";
  }
  std::cout << "\n";
}




void Parameters::check_parameters() const
{
  require(file_exists(mesh_filename), "File '" + mesh_filename + "' doesn't "
          "exist");
  require(file_exists(properties_filename), "File '" + properties_filename +
          "' doesn't exist");

  require(n_rect_elements_x > 0, "nnx parameter is wrong: " +
          d2s(n_rect_elements_x));
  require(n_rect_elements_z > 0, "nnz parameter is wrong: " +
          d2s(n_rect_elements_z));
}




void Parameters::add_option(const std::string &key, ParamBase *parameter)
{
  std::pair<ParaMap::iterator, bool> result =
      _parameters.insert(std::pair<std::string, ParamBase*>(key, parameter));

  require(result.second, "The insertion of the option with the key word '" +
          key + "' failed.");
}




void Parameters::update_longest_string_key_len() const
{
  _longest_string_key_len = 0;

  ParaMap::const_iterator iter = _parameters.begin();

  for (; iter != _parameters.end(); ++iter)
  {
    const int len_key_string = iter->first.size();

    if (len_key_string > _longest_string_key_len)
      _longest_string_key_len = len_key_string;
  }
}




void Parameters::update_longest_string_value_len() const
{
  _longest_string_value_len = 0;

  ParaMap::const_iterator iter = _parameters.begin();

  for (; iter != _parameters.end(); ++iter)
  {
    const ParamBase *par = iter->second;
    const int len_value_string = par->str().size();

    if (len_value_string > _longest_string_value_len)
      _longest_string_value_len = len_value_string;
  }
}




//==============================================================================
//
// Compare two pairs containing info about parameters by the priority of the
// parameters
//
//==============================================================================
bool compare_by_parameter_priority(const std::pair<std::string, ParamBase*> &a,
                                   const std::pair<std::string, ParamBase*> &b)
{
  return (a.second->get_priority() < b.second->get_priority());
}
