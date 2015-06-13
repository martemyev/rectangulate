#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include <climits>
#include <sstream>
#include <map>

/**
 * Default file name for the parameters representing file names.
 */
const std::string DEFAULT_FILE_NAME = "no-file";

/**
 * Default length of strings for printing aligned key words and values of the
 * parameters.
 */
const int DEFAULT_PRINT_LEN = 15;

/**
 * To unify the representation of the parameters when they are printed to
 * output stream, they all have this length (space is added to those which
 * are shorter).
 */
const int PARAM_OUT_LENGTH = 20;

/**
 * The space between the name of the option (to define a parameter) and its
 * description.
 */
const int SPACE_BETWEEN = 5;




/**
 * Base abstract class for a one parameter.
 */
class ParamBase
{
public:

  ParamBase()
    : _description(""),
      _priority(0)
  { }

  ParamBase(const std::string &desc, int priority)
    : _description(desc),
      _priority(priority)
  { }

  ParamBase(const ParamBase &pb)
    : _description(pb._description),
      _priority(pb._priority)
  { }

  ParamBase& operator=(const ParamBase &pb)
  {
    _description = pb._description;
    _priority = pb._priority;
    return *this;
  }

  virtual ~ParamBase() { }

  std::string get_description() const { return _description; }

  int get_priority() const { return _priority; }

  /// Read the value of the parameter from a string
  /// @param from a string from which the value is read
  virtual void read(const std::string &from) = 0;

  /// Convert the value of the parameter into a string
  virtual std::string str() const = 0;

protected:

  /// Description of the parameter (appears when help is invoked, for example)
  std::string _description;

  /// This attribute controls the order of outputting the parameters (when help
  /// is invoked, for example). Without this parameter (or, when there are
  /// several parameters with the same priority) the parameters appear in
  /// alphabetical order of keys.
  int _priority;
};




/**
 * Template class for a one parameter.
 */
template <typename T>
class OneParam : public ParamBase
{
public:

  /// Constructor takes a description of the parameter and a pointer to a value.
  /// It also may take a priority value of the parameter, however it's optional
  OneParam(const std::string &desc,
           T* val,
           int priority = INT_MAX)
    : ParamBase(desc, priority),
      _value(val)
  { }

  OneParam(const OneParam &op)
    : ParamBase(op),
      _value(op._value)
  { }

  OneParam& operator=(const OneParam &op)
  {
    ParamBase::operator=(op);
    _value = op._value;
    return *this;
  }

  virtual ~OneParam() { }

  /// Value of the parameter is saved somewhere else, and here we keep the
  /// pointer to it. This member is public, yes.
  T* _value;

  /// Read the value of the parameter from a string
  /// @param a string from which the value is read
  virtual void read(const std::string &from)
  {
    std::istringstream is(from);
    is >> *(_value);
  }

  /// Convert the value of the parameter into a string
  virtual std::string str() const
  {
    std::ostringstream os;
    os << *(_value);
    return os.str();
  }
};




/**
 * This class takes care of the parameters of the current program.
 */
class Parameters
{
public:

  Parameters();

  ~Parameters() { }

  void read_command_line(int argc, char **argv);

  std::string mesh_filename;
  std::string properties_filename;
  int         n_rect_elements_x;
  int         n_rect_elements_z;

private:
  Parameters(const Parameters&);
  Parameters& operator =(const Parameters&);

  void add_option(const std::string &key, ParamBase *parameter);

  void show_options() const;
  void print_parameters() const;
  void check_parameters() const;

  void update_longest_string_key_len() const;
  void update_longest_string_value_len() const;

  /// The map between the key word representing a parameter, and its value (and
  /// maybe other attributes such as description).
  std::map<std::string, ParamBase*> _parameters;

  /// To allow the shorter name for the map.
  typedef std::map<std::string, ParamBase*> ParaMap;

  /// Length of the longest string representing the key words of the parameters.
  mutable int _longest_string_key_len;

  /// Length of the longest string representing the values of the parameters.
  mutable int _longest_string_value_len;
};


//==============================================================================
//
// Compare two pairs containing info about parameters by the priority of the
// parameters
//
//==============================================================================
bool compare_by_parameter_priority(const std::pair<std::string, ParamBase*> &a,
                                   const std::pair<std::string, ParamBase*> &b);

#endif // PARAMETERS_HPP
