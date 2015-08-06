#include "rectangular_mesh.hpp"
#include "utilities.hpp"
#include "rectangle.hpp"
#include "triangular_mesh.hpp"
#include "triangle.hpp"

#include <fstream>
#include <random>

const int    VOID_MATERIAL_ID       = -1;
const double VOID_MATERIAL_PROPERTY = -1;
typedef float OUT_FLOAT_TYPE;



RectangularMesh::RectangularMesh(const Point2 &beg,
                                 const Point2 &end,
                                 int nx_,
                                 int nz_)
  : _min_coord(beg),
    _max_coord(end),
    _n_elements_x(nx_),
    _n_elements_z(nz_),
    _n_vertices(0),
    _vertices(nullptr),
    _n_elements(0),
    _elements(nullptr)
{
  require(_min_coord < _max_coord, "The max point of the mesh (" +
          d2s(_max_coord) + ") must be bigger than the min point (" +
          d2s(_min_coord) + ")");
  require(_n_elements_x > 0, "The number of elements in x-direction (" +
          d2s(_n_elements_x) + ") must be >0");
  require(_n_elements_z > 0, "The number of elements in z-direction (" +
          d2s(_n_elements_z) + ") must be >0");
}




RectangularMesh::~RectangularMesh()
{
  delete[] _elements;
  delete[] _vertices;
}




void RectangularMesh::build()
{
  std::cout << "Building a rectangular mesh..." << std::endl;
  double t_begin = get_wall_time();

  build_vertices();
  build_elements();

  std::cout << "Building a rectangular mesh is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}




void RectangularMesh::build_vertices()
{
  const double x0 = _min_coord.x(); // limits in x-direction
  const double x1 = _max_coord.x();
  const double z0 = _min_coord.z(); // limits in z-direction
  const double z1 = _max_coord.z();

  require(x1 > x0 && z1 > z0, "Incorrect limits of the rectangular domain");

  const double hx = (x1 - x0) / _n_elements_x; // step in x-direction
  const double hz = (z1 - z0) / _n_elements_z; // step in z-direction

  double *coord_x = new double[_n_elements_x + 1];
  for (int i = 0; i < _n_elements_x; ++i)
    coord_x[i] = x0 + i * hx;
  coord_x[_n_elements_x] = x1;

  double *coord_z = new double[_n_elements_z + 1];
  for (int i = 0; i < _n_elements_z; ++i)
    coord_z[i] = z0 + i * hz;
  coord_z[_n_elements_z] = z1;

  _n_vertices = (_n_elements_x + 1) * (_n_elements_z + 1);
  _vertices = new Point2[_n_vertices];

  int ver = 0; // the number of a current vertex
  for (int i = 0; i < _n_elements_z + 1; ++i)
  {
    const double z = coord_z[i];
    for (int j = 0; j < _n_elements_x + 1; ++j)
    {
      const double x = coord_x[j];
      _vertices[ver] = Point2(x, z);
      ++ver;
    }
  }

  delete[] coord_z;
  delete[] coord_x;
}




void RectangularMesh::build_elements()
{
  _n_elements = _n_elements_x * _n_elements_z; // the total number of rectangles
  _elements = new Rectangle[_n_elements]; // allocate the memory for all rectangles

  // the numbers of vertices describing a rectangle
  int vert[Rectangle::N_VERTICES];

  int rect = 0; // the number of a current rectangle
  for (int i = 0; i < _n_elements_z; ++i)
  {
    for (int j = 0; j < _n_elements_x; ++j)
    {
      // serial number of the cell is its index
      const int index = _n_elements_x*i + j;
      // numeration of rectangle's vertices is the following
      // 0 --- 1  ------> X
      // |     |  |
      // |     |  \/ Z
      // 2 --- 3
      vert[0] = i * (_n_elements_x + 1) + j;
      vert[1] = i * (_n_elements_x + 1) + j + 1;
      vert[2] = (i + 1) * (_n_elements_x + 1) + j;
      vert[3] = (i + 1) * (_n_elements_x + 1) + j + 1;
      _elements[rect] = Rectangle(vert, index);
      ++rect;
    }
  }
}



void RectangularMesh::assign_material_id(const TriangularMesh &tri_mesh,
                                         int n_rand_points)
{
  // First, we pass through the elements of the 'tri_mesh', consider each
  // triangle, generate n_rand_points random points in it, and assign the
  // material IDs to the coinciding elements of this rectangular grid. This goes
  // very fast. Then we consider the elements of the rectangular grid that are
  // still have no assigned properties. For them we go in a standard way, and
  // look for the elements of the 'tri_mesh' containing centers of the
  // rectangular cells. This is slow, but hopefully the number of unassigned
  // rectangular elements much smaller, than if we go standard way.

  double t0 = get_wall_time();
  std::cout << "Assigning material IDs..." << std::endl;

  // all unassigned from the beginning
  std::vector<bool> assigned(_n_elements, false);
  // there is a chance that several cells from the 'tri_mesh' could have their
  // centers in one rectangular element. To make the most correct assignment
  // we consider the distance between them to the center of the rectangle.
  // Assignment happens for the 'tri_cell' whose center is the closest to
  // center of the rectangle.
  std::vector<double> dist_to_center(_n_elements);
  // For the reason mentioned above let's compute and keep the centers of the
  // rectangle in advance.
  std::vector<Point2> rect_centers(_n_elements);
  for (int el = 0; el < _n_elements; ++el)
    rect_centers[el] = _elements[el].center(_vertices);

  std::mt19937 rand_engine; // random generator engine

  Point2 tri_vertices[Triangle::N_VERTICES]; // vertices of a triangle

  bool throw_exc = false; // we do not throw an exception if we can't find
                          // a proper rect cell for the one from the 'tri_mesh'
  for (int el = 0; el < tri_mesh.n_elements(); ++el)
  {
    const Triangle& tri_cell = tri_mesh.element(el);
    tri_cell.get_vertices(tri_mesh.get_vertices(), tri_vertices);
    const Point2& A = tri_vertices[0];
    const Point2& B = tri_vertices[1];
    const Point2& C = tri_vertices[2];
    for (int r = 0; r < n_rand_points; ++r)
    {
      std::uniform_real_distribution<double> distribution(0.0, 1.0);
      const double r1 = distribution(rand_engine);
      const double r2 = distribution(rand_engine);
      const double r3 = sqrt(r1);
      const Point2 tri_point = (1.-r3)*A + r3*(1.-r2)*B + r3*r2*C;

      const int rect_cell_number = this->find_element(tri_point, throw_exc);
      if (rect_cell_number != -1) // if the element is found
      {
        // if it wasn't assigned, assign properties and save the distance
        // between the point and the center of the rectangle
        if (!assigned[rect_cell_number])
        {
          _elements[rect_cell_number].set_material_id(tri_cell.material_id());
          assigned[rect_cell_number] = true;
          dist_to_center[rect_cell_number] = (rect_centers[rect_cell_number] -
                                              tri_point).L2_norm();
        } else { // if the rect cell has been already assigned, compare the new
                 // distance with the previously computed one
          double dist = (rect_centers[rect_cell_number]-tri_point).L2_norm();
          if (dist < dist_to_center[rect_cell_number])
          {
            _elements[rect_cell_number].set_material_id(tri_cell.material_id());
            dist_to_center[rect_cell_number] = dist;
          }
        }
      }
    }
  }

  int n_first_assigned = 0;
  for (size_t i = 0; i < assigned.size(); ++i)
    n_first_assigned += (assigned[i] ? 1 : 0);
  std::cout << "assign properties from triangular mesh to part of the rect one "
               "time = " << get_wall_time() - t0 << "\nN assigned rect cells "
            << n_first_assigned << " what is "
            << static_cast<int>(n_first_assigned*100/_n_elements)
            << "% of total number" << std::endl;

  // Here - there is a chance, that all rectangular elements are already
  // assigned. This is possible when the 'tri_mesh' is much finer than the
  // rectangular one. However, we check the rectangular elements to be sure.
  throw_exc = false; // again, we don't throw an exception if an unassigned
                     // rectangular element was not found in the 'tri_mesh'

  t0 = get_wall_time();
  int n_assigned = n_first_assigned;
  for (int el = 0; el < _n_elements; ++el)
  {
    if (!assigned[el])
    {
      const Point2 center = _elements[el].center(_vertices);
      const int tri_cell_number = tri_mesh.find_element(center, throw_exc);
      if (tri_cell_number != -1)
      {
        const Triangle& tri_cell = tri_mesh.element(tri_cell_number);
        _elements[el].set_material_id(tri_cell.material_id());
      }
      else
        _elements[el].set_material_id(VOID_MATERIAL_ID);

      assigned[el] = true;
      ++n_assigned;
      { // write down the statistics
        const int elements_left = _n_elements - n_assigned;
        const int percent_left = int(elements_left*100/_n_elements);
        if (elements_left % 1000 == 0)
          std::cout << elements_left << "(" << percent_left << "%) left to "
                       "assign, time elapsed = "
                    << get_wall_time() - t0 << std::endl;
      }
    }
  }

  std::cout << "Assigning material IDs is done. Time = "
            << get_wall_time() - t0 << std::endl;
}



int RectangularMesh::find_element(const Point2 &point,
                                  bool throw_exception) const
{
  const double px = point.x(); // coordinates of the point of interest
  const double pz = point.z();

  const double x0 = _min_coord.x(); // limits of the rect mesh
  const double x1 = _max_coord.x();
  const double z0 = _min_coord.z();
  const double z1 = _max_coord.z();

  // check that the point is within the mesh
  const double tol = FIND_CELL_TOLERANCE;
  if (px < x0 - tol || px > x1 + tol ||
      pz < z0 - tol || pz > z1 + tol)
  {
    if (throw_exception)
      require(false, "The given point " + d2s(point) + " doesn't belong "
              "to the rectangular mesh");

    return -1; // to show that the point in not here
  }

  // since the elements of the rectangular mesh are numerated in the following
  // way:
  // -----------
  // | 0  | 1  |
  // -----------
  // | 2  | 3  |
  // -----------
  // we can simplify the search of the element containing the given point:

  const double hx = (x1 - x0) / _n_elements_x;
  const double hz = (z1 - z0) / _n_elements_z;

  const int nx = std::min(static_cast<int>((px-x0)/hx), _n_elements_x-1);
  const int nz = std::min(static_cast<int>((pz-z0)/hz), _n_elements_z-1);

  if (nx < 0 || nx >= _n_elements_x ||
      nz < 0 || nz >= _n_elements_z)
  {
    if (throw_exception)
      require(false, "The rectangular element for the point " + d2s(point) +
              " wasn't found");

    return -1; // to show that the point in not here
  }

  return (nz*_n_elements_x + nx);
}



void RectangularMesh::
write_binary_files(const std::string &prop_filename,
                   std::vector<std::string> &out_filenames) const
{
  std::cout << "Writing binary files..." << std::endl;
  double t_begin = get_wall_time();

  // Create the map between the material IDs and the media properties. The
  // material IDs must exactly the same as in the given triangular mesh.
  std::map<int, std::vector<double> > properties;

  get_properties(prop_filename, properties);

  const int n_properties = properties.begin()->second.size();

  out_filenames.clear();
  out_filenames.resize(n_properties);

  std::ofstream *out = new std::ofstream[n_properties];
  for (int i = 0; i < n_properties; ++i)
  {
    out_filenames[i] = "property" + d2s(i) + ".bin";
    out[i].open(out_filenames[i].c_str(), std::ios::binary);
    require(out[i], "File '" + out_filenames[i] + "' can't be opened");
  }

  for (int el = 0; el < _n_elements; ++el)
  {
    const int matID = _elements[el].material_id();
    std::vector<double> values;
    if (matID == VOID_MATERIAL_ID)
    {
      values.resize(n_properties, VOID_MATERIAL_PROPERTY);
    }
    else
    {
      std::map<int, std::vector<double> >::const_iterator iter =
          properties.find(matID);
      require(iter != properties.end(), "The material ID " + d2s(matID) +
              " wasn't found in the properties file");
      values = iter->second;
    }

    for (int i = 0; i < n_properties; ++i)
    {
      OUT_FLOAT_TYPE val = values[i];
      out[i].write(reinterpret_cast<char*>(&val), sizeof(OUT_FLOAT_TYPE));
    }
  }

  for (int i = 0; i < n_properties; ++i)
    out[i].close();

  delete[] out;

  std::cout << "Writing binary files is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}



void get_properties(const std::string &filename,
                    std::map<int, std::vector<double> > &properties)
{
  std::cout << "Getting properties..." << std::endl;
  double t_begin = get_wall_time();

  std::ifstream in(filename.c_str());
  require(in, "File " + filename + " can't be opened");

  int materialID;
  std::vector<double> values;
  int n_values = 0; // default value

  std::string line;

  while (getline(in, line)) // read the file line-by-line
  {
    // if the line is empty or starts with '#' (a comment), we skip it
    if (line.empty() || line[0] == '#') continue;

    values.clear();
    if (n_values == 0) // first run
      values.reserve(10); // we don't expect more than 10 parameters by default
    else
      values.reserve(n_values);

    std::istringstream instr(line);
    instr >> materialID;
    double val;
    while (instr >> val)
      values.push_back(val);

    if (n_values == 0)
      n_values = values.size();
    else
      require(n_values == (int)values.size(), "The number of parameters in the "
              "first appearance was " + d2s(n_values) + ", but in some line "
              "there are " + d2s(values.size()) + " of them");

    // insert the vector into the map
    std::pair<std::map<int, std::vector<double> >::const_iterator, bool> res =
        properties.insert(std::pair<int, std::vector<double> >(materialID, values));

    // check that the insertion was successfull
    require(res.second, "Insertion of the values for the material with ID =" +
            d2s(materialID) + " failed. Check the data (file = " + filename +
            ")");
  }

  in.close();

  std::cout << "Getting properties is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}



void convert_to_xz(const std::vector<std::string>& out_filenames,
                   int nnx,
                   int nnz)
{
  std::cout << "Converting to xz plane..." << std::endl;
  double t_begin = get_wall_time();

  const int n_files = out_filenames.size();

  for (int f = 0; f < n_files; ++f)
  {
    const std::string filename = out_filenames[f];

    std::ifstream in(filename.c_str(), std::ios::binary);
    require(in, "File '" + filename + "' can't be opened");

    in.seekg(0, in.end); // jump to the end of the file
    int length = in.tellg(); // total length of the file in bytes
    int size_value = length / (nnx*nnz); // size (in bytes) of one value

    require(length % (nnx*nnz) == 0, "The number of bytes in the file " +
            filename + " is not divisible by the number of elements");

    in.seekg(0, in.beg); // jump to the beginning of the file

    double *values = new double[nnx*nnz];

    if (size_value == sizeof(double))
    {
      in.read((char*)values, nnx*nnz*size_value); // read all at once

      require(nnx*nnz == (int)in.gcount(), "The number of successfully read "
              "elements is different from the expected one");
    }
    else if (size_value == sizeof(float))
    {
      float val = 0;
      for (int i = 0; i < nnx*nnz; ++i)  // read element-by-element
      {
        in.read((char*)&val, size_value); // read a 'float' value
        values[i] = val;                  // convert it to a 'double' value
      }
    }
    else require(false, "Unknown size of an element in bytes");

    in.close();

    const std::string fname_out = file_stem(filename) + "_xz" +
                                  file_extension(filename);
    std::ofstream out(fname_out.c_str(), std::ios::binary);
    require(out, "File '" + fname_out + "' can't be opened");

    for (int i = nnz-1; i >= 0; --i)
    {
      for (int j = 0; j < nnx; ++j)
        if (size_value == sizeof(double))
          out.write((char*)&values[i*nnx+j], size_value);
        else
        {
          float val = values[i*nnx+j];
          out.write((char*)&val, size_value);
        }
    }

    out.close();

    delete[] values;
  }

  std::cout << "Converting to xz plane is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}



void convert_to_node_values(const std::vector<std::string>& out_filenames,
                            int nnx,
                            int nnz)
{
  std::cout << "Converting to node values..." << std::endl;
  double t_begin = get_wall_time();

  const int n_files = out_filenames.size();

  for (int f = 0; f < n_files; ++f)
  {
    const std::string filename = out_filenames[f];

    std::ifstream in(filename.c_str(), std::ios::binary);
    require(in, "File '" + filename + "' can't be opened");

    in.seekg(0, in.end); // jump to the end of the file
    int length = in.tellg(); // total length of the file in bytes
    int size_value = length / (nnx*nnz); // size (in bytes) of one value

    require(length % (nnx*nnz) == 0, "The number of bytes in the file " +
            filename + " is not divisible by the number of elements");

    in.seekg(0, in.beg); // jump to the beginning of the file

    double *values = new double[nnx*nnz];

    if (size_value == sizeof(double))
    {
      in.read((char*)values, nnx*nnz*size_value); // read all at once

      require(nnx*nnz == (int)in.gcount(), "The number of successfully read "
              "elements is different from the expected one");
    }
    else if (size_value == sizeof(float))
    {
      float val = 0;
      for (int i = 0; i < nnx*nnz; ++i)  // read element-by-element
      {
        in.read((char*)&val, size_value); // read a 'float' value
        values[i] = val;                  // convert it to a 'double' value
      }
    }
    else require(false, "Unknown size of an element in bytes");

    in.close();

    double *values_at_nodes = new double[(nnx+1)*(nnz+1)];

    for (int i = 0; i < nnz+1; ++i)
    {
      for (int j = 0; j < nnx+1; ++j)
      {
        double sum = 0.0;
        int count = 0;

        if (i-1 >= 0)  { sum += values[i-1]; ++count; }
        if (i+1 < nnz) { sum += values[i+1]; ++count; }
        if (j-1 >= 0)  { sum += values[j-1]; ++count; }
        if (j+1 < nnx) { sum += values[j+1]; ++count; }

        const int el = i*(nnx+1) + j;
        values_at_nodes[el] = sum / count;
      }
    }

    const std::string fname_out = file_stem(filename) + "_nodes" +
                                  file_extension(filename);
    std::ofstream out(fname_out.c_str(), std::ios::binary);
    require(out, "File '" + fname_out + "' can't be opened");

    if (size_value == sizeof(double))
    {
      out.write((char*)values_at_nodes, (nnx+1)*(nnz+1)*size_value);
    }
    else if (size_value == sizeof(float))
    {
      for (int i = 0; i < nnz+1; ++i)
      {
        for (int j = 0; j < nnx+1; ++j)
        {
          float val = values_at_nodes[i*(nnx+1)+j];
          out.write((char*)&val, size_value);
        }
      }
    }
    else require(false, "Unknown size of an element in bytes");

    out.close();

    delete[] values_at_nodes;
    delete[] values;
  }

  std::cout << "Converting to node values is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}


