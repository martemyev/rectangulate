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
    _verts_ID(nullptr),
    _cells_ID(nullptr)
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
  delete[] _cells_ID;
  delete[] _verts_ID;
}




void RectangularMesh::
assign_material_id_in_cells(const TriangularMesh &tri_mesh)
{
  double t0 = get_wall_time();
  std::cout << "Assigning material IDs in cells..." << std::endl;

  _cells_ID = new int[_n_elements_x*_n_elements_z];

  const double x0 = _min_coord.x(); // limits in x-direction
  const double x1 = _max_coord.x();
  const double z0 = _min_coord.z(); // limits in z-direction
  const double z1 = _max_coord.z();
  const double hx = (x1 - x0) / _n_elements_x; // step in x-direction
  const double hz = (z1 - z0) / _n_elements_z; // step in z-direction

  const bool throw_exception = false;

  for (int iz = 0; iz < _n_elements_z; ++iz)
  {
    const double zcen = z0 + (iz+0.5)*hz; // z-coord of a cell center
    int triangle_index = 0; // reset the index from the previous search
    for (int ix = 0; ix < _n_elements_x; ++ix)
    {
      const double xcen = x0 + (ix+0.5)*hx; // x-coord of a cell center
      Point2 cell_center(xcen, zcen);
      require(tri_mesh.contains_point(cell_center), "The triangular mesh "
              "doesn't contain the point: " + d2s(cell_center));

      // every search along the x-line at the same z-coordinate starts with the
      // previously found triangle
      triangle_index = tri_mesh.find_element(cell_center,
                                             triangle_index,
                                             throw_exception);

      if (triangle_index < 0) // if it wasn't found
      {
        std::cout << "  full search for point " << cell_center << std::endl;
        triangle_index = tri_mesh.full_search(cell_center);
      }

      const int mat_ID = tri_mesh.element(triangle_index).material_id();
      _cells_ID[iz*_n_elements_x + ix] = mat_ID;
    }
  }

  std::cout << "Assigning material IDs in cells is done. Time = "
            << get_wall_time() - t0 << std::endl;
}


void RectangularMesh::
assign_material_id_at_nodes(const TriangularMesh &tri_mesh)
{
  double t0 = get_wall_time();
  std::cout << "Assigning material IDs at nodes..." << std::endl;

  _verts_ID = new int[(_n_elements_x+1)*(_n_elements_z+1)];

  const double x0 = _min_coord.x(); // limits in x-direction
  const double x1 = _max_coord.x();
  const double z0 = _min_coord.z(); // limits in z-direction
  const double z1 = _max_coord.z();
  const double hx = (x1 - x0) / _n_elements_x; // step in x-direction
  const double hz = (z1 - z0) / _n_elements_z; // step in z-direction

  const bool throw_exception = false;

  for (int iz = 0; iz < _n_elements_z+1; ++iz)
  {
    const double z = (iz == _n_elements_z ? z1 : z0 + iz*hz);
    int triangle_index = 0; // reset the index from the previous search
    for (int ix = 0; ix < _n_elements_x+1; ++ix)
    {
      const double x = (ix == _n_elements_x ? x1 : x0 + ix*hx);
      Point2 vertex(x, z);
      require(tri_mesh.contains_point(vertex), "The triangular mesh doesn't "
              "contain the point: " + d2s(vertex));

      // every search along the x-line at the same z-coordinate starts with the
      // previously found triangle
      triangle_index = tri_mesh.find_element(vertex,
                                             triangle_index,
                                             throw_exception);

      if (triangle_index < 0) // if it wasn't found
      {
        std::cout << "  full search for point " << vertex << std::endl;
        triangle_index = tri_mesh.full_search(vertex);
      }

      const int mat_ID = tri_mesh.element(triangle_index).material_id();
      _verts_ID[iz*(_n_elements_x+1)+ix] = mat_ID;
    }
  }

  std::cout << "Assigning material IDs at nodes is done. Time = "
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
write_binary_files_in_cells(const std::string &prop_filename,
                            std::vector<std::string> &filenames_out) const
{
  std::cout << "Writing binary files in cells..." << std::endl;
  double t_begin = get_wall_time();

  // Create the map between the material IDs and the media properties. The
  // material IDs must exactly the same as in the given triangular mesh.
  std::map<int, std::vector<double> > properties;

  get_properties(prop_filename, properties);

  const int n_properties = properties.begin()->second.size();

  filenames_out.resize(n_properties);

  std::ofstream *out = new std::ofstream[n_properties];
  for (int i = 0; i < n_properties; ++i)
  {
    filenames_out[i] = "property" + d2s(i) + "_cells.bin";
    out[i].open(filenames_out[i].c_str(), std::ios::binary);
    require(out[i], "File '" + filenames_out[i] + "' can't be opened");
  }

  for (int i = 0; i < _n_elements_x*_n_elements_z; ++i)
  {
    const int matID = _cells_ID[i];
    std::map<int, std::vector<double> >::const_iterator iter =
        properties.find(matID);
    require(iter != properties.end(), "The material ID " + d2s(matID) +
            " wasn't found in the properties file");
    const std::vector<double> values = iter->second;

    for (int p = 0; p < n_properties; ++p)
    {
      OUT_FLOAT_TYPE val = values[p];
      out[p].write(reinterpret_cast<char*>(&val), sizeof(OUT_FLOAT_TYPE));
    }
  }

  for (int i = 0; i < n_properties; ++i)
    out[i].close();

  delete[] out;

  std::cout << "Writing binary files in cells is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}




void RectangularMesh::
write_ASCII_files_at_nodes(const std::string &prop_filename) const
{
  std::cout << "Writing ASCII files at nodes..." << std::endl;
  double t_begin = get_wall_time();

  // Create the map between the material IDs and the media properties. The
  // material IDs must exactly the same as in the given triangular mesh.
  std::map<int, std::vector<double> > properties;

  get_properties(prop_filename, properties);

  const int n_properties = properties.begin()->second.size();

  std::vector<std::string> filenames_out(n_properties);

  std::ofstream *out = new std::ofstream[n_properties];
  for (int i = 0; i < n_properties; ++i)
  {
    filenames_out[i] = "property" + d2s(i) + "_nodes.txt";
    out[i].open(filenames_out[i].c_str(), std::ios::binary);
    require(out[i], "File '" + filenames_out[i] + "' can't be opened");
  }

  for (int i = 0; i < (_n_elements_x+1)*(_n_elements_z+1); ++i)
  {
    const int matID = _verts_ID[i];
    std::map<int, std::vector<double> >::const_iterator iter =
        properties.find(matID);
    require(iter != properties.end(), "The material ID " + d2s(matID) +
            " wasn't found in the properties file");
    const std::vector<double> values = iter->second;

    for (int p = 0; p < n_properties; ++p)
      out[p] << values[p] << "\n";
  }

  for (int i = 0; i < n_properties; ++i)
    out[i].close();

  delete[] out;

  std::cout << "Writing ASCII files at nodes is done. Time = "
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



void convert_in_cells_to_xz(const std::vector<std::string>& filenames,
                            int nnx,
                            int nnz)
{
  std::cout << "Converting to xz plane..." << std::endl;
  double t_begin = get_wall_time();

  const int n_files = filenames.size();

  for (int f = 0; f < n_files; ++f)
  {
    const std::string filename = filenames[f];

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
