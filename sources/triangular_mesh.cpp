#include "triangular_mesh.hpp"
#include "utilities.hpp"
#include "triangle.hpp"

#include <fstream>
#include <map>
#include <set>




TriangularMesh::TriangularMesh()
  : _min_coord(),
    _max_coord(),
    _vertices(),
    _elements()
{ }



TriangularMesh::~TriangularMesh()
{
  clear();
}



void TriangularMesh::clear()
{
  _min_coord = Point2();
  _max_coord = Point2();

  _vertices.clear();

  auto elem = _elements.begin();
  auto end  = _elements.end();
  for (; elem != end; ++elem)
    delete *elem;
  _elements.clear();
}



void TriangularMesh::read_msh(const std::string &meshfile)
{
  std::cout << "Reading a triangular mesh..." << std::endl;
  double t_begin = get_wall_time();

  clear(); // in case the mesh was initilized before

  // open the mesh file
  std::ifstream in(meshfile.c_str());
  require(in, "File " + meshfile + " cannot be opened!");

  std::string str;
  in >> str; // the first string of Gmsh file is "$MeshFormat"
  expect(str == "$MeshFormat",
         "The first string of the Gmsh file " + meshfile + " doesn't equal to "
         "\"$MeshFormat\". The actual string is \"" + str + "\"");

  // Read the information about the mesh
  double version;
  int binary, dsize;
  in >> version >> binary >> dsize;
  expect(version >= 2.2,
         "The version of Gmsh's mesh is unexpected (" + d2s(version) + ").");
  expect(dsize == sizeof(double),
         "The size of Gmsh's double (" + d2s(dsize) + ") doesn't equal to size "
         "of double type (" + d2s(sizeof(double)) + ")");

  getline(in, str); // read some empty string

  // There is additional 1 (the number - one) in binary format
  if (binary)
  {
    int one;
    in.read(reinterpret_cast<char*>(&one), sizeof(int));
    require(one == 1, "The binary one (" + d2s(one) + ") doesn't equal to 1!");
  }

  // We make a map between a serial number of the vertex and its number in the
  // file. It will help us when we create mesh elements.
  std::map<int, int> vertices_map;

  const int n_coords = 3; // in Gmsh the points are all 3D

  // Read lines of mesh file.
  // If we face specific keyword, we'll treat the section.
  while (in >> str)
  {
    if (str == "$Nodes") // read the mesh vertices
    {
      // the number of all mesh vertices (that are saved in the file)
      int n_vertices;
      in >> n_vertices; // read that number
      _vertices.resize(n_vertices); // allocate the memory for mesh vertices
      getline(in, str); // read some empty string

      int number; // the global number of the vertex
      // Cartesian coordinates of the vertex (Gmsh produces 3D mesh regardless
      // its real dimension)
      double coord[n_coords];
      // limits of the computational domain (maximal and minimal values of the
      // nodal coordinates)
      std::vector<double> maxcoord(n_coords);
      std::vector<double> mincoord(n_coords);

      // read vertices
      for (int ver = 0; ver < n_vertices; ++ver)
      {
        if (binary) // binary format
        {
          // a global number of each node
          in.read(reinterpret_cast<char*>(&number), sizeof(unsigned int));
          // node coordinates
          in.read(reinterpret_cast<char*>(coord), n_coords*sizeof(double));
        }
        else // ASCII format
        {
          in >> number; // read a global number of a node
          for (int i = 0; i < n_coords; ++i)
            in >> coord[i]; // read the coordinates of the node
        }

        if (ver == 0) // for the first vertex
        {
          // initialization of the max and min coordinates
          for (int i = 0; i < n_coords; ++i)
            maxcoord[i] = mincoord[i] = coord[i];
        }
        else // for the other vertices
        {
          for (int i = 0; i < n_coords; ++i)
          {
            // searching max and min coordinates
            maxcoord[i] = std::max(maxcoord[i], coord[i]);
            mincoord[i] = std::min(mincoord[i], coord[i]);
          }
        }

        _vertices[ver] = Point2(coord[0], coord[1]); // save the vertex
        vertices_map[number] = ver; // add the number of vertex to the map
      }
      // these points may or may not be one of the mesh vertices if the domain
      // has curvilinear boundaries
      _max_coord = Point2(maxcoord[0], maxcoord[1]);
      _min_coord = Point2(mincoord[0], mincoord[1]);

      expect(n_vertices == (int)vertices_map.size(),
             "Vertices numbers are not unique: n_vertices = " +
             d2s<int>(n_vertices) + " vertices_map.size() = " +
             d2s<unsigned>(vertices_map.size()));

    } // read the vertices

    else if (str == "$Elements") // read the mesh elements
    {
      int nelements; // the number of mesh elements
      in >> nelements; // read that number
      getline(in, str); // empty string

      int number; // the serial number of the element [1, nElements]
      int el_type; // the type of the element (1 - line, 2 - triangle, etc)
      int n_tags; // the number of tags describing the element
      int phys_domain; // the physical domain where the element takes place
      int elem_domain; // the elementary domain where the element takes place
      int n_partitions = 0; // number of partitions in which the element takes place
      int partition; // the partition which the element belongs to
      // "ghost cells" are other partitions which this element is connected with
      std::vector<int> ghost_cells;
      std::set<int> partitions; // a list of partitions

      // The map between the type of the element, and the number of nodes that
      // describes it. We mention here only expected elements that can appear
      // in a triangular (only triangular - other classes representing other
      // meshes can implement it in different way) mesh
      std::map<int, int> type_nodes;
      type_nodes[1] = 2; // 2-nodes line
      type_nodes[2] = 3; // 3-nodes triangle
      type_nodes[3] = 4; // 4-nodes quadrangle
      type_nodes[4] = 4; // 4-nodes tetrahedron
      type_nodes[5] = 8; // 8-nodes hexahedron
      type_nodes[15]= 1; // 1-node point

      if (binary) // binary format
      {
        int n_elem_part = 0; // part of all elements
        const int header_size = 3; // header consists of 3 numbers: type of the
                                   // element, number of elements of this type,
                                   // and number of tags
        int header[header_size];
        int n_elem_type; // number of elements of a specific type

        while (n_elem_part < nelements)
        {
          in.read(reinterpret_cast<char*>(header), header_size*sizeof(int));
          el_type     = header[0];
          n_elem_type = header[1];
          n_tags      = header[2];

          n_elem_part += n_elem_type;

          // how many vertices (nodes) describe the element
          std::map<int, int>::const_iterator el_type_iter =
              type_nodes.find(el_type);

          // check that the element has known (and expected) type
          require(el_type_iter != type_nodes.end(), "This type of the Gmsh's "
                  "element (" + d2s(el_type) + ") in the mesh file \"" +
                  meshfile + "\" is unknown");

          const int n_elem_nodes = el_type_iter->second; // the number of nodes
          std::vector<int> data(1 + n_tags + n_elem_nodes); // data for each element

          if (el_type == 2)                 // 3-nodes triangle
            _elements.reserve(n_elem_type); // allocate the memory

          for (int el = 0; el < n_elem_type; ++el)
          {
            in.read(reinterpret_cast<char*>(&data[0]), data.size()*sizeof(int));
            int dd = 0; // index for data array
            number = data[dd++];
            // physical domain - the most important value (to distinguish
            // materials with different properties)
            phys_domain = (n_tags > 0) ? data[dd++] : 0;
            // elementary domain - to distinguish different geometrical domains
            // (usually - rarely used)
            elem_domain = (n_tags > 1) ? data[dd++] : 0;
            // the number of tags is bigger than 2 if there are some partitions
            if (n_tags > 2)
            {
              // the number of partitions where this elements takes place
              n_partitions = data[dd++];
              expect(n_partitions >= 1, "The number of tags is more than 2. "
                     "That means that we have partitions. But the number of "
                     "partitions is " + d2s<int>(n_partitions));
              // the partition which the element belongs to
              partition = data[dd++] - 1; // we do (-1) since we associate the number
                                          // of partition (which is numerated from 1)
                                          // with the number of coarse element (which
                                          // is numerated from 0)
              partitions.insert(partition);
              // "ghost cells"
              if (n_partitions > 1) // if the element is on the boundary between
                                    // the partitions, it is described by "ghost
                                    // cells" as well
              {
                ghost_cells.resize(n_partitions - 1);
                for (int gc = 0; gc < n_partitions - 1; ++gc)
                {
                  // 'minus' since ghost cells are described by number of
                  // partition with the negative sign
                  ghost_cells[gc] = -data[dd++];
                  expect(ghost_cells[gc] > 0,
                         "The number of the ghost cell (positive one) is "
                         "unexpected (" + d2s<int>(ghost_cells[gc]) + ")");
                  // we decrease by 1 for the same reason as in case of the
                  // number of partition
                  --ghost_cells[gc];
                }
              }
            } // in case of partitions

            std::vector<int> nodes(n_elem_nodes); // allocate memory for nodes
            for (int i = 0; i < n_elem_nodes; ++i)
            {
              // vertices can be numerated not sequentially (or at least not
              // from 0)
              nodes[i] = vertices_map.find(data[dd++])->second;
            }

            data.clear();

            const int extra_id = (n_partitions >= 1 ? partition : elem_domain);

            // add the new element in the list
            if (el_type == 2) // 3-nodes triangle
            {
              _elements.push_back(new Triangle(nodes,
                                               _elements.size(), // serial number
                                               phys_domain,
                                               extra_id));
            }
          } // pass through all elements of one type
        }

        expect(n_elem_part == nelements, "We read " + d2s(n_elem_part) +
               " elements, whereas there are " + d2s(nelements) + " of them");

      } // binary format

      else // ASCII format
      {
        for (int el = 0; el < nelements; ++el)
        {
          // read serial number, type of an element, and number of tags
          in >> number >> el_type >> n_tags;
          std::vector<int> data(n_tags); // allocate the memory for some data
          for (int i = 0; i < n_tags; ++i) // read this information
            in >> data[i];
          // physical domain - the most important value (to distinguish
          // materials with different properties)
          phys_domain = (n_tags > 0) ? data[0] : 0;
          // elementary domain - to distinguish different geometrical domains
          // (usually - rarely used)
          elem_domain = (n_tags > 1) ? data[1] : 0;
          // the number of tags is bigger than 2 if there are some partitions
          if (n_tags > 2)
          {
            // the number of partitions where this elements takes place
            n_partitions = data[2];
            expect(n_partitions >= 1, "The number of tags is more than 2. That "
                   "means that we have partitions. But the number of "
                   "partitions is " + d2s<int>(n_partitions));
            // the partition which the element belongs to
            partition = data[3] - 1; // we do (-1) since we associate the number
                                     // of partition (which is numerated from 1)
                                     // with the number of coarse element (which
                                     // is numerated from 0)
            partitions.insert(partition);
            // "ghost cells"
            if (n_partitions > 1) // if the element is on the boundary between
                                  // the partitions, it is described by "ghost
                                  // cells" as well
            {
              ghost_cells.resize(n_partitions - 1);
              for (int gc = 0; gc < n_partitions - 1; ++gc)
              {
                // 'minus' since ghost cells are described by number of
                // partition with the negative sing
                ghost_cells[gc] = -data[4 + gc];
                expect(ghost_cells[gc] > 0,
                       "The number of the ghost cell (positive one) is "
                       "unexpected (" + d2s<int>(ghost_cells[gc]) + ")");
                // we decrease by 1 for the same reason as in case of the
                // number of partition
                --ghost_cells[gc];
              }
            }
          }

          data.clear(); // clear the memory

          // how many vertices (nodes) describe the element
          std::map<int, int>::const_iterator el_type_iter =
              type_nodes.find(el_type);

          // check that the element has known (and expected) type
          require(el_type_iter != type_nodes.end(),
                  "This type of the Gmsh's element (" + d2s(el_type) +
                  ") in the mesh file \"" + meshfile + "\" is unknown");

          const int n_elem_nodes = el_type_iter->second; // the number of nodes
          std::vector<int> nodes(n_elem_nodes); // allocate memory for nodes
          for (int i = 0; i < n_elem_nodes; ++i)
          {
            in >> nodes[i]; // read the numbers of nodes
            // vertices can be numerated not sequentially
            // (or at least not from 0)
            nodes[i] = vertices_map.find(nodes[i])->second;
          }

          const int extra_id = (n_partitions >= 1 ? partition : elem_domain);

          // add the new element in the list
          if (el_type == 2) // 3-nodes triangle
          {
            _elements.push_back(new Triangle(nodes,
                                             _elements.size(), // serial number
                                             phys_domain,
                                             extra_id));
          }
        } // loop over mesh elements

        // check some expectations
        expect(number == nelements, "The number of the last read Gmsh's "
               "element (" + d2s<int>(number) + ") is not equal to the amount "
               "of all elements in the mesh (" + d2s<int>(nelements) + ")");

      } // ASCII format

      // requirements after reading elements
      require(!_vertices.empty(), "There are no vertices in a mesh described in"
              " the file " + meshfile);
      require(!_elements.empty(), "There are no triangles in the mesh described"
              " in the file " + meshfile);

      // check that the numbers of partitions are sequential
      for (int par = 0; par < (int)partitions.size(); ++par)
        require(*partitions.find(par) == par, "The numeration of the partitions"
                " is not dense");

    } // read the elements
  } // read the mesh file

  in.close(); // close the file

  std::cout << "Reading a triangular mesh is done. Time = "
            << get_wall_time() - t_begin << std::endl;
}



const Triangle& TriangularMesh::element(int number) const
{
  expect(number >= 0 && number < (int)_elements.size(), "The requested element "
         "(number " + d2s(number) + ") is out of range [0, " +
         d2s(_elements.size()) + ")");

  return *_elements[number];
}



int TriangularMesh::find_element(const Point2 &point,
                                 bool throw_exception) const
{
  const double px = point.x();
  const double pz = point.z();

  // if the point is out of the region bounded by the rectangle connecting the
  // min_coord and max_coord (i.e. min and max points) of the mesh, that means
  // that the point clearly doesn't belong to the mesh
  const double x0 = _min_coord.x();
  const double x1 = _max_coord.x();
  const double z0 = _min_coord.z();
  const double z1 = _max_coord.z();
  if (px < x0 - FIND_CELL_TOLERANCE ||
      px > x1 + FIND_CELL_TOLERANCE ||
      pz < z0 - FIND_CELL_TOLERANCE ||
      pz > z1 + FIND_CELL_TOLERANCE)
  {
    if (throw_exception)
      require(false, "A triangle containing point " + d2s(point) + " was not "
              "found");

    return -1; // an element wan't found
  }

  // If the point is within the region bounded by min and max possible points,
  // that also doesn't mean that the point belong to the mesh, but in this case
  // we have to check all elements.

  for (size_t el = 0; el < _elements.size(); ++el)
  {
    if(_elements[el]->contains_point(point, _vertices))
      return el;
  }

  // We can be here only if we didn't find a cell
  if (throw_exception)
    require(false, "A triangle containing point " + d2s(point) + " was not "
            "found");

  return -1; // to show that the element wasn't found
}
