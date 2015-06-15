#ifndef RECTANGULAR_MESH_HPP
#define RECTANGULAR_MESH_HPP

#include "config.hpp"
#include "point.hpp"

#include <vector>
#include <map>

class Rectangle;
class TriangularMesh;



class RectangularMesh
{
public:

  /**
   * @param beg - the starting point (left low corner) of the domain.
   *              this point will initialize _min_coord.
   * @param end - the end point (right high corner) of the domain
   *              this point will initialize _max_coord.
   * @param nx - the number of the rectangular elements in x-direction
   * @param nz - the number of the rectangular elements in z-direction
   */
  RectangularMesh(const Point2 &beg,
                  const Point2 &end,
                  int nx,
                  int nz);

  ~RectangularMesh();

  void build();

  void assign_material_id(const TriangularMesh &tri_mesh);

  void write_binary_files(const std::string &filename) const;


protected: // ========================== PROTECTED =============================

  Point2 _min_coord; ///< This point is the left upper corner of the domain.
  Point2 _max_coord; ///< This point is the right lower corner of the domain.

  int _n_elements_x; ///< The number of the mesh elements in x-direction.
  int _n_elements_z; ///< The number of the mesh elements in z-direction

  int _n_vertices;   ///< The number of the mesh vertices.
  Point2 *_vertices; ///< Vector of vertices.

  int _n_elements;      ///< The number of the mesh cells.
  Rectangle *_elements; ///< Mesh elements

  RectangularMesh(const RectangularMesh &mesh);
  RectangularMesh& operator =(const RectangularMesh &mesh);

  void build_vertices();
  void build_elements();

  int find_element(const Point2 &point, bool throw_exception) const;
};


void get_properties(const std::string &filename,
                    std::map<int, std::vector<double> > &properties);


#endif // RECTANGULAR_MESH_HPP
