#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include "config.hpp"
#include "point.hpp"



class Rectangle
{
public:

  /**
   * The number of vertices of a rectangle
   */
  static const int N_VERTICES = 4;

  Rectangle();

  /**
   * @param vert - indices of the vertices
   * @param points - vertices of the mesh
   * @param number - a (serial) number (index) of the rectangle
   * @param mat_id - material ID
   */
  Rectangle(int *vert,
            Point2 *points,
            int number = 0,
            int mat_id = 0);

  /**
   * @param v1 - first vertex
   * @param v2 - second vertex
   * @param v3 - third vertex
   * @param v4 - fourth vertex
   * @param points - vertices of the mesh
   * @param number - a (serial) number (index) of the rectangle
   * @param mat_id - material ID
   */
  Rectangle(int v1, int v2, int v3, int v4,
            Point2 *points,
            int number = 0,
            int mat_id = 0);

  Rectangle(const Rectangle &rect);
  Rectangle& operator =(const Rectangle &rect);

  ~Rectangle() { }

  int number() const { return _number; }
  int material_id() const { return _material_id; }
  void set_material_id(int mat_id) { _material_id = mat_id; }

  const int* get_vertices() const { return _vertices; }
  int vertex(int number) const;

  Point2 center() const { return _center; }



protected: // =========================== PROTECTED ============================

  int _vertices[N_VERTICES]; ///< Indices of the vertices
  int _number;               ///< Index
  int _material_id;          ///< ID representing some domain (material)
  Point2 _center;            ///< Center of the rectangle
};



#endif // RECTANGLE_HPP
