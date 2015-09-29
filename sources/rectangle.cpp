#include "rectangle.hpp"
#include "point.hpp"
#include "utilities.hpp"

#include <cmath>



Rectangle::Rectangle()
  : _number(0)
  , _material_id(0)
  , _center()
{ }



Rectangle::Rectangle(int *vert, Point2 *points, int num, int mat_id)
  : _number(num)
  , _material_id(mat_id)
  , _center()
{
  for (int i = 0; i < N_VERTICES; ++i)
    _vertices[i] = vert[i];
  const double x0 = points[_vertices[0]].x();
  const double z0 = points[_vertices[0]].z();
  const double x1 = points[_vertices[3]].x();
  const double z1 = points[_vertices[3]].z();
  _center = Point2(0.5*(x0+x1), 0.5*(z0+z1));
}



Rectangle::Rectangle(int v1, int v2, int v3, int v4,
                     Point2 *points,
                     int num,
                     int mat_id)
  : _number(num)
  , _material_id(mat_id)
  , _center()
{
  _vertices[0] = v1;
  _vertices[1] = v2;
  _vertices[2] = v3;
  _vertices[3] = v4;
  const double x0 = points[_vertices[0]].x();
  const double z0 = points[_vertices[0]].z();
  const double x1 = points[_vertices[3]].x();
  const double z1 = points[_vertices[3]].z();
  _center = Point2(0.5*(x0+x1), 0.5*(z0+z1));
}



Rectangle::Rectangle(const Rectangle &rect)
  : _number(rect._number)
  , _material_id(rect._material_id)
  , _center(rect._center)
{
  for (int i = 0; i < N_VERTICES; ++i)
    _vertices[i] = rect._vertices[i];
}



Rectangle& Rectangle::operator =(const Rectangle &rect)
{
  for (int i = 0; i < N_VERTICES; ++i)
    _vertices[i] = rect._vertices[i];
  _number = rect._number;
  _material_id = rect._material_id;
  _center = rect._center;
  return *this;
}
