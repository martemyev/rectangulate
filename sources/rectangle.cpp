#include "rectangle.hpp"
#include "point.hpp"
#include "utilities.hpp"

#include <cmath>




Rectangle::Rectangle()
  : _number(0),
    _material_id(0)
{ }



Rectangle::Rectangle(int *vert, int number_, int mat_id)
  : _number(number_),
    _material_id(mat_id)
{
  for (int i = 0; i < N_VERTICES; ++i)
    _vertices[i] = vert[i];
}



Rectangle::Rectangle(int v1,
                     int v2,
                     int v3,
                     int v4,
                     int number_,
                     int mat_id)
  : _number(number_),
    _material_id(mat_id)
{
  _vertices[0] = v1;
  _vertices[1] = v2;
  _vertices[2] = v3;
  _vertices[3] = v4;
}



Rectangle::Rectangle(const Rectangle &rect)
  : _number(rect._number),
    _material_id(rect._material_id)
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
  return *this;
}



Point2 Rectangle::center(const Point2 *points) const
{
  // the numeration of rectangle's vertices is the following
  // 0 --- 1
  // |     |
  // 2 --- 3

  const Point2& v0 = points[_vertices[0]];
  const Point2& v3 = points[_vertices[3]];

  const double xcen = 0.5 * (v0.x() + v3.x());
  const double zcen = 0.5 * (v0.z() + v0.z());

  return Point2(xcen, zcen);
}



