#include "triangle.hpp"
#include "utilities.hpp"
#include "point.hpp"

#include <cmath>



Triangle::Triangle()
  : _vertices(),
    _number(0),
    _material_id(0),
    _partition_id(0)
{ }



Triangle::Triangle(const std::vector<int> &ver,
                   int num,
                   int mat_id,
                   int part_id)
  : _vertices(ver),
    _number(num),
    _material_id(mat_id),
    _partition_id(part_id)
{ }



Triangle::Triangle(const Triangle &tri)
  : _vertices(tri._vertices),
    _number(tri._number),
    _material_id(tri._material_id),
    _partition_id(tri._partition_id)
{ }



Triangle& Triangle::operator =(const Triangle &tri)
{
  _vertices     = tri._vertices;
  _number       = tri._number;
  _material_id  = tri._material_id;
  _partition_id = tri._partition_id;

  return *this;
}



Point2 Triangle::center(const std::vector<Point2> &mesh_points) const
{
  Point2 vert[Triangle::N_VERTICES];
  for (int i = 0; i < Triangle::N_VERTICES; ++i)
    vert[i] = mesh_points[_vertices[i]];

  double xcen = 0.0;
  double zcen = 0.0;

  for (int i = 0; i < Triangle::N_VERTICES; ++i)
  {
    xcen += vert[i].x();
    zcen += vert[i].z();
  }

  xcen /= Triangle::N_VERTICES;
  zcen /= Triangle::N_VERTICES;

  return Point2(xcen, zcen);
}



bool Triangle::contains_point(const Point2 &point,
                              const std::vector<Point2> &mesh_points) const
{
  // we create 3 triangles defined by the vertices of the triangle and by the
  // point of interest. if the areas of these 3 triangles coincide with the
  // area of the original triangle - the point is inside the triangle

  const Point2& v0 = mesh_points[_vertices[0]];
  const Point2& v1 = mesh_points[_vertices[1]];
  const Point2& v2 = mesh_points[_vertices[2]];

  const double area_orig = area(v0, v1, v2);

  expect(area_orig > VERY_SMALL_NUMBER, "The area of the triangle is too "
         "small: " + d2s<double>(area_orig, true, 12));

  const double area_0 = area(v0, v1, point);
  const double area_1 = area(v0, point, v2);
  const double area_2 = area(point, v1, v2);

  const double diff = fabs((area_orig - area_0 - area_1 - area_2) / area_orig);

  if (diff < FIND_CELL_TOLERANCE)
    return true;

  return false;
}



void Triangle::get_vertices(const std::vector<Point2> &points,
                            Point2 *vertices) const
{
  for (int v = 0; v < N_VERTICES; ++v)
  {
    const int vert = _vertices[v];
    expect(vert >= 0 && vert < static_cast<int>(points.size()), "The " + d2s(v)+
           "-th vertex (" + d2s(vert) + " is out of range [0, n_points), where "
           "n_points = points.size(), where points is a given vector");

    vertices[v] = points[vert];
  }
}



double Triangle::area(const Point2& v0,
                      const Point2& v1,
                      const Point2& v2) const
{
  const double x0 = v0.x();
  const double x1 = v1.x();
  const double x2 = v2.x();

  const double z0 = v0.z();
  const double z1 = v1.z();
  const double z2 = v2.z();

  return 0.5 * fabs((x1 - x0) * (z2 - z0) -
                    (x2 - x0) * (z1 - z0));
}

