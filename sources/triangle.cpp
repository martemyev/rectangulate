#include "triangle.hpp"



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


