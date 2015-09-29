#include "point.hpp"
#include "utilities.hpp"

#include <cmath>



Point2::Point2(const double coordinates[])
  : _x(coordinates[0]),
    _z(coordinates[1])
{ }



Point2::Point2(const std::vector<double> &coordinates)
  : _x(coordinates[0]),
    _z(coordinates[1])
{
  expect(N_COORD == (int)coordinates.size(), "The size of the provided vector "
         "of coordinates (" + d2s(coordinates.size()) + ") doesn't coincide "
         "with the required number of coordinates (" + d2s(N_COORD) + ")");
}



Point2::Point2(double x_coord,
               double z_coord)
  : _x(x_coord),
    _z(z_coord)
{ }



Point2::Point2(const Point2 &p)
  : _x(p._x),
    _z(p._z)
{ }



Point2& Point2::operator =(const Point2 &p)
{
  _x = p._x;
  _z = p._z;
  return *this;
}



Point2& Point2::operator *=(double d)
{
  _x *= d;
  _z *= d;
  return *this;
}



Point2& Point2::operator /=(double d)
{
  expect(fabs(d) > VERY_SMALL_NUMBER, "Divide by very small number (" +
         d2s<double>(d, true) + ") may cause big problems");

  _x /= d;
  _z /= d;
  return *this;
}



Point2& Point2::operator +=(const Point2 &p)
{
  _x += p._x;
  _z += p._z;
  return *this;
}



Point2& Point2::operator -=(const Point2 &p)
{
  _x -= p._x;
  _z -= p._z;
  return *this;
}



Point2 operator *(double d, const Point2 &point)
{
  Point2 p = point;
  p._x *= d;
  p._z *= d;
  return p;
}



Point2 operator +(const Point2 &p1, const Point2 &p2)
{
  Point2 res = p1;
  res += p2;
  return res;
}



Point2 operator -(const Point2 &p1, const Point2 &p2)
{
  Point2 res = p1;
  res -= p2;
  return res;
}



std::ostream& operator <<(std::ostream &os, const Point2 &p)
{
  os << "(" << p._x << ", " << p._z << ")";
  return os;
}



bool Point2::compare_by_x(const Point2 &a, const Point2 &b)
{
  return (a._x < b._x);
}



bool Point2::compare_by_z(const Point2 &a, const Point2 &b)
{
  return (a._z < b._z);
}



bool Point2::operator ()(const Point2 &a, const Point2 &b) const
{
  return compare_by_norm(a, b);
}



bool Point2::operator <(const Point2 &p) const
{
  return compare_by_norm(*this, p);
}



bool Point2::compare_by_norm(const Point2 &a, const Point2 &b)
{
  const Point2 c = a - b;

  return (c.L2_norm() > FLOAT_NUMBERS_EQUALITY_TOLERANCE &&
          a.L2_norm() < b.L2_norm());
}



double dot_product(const Point2 &a, const Point2 &b)
{
  return (a._x * b._x +
          a._z * b._z);
}



double Point2::L2_norm() const
{
  return sqrt(_x*_x + _z*_z);
}

