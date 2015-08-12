#ifndef POINT_HPP
#define POINT_HPP

#include "config.hpp"
#include <iostream>
#include <vector>


/**
 * @brief 2D point
 *
 * A point in 2-dimensional space: (x,z). Because of the geophysical nature of
 * possible applications we adopt this notation: the axis are 'x' and 'z'.
 */
class Point2
{
public:

  /**
   * @brief Number of coordinates describing a point
   *
   * The number of Cartesian coordinates, that describe the point. Here we
   * always use 2 coordinates to describe a point.
   */
  static const int N_COORD = 2;

  Point2() : _x(0.), _z(0.) { }

  Point2(const double coordinates[]);

  Point2(const std::vector<double> &coordinates);

  Point2(double x_coord,
         double z_coord);

  Point2(const Point2 &p);

  ~Point2() { }

  Point2& operator =(const Point2 &p);

  /**
   * Get an x-coordinate of the point
   */
  double x() const { return _x; }

  /**
   * Get a z-coordinate of the point
   */
  double z() const { return _z; }

  /**
   * Get/set an x-coordinate of the point
   */
  double& x() { return _x; }

  /**
   * Get/set a z-coordinate of the point
   */
  double& z() { return _z; }

  /**
   * @brief point *= scalar
   *
   * Multiply the point by a scalar. In this case all coordinates are multiplied
   * by this number. It is a scaling.
   * @param d - scaling factor
   */
  Point2& operator *=(double d);

  /**
   * @brief point /= scalar
   *
   * Divide the point by some number. In this case all coordinates are divided
   * by this number. It is a scaling.
   * @param d - denominator
   */
  Point2& operator /=(double d);

  /**
   * @brief point_a += point_b
   *
   * Add a point to this one: point_a = point_a + point_b. The points are added
   * by coordinates.
   * @param p - a point to be added
   */
  Point2& operator +=(const Point2 &p);

  /**
   * @brief point_a -= point_b
   *
   * Subtract a point from this one: point_a = point_a - point_b. The points
   * are subtracted by coordinates.
   * @param p - a point to be subtracted
   */
  Point2& operator -=(const Point2 &p);

  /**
   * @brief Scaling of the point
   *
   * Multiplication of all the coordinates of the given point by a scalar.
   *
   * @param d - a scalar
   * @param point - a point to be scaled
   * @return a newly created scaled point
   */
  friend Point2 operator *(double d, const Point2 &point);

  /**
   * @brief point3 = point1 + point2
   *
   * Add one point to another one: p = p1 + p2. The points are added by
   * coordinates.
   * @param p1 - one point (left operand)
   * @param p2 - another point (right operand)
   * @return p = p1 + p2
   */
  friend Point2 operator +(const Point2 &p1, const Point2 &p2);

  /**
   * @brief point3 = point1 - point2
   *
   * Subtract one point from another one: p = p1 - p2. The points are
   * substracted by coordinates.
   * @param p1 - one point (left operand)
   * @param p2 - another point (right operand)
   * @return p = p1 - p2
   */
  friend Point2 operator -(const Point2 &p1, const Point2 &p2);

  /**
   * @brief |point1| < |point2|
   *
   * Comparison between two points based on their norms (i.e. distance from the
   * origin). This function simply calls
   * compare_by_norm(const Point &a, const Point &b).
   * @param a - one point for comparison
   * @param b - another point for comparison
   * @return true if 'b' point is farther from the origin than 'a' point
   */
  bool operator ()(const Point2 &a, const Point2 &b) const;

  /**
   * @brief |this_point| < |another_point|
   *
   * Comparison between this and another points based on their norms (i.e.
   * distance from the origin). This function simply calls
   * compare_by_norm(const Point &a, const Point &b).
   * @param p - another point for comparison
   * @return true if 'p' point is farther from the origin than 'this' point
   */
  bool operator <(const Point2 &p) const;

  /**
   * @brief Output to a stream
   *
   * Output by '<<' to some stream
   * @param os - output stream (e.g. std::cout)
   * @param p - the point itself
   * @return changed output stream
   */
  friend std::ostream& operator <<(std::ostream &os, const Point2 &p);

  /**
   * @brief point1.x < point2.x
   *
   * The comparison between points which is based on their x-coordinates.
   * The bigger x-coordinates, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their x-coordinates.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger x-coordinate
   */
  static bool compare_by_x(const Point2 &a, const Point2 &b);

  /**
   * @brief point1.z < point2.z
   *
   * The comparison between points which is based on their z-coordinates.
   * The bigger z-coordinates, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their z-coordinates.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger z-coordinate
   */
  static bool compare_by_z(const Point2 &a, const Point2 &b);

  /**
   * @brief |point1| < |point2|
   *
   * The comparison between points which is based on the norm of the
   * vector between the point and the origin (0, 0).
   * The farther from origin, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their distance from the origin.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger norm
   */
  static bool compare_by_norm(const Point2 &a, const Point2 &b);

  /**
   * Dot (inner) product of the two given points (which are 2D vectors).
   */
  friend double dot_product(const Point2 &a, const Point2 &b);

  /**
   * L2 norm of 'this' point (which is a 2D vector).
   */
  double L2_norm() const;


private: // ================= PRIVATE =======================

  /**
   * Cartesian coordinates of the point.
   */
  double _x, _z;
};






class PhysicalPoint2
{
public:

  PhysicalPoint2();
  PhysicalPoint2(const PhysicalPoint2& pp);
  PhysicalPoint2& operator =(const PhysicalPoint2& pp);
  PhysicalPoint2(double x_, double z_, int mat_ID);
  ~PhysicalPoint2() { }

  Point2 point() const { return _point; }
  int material_ID() const { return _material_ID; }
  void material_ID(int matID) { _material_ID = matID; }

private:

  Point2 _point; ///< Coordinates of the physical point
  int _material_ID; ///< Material ID of the physical point (to define physical
                    ///< properties)
};



#endif // POINT_HPP
