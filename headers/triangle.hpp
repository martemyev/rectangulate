#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include "config.hpp"

#include <vector>

class Point2;



class Triangle
{
public:

  static const int N_VERTICES = 3;

  Triangle();

  Triangle(const std::vector<int> &ver,
           int number,
           int mat_id,
           int part_id = 0);

  Triangle(const Triangle &tri);

  Triangle& operator =(const Triangle &tri);

  ~Triangle() { }

  int material_id() const { return _material_id; }

  Point2 center(const std::vector<Point2> &points) const;

  bool contains_point(const Point2 &point,
                      const std::vector<Point2> &mesh_points,
                      double *areas_ratio = nullptr) const;

  /**
   * Get triangle's vertices with the coordinates.
   */
  void get_vertices(const std::vector<Point2> &points,
                    Point2 *vertices) const;

  int vertex(int num) const;

private:

  std::vector<int> _vertices; ///< Vertices indices

  int _number; ///< Index of the element

  int _material_id; ///< ID of the physical domain where the element takes place

  int _partition_id; ///< ID of partition (in case of domain decomposition)

  double area(const Point2& v0, const Point2& v1, const Point2& v2) const;
};



#endif // TRIANGLE_HPP
