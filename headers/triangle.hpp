#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include "config.hpp"

#include <vector>



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



private:

  std::vector<int> _vertices; ///< Vertices indices

  int _number; ///< Index of the element

  int _material_id; ///< ID of the physical domain where the element takes place

  int _partition_id; ///< ID of partition (in case of domain decomposition)

};



#endif // TRIANGLE_HPP
