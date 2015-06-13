#include "parameters.hpp"
#include "triangular_mesh.hpp"
//#include "rectangular_mesh.hpp"
//#include "point.hpp"

#include <iostream>



int main(int argc, char **argv)
{
  try
  {
    Parameters param;
    param.read_command_line(argc, argv);

    TriangularMesh tri_mesh;
    tri_mesh.read_msh(param.mesh_filename);

//    const Point2& min = tri_mesh.min_coord();
//    const Point2& max = tri_mesh.max_coord();

//    RectangularMesh rect_mesh;
//    rect_mesh.build(min, max, param.n_rect_elements_x, param.n_rect_elements_z);
//    rect_mesh.assign_material_id(tri_mesh);
//    rect_mesh.write_binary_files(param.properties_filename);
  }
  catch(int)
  {
    return 1;
  }
  catch(const std::exception &e)
  {
    std::cerr << "\n\n" << e.what() << std::endl;
    return 2;
  }
  catch(...)
  {
    std::cerr << "\n\nUnknown exception has been thrown!" << std::endl;
    return 3;
  }
  
  return 0;
}
