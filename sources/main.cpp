#include "parameters.hpp"
#include "triangular_mesh.hpp"
#include "rectangular_mesh.hpp"

#include <iostream>



int main(int argc, char **argv)
{
  try
  {
    Parameters param;
    param.read_command_line(argc, argv);

    TriangularMesh tri_mesh;
    tri_mesh.read_msh(param.mesh_filename);

    const Point2 min = tri_mesh.min_coord();
    const Point2 max = tri_mesh.max_coord();

    const int nnx = (max.x() - min.x()) / param.h_rect_x;
    const int nnz = (max.z() - min.z()) / param.h_rect_z;

    std::cout << "nnx = " << nnx << "\nnnz = " << nnz << std::endl;

    RectangularMesh rect_mesh(min, max, nnx, nnz);
    rect_mesh.build();
    rect_mesh.assign_material_id(tri_mesh);

    std::vector<std::string> out_filenames;
    rect_mesh.write_binary_files(param.properties_filename, out_filenames);

    convert_to_xz(out_filenames, nnx, nnz);
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
