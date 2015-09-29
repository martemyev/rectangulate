#include "parameters.hpp"
#include "rectangular_mesh.hpp"
#include "triangular_mesh.hpp"
#include "utilities.hpp"

#include <iostream>



int main(int argc, char **argv)
{
  auto t_begin = get_wall_time();

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

    if (param.assign_cells)
    {
      rect_mesh.assign_material_id_in_cells(tri_mesh);
      std::vector<std::string> filenames_out;
      rect_mesh.write_binary_files_in_cells(param.properties_filename,
                                            filenames_out);
      convert_in_cells_to_xz(filenames_out, nnx, nnz);
    }

    if (param.assign_nodes)
    {
      rect_mesh.assign_material_id_at_nodes(tri_mesh);
      rect_mesh.write_ASCII_files_at_nodes(param.properties_filename);
    }
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

  std::cout.precision(8);
  std::cout << "\nTOTAL TIME\n";
  std::cout << "wall time = " << get_wall_time() - t_begin << " seconds"
            << std::endl;
}

