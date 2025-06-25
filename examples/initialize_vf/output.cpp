#include <iostream>

#include "examples/initialize_vf/functions.h"

#include <fstream>
#include <iomanip>
#include <string>

// volume fraction data
void write_vtr(const std::string& filepath, const Data<double>& vf, const BasicMesh& mesh) {

    std::ofstream out(filepath);
    if (!out.is_open()) {
    std::cerr << "Failed to open " << filepath << " for writing.\n";
    return;
}
    out << std::fixed << std::setprecision(6);

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <RectilinearGrid WholeExtent=\"0 " << mesh.imax() + 1 << " 0 " << mesh.jmax() + 1 << " 0 0\">\n";
    out << "    <Piece Extent=\"0 " << mesh.imax() + 1 << " 0 " << mesh.jmax() + 1 << " 0 0\">\n";

    out << "      <CellData Scalars=\"volumeFraction\">\n";
    out << "        <DataArray type=\"Float32\" Name=\"volumeFraction\" format=\"ascii\">\n";
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
            out << "          " << vf(i, j) << "\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </CellData>\n";
    
    out << "      <Coordinates>\n";
    out << "        <DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">\n";
    for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
        out << "          " << mesh.x(i) << "\n";
    }
    out << "        </DataArray>\n";

    out << "        <DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
        out << "          " << mesh.y(j) << "\n";
    }
    out << "        </DataArray>\n";

    out << "        <DataArray type=\"Float32\" Name=\"Z\" format=\"ascii\">\n";
    out << "          0.0\n"; 
    out << "        </DataArray>\n";
    out << "      </Coordinates>\n";

    out << "    </Piece>\n";
    out << "  </RectilinearGrid>\n";
    out << "</VTKFile>\n";
    out.close();
}

// nested grid
void write_vtu(const std::string& filepath, const std::vector<const Cell*>& cells) {

  std::ofstream out(filepath);
  if (!out.is_open()) {
    std::cerr << "Failed to open " << filepath << " for writing.\n";
    return;
}

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << 4 * cells.size()
      << "\" NumberOfCells=\"" << cells.size() << "\">\n";

  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const Cell* c : cells) {
    double x = c->x, y = c->y, dx = c->dx;
    out << x        << " " << y        << " 0\n";
    out << x + dx   << " " << y        << " 0\n";
    out << x + dx   << " " << y + dx   << " 0\n";
    out << x        << " " << y + dx   << " 0\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  out << "      <Cells>\n";

  // Connectivity
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int i = 0; i < cells.size(); ++i) {
    int base = 4 * i;
    out << base << " " << base + 1 << " " << base + 2 << " " << base + 3 << "\n";
  }
  out << "        </DataArray>\n";

  // Offsets
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int i = 1; i <= cells.size(); ++i)
    out << 4 * i << "\n";
  out << "        </DataArray>\n";

  // Types (VTK_QUAD = 9)
  out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < cells.size(); ++i)
    out << "9\n";
  out << "        </DataArray>\n";

  out << "      </Cells>\n";

  // === No scalar data ===
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}

