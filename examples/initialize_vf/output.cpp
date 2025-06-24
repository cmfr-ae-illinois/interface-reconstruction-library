#include <iostream>

#include "examples/initialize_vf/functions.h"

#include <fstream>
#include <iomanip>
#include <string>

void write_vtr(const std::string& filepath, const Data<double>& vf, const BasicMesh& mesh) {
    int imin = mesh.imin(), imax = mesh.imax();
    int jmin = mesh.jmin(), jmax = mesh.jmax();
    int Nx = imax - imin;
    int Ny = jmax - jmin;

    std::ofstream out(filepath);
    out << std::fixed << std::setprecision(6);

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <RectilinearGrid WholeExtent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n";
    out << "    <Piece Extent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n";

    out << "      <CellData Scalars=\"volumeFraction\">\n";
    out << "        <DataArray type=\"Float32\" Name=\"volumeFraction\" format=\"ascii\">\n";
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            out << "          " << vf(i, j) << "\n";
        }
    }

    out << "        </DataArray>\n";
    out << "      </CellData>\n";

    
    out << "      <Coordinates>\n";

    
    out << "        <DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">\n";
    for (int i = imin; i <= imax; ++i) {
        out << "          " << mesh.x(i) << "\n";
    }
    out << "        </DataArray>\n";

    
    out << "        <DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
    for (int j = jmin; j <= jmax; ++j) {
        out << "          " << mesh.y(j) << "\n";
    }
    out << "        </DataArray>\n";

    
    out << "        <DataArray type=\"Float32\" Name=\"Z\" format=\"ascii\">\n";
    out << "          0.0\n";  // only one Z-slice
    out << "        </DataArray>\n";

    out << "      </Coordinates>\n";
    out << "    </Piece>\n";
    out << "  </RectilinearGrid>\n";
    out << "</VTKFile>\n";

    out.close();
}