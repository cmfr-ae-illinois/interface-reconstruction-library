#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include "stencil_rotator.h"

#include <vtkSmartPointer.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkAppendFilter.h>
#include <vtkCellLocator.h>
#include <vtkPolygon.h>
#include <vtkCellCenters.h>

#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <stdexcept>

namespace IRL {

void classify_simulation(IRL::Classifier& classifier, const std::string& filenameNGA, const std::string& filenamePlic,
    int cannonicalize_symmetries = 0, int include_Moments = 1, bool include_Surface_Area = false, bool include_Eigenvalues = false, float noise_stddev = 0.0f, float epsilon_connectivity = 1e-12f, std::vector<int>* savedClasses = nullptr, int downsample_factor = 2) 
    {
    auto reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(filenameNGA.c_str());
    reader->Update();

    vtkRectilinearGrid* grid = reader->GetOutput();
    if (!grid) {
        throw std::runtime_error("Failed to read rectilinear grid from file: " + filenameNGA);
    }

    vtkCellData* cellData = grid->GetCellData();
    if (!cellData) {
        throw std::runtime_error("Rectilinear grid has no cell data: " + filenameNGA);
    }

    vtkDataArray* vofArray =
        vtkDataArray::SafeDownCast(cellData->GetArray("VOF"));

    vtkDataArray* liquidBaryArray =
        vtkDataArray::SafeDownCast(cellData->GetArray("liquid_bary"));

    vtkDataArray* surfAreaArray =
        vtkDataArray::SafeDownCast(cellData->GetArray("surf_area"));

    if (!vofArray) {
        throw std::runtime_error("Missing required cell array 'VOF' in file: " + filenameNGA);
    }

    if (!liquidBaryArray) {
        throw std::runtime_error("Missing required cell array 'liquid_bary' in file: " + filenameNGA);
    }

    if (liquidBaryArray->GetNumberOfComponents() != 3) {
        throw std::runtime_error("'liquid_bary' must have 3 components.");
    }

    int dims[3];
    grid->GetDimensions(dims);  // dims are for points, so cells are dims-1
    int cellDims[3] = {dims[0] - 1, dims[1] - 1, dims[2] - 1};


    //DEbug, print first values of the barycenter
    int counter = 0;
    for (int i = 0; i < cellDims[0] * cellDims[1] * cellDims[2]; ++i) {
        double alpha = vofArray->GetComponent(i, 0);
        if (alpha > 1e-12 && alpha < 1.0 - 1e-12) {
            std::cout << "First 10 relevant liquid barycenter values:" << std::endl;
            double bary[3];
            liquidBaryArray->GetTuple(i, bary);
            std::cout << "Cell " << i << ": vfrac: " << alpha << " bary: (" << bary[0] << ", " << bary[1] << ", " << bary[2] << ")" << std::endl;
            counter ++;
        }
        if (counter >= 20) {
            break;
        }
    }
    counter = 0;
    // loop over all barycenters and see if there are any not equal to 0
    int non_zero_barycenters = 0;
    for (int i = 0; i < liquidBaryArray->GetNumberOfTuples(); ++i) {
        double alpha = vofArray->GetComponent(i, 0);
        double bary[3];
        liquidBaryArray->GetTuple(i, bary);
        if ((std::abs(bary[0]) > 1e-12 || std::abs(bary[1]) > 1e-12 || std::abs(bary[2]) > 1e-12) && alpha > 1e-12 && alpha < 1.0 - 1e-12) {
            non_zero_barycenters++;
            if (counter < 20) {
                std::cout << "Non-zero barycenter found at cell " << i << ": vfrac: " << alpha << " bary: (" << bary[0] << ", " << bary[1] << ", " << bary[2] << ")" << std::endl;
            }
            counter++;
        }
    }

    std::cout << "Number of non-zero barycenters in grid: " << non_zero_barycenters << std::endl;

    int downsample[3] = {downsample_factor, downsample_factor, downsample_factor};
    int newCellDims[3] = {
        cellDims[0] / downsample[0],
        cellDims[1] / downsample[1],
        cellDims[2] / downsample[2]
    };
    int newDims[3] = {
        newCellDims[0] + 1,
        newCellDims[1] + 1,
        newCellDims[2] + 1
    };

    // Read PLIC file
    auto plic_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    plic_reader->SetFileName(filenamePlic.c_str());
    plic_reader->Update();

    vtkUnstructuredGrid* plic_grid = plic_reader->GetOutput();
    if (!plic_grid) {
        throw std::runtime_error("Failed to read PLIC unstructured grid from file: " + filenamePlic);
    }
    
    // Cell cell locator and plic centroids
    auto locator_finegrid = vtkSmartPointer<vtkCellLocator>::New();
    locator_finegrid->SetDataSet(grid);
    locator_finegrid->BuildLocator();
    auto grid_cell_center_filter = vtkSmartPointer<vtkCellCenters>::New();
    grid_cell_center_filter->SetInputData(grid);
    grid_cell_center_filter->Update();
    auto grid_cell_centers = grid_cell_center_filter->GetOutput()->GetPoints();
    auto cell_center_filter = vtkSmartPointer<vtkCellCenters>::New();
    cell_center_filter->SetInputData(plic_grid);
    cell_center_filter->Update();
    auto cell_centers = cell_center_filter->GetOutput()->GetPoints();
    
    // Use liquid barycenter directly from the VTR file
    vtkDataArray* liquid_barycenter = liquidBaryArray;

    // Create new coordinate arrays
    auto newXCoords = vtkSmartPointer<vtkDoubleArray>::New();
    auto newYCoords = vtkSmartPointer<vtkDoubleArray>::New();
    auto newZCoords = vtkSmartPointer<vtkDoubleArray>::New();

    auto oldX = grid->GetXCoordinates();
    auto oldY = grid->GetYCoordinates();
    auto oldZ = grid->GetZCoordinates();

    for (int i = 0; i < newDims[0]; ++i)
        newXCoords->InsertNextValue(oldX->GetComponent(i * downsample[0], 0));
    for (int j = 0; j < newDims[1]; ++j)
        newYCoords->InsertNextValue(oldY->GetComponent(j * downsample[1], 0));
    for (int k = 0; k < newDims[2]; ++k)
        newZCoords->InsertNextValue(oldZ->GetComponent(k * downsample[2], 0));

    // Create new grid
    auto downsampledGrid = vtkSmartPointer<vtkRectilinearGrid>::New();
    downsampledGrid->SetDimensions(newDims);
    downsampledGrid->SetXCoordinates(newXCoords);
    downsampledGrid->SetYCoordinates(newYCoords);
    downsampledGrid->SetZCoordinates(newZCoords);

    // Create new VOF array
    auto newVOF = vtkSmartPointer<vtkDoubleArray>::New();
    newVOF->SetName("VOF");
    newVOF->SetNumberOfComponents(1);
    newVOF->SetNumberOfTuples(newCellDims[0] * newCellDims[1] * newCellDims[2]);

    // Create new Liquid_Barycenter array
    auto newLiquidBarycenter = vtkSmartPointer<vtkDoubleArray>::New();
    newLiquidBarycenter->SetName("LiquidBarycenter");
    newLiquidBarycenter->SetNumberOfComponents(3);
    newLiquidBarycenter->SetNumberOfTuples(newCellDims[0] * newCellDims[1] * newCellDims[2]);
    /*
    auto newFirstMoment = vtkSmartPointer<vtkDoubleArray>::New();
    newFirstMoment->SetName("M1");
    newFirstMoment->SetNumberOfComponents(3);
    newFirstMoment->SetNumberOfTuples(newCellDims[0] * newCellDims[1] * newCellDims[2]);
    */
    
    // Downsample cell data
    int idx = 0;
    for (int K = 0; K < newCellDims[2]; ++K) {
        for (int J = 0; J < newCellDims[1]; ++J) {
            for (int I = 0; I < newCellDims[0]; ++I) {

                double sum_alpha = 0.0;
                double sum_m1x = 0.0, sum_m1y = 0.0, sum_m1z = 0.0;
                int count = 0;

                // Aggregate fine cells within this coarse block
                for (int dK = 0; dK < downsample[2]; ++dK) {
                    for (int dJ = 0; dJ < downsample[1]; ++dJ) {
                        for (int dI = 0; dI < downsample[0]; ++dI) {
                            int oldI = I * downsample[0] + dI;
                            int oldJ = J * downsample[1] + dJ;
                            int oldK = K * downsample[2] + dK;

                            int oldIdx = oldK * (cellDims[0] * cellDims[1]) + oldJ * cellDims[0] + oldI;

                            // Fine-grid volume fraction and barycenter
                            double alpha = vofArray->GetComponent(oldIdx, 0);
                            double bary[3];
                            liquid_barycenter->GetTuple(oldIdx, bary);

                            // Accumulate weighted first moments
                            sum_alpha += alpha;
                            sum_m1x += alpha * bary[0];
                            sum_m1y += alpha * bary[1];
                            sum_m1z += alpha * bary[2];

                            ++count;
                        }
                    }
                }

                // Compute coarse cell barycenter (physical coordinates)
                double cx = 0.0, cy = 0.0, cz = 0.0;
                if (sum_alpha > 1e-16) {
                    cx = sum_m1x / sum_alpha;
                    cy = sum_m1y / sum_alpha;
                    cz = sum_m1z / sum_alpha;
                }

                newVOF->SetValue(idx, sum_alpha / count);
                newLiquidBarycenter->SetTuple3(idx, cx, cy, cz);

                ++idx;
            }
        }
    }
    // Downsample cell data and convert to coarse-cell unit coordinates
    /*
    int idx = 0;
    for (int K = 0; K < newCellDims[2]; ++K) {
        for (int J = 0; J < newCellDims[1]; ++J) {
            for (int I = 0; I < newCellDims[0]; ++I) {

                double sum_liquid_volume = 0.0;
                double sum_cell_volume = 0.0;

                // Physical first moment over the whole coarse cell
                double sum_m1x_phys = 0.0;
                double sum_m1y_phys = 0.0;
                double sum_m1z_phys = 0.0;

                // Physical surface area over the whole coarse cell
                //double sum_area_phys = 0.0;

                // Coarse-cell bounds
                const double x0c = newXCoords->GetValue(I);
                const double x1c = newXCoords->GetValue(I + 1);
                const double y0c = newYCoords->GetValue(J);
                const double y1c = newYCoords->GetValue(J + 1);
                const double z0c = newZCoords->GetValue(K);
                const double z1c = newZCoords->GetValue(K + 1);

                const double dxc = x1c - x0c;
                const double dyc = y1c - y0c;
                const double dzc = z1c - z0c;

                const double cxc = 0.5 * (x0c + x1c);
                const double cyc = 0.5 * (y0c + y1c);
                const double czc = 0.5 * (z0c + z1c);

                const double coarse_cell_volume = dxc * dyc * dzc;

                // Aggregate fine cells within this coarse block
                for (int dK = 0; dK < downsample[2]; ++dK) {
                    for (int dJ = 0; dJ < downsample[1]; ++dJ) {
                        for (int dI = 0; dI < downsample[0]; ++dI) {
                            const int oldI = I * downsample[0] + dI;
                            const int oldJ = J * downsample[1] + dJ;
                            const int oldK = K * downsample[2] + dK;

                            const int oldIdx =
                                oldK * (cellDims[0] * cellDims[1])
                            + oldJ * cellDims[0]
                            + oldI;

                            // Fine-cell bounds
                            const double x0f = oldX->GetComponent(oldI,     0);
                            const double x1f = oldX->GetComponent(oldI + 1, 0);
                            const double y0f = oldY->GetComponent(oldJ,     0);
                            const double y1f = oldY->GetComponent(oldJ + 1, 0);
                            const double z0f = oldZ->GetComponent(oldK,     0);
                            const double z1f = oldZ->GetComponent(oldK + 1, 0);

                            const double fine_cell_volume =
                                (x1f - x0f) * (y1f - y0f) * (z1f - z0f);

                            const double alpha = vofArray->GetComponent(oldIdx, 0);
                            const double fine_liquid_volume = alpha * fine_cell_volume;

                            double bary[3];
                            liquid_barycenter->GetTuple(oldIdx, bary);

                            double m1_phys[3] = {
                                alpha * bary[0],
                                alpha * bary[1],
                                alpha * bary[2]
                            };

                            //const double area_phys = surface_area->GetComponent(oldIdx, 0);

                            sum_cell_volume += fine_cell_volume;
                            sum_liquid_volume += fine_liquid_volume;

                            // Raw physical first moments are additive.
                            sum_m1x_phys += m1_phys[0];
                            sum_m1y_phys += m1_phys[1];
                            sum_m1z_phys += m1_phys[2];

                            // Physical surface areas are additive.
                            //sum_area_phys += area_phys;
                        }
                    }
                }

                double alpha_coarse = 0.0;
                if (coarse_cell_volume > 1e-16) {
                    alpha_coarse = sum_liquid_volume / coarse_cell_volume;
                }
                // Convert physical first moment to moment in the coarse unit cell
                double m1x_unit = 0.0;
                double m1y_unit = 0.0;
                double m1z_unit = 0.0;

                if (sum_liquid_volume > 1e-16 && coarse_cell_volume > 1e-16) {
                    m1x_unit =
                        (sum_m1x_phys - cxc * sum_liquid_volume)
                        / (dxc * coarse_cell_volume);

                    m1y_unit =
                        (sum_m1y_phys - cyc * sum_liquid_volume)
                        / (dyc * coarse_cell_volume);

                    m1z_unit =
                        (sum_m1z_phys - czc * sum_liquid_volume)
                        / (dzc * coarse_cell_volume);
                }

                // Convert physical surface area to unit-cell surface area
                
                double area_unit = 0.0;
                if (coarse_cell_volume > 1e-16) {
                    area_unit =
                        sum_area_phys / std::pow(coarse_cell_volume, 2.0 / 3.0);
                }

                newVOF->SetValue(idx, alpha_coarse);
                newFirstMoment->SetTuple3(idx, m1x_unit, m1y_unit, m1z_unit);
                //newSurfaceArea->SetValue(idx, area_unit);

                ++idx;
            }
        }
    }*/

    // Add arrays to grid
    downsampledGrid->GetCellData()->AddArray(newVOF);
    downsampledGrid->GetCellData()->AddArray(newLiquidBarycenter);
    //downsampledGrid->GetCellData()->AddArray(newFirstMoment);

    // Convert point dims to cell dims (cells = points - 1)
    int nx = newCellDims[0];
    int ny = newCellDims[1];
    int nz = newCellDims[2];

    std::cout << "Downsampled grid cell dimensions:" << std::endl;
    std::cout << "  nx = " << nx << std::endl;
    std::cout << "  ny = " << ny << std::endl;
    std::cout << "  nz = " << nz << std::endl;

    // Create new cell data array for interface type
    auto interface_type = vtkSmartPointer<vtkIntArray>::New();
    interface_type->SetName("InterfaceType");
    interface_type->SetNumberOfComponents(1);
    interface_type->SetNumberOfTuples(grid->GetNumberOfCells());
    for (int i = 0; i < grid->GetNumberOfCells(); i++) {
        interface_type->SetValue(i, -1);
    }
    // Create certainty array
    auto certainty = vtkSmartPointer<vtkFloatArray>::New();
    certainty->SetName("Certainty");
    certainty->SetNumberOfComponents(1);
    certainty->SetNumberOfTuples(grid->GetNumberOfCells());
    for (int i = 0; i < grid->GetNumberOfCells(); i++) {
        certainty->SetValue(i, 0.0);
    }

    // Create new cell data array for interface type on PLIC surface
    auto plic_interface_type = vtkSmartPointer<vtkIntArray>::New();
    plic_interface_type->SetName("InterfaceType");
    plic_interface_type->SetNumberOfComponents(1);
    plic_interface_type->SetNumberOfTuples(plic_grid->GetNumberOfCells());
    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        plic_interface_type->SetValue(i, -1);
    }
    // Create certainty array on PLIC surface
    auto plic_certainty = vtkSmartPointer<vtkFloatArray>::New();
    plic_certainty->SetName("Certainty");
    plic_certainty->SetNumberOfComponents(1);
    plic_certainty->SetNumberOfTuples(plic_grid->GetNumberOfCells());
    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        plic_certainty->SetValue(i, 0.0);
    }

    int stencil_size_reader = classifier.getStencilSize();
    int half = stencil_size_reader / 2;

    //double epsilon = 1e-10; // threshold for considering a cell as "interface" (0 < vfrac < 1) vs "empty" (vfrac=0) or "full" (vfrac=1)

    int no_filled_cells = 0;
    int no_paraboloids = 0;
    int no_cylinders = 0;
    int no_spheres = 0;
    int no_sheets = 0;
    int no_cut_sheets = 0;
    int no_ligament_tip = 0;
    int no_isolated = 0; // for pre-processing

    double start_time = static_cast<double>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    // Loop through interior cells (excluding half-cell boundary)
    for (int i = half; i < nx - half; ++i) {
        for (int j = half; j < ny - half; ++j) {
            for (int k = half; k < nz - half; ++k) {

                vtkIdType centerCellId = i + j * nx + k * nx * ny;

                // Skip if no interface
                double vof_center = newVOF->GetComponent(centerCellId, 0);
                if (vof_center <= epsilon_connectivity || vof_center >= 1.0 - epsilon_connectivity) {
                    continue;
                }

                //int center_sid = structure_id[centerCellId]; // for pre-processing Herrmann style

                no_filled_cells++;

                // Compute physical center and size of central cell
                double x0_c = newXCoords->GetValue(i);
                double x1_c = newXCoords->GetValue(i + 1);
                double y0_c = newYCoords->GetValue(j);
                double y1_c = newYCoords->GetValue(j + 1);
                double z0_c = newZCoords->GetValue(k);
                double z1_c = newZCoords->GetValue(k + 1);

                double dx_c = x1_c - x0_c;
                double dy_c = y1_c - y0_c;
                double dz_c = z1_c - z0_c;
                double cx_c = 0.5 * (x0_c + x1_c);
                double cy_c = 0.5 * (y0_c + y1_c);
                double cz_c = 0.5 * (z0_c + z1_c);

                // Allocate stencil containers
                std::vector<std::vector<std::vector<double>>> vfrac(
                    stencil_size_reader,
                    std::vector<std::vector<double>>(
                        stencil_size_reader,
                        std::vector<double>(stencil_size_reader, 0.0)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoment(
                    stencil_size_reader,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size_reader,
                        std::vector<Eigen::Vector3d>(stencil_size_reader, Eigen::Vector3d::Zero())));

                // Loop over neighborhood cells in stencil
                
                for (int di = -half; di <= half; ++di) {
                    for (int dj = -half; dj <= half; ++dj) {
                        for (int dk = -half; dk <= half; ++dk) {

                            int ni = i + di;
                            int nj = j + dj;
                            int nk = k + dk;

                            vtkIdType cellId = ni + nj * nx + nk * nx * ny;
                            if (cellId < 0 || cellId >= newVOF->GetNumberOfTuples())
                                continue;

                            // Get local cell size and center
                            double x0 = newXCoords->GetValue(ni);
                            double x1 = newXCoords->GetValue(ni + 1);
                            double y0 = newYCoords->GetValue(nj);
                            double y1 = newYCoords->GetValue(nj + 1);
                            double z0 = newZCoords->GetValue(nk);
                            double z1 = newZCoords->GetValue(nk + 1);

                            double dx = x1 - x0;
                            double dy = y1 - y0;
                            double dz = z1 - z0;
                            double cx = 0.5 * (x0 + x1);
                            double cy = 0.5 * (y0 + y1);
                            double cz = 0.5 * (z0 + z1);

                            // If center has no structure id (shouldn't happen if vof_center is interface), skip cell
                            // And if neighbor is not part of the same structure, ignore it for classification.
                            double alpha = newVOF->GetComponent(cellId, 0);
                            double bary[3];
                            newLiquidBarycenter->GetTuple(cellId, bary);

                            // Step 1: normalize position of cell center relative to central cell
                            // The cell centers form a grid in the training coordinate system
                            double center_rel_x = (cx - cx_c) / dx_c;
                            double center_rel_y = (cy - cy_c) / dy_c;
                            double center_rel_z = (cz - cz_c) / dz_c;

                            // Step 2: normalize barycenter relative to its own cell center
                            double local_rel_x = (bary[0] - cx_c) / dx;
                            double local_rel_y = (bary[1] - cy_c) / dy;
                            double local_rel_z = (bary[2] - cz_c) / dz;

                            // Step 3: combine both transforms
                            // The total normalized barycenter position (in training coordinate frame)
                            double bx_unit = local_rel_x;
                            double by_unit = local_rel_y;
                            double bz_unit = local_rel_z;

                            // --- Compute first moment in this normalized unit frame ---
                            Eigen::Vector3d m1(alpha * bx_unit, alpha * by_unit, alpha * bz_unit);

                            int si = di + half;
                            int sj = dj + half;
                            int sk = dk + half;

                            vfrac[si][sj][sk] = alpha;
                            firstMoment[si][sj][sk] = m1;
                        }
                    }
                }/*
                for (int di = -half; di <= half; ++di) {
                    for (int dj = -half; dj <= half; ++dj) {
                        for (int dk = -half; dk <= half; ++dk) {

                            int ni = i + di;
                            int nj = j + dj;
                            int nk = k + dk;

                            vtkIdType cellId = ni + nj * nx + nk * nx * ny;
                            if (cellId < 0 || cellId >= newVOF->GetNumberOfTuples())
                                continue;

                            double alpha = newVOF->GetComponent(cellId, 0);
                            double m1[3] = {0.0, 0.0, 0.0};
                            newFirstMoment->GetTuple(cellId, m1);
                            //double area = newSurfaceArea->GetComponent(cellId, 0);

                            int si = di + half;
                            int sj = dj + half;
                            int sk = dk + half;

                            vfrac[si][sj][sk] = alpha;
                            firstMoment[si][sj][sk] = Eigen::Vector3d(m1[0], m1[1], m1[2]);
                            //surfaceArea[si][sj][sk] = area;
                        }
                    }
                }*/

                // Flatten stencil into 1D vector
                std::vector<double> flattened_state;
                for (int si = 0; si < stencil_size_reader; ++si) {
                    for (int sj = 0; sj < stencil_size_reader; ++sj) {
                        for (int sk = 0; sk < stencil_size_reader; ++sk) {
                            if (include_Moments >= 0) {
                                flattened_state.push_back(vfrac[si][sj][sk]);
                            }
                            if (include_Moments >= 1) {
                                flattened_state.push_back(firstMoment[si][sj][sk].x());
                                flattened_state.push_back(firstMoment[si][sj][sk].y());
                                flattened_state.push_back(firstMoment[si][sj][sk].z());
                            }
                        }
                    }
                }

                // Make flattened_state a float vector
                std::vector<float> flattened_state_float(flattened_state.begin(), flattened_state.end());

                if (include_Moments >= 2) {
                    Eigen::Matrix3d secondMoment = IRL::compute2ndMoment(flattened_state_float, stencil_size_reader, 1); // from_ith_moment=1
                    flattened_state_float.push_back(static_cast<float>(secondMoment(0, 0))); // xx
                    flattened_state_float.push_back(static_cast<float>(secondMoment(1, 1))); // yy
                    flattened_state_float.push_back(static_cast<float>(secondMoment(2, 2))); // zz
                    flattened_state_float.push_back(static_cast<float>(secondMoment(0, 1))); // xy
                    flattened_state_float.push_back(static_cast<float>(secondMoment(0, 2))); // xz
                    flattened_state_float.push_back(static_cast<float>(secondMoment(1, 2))); // yz
                }

                if (include_Eigenvalues) {
                    IRL::appendInertiaEigenvalues(flattened_state_float, stencil_size_reader, 1);
                }

                //Preprocess stencil
                
                IRL::preprocess_stencil(flattened_state_float,
                        stencil_size_reader,
                        cannonicalize_symmetries,
                        include_Moments,
                        include_Surface_Area,
                        include_Eigenvalues,
                        noise_stddev,
                        epsilon_connectivity);

                // Classify
                std::vector<float> out_probs;
                int predicted_class = classifier.classify(flattened_state_float, &out_probs);
                float max_prob = *std::max_element(out_probs.begin(), out_probs.end());
                certainty->SetValue(centerCellId, max_prob);

                switch (predicted_class) {
                    case 0:
                        no_paraboloids++;
                        interface_type->SetValue(centerCellId, 0);
                        break;
                    case 1:
                        no_cylinders++;
                        interface_type->SetValue(centerCellId, 1);
                        break;
                    case 2:
                        no_spheres++;
                        interface_type->SetValue(centerCellId, 2);
                        break;
                    case 3:
                        no_sheets++;
                        interface_type->SetValue(centerCellId, 3);
                        break;
                    case 4:
                        no_ligament_tip++;
                        interface_type->SetValue(centerCellId, 4);
                        break;
                    case 5:
                        no_cut_sheets++;
                        interface_type->SetValue(centerCellId, 5);

                        // Debugging: print stencil and moments for first few class 5 predictions
                        static int printed_class5 = 0;
                        if (printed_class5 < 5) {
                            ++printed_class5;

                            std::vector<CellData> dbg_stencil;
                            SecondMoments dbg_I{};
                            SecondMoments* dbg_Ip = (include_Moments >= 2) ? &dbg_I : nullptr;
                            Eigenvalues dbg_eig{};
                            Eigenvalues* dbg_eigp = include_Eigenvalues ? &dbg_eig : nullptr;

                            IRL::unpackStencil(
                                flattened_state_float,
                                dbg_stencil,
                                dbg_Ip,
                                dbg_eigp,
                                stencil_size_reader,
                                include_Moments,
                                include_Surface_Area,
                                include_Eigenvalues
                            );

                            std::cout << "\n============================================================\n";
                            std::cout << "CLASS 5 STENCIL DEBUG #" << printed_class5 << "\n";
                            std::cout << "domain cell = [" << i << "," << j << "," << k << "]\n";
                            std::cout << "centerCellId = " << centerCellId << "\n";
                            std::cout << "vof_center = " << vof_center << "\n";
                            std::cout << "certainty = " << max_prob << "\n";
                            std::cout << "stencil_size = " << stencil_size_reader << "\n";
                            std::cout << "Include moments: " << include_Moments << "\n";
                            std::cout << "Include surface area: " << include_Surface_Area << "\n";
                            std::cout << "Include eigenvalues: " << include_Eigenvalues << "\n";

                            for (int sk = 0; sk < stencil_size_reader; ++sk) {
                                std::cout << "\n----- z-slice k = " << sk << " -----\n";

                                for (int sj = 0; sj < stencil_size_reader; ++sj) {
                                    for (int si = 0; si < stencil_size_reader; ++si) {
                                        const CellData& c = dbg_stencil[cellIndex(si, sj, sk, stencil_size_reader)];

                                        float cx = 0.0f, cy = 0.0f, cz = 0.0f;
                                        float surfArea = 0.0f;
                                        if (include_Moments >= 1 && c.vfrac > 1.0e-12f) {
                                            cx = c.mx / c.vfrac;
                                            cy = c.my / c.vfrac;
                                            cz = c.mz / c.vfrac;
                                            if (include_Surface_Area) {
                                                surfArea = c.area;
                                            }
                                        }

                                        std::cout << "[" << si << "," << sj << "," << sk << "] "
                                                << c.vfrac << " "
                                                << "(" << cx << ", " << cy << ", " << cz << ")" << " A=" << surfArea;

                                        if (si < stencil_size_reader - 1) {
                                            std::cout << "    ";
                                        }
                                    }
                                    std::cout << "\n";
                                }
                            }

                            if (include_Moments >= 2) {
                                std::cout << "Second moments: "
                                        << dbg_I.Ixx << ", "
                                        << dbg_I.Iyy << ", "
                                        << dbg_I.Izz << ", "
                                        << dbg_I.Ixy << ", "
                                        << dbg_I.Ixz << ", "
                                        << dbg_I.Iyz << "\n";
                            }

                            if (include_Eigenvalues) {
                                std::cout << "Eigenvalues: "
                                        << dbg_eig.lambda1 << ", "
                                        << dbg_eig.lambda2 << ", "
                                        << dbg_eig.lambda3 << "\n";
                            }

                            std::cout << "============================================================\n";
                        }

                        break;
                    default:
                        std::cerr << "Warning: unknown predicted_class = " << predicted_class << std::endl;
                        break;
                }

                if (savedClasses) {
                    savedClasses->push_back(predicted_class);
                }
            }
        }
    }

    double end_time = static_cast<double>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::cout << "Classification time (s): " << (end_time - start_time) / 1e9 << std::endl;
  

    std::cout << "\n=== Classification Summary ===" << std::endl;
    std::cout << "Total cells with interface:   " << no_filled_cells << std::endl;
    std::cout << "Paraboloids:          " << no_paraboloids << std::endl;
    std::cout << "Cylinders:            " << no_cylinders << std::endl;
    std::cout << "Spheres:              " << no_spheres << std::endl;
    std::cout << "Sheets:               " << no_sheets << std::endl;
    std::cout << "Ligament tips:        " << no_ligament_tip << std::endl;
    std::cout << "Sheet Edges:          " << no_cut_sheets << std::endl;

    // Write Grid file
    grid->GetCellData()->AddArray(interface_type);
    grid->GetCellData()->AddArray(certainty);
    auto grid_writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
    grid_writer->SetFileName("grid.vtr");
    grid_writer->SetInputData(downsampledGrid);
    grid_writer->Write();

    // Cell locator
    auto locator = vtkSmartPointer<vtkCellLocator>::New();
    locator->SetDataSet(downsampledGrid);
    locator->BuildLocator();
    
    // Convert interface_type to PLIC array
    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        auto type = interface_type->GetValue(locator->FindCell(cell_centers->GetPoint(i)));
        plic_interface_type->SetValue(i, type);
        auto cert = certainty->GetValue(locator->FindCell(cell_centers->GetPoint(i)));
        plic_certainty->SetValue(i, cert);
    }

    // Write PLIC file
    plic_grid->GetCellData()->AddArray(plic_interface_type);
    plic_grid->GetCellData()->AddArray(plic_certainty);
    auto plic_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    plic_writer->SetFileName("plic.vtu");
    plic_writer->SetInputData(plic_grid);
    plic_writer->Write();
}
} // namespace IRL