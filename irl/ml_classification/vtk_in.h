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

namespace IRL {

void classify_simulation(IRL::Classifier& classifier, const std::string& filenameNGA, const std::string& filenamePlic,
    int cannonicalize_symmetries = 0, int include_Moments = 1, bool include_Eigenvalues = false, float noise_stddev = 0.0f, float epsilon_connectivity = 1e-10f, std::vector<int>* savedClasses = nullptr, int downsample_factor = 2) 
    {
    auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    reader->SetCaseFileName(filenameNGA.c_str());
    //reader->SetCaseFileName("/home/quirin/mlcfd/Repositories/jet/nga.case");
    //reader->SetCaseFileName("/home/quirin/jet/nga.case");
    reader->UpdateInformation();

    vtkInformation* info = reader->GetOutputInformation(0);
    int nSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> timeSteps(nSteps);
    info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeSteps.data());
    reader->SetTimeValue(timeSteps[nSteps - 1]);
    reader->Update();
    auto mb = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    auto block = mb->GetBlock(0);
    auto dataSet = vtkDataSet::SafeDownCast(block);
    vtkCellData* cellData = dataSet->GetCellData();
    vtkDataArray* vofArray = cellData->GetArray("VOF");
    auto grid = vtkRectilinearGrid::SafeDownCast(dataSet);

    int dims[3];
    grid->GetDimensions(dims);  // dims are for points, so cells are dims-1
    int cellDims[3] = {dims[0] - 1, dims[1] - 1, dims[2] - 1};

    int downsample[3] = {downsample_factor, downsample_factor, downsample_factor}; //Change this to {1,1,1} to disable downsampling
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
    auto plic_reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    // plic_reader->SetCaseFileName("/home/quirin/jet/plic.case");
    plic_reader->SetCaseFileName(filenamePlic.c_str());
    plic_reader->UpdateInformation();

    vtkInformation* plic_info = plic_reader->GetOutputInformation(0);
    int plic_nSteps = plic_info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> plic_timeSteps(plic_nSteps);
    plic_info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), plic_timeSteps.data());
    plic_reader->SetTimeValue(plic_timeSteps[plic_nSteps - 1]);
    plic_reader->Update();
    auto plic_mb = vtkMultiBlockDataSet::SafeDownCast(plic_reader->GetOutput());
    auto plic_block_merger = vtkSmartPointer<vtkAppendFilter>::New();
    for (int i = 0; i < plic_mb->GetNumberOfBlocks(); i++) {
        auto plic_block = plic_mb->GetBlock(i);
        if (plic_block) { plic_block_merger->AddInputData(plic_block);}
    }
    plic_block_merger->Update();
    auto plic_grid = vtkUnstructuredGrid::SafeDownCast(plic_block_merger->GetOutput());
    
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
    
    // Regenerate first moments and surface area from PLIC
    auto liquid_barycenter = vtkSmartPointer<vtkDoubleArray>::New();
    liquid_barycenter->SetName("M1"); //WRONG? "Liquid_Barycenter"
    liquid_barycenter->SetNumberOfComponents(3);
    liquid_barycenter->SetNumberOfTuples(grid->GetNumberOfCells());
    //auto surface_area = vtkSmartPointer<vtkDoubleArray>::New();
    //surface_area->SetName("Surface Area");
    //surface_area->SetNumberOfComponents(1);
    //surface_area->SetNumberOfTuples(grid->GetNumberOfCells());
    for (int i = 0; i < grid->GetNumberOfCells(); i++) {
        double center[3];
        grid_cell_centers->GetPoint(i, center);
        liquid_barycenter->SetTuple3(i, center[0], center[1], center[2]);
        //surface_area->SetValue(i, 0.0);
    }
    double max_vfrac_error = 0.0, rms_vfrac_error = 0.0;
    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        auto plic_cell = plic_grid->GetCell(i);
        if (plic_cell->GetCellType() == VTK_POLYGON) {
            // Cast PLIC cell into vtkPolygon
            auto polygon = vtkPolygon::SafeDownCast(plic_cell);
            // Compute PLIC normal and first point coordinates from polygon
            double normal[3], firstpoint[3];
            polygon->ComputeNormal(polygon->GetPoints(), normal);
            polygon->GetPoints()->GetPoint(0, firstpoint);
            // Convert to IRL normal
            const auto irl_normal = IRL::Normal(normal[0], normal[1], normal[2]);
            // Compute plane distance
            const double irl_distance = normal[0] * firstpoint[0] + normal[1] * firstpoint[1] + normal[2] * firstpoint[2];
            // Construct IRL plane matching PLIC
            const auto irl_separator = IRL::PlanarSeparator::fromOnePlane(IRL::Plane(irl_normal, irl_distance));
            // Get local cell ID
            auto cell_id = locator_finegrid->FindCell(cell_centers->GetPoint(i));
            if (cell_id < 0 || cell_id >= grid->GetNumberOfCells()) {
                std::cerr << "Warning: could not find valid cell for PLIC polygon " << i << ", skipping..." << std::endl;
                continue; // skip invalid polygons (outside grid)
            }
            // Construct IRL rectangular cuboid matching local cell
            double bds[6];
            grid->GetCellBounds(cell_id, bds);
            const auto irl_cell = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(bds[0], bds[2], bds[4]), IRL::Pt(bds[1], bds[3], bds[5]));
            // Get volume moments for that cell
            auto moments = IRL::getVolumeMoments<IRL::VolumeMoments>(irl_cell, irl_separator);
            moments.centroid() /= IRL::safelyEpsilon(moments.volume());
            liquid_barycenter->SetTuple3(cell_id, moments.centroid()[0], moments.centroid()[1], moments.centroid()[2]); // was i instead of cell_id
            // Compute vfrac to monitor error
            const double vfrac = moments.volume() / irl_cell.calculateVolume();
            const double vfrac_error = std::fabs(vfrac - vofArray->GetComponent(cell_id, 0));
            max_vfrac_error = std::max(max_vfrac_error, vfrac_error);
            rms_vfrac_error += vfrac_error * vfrac_error;

            // Get surface area from polygon
            const double area = polygon->ComputeArea(); //Actually returns area
            //surface_area->SetValue(cell_id, area);
        }
    }
    rms_vfrac_error = std::sqrt(rms_vfrac_error / static_cast<double>(plic_grid->GetNumberOfCells()));
    std::cout << "MAX reconstruted vfrac error = " << max_vfrac_error << std::endl;
    std::cout << "RMS reconstruted vfrac error = " << rms_vfrac_error << " (should be O(FLOAT_EPSILON) = 1.2e-7)" << std::endl;

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

    //Create new surface area array
    

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

    // Add arrays to grid
    downsampledGrid->GetCellData()->AddArray(newVOF);
    downsampledGrid->GetCellData()->AddArray(newLiquidBarycenter);


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

    double epsilon = 1e-12; // threshold for considering a cell as "interface" (0 < vfrac < 1) vs "empty" (vfrac=0) or "full" (vfrac=1)

    int no_filled_cells = 0;
    int no_paraboloids = 0;
    int no_cylinders = 0;
    int no_spheres = 0;
    int no_sheets = 0;
    int no_cut_sheets = 0;
    int no_ligament_tip = 0;
    int no_isolated = 0; // for pre-processing

    double start_time = static_cast<double>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    // Preprocessing Herrmann2010a below
    /*
    std::vector<int> structure_id(newVOF->GetNumberOfTuples(), -1);

    auto coarseId = [&](int ii, int jj, int kk) -> vtkIdType {
        return ii + jj * nx + kk * nx * ny;
    };

    static const int n6[6][3] = {
        {+1, 0, 0}, {-1, 0, 0},
        { 0,+1, 0}, { 0,-1, 0},
        { 0, 0,+1}, { 0, 0,-1}
    };

    int next_id = 0;

    // Loop over all coarse cells and flood-fill each untagged liquid component
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {

                vtkIdType seed = coarseId(i,j,k);
                if (structure_id[seed] != -1) continue;

                double a_seed = newVOF->GetComponent(seed, 0);
                if (a_seed <= epsilon) continue;   // not in "liquid region" => no structure id

                // Start BFS for this structure
                const int my_id = next_id++;
                structure_id[seed] = my_id;

                std::vector<vtkIdType> q;
                q.reserve(1024);
                q.push_back(seed);

                for (size_t qi = 0; qi < q.size(); ++qi) {
                    vtkIdType cid = q[qi];

                    // Convert cid -> (ci,cj,ck)
                    int ck = static_cast<int>(cid / (nx * ny));
                    int rem = static_cast<int>(cid - ck * nx * ny);
                    int cj = rem / nx;
                    int ci = rem - cj * nx;

                    for (int n = 0; n < 6; ++n) {
                        int ni = ci + n6[n][0];
                        int nj = cj + n6[n][1];
                        int nk = ck + n6[n][2];
                        if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz) continue;

                        vtkIdType nid = coarseId(ni,nj,nk);
                        if (structure_id[nid] != -1) continue;

                        if (newVOF->GetComponent(nid, 0) > epsilon) {
                            structure_id[nid] = my_id;
                            q.push_back(nid);
                        }
                    }
                }
            }
        }
    }

    std::cout << "Global structure identification (coarse grid): "
            << next_id << " structures tagged." << std::endl;
    */

    // Loop through interior cells (excluding half-cell boundary)
    for (int i = half; i < nx - half; ++i) {
        for (int j = half; j < ny - half; ++j) {
            for (int k = half; k < nz - half; ++k) {

                vtkIdType centerCellId = i + j * nx + k * nx * ny;

                // Skip if no interface
                double vof_center = newVOF->GetComponent(centerCellId, 0);
                if (vof_center <= epsilon || vof_center >= 1.0 - epsilon) {
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

                            // Retrieve data for Herrmann2010a preprocessing
                            //int sid = structure_id[cellId];

                            // If center has no structure id (shouldn't happen if vof_center is interface), skip cell
                            // And if neighbor is not part of the same structure, ignore it for classification.
                            double alpha = newVOF->GetComponent(cellId, 0);
                            double bary[3];
                            newLiquidBarycenter->GetTuple(cellId, bary);
                            
                            // Below is Herrmann2010a preprocessing
                            /*
                            if (preProcess) {
                                if (center_sid < 0) {
                                    // no liquid structure attached to center; safest is to skip classification. Should not happen.
                                    continue; 
                                }
                                if (sid != center_sid) {
                                    alpha = 0.0;
                                    bary[0] = cx; bary[1] = cy; bary[2] = cz; // doesn't matter since alpha=0, but keeps physics correct
                                }
                            }
                            */

                            // Step 1: normalize position of cell center relative to central cell
                            // The cell centers form a grid in the training coordinate system
                            double center_rel_x = (cx - cx_c) / dx_c;
                            double center_rel_y = (cy - cy_c) / dy_c;
                            double center_rel_z = (cz - cz_c) / dz_c;

                            // Step 2: normalize barycenter relative to its own cell center
                            double local_rel_x = (bary[0] - cx) / dx;
                            double local_rel_y = (bary[1] - cy) / dy;
                            double local_rel_z = (bary[2] - cz) / dz;

                            // Step 3: combine both transforms
                            // The total normalized barycenter position (in training coordinate frame)
                            double bx_unit = center_rel_x + local_rel_x;
                            double by_unit = center_rel_y + local_rel_y;
                            double bz_unit = center_rel_z + local_rel_z;

                            // --- Compute first moment in this normalized unit frame ---
                            Eigen::Vector3d m1(alpha * bx_unit, alpha * by_unit, alpha * bz_unit);

                            int si = di + half;
                            int sj = dj + half;
                            int sk = dk + half;

                            vfrac[si][sj][sk] = alpha;
                            firstMoment[si][sj][sk] = m1;
                        }
                    }
                }

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
                if (include_Moments >= 2) {
                    Eigen::Matrix3d secondMoment = IRL::compute2ndMoment(flattened_state, stencil_size_reader, /*from_ith_moment=*/1);
                    flattened_state.push_back(secondMoment(0, 0)); // xx
                    flattened_state.push_back(secondMoment(1, 1)); // yy
                    flattened_state.push_back(secondMoment(2, 2)); // zz
                    flattened_state.push_back(secondMoment(0, 1)); // xy
                    flattened_state.push_back(secondMoment(0, 2)); // xz
                    flattened_state.push_back(secondMoment(1, 2)); // yz
                }

                if (include_Eigenvalues) {
                    IRL::appendInertiaEigenvalues(flattened_state, stencil_size_reader, 1);
                }

                //print length of flattened state for debugging
                //std::cout << "Flattened state length: " << flattened_state.size() << std::endl;

                // Make flattened_state a float vector
                std::vector<float> flattened_state_float(flattened_state.begin(), flattened_state.end());

                //Preprocess stencil
                IRL::preprocess_stencil(flattened_state_float, stencil_size_reader, cannonicalize_symmetries, include_Moments, include_Eigenvalues, noise_stddev, epsilon_connectivity);
                

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
                    case 5: {
                        no_cut_sheets++;
                        interface_type->SetValue(centerCellId, 5);

                        /*
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
                                include_Eigenvalues
                            );

                            std::cout << "\n============================================================\n";
                            std::cout << "CLASS 5 STENCIL DEBUG #" << printed_class5 << "\n";
                            std::cout << "domain cell = [" << i << "," << j << "," << k << "]\n";
                            std::cout << "centerCellId = " << centerCellId << "\n";
                            std::cout << "vof_center = " << vof_center << "\n";
                            std::cout << "certainty = " << max_prob << "\n";
                            std::cout << "stencil_size = " << stencil_size_reader << "\n";

                            for (int sk = 0; sk < stencil_size_reader; ++sk) {
                                std::cout << "\n----- z-slice k = " << sk << " -----\n";

                                for (int sj = 0; sj < stencil_size_reader; ++sj) {
                                    for (int si = 0; si < stencil_size_reader; ++si) {
                                        const CellData& c = dbg_stencil[cellIndex(si, sj, sk, stencil_size_reader)];

                                        float cx = 0.0f, cy = 0.0f, cz = 0.0f;
                                        if (include_Moments >= 1 && c.vfrac > 1.0e-12f) {
                                            cx = c.mx / c.vfrac;
                                            cy = c.my / c.vfrac;
                                            cz = c.mz / c.vfrac;
                                        }

                                        std::cout << "[" << si << "," << sj << "," << sk << "] "
                                                << c.vfrac << " "
                                                << "(" << cx << ", " << cy << ", " << cz << ")";

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
                        */
                        /*
                        if (certainty->GetValue(centerCellId) < 0.95f) {
                            predicted_class = 3; // re-assign to cut sheet if certainty is low, since class 3 is "catch-all" for hard cases and should be most robust
                            no_cut_sheets--;
                            no_sheets++;
                            interface_type->SetValue(centerCellId, 3);
                        }
                        */
                        break;

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