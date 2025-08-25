// #include <torch/torch.h>
#include <iostream>
#include <vector>
#include <string>

// #include "irl/ml_classification/data_gen.h"
// #include "irl/ml_classification/data_set.h"
// #include "irl/ml_classification/net.h"

#include <random> //REMOVE, was just for testing

#include <vtkSmartPointer.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkRectilinearGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>
#include <vtkUnstructuredGridWriter.h>

template <class CellData, class CoordData>
void WriteField(const CoordData& coordsx, const CoordData& coordsy, const CoordData& coordsz,
                const CellData& field, const std::string& field_name) {
  const auto file_name = "./" + field_name + ".vtr"; 

  double start = 0.0;
  double end = 1.0;

  const int ncellsx = coordsx.size() - 1;
  const int ncellsy = coordsy.size() - 1;
  const int ncellsz = coordsz.size() - 1;

  FILE* file;
  file = fopen(file_name.c_str(), "w");

  fprintf(file, "<VTKFile type=\"RectilinearGrid\">\n");
  fprintf(file, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0,
          ncellsx, 0, ncellsy, 0, ncellsz);
  fprintf(file, "<Piece Extent=\"%d %d %d %d %d %d\">\n", 0, ncellsx, 0, ncellsy,
          0, ncellsz);

  fprintf(file, "<Coordinates>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncellsx + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coordsx[i])); //QUESTON: coords is +1 longer than ncells...?
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncellsy + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coordsy[i]));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncellsz + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coordsz[i]));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file, "</Coordinates>\n");

  fprintf(file, "<PointData>\n</PointData>\n");

  fprintf(file, "<CellData Scalars=\"");
  fprintf(file, "%s ", field_name.c_str());
  fprintf(file, "\" >\n");
  fprintf(file,
          "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n",
          field_name.c_str());
  for (int k = 0; k < ncellsz; ++k) {
    for (int j = 0; j < ncellsy; ++j) {
      for (int i = 0; i < ncellsx; ++i) {
        fprintf(file, "%15.8E ", static_cast<double>(field[i][j][k]));
      }
    }
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
  fclose(file);
}

int main (int argc, char* argv[]) {

    auto reader_fields = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    reader_fields->SetCaseFileName("/home/quirin/jet/nga.case");  // Adjust the path as needed
    reader_fields->UpdateInformation();  // Populate meta info like time steps
    vtkInformation* info_fields = reader_fields->GetOutputInformation(0);

    auto reader_surface = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    reader_surface->SetCaseFileName("/home/quirin/jet/plic.case");  // Adjust the path as needed
    reader_surface->UpdateInformation();  // Populate meta info like time steps
    vtkInformation* info_surface = reader_surface->GetOutputInformation(0);

    if (info_fields->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()) && info_surface->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        int nSteps = info_fields->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        if (nSteps != info_surface->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
            std::cerr << "PLIC file does not have same number of steps as mesh file." << std::endl;
            return 1;
        }

        std::vector<double> timeSteps(nSteps);
        info_fields->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0]);


        std::cout << "Number of timesteps: " << nSteps << std::endl;
        nSteps = 50;

        for (int t = 49; t < nSteps; ++t) {
            reader_fields->SetTimeValue(timeSteps[t]);
            reader_fields->Update();
            reader_surface->SetTimeValue(timeSteps[t]);
            reader_surface->Update();

            auto output_fields = vtkMultiBlockDataSet::SafeDownCast(reader_fields->GetOutput());
            if (!output_fields) {
                std::cerr << "Mesh reader did not return a vtkMultiBlockDataSet." << std::endl;
                return 1;
            }

            auto output_surface = vtkMultiBlockDataSet::SafeDownCast(reader_surface->GetOutput());
            if (!output_surface) {
                std::cerr << "PLIC reader did not return a vtkMultiBlockDataSet." << std::endl;
                return 1;
            }

            std::cout << "Number of mesh blocks: " << output_fields->GetNumberOfBlocks() << std::endl;
            for (unsigned int i = 0; i < output_fields->GetNumberOfBlocks(); ++i) {
                auto block = output_fields->GetBlock(i);
                if (block) {
                    std::cout << "Block " << i << " is type: " << block->GetClassName() << std::endl;
                } else {
                    std::cout << "Block " << i << " is nullptr." << std::endl;
                }
            }

            std::cout << "Number of surface blocks: " << output_surface->GetNumberOfBlocks() << std::endl;
            for (unsigned int i = 0; i < output_surface->GetNumberOfBlocks(); ++i) {
                auto block = output_surface->GetBlock(i);
                if (block) {
                    std::cout << "Block " << i << " is type: " << block->GetClassName() << std::endl;
                } else {
                    std::cout << "Block " << i << " is nullptr." << std::endl;
                }
            }

            // Merge all PLIC blocks
            vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
            for (unsigned int i = 0; i < output_surface->GetNumberOfBlocks(); ++i) {
                appendFilter->AddInputData(output_surface->GetBlock(i));
            }
            appendFilter->Update();
            vtkSmartPointer<vtkUnstructuredGrid> mergedData = vtkSmartPointer<vtkUnstructuredGrid>::New();
            mergedData->DeepCopy(appendFilter->GetOutput());

            auto mesh = vtkRectilinearGrid::SafeDownCast(reader_fields->GetOutput()->GetBlock(0));
            if (!mesh) {
                std::cerr << "Could not access block 0 of the dataset." << std::endl;
                return 1;
            }

            const auto dimensions = mesh->GetDimensions();
            const int ncellsx = dimensions[0] - 1;
            const int ncellsy = dimensions[1] - 1;
            const int ncellsz = dimensions[2] - 1;

            std::cout << "Dimensions of grid are: " << ncellsx << ", " << ncellsy << ", " << ncellsz << std::endl;

            auto px = mesh->GetXCoordinates();
            auto py = mesh->GetYCoordinates();
            auto pz = mesh->GetZCoordinates();
            std::vector<double> coordsx(ncellsx + 1, 0.0);
            std::vector<double> coordsy(ncellsy + 1, 0.0);
            std::vector<double> coordsz(ncellsz + 1, 0.0);
            for (int i = 0; i < ncellsx + 1; i++) coordsx[i] = px->GetComponent(i, 0);
            for (int j = 0; j < ncellsy + 1; j++) coordsy[j] = py->GetComponent(j, 0);
            for (int k = 0; k < ncellsz + 1; k++) coordsz[k] = pz->GetComponent(k, 0);

            auto cellData = mesh->GetCellData();
            int numArrays = cellData->GetNumberOfArrays();
            std::cout << "Available cell data arrays: " << numArrays << std::endl;
            for (int i = 0; i < numArrays; ++i) {
                std::cout << "Array " << i << ": " << cellData->GetArrayName(i) << std::endl;
            }

            vtkDataArray* vofArray = mesh->GetCellData()->GetArray("VOF");
            if (!vofArray) {
                std::cerr << "VOF array not found in the dataset." << std::endl;
                return 1;
            }
            std::cout << "Number of tuples = " << vofArray->GetNumberOfTuples() << std::endl;
            std::cout << "VOF values (first 10):" << std::endl;
            for (int i = 0; i < 10; ++i) {
                std::cout << vofArray->GetComponent(i, 0) << " ";
            }
            std::cout << std::endl;
        
            auto vfrac = std::vector<std::vector<std::vector<double>>>(ncellsx, std::vector<std::vector<double>>(ncellsy, std::vector<double>(ncellsz, 0.0)));

            for (int i = 0; i < ncellsx; ++i) {
                for (int j = 0; j < ncellsy; ++j) {
                    for (int k = 0; k < ncellsz; ++k) {
                        int ids[3] = {i,j,k};
                        vfrac[i][j][k] = static_cast<double>(vofArray->GetComponent(mesh->ComputeCellId(ids), 0));
                    }
                }
            }

            WriteField(coordsx, coordsy, coordsz, vfrac, std::string("vfrac_")+std::to_string(t));

            // Writing merged PLIC
            vtkSmartPointer<vtkUnstructuredGridWriter> writer_surface = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
            const auto file_name = "plic_" + std::to_string(t) + ".vtk"; 
            writer_surface->SetFileName(file_name.c_str());
            writer_surface->SetInputData(mergedData);
            writer_surface->Write();
        }
    } else {
        std::cerr << "No time steps found in dataset." << std::endl;
    }
    return 0;
}