#include <torch/torch.h>
#include <iostream>
#include <vector>
#include <string>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/data_set.h"
#include "irl/ml_classification/net.h"

#include <random> //REMOVE, was just for testing

#include <vtkSmartPointer.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkRectilinearGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>

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

    auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    reader->SetCaseFileName("/home/quirin/jet/nga.case");  // Adjust the path as needed
    reader->UpdateInformation();  // Populate meta info like time steps
    vtkInformation* info = reader->GetOutputInformation(0);
    if (info->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        int nSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        std::vector<double> timeSteps(nSteps);
        info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0]);


        std::cout << "Number of timesteps: " << nSteps << std::endl;
        nSteps = 50;

        for (int t = 49; t < nSteps; ++t) {
            reader->SetTimeValue(timeSteps[t]);
            reader->Update();

            auto output = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
            if (!output) {
                std::cerr << "Reader did not return a vtkMultiBlockDataSet." << std::endl;
                return 1;
            }

            std::cout << "Number of blocks: " << output->GetNumberOfBlocks() << std::endl;

            for (unsigned int i = 0; i < output->GetNumberOfBlocks(); ++i) {
                auto block = output->GetBlock(i);
                if (block) {
                    std::cout << "Block " << i << " is type: " << block->GetClassName() << std::endl;
                } else {
                    std::cout << "Block " << i << " is nullptr." << std::endl;
                }
            }

            auto block0 = vtkRectilinearGrid::SafeDownCast(reader->GetOutput()->GetBlock(0));
            if (!block0) {
                std::cerr << "Could not access block 0 of the dataset." << std::endl;
                return 1;
            }

            const auto dimensions = block0->GetDimensions();
            const int ncellsx = dimensions[0] - 1;
            const int ncellsy = dimensions[1] - 1;
            const int ncellsz = dimensions[2] - 1;

            std::cout << "Dimensions of grid are: " << ncellsx << ", " << ncellsy << ", " << ncellsz << std::endl;

            auto px = block0->GetXCoordinates();
            auto py = block0->GetYCoordinates();
            auto pz = block0->GetZCoordinates();
            std::vector<double> coordsx(ncellsx + 1, 0.0);
            std::vector<double> coordsy(ncellsy + 1, 0.0);
            std::vector<double> coordsz(ncellsz + 1, 0.0);
            for (int i = 0; i < ncellsx + 1; i++) coordsx[i] = px->GetComponent(i, 0);
            for (int j = 0; j < ncellsy + 1; j++) coordsy[j] = py->GetComponent(j, 0);
            for (int k = 0; k < ncellsz + 1; k++) coordsz[k] = pz->GetComponent(k, 0);

            auto cellData = block0->GetCellData();
            int numArrays = cellData->GetNumberOfArrays();
            std::cout << "Available cell data arrays: " << numArrays << std::endl;
            for (int i = 0; i < numArrays; ++i) {
                std::cout << "Array " << i << ": " << cellData->GetArrayName(i) << std::endl;
            }

            vtkDataArray* vofArray = block0->GetCellData()->GetArray("VOF");
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
                        const int index = k * (ncellsy * ncellsx) + j * ncellsx + i;
                        vfrac[i][j][k] = static_cast<double>(vofArray->GetComponent(index, 0));
                        // Each z-slice has ncellsy * ncellsx elements, move forward k slices.
                        // Each row has ncellsx elements, so this skips over j rows
                        // Move to the correct column in a row
                    }
                }
            }

            WriteField(coordsx, coordsy, coordsz, vfrac, std::string("vfrac_")+std::to_string(t));
        }
    } else {
        std::cerr << "No time steps found in dataset." << std::endl;
    }
    return 0;
}