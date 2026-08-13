#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_TPP_

#include "examples/implicit_surface_reconstruction/binary.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace sparse_moment_io {

constexpr std::uint32_t kVersion = 1;

struct Header {
  std::uint32_t nx = 0;
  std::uint32_t ny = 0;
  std::uint32_t nz = 0;
  std::uint64_t mask_word_count = 0;
  std::uint64_t mixed_cell_count = 0;
};

template <class T>
void writeValue(std::ofstream* stream, const T& value) {
  static_assert(std::is_trivially_copyable<T>::value,
                "Binary values must be trivially copyable");
  stream->write(reinterpret_cast<const char*>(&value), sizeof(T));
  if (!*stream) throw std::runtime_error("Failed while writing moment file");
}

template <class T>
void readValue(std::ifstream* stream, T* value) {
  static_assert(std::is_trivially_copyable<T>::value,
                "Binary values must be trivially copyable");
  stream->read(reinterpret_cast<char*>(value), sizeof(T));
  if (!*stream) throw std::runtime_error("Truncated sparse moment file");
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void writeHeader(std::ofstream* stream, const BasicMesh& mesh,
                 const InsideCellMask& mask,
                 const std::uint64_t mixed_cell_count) {
  writeValue(stream, kVersion);
  writeValue(stream, static_cast<std::uint32_t>(VM_ORDER));
  writeValue(stream, static_cast<std::uint32_t>(SM_ORDER));
  writeValue(stream, static_cast<std::uint32_t>(mesh.getNx()));
  writeValue(stream, static_cast<std::uint32_t>(mesh.getNy()));
  writeValue(stream, static_cast<std::uint32_t>(mesh.getNz()));
  writeValue(stream, static_cast<std::uint64_t>(mask.words().size()));
  writeValue(stream, mixed_cell_count);
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
Header readHeader(std::ifstream* stream) {
  std::uint32_t version = 0, vm_order = 0, sm_order = 0;
  Header header;
  readValue(stream, &version);
  readValue(stream, &vm_order);
  readValue(stream, &sm_order);
  readValue(stream, &header.nx);
  readValue(stream, &header.ny);
  readValue(stream, &header.nz);
  readValue(stream, &header.mask_word_count);
  readValue(stream, &header.mixed_cell_count);

  if (version != kVersion) {
    throw std::runtime_error("Unsupported sparse moment file version");
  }
  if (vm_order != VM_ORDER || sm_order != SM_ORDER) {
    throw std::runtime_error("Moment orders do not match sparse moment file");
  }

  const std::uint64_t cell_count =
      static_cast<std::uint64_t>(header.nx) * header.ny * header.nz;
  if (header.mask_word_count != (cell_count + 63) / 64) {
    throw std::runtime_error("Invalid inside-mask size in sparse moment file");
  }
  return header;
}

inline InsideCellMask readMask(std::ifstream* stream, const Header& header) {
  const std::size_t cell_count =
      static_cast<std::size_t>(header.nx) * header.ny * header.nz;
  InsideCellMask mask(cell_count);
  for (std::uint64_t n = 0; n < header.mask_word_count; ++n) {
    readValue(stream, &mask.words()[static_cast<std::size_t>(n)]);
  }
  return mask;
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
std::vector<SparseMixedCellMoments<VM_ORDER, SM_ORDER>> readMixedCells(
    std::ifstream* stream, const Header& header) {
  using Record = SparseMixedCellMoments<VM_ORDER, SM_ORDER>;
  if (header.mixed_cell_count >
      static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error("Too many mixed records in sparse moment file");
  }

  std::vector<Record> records(
      static_cast<std::size_t>(header.mixed_cell_count));
  std::uint32_t previous_index = 0;
  bool have_previous = false;
  const std::uint64_t cell_count =
      static_cast<std::uint64_t>(header.nx) * header.ny * header.nz;
  for (auto& record : records) {
    readValue(stream, &record.linear_index);
    if (record.linear_index >= cell_count ||
        (have_previous && record.linear_index <= previous_index)) {
      throw std::runtime_error(
          "Mixed-cell indices are invalid or not strictly ordered");
    }
    previous_index = record.linear_index;
    have_previous = true;
    for (std::size_t n = 0; n < Record::NV; ++n)
      readValue(stream, &record.volume[n]);
    for (std::size_t n = 0; n < Record::NS; ++n)
      readValue(stream, &record.surface[n]);
  }
  return records;
}

template <class Function>
void forEachSetBit(const InsideCellMask& mask, Function&& function) {
  const auto& words = mask.words();
  for (std::size_t word_index = 0; word_index < words.size(); ++word_index) {
    std::uint64_t word = words[word_index];
    while (word != 0) {
#if defined(__GNUC__) || defined(__clang__)
      const unsigned int bit = static_cast<unsigned int>(__builtin_ctzll(word));
#else
      unsigned int bit = 0;
      while (((word >> bit) & std::uint64_t{1}) == 0) ++bit;
#endif
      const std::size_t index = word_index * 64 + bit;
      if (index < mask.cellCount()) function(index);
      word &= word - 1;
    }
  }
}

inline void kahanAdd(double* sum, double* compensation, const double value) {
  const double corrected = value - *compensation;
  const double updated = *sum + corrected;
  *compensation = (updated - *sum) - corrected;
  *sum = updated;
}

}  // namespace sparse_moment_io

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void writeMomentsToBinary(
    const BasicMesh& mesh, const InsideCellMask& inside_cells,
    const std::vector<SparseMixedCellMoments<VM_ORDER, SM_ORDER>>& mixed_cells,
    const std::string& filename) {
  using Record = SparseMixedCellMoments<VM_ORDER, SM_ORDER>;
  const std::size_t expected_cells =
      static_cast<std::size_t>(mesh.getNx()) * mesh.getNy() * mesh.getNz();
  if (inside_cells.cellCount() != expected_cells) {
    throw std::runtime_error("Inside mask does not match the mesh");
  }

  std::ofstream stream(filename, std::ios::binary);
  if (!stream)
    throw std::runtime_error("Cannot open file for writing: " + filename);

  sparse_moment_io::writeHeader<VM_ORDER, SM_ORDER>(
      &stream, mesh, inside_cells,
      static_cast<std::uint64_t>(mixed_cells.size()));
  for (const std::uint64_t word : inside_cells.words())
    sparse_moment_io::writeValue(&stream, word);
  for (const Record& record : mixed_cells) {
    sparse_moment_io::writeValue(&stream, record.linear_index);
    for (std::size_t n = 0; n < Record::NV; ++n)
      sparse_moment_io::writeValue(&stream, record.volume[n]);
    for (std::size_t n = 0; n < Record::NS; ++n)
      sparse_moment_io::writeValue(&stream, record.surface[n]);
  }
  std::cout << "✅ Sparse moments written to binary file: " << filename
            << std::endl;
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void readMomentsFromBinary(
    const std::string& filename,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments) {
  if (moments == nullptr) throw std::invalid_argument("moments is null");
  using VM = IRL::GeneralMoments3D<VM_ORDER>;
  using Record = SparseMixedCellMoments<VM_ORDER, SM_ORDER>;

  std::ifstream stream(filename, std::ios::binary);
  if (!stream)
    throw std::runtime_error("Cannot open file for reading: " + filename);
  const auto header = sparse_moment_io::readHeader<VM_ORDER, SM_ORDER>(&stream);
  const BasicMesh& mesh = moments->getMesh();
  if (header.nx != static_cast<std::uint32_t>(mesh.getNx()) ||
      header.ny != static_cast<std::uint32_t>(mesh.getNy()) ||
      header.nz != static_cast<std::uint32_t>(mesh.getNz())) {
    throw std::runtime_error("Sparse moment file dimensions do not match mesh");
  }
  const InsideCellMask mask = sparse_moment_io::readMask(&stream, header);
  const auto records =
      sparse_moment_io::readMixedCells<VM_ORDER, SM_ORDER>(&stream, header);

  for (int i = mesh.imin(); i <= mesh.imax(); ++i)
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j)
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) (*moments)(i, j, k) = {};

  sparse_moment_io::forEachSetBit(mask, [&](const std::size_t index) {
    int i, j, k;
    getCellIndicesFromLinearIndex(index, mesh.getNy(), mesh.getNz(), &i, &j,
                                  &k);
    i += mesh.imin();
    j += mesh.jmin();
    k += mesh.kmin();
    const IRL::RectangularCuboid cell = IRL::RectangularCuboid::fromBoundingPts(
        IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
        IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
    (*moments)(i, j, k).first = IRL::getVolumeMoments<VM>(cell);
  });
  for (const Record& record : records) {
    int i, j, k;
    getCellIndicesFromLinearIndex(record.linear_index, mesh.getNy(),
                                  mesh.getNz(), &i, &j, &k);
    i += mesh.imin();
    j += mesh.jmin();
    k += mesh.kmin();
    auto& destination = (*moments)(i, j, k);
    for (std::size_t n = 0; n < Record::NV; ++n)
      destination.first[n] = record.volume[n];
    for (std::size_t n = 0; n < Record::NS; ++n)
      destination.second[n] = record.surface[n];
  }
  std::cout << "✅ Sparse moments loaded from binary file: " << filename
            << std::endl;
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void coarsenMomentsFromBinary(
    const std::string& fine_filename, int factor,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* coarse_moments) {
  if (factor <= 0) throw std::invalid_argument("factor must be > 0");
  if (coarse_moments == nullptr)
    throw std::invalid_argument("coarse_moments is null");

  using VM = IRL::GeneralMoments3D<VM_ORDER>;
  using Record = SparseMixedCellMoments<VM_ORDER, SM_ORDER>;
  std::ifstream stream(fine_filename, std::ios::binary);
  if (!stream)
    throw std::runtime_error("Cannot open fine file: " + fine_filename);
  const auto header = sparse_moment_io::readHeader<VM_ORDER, SM_ORDER>(&stream);
  const InsideCellMask mask = sparse_moment_io::readMask(&stream, header);
  const auto records =
      sparse_moment_io::readMixedCells<VM_ORDER, SM_ORDER>(&stream, header);

  const BasicMesh& mesh = coarse_moments->getMesh();
  if (header.nx != static_cast<std::uint32_t>(factor * mesh.getNx()) ||
      header.ny != static_cast<std::uint32_t>(factor * mesh.getNy()) ||
      header.nz != static_cast<std::uint32_t>(factor * mesh.getNz())) {
    throw std::runtime_error(
        "Fine sparse dimensions must equal factor times coarse dimensions");
  }

  Data<VM> volume_compensation(&mesh);
  Data<IRL::GeneralSurfaceMoments3D<SM_ORDER>> surface_compensation(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i)
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j)
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        (*coarse_moments)(i, j, k) = {};
        volume_compensation(i, j, k) = {};
        surface_compensation(i, j, k) = {};
      }

  const double fine_dx = mesh.lx() / static_cast<double>(header.nx);
  const double fine_dy = mesh.ly() / static_cast<double>(header.ny);
  const double fine_dz = mesh.lz() / static_cast<double>(header.nz);
  const double x_lower = mesh.x(mesh.imin());
  const double y_lower = mesh.y(mesh.jmin());
  const double z_lower = mesh.z(mesh.kmin());

  auto add_volume = [&](const std::size_t index, const VM& volume) {
    int i, j, k;
    getCellIndicesFromLinearIndex(index, header.ny, header.nz, &i, &j, &k);
    const int I = mesh.imin() + i / factor;
    const int J = mesh.jmin() + j / factor;
    const int K = mesh.kmin() + k / factor;
    auto& sum = (*coarse_moments)(I, J, K).first;
    auto& compensation = volume_compensation(I, J, K);
    for (std::size_t n = 0; n < Record::NV; ++n)
      sparse_moment_io::kahanAdd(&sum[n], &compensation[n], volume[n]);
  };

  sparse_moment_io::forEachSetBit(mask, [&](const std::size_t index) {
    int i, j, k;
    getCellIndicesFromLinearIndex(index, header.ny, header.nz, &i, &j, &k);
    const IRL::Pt lower(x_lower + i * fine_dx, y_lower + j * fine_dy,
                        z_lower + k * fine_dz);
    const IRL::Pt upper(lower.x() + fine_dx, lower.y() + fine_dy,
                        lower.z() + fine_dz);
    const auto cell = IRL::RectangularCuboid::fromBoundingPts(lower, upper);
    add_volume(index, IRL::getVolumeMoments<VM>(cell));
  });

  for (const Record& record : records) {
    int i, j, k;
    getCellIndicesFromLinearIndex(record.linear_index, header.ny, header.nz, &i,
                                  &j, &k);
    const int I = mesh.imin() + i / factor;
    const int J = mesh.jmin() + j / factor;
    const int K = mesh.kmin() + k / factor;
    auto& sum = (*coarse_moments)(I, J, K);
    auto& volume_c = volume_compensation(I, J, K);
    auto& surface_c = surface_compensation(I, J, K);
    for (std::size_t n = 0; n < Record::NV; ++n)
      sparse_moment_io::kahanAdd(&sum.first[n], &volume_c[n], record.volume[n]);
    for (std::size_t n = 0; n < Record::NS; ++n)
      sparse_moment_io::kahanAdd(&sum.second[n], &surface_c[n],
                                 record.surface[n]);
  }
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_TPP_
