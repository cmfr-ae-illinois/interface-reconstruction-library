// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_H_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_H_

#include <utility>
#include <iomanip>
#include <typeinfo>

#include "irl/data_structures/chained_block_storage.h"
#include "irl/data_structures/small_vector.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

// add log for memory operations
// #define DEBUG_MEMORY

namespace IRL {

template <typename ObjectType>
constexpr size_t sizeof_round() {
    size_t the_size = sizeof(ObjectType);
    return (the_size + 15) & ~size_t(15); // Rounds up to the nearest multiple of 16
}

// Overload for objects
template <typename ObjectType>
constexpr size_t sizeof_round(const ObjectType&) {
    return sizeof_round<ObjectType>();
}

// Class that has storage in raw bytes (char)
class HalfEdgeStorage {
 public:
  HalfEdgeStorage(void)
      : data_blocks_start_m(),
        block_size_m(),
        open_block_m(static_cast<std::size_t>(-1)),
        free_location_m(nullptr),
        open_block_end_m(nullptr) {}

  char* operator[](const std::size_t a_index) {
    assert(!data_blocks_start_m.empty());
    assert(a_index < std::min(this->size(), this->blockSize(0)));
    return data_blocks_start_m[0] + a_index;
  }
  const char* operator[](const std::size_t a_index) const {
    assert(!data_blocks_start_m.empty());
    assert(a_index < std::min(this->size(), this->blockSize(0)));
    return data_blocks_start_m[0] + a_index;
  }

  std::size_t blockSize(std::size_t a_index) const {
    assert(data_blocks_start_m.size() == block_size_m.size());
    assert(a_index < block_size_m.size());
    return block_size_m[a_index];
  }

  void resize(const std::size_t a_size) {
    #ifdef DEBUG_MEMORY
    std::cout << "HalfEdge storage " << static_cast<void*>(this) << " is resized to " << a_size << std::endl;
    #endif
    this->resetToSize(a_size);
    open_block_m = 0;
    free_location_m = data_blocks_start_m[0] + a_size;
    open_block_end_m = data_blocks_start_m[0] + block_size_m[0];
    assert(this->size() == a_size);
  }

  void resizeFor(const std::size_t a_size) {
    #ifdef DEBUG_MEMORY
    std::cout << "HalfEdge storage " << static_cast<void*>(this) << " is resized for " << a_size << std::endl;
    #endif
    this->resetToSize(a_size);
    open_block_m = 0;
    free_location_m = data_blocks_start_m[0];
    open_block_end_m = data_blocks_start_m[0] + block_size_m[0];
    assert(this->size() == 0);
  }

  std::size_t size(void) const {
    if (data_blocks_start_m.empty()) {
      return 0;
    } else {
      std::size_t size = 0;
      for (std::size_t n = 0; n < open_block_m; ++n) {
        size += this->blockSize(n);
      }
      assert(data_blocks_start_m[open_block_m] <= free_location_m);
      assert(data_blocks_start_m[open_block_m] + block_size_m[open_block_m] >=
             free_location_m);
      size += static_cast<std::size_t>(free_location_m -
                                       data_blocks_start_m[open_block_m]);
      return size;
    }
  }

  template <class ObjectType>
  ObjectType* getNewObject(void) {
    #ifdef DEBUG_MEMORY
    std::cout << "HalfEdge storage " << static_cast<void*>(this) << " create a new object of type " << typeid(ObjectType).name() << ", the size is " << sizeof_round<ObjectType>() << " bytes or 0x" << std::hex << sizeof_round<ObjectType>() << std::dec << std::endl;
    #endif
    auto OBJECT_SIZE = sizeof_round<ObjectType>(); // float128 require to be memory aligned with adress that are multiple of 16, so we restrict object to have a size that is also a multiple of 16.
    if (free_location_m + OBJECT_SIZE > open_block_end_m) {
      #ifdef DEBUG_MEMORY
      std::cout << "not enough place in the current block, switching to the next one" << std::endl;
      #endif
      if (open_block_m + 1 == data_blocks_start_m.size()) {
        this->allocateAndAppendNewBlock();
      } else {
        ++open_block_m;
        free_location_m = data_blocks_start_m[open_block_m];
        open_block_end_m = free_location_m + block_size_m[open_block_m];
      }
    }
    assert(free_location_m + OBJECT_SIZE <= open_block_end_m);
    ObjectType* const new_object =
        reinterpret_cast<ObjectType*>(free_location_m);
    free_location_m += OBJECT_SIZE;
    #ifdef DEBUG_MEMORY
    std::cout << "new object allocated at " << static_cast<void*>(new_object) << ", next allocation should be at " << static_cast<void*>(free_location_m) << std::endl;
    #endif
    return new_object;
  }

  std::size_t capacity(void) const {
    std::size_t capacity = 0;
    for (std::size_t n = 0; n < data_blocks_start_m.size(); ++n) {
      capacity += this->blockSize(n);
    }
    return capacity;
  }

  HalfEdgeStorage(const HalfEdgeStorage& a_rhs) noexcept {
    // this->resize(a_rhs.size());
    // assert(this->size() == a_rhs.size());
    // if (a_rhs.size() > 0) {
    //   assert(a_rhs.size() <= a_rhs.blockSize(0));
    //   this->copyData(a_rhs.data_blocks_start_m[0], data_blocks_start_m[0],
    //                  a_rhs.size());
    // }
    #ifdef DEBUG_MEMORY
    std::cout << "initializing HalfEdge storage " << static_cast<void*>(this) << " based on another HalfEdge storage" << std::endl;
    #endif
    for (auto& block : data_blocks_start_m) {
        #ifdef DEBUG_MEMORY
        std::cout << "deleting block " << static_cast<void*>(block) << std::endl;
        #endif
        ::operator delete(block);
        block = nullptr;
    }
    this->data_blocks_start_m.resize(0);
    this->block_size_m.resize(0);
    for (std::size_t n = 0; n < a_rhs.data_blocks_start_m.size(); ++n) {
      std::size_t new_block_size = a_rhs.block_size_m[n];
      this->data_blocks_start_m.push_back(
          static_cast<char*>(::operator new(new_block_size)));
      this->block_size_m.push_back(new_block_size);
      #ifdef DEBUG_MEMORY
      std::cout << "creating new block " << static_cast<void*>(data_blocks_start_m[n]) << " of size " << new_block_size << " bytes or 0x" << std::hex << new_block_size << std::dec << std::endl;
      #endif
      if (n != a_rhs.open_block_m) {
        this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
                        block_size_m[n]);
        #ifdef DEBUG_MEMORY
        std::cout << "and copying " << block_size_m[n] << "bytes of data or 0x" << std::hex << block_size_m[n] << std::dec << std::endl;
        #endif
      } else {
        this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
                        a_rhs.free_location_m - a_rhs.data_blocks_start_m[n]);
        #ifdef DEBUG_MEMORY
        std::cout << "and copying " << a_rhs.free_location_m - a_rhs.data_blocks_start_m[n] << "bytes of data or 0x" << std::hex << a_rhs.free_location_m - a_rhs.data_blocks_start_m[n] << std::dec << std::endl;
        #endif
      }
    }
    if (a_rhs.open_block_m != -1) {
      this->open_block_m = a_rhs.open_block_m;
      this->open_block_end_m = a_rhs.open_block_end_m - 
              a_rhs.data_blocks_start_m[a_rhs.open_block_m] + this->data_blocks_start_m[this->open_block_m];
      this->free_location_m = a_rhs.free_location_m - 
              a_rhs.data_blocks_start_m[a_rhs.open_block_m] + this->data_blocks_start_m[this->open_block_m];
    } else {
      this->open_block_m = static_cast<std::size_t>(-1);
      this->open_block_end_m = nullptr;
      this->free_location_m = nullptr;
    }
    assert(this->size() == a_rhs.size());
    #ifdef DEBUG_MEMORY
    std::cout << "done initializing HalfEdge storage " << static_cast<void*>(this) << std::endl;
    #endif
  }

  // NOTE: ONLY EVER COPIES FIRST BLOCK FROM RHS.
  HalfEdgeStorage(const HalfEdgeStorage&& a_rhs) noexcept = delete;

  HalfEdgeStorage& operator=(const HalfEdgeStorage& a_rhs) noexcept {
    #ifdef DEBUG_MEMORY
    std::cout << "copying a HalfEdge storage to HalfEdge storage " << static_cast<void*>(this) << std::endl;
    #endif
    if (this != &a_rhs) {
      // this->resize(a_rhs.size());
      // assert(this->size() == a_rhs.size());
      // if (a_rhs.size() > 0) {
      //   assert(a_rhs.size() <= a_rhs.blockSize(0));
      //   for (std::size_t n = 0; n < data_blocks_start_m.size(); ++n) {
      //     this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
      //                   block_size_m[n]);
      //   }
      // }
      std::size_t old_size = this->data_blocks_start_m.size();
      std::size_t smaller_value = old_size <= a_rhs.data_blocks_start_m.size() ? old_size : a_rhs.data_blocks_start_m.size();
      for (std::size_t n = 0; n < smaller_value; n++) {
        if (a_rhs.block_size_m[n] != this->block_size_m[n]) {
          #ifdef DEBUG_MEMORY
          std::cout << "deleting block " << static_cast<void*>(data_blocks_start_m[n]) << std::endl;
          #endif
          ::operator delete(this->data_blocks_start_m[n]);
          std::size_t new_block_size = a_rhs.block_size_m[n];
          this->data_blocks_start_m[n] = static_cast<char*>(::operator new(new_block_size));
          this->block_size_m[n] = new_block_size;
          #ifdef DEBUG_MEMORY
          std::cout << "creating new block " << static_cast<void*>(data_blocks_start_m[n]) << " of size " << new_block_size << " bytes or 0x" << std::hex << new_block_size << std::dec << std::endl;
          #endif
        }
        if (n != a_rhs.open_block_m) {
          this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
                          a_rhs.block_size_m[n]);
        #ifdef DEBUG_MEMORY
        std::cout << "and copying " << block_size_m[n] << "bytes of data or 0x" << std::hex << block_size_m[n] << std::dec << std::endl;
        #endif
        } else {
          this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
                          a_rhs.free_location_m - a_rhs.data_blocks_start_m[n]);
          #ifdef DEBUG_MEMORY
          std::cout << "and copying " << a_rhs.free_location_m - a_rhs.data_blocks_start_m[n] << "bytes of data or 0x" << std::hex << a_rhs.free_location_m - a_rhs.data_blocks_start_m[n] << std::dec << std::endl;
          #endif
        }
      }
      if (old_size <= a_rhs.data_blocks_start_m.size()) {
        for (std::size_t n = old_size; n < a_rhs.data_blocks_start_m.size(); ++n) {
          std::size_t new_block_size = a_rhs.block_size_m[n];
          this->data_blocks_start_m.push_back(
              static_cast<char*>(::operator new(new_block_size)));
          this->block_size_m.push_back(new_block_size);
          #ifdef DEBUG_MEMORY
          std::cout << "creating new block " << static_cast<void*>(data_blocks_start_m[n]) << " of size " << new_block_size << " bytes or 0x" << std::hex << new_block_size << std::dec  << std::endl;
          #endif
          if (n != a_rhs.open_block_m) {
            this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
                            block_size_m[n]);
        #ifdef DEBUG_MEMORY
        std::cout << "and copying " << block_size_m[n] << "bytes of data or 0x" << std::hex << block_size_m[n] << std::dec << std::endl;
        #endif
          } else {
            this->copyData(a_rhs.data_blocks_start_m[n], data_blocks_start_m[n],
                            a_rhs.free_location_m - a_rhs.data_blocks_start_m[n]);
          #ifdef DEBUG_MEMORY
          std::cout << "and copying " << a_rhs.free_location_m - a_rhs.data_blocks_start_m[n] << "bytes of data or 0x" << std::hex << a_rhs.free_location_m - a_rhs.data_blocks_start_m[n] << std::dec << std::endl;
          #endif
          }
        }
      } else {
        for (std::size_t n = smaller_value; n < old_size; ++n) {
          #ifdef DEBUG_MEMORY
          std::cout << "deleting block " << static_cast<void*>(data_blocks_start_m[n]) << std::endl;
          #endif
          ::operator delete(this->data_blocks_start_m[n]);
        }
        this->data_blocks_start_m.resize(a_rhs.data_blocks_start_m.size());
        this->block_size_m.resize(a_rhs.data_blocks_start_m.size());
      }
      if (a_rhs.open_block_m != -1) {
        this->open_block_m = a_rhs.open_block_m;
        this->open_block_end_m = a_rhs.open_block_end_m - 
                a_rhs.data_blocks_start_m[a_rhs.open_block_m] + this->data_blocks_start_m[this->open_block_m];
        this->free_location_m = a_rhs.free_location_m - 
                a_rhs.data_blocks_start_m[a_rhs.open_block_m] + this->data_blocks_start_m[this->open_block_m];
      } else {
        this->open_block_m = static_cast<std::size_t>(-1);
        this->open_block_end_m = nullptr;
        this->free_location_m = nullptr;
      }
      assert(this->size() == a_rhs.size());
    }
    #ifdef DEBUG_MEMORY
    std::cout << "done copying to HalfEdge storage " << static_cast<void*>(this) << std::endl;
    #endif
    return (*this);
  }

  HalfEdgeStorage& operator=(const HalfEdgeStorage&& a_rhs) noexcept = delete;

  ~HalfEdgeStorage(void) {
    #ifdef DEBUG_MEMORY
    std::cout << "deleting HalfEdge storage " << static_cast<void*>(this) << std::endl;
    #endif
    for (auto& block : data_blocks_start_m) {
      #ifdef DEBUG_MEMORY
      std::cout << "deleting block " << static_cast<void*>(block) << std::endl;
      #endif
      ::operator delete(block);
    }
    open_block_m = 0;
    free_location_m = nullptr;
  }

 private:
  static void copyData(const char* __restrict__ a_other_data,
                       char* const __restrict__ a_my_data,
                       const std::size_t a_size) {
    std::memcpy(a_my_data, a_other_data, a_size);
  }

  void resetToSize(const std::size_t a_size) {
    #ifdef DEBUG_MEMORY
    std::cout << "HalfEdge storage " << static_cast<void*>(this) << " : reseting to size : " << a_size << std::endl;
    #endif
    if (data_blocks_start_m.empty() || this->blockSize(0) < a_size) {
      for (auto& block : data_blocks_start_m) {
        #ifdef DEBUG_MEMORY
        std::cout << "deleting block " << static_cast<void*>(block) << std::endl;
        #endif
        ::operator delete(block);
        block = nullptr;
      }
      data_blocks_start_m.resize(1);
      block_size_m.resize(1);
      data_blocks_start_m[0] = static_cast<char*>(::operator new(a_size));
      block_size_m[0] = a_size;
      #ifdef DEBUG_MEMORY
      std::cout << "creating new block " << static_cast<void*>(data_blocks_start_m[0]) << " of size " << a_size << " bytes or 0x" << std::hex << a_size << std::dec << std::endl;
      #endif
    }
  }

  void allocateAndAppendNewBlock(void) {
    // GROWTH AS FRACTION OF CURRENT SIZE
    static constexpr double growth = 1.0;
    const std::size_t new_block_size =
        static_cast<std::size_t>(growth * static_cast<double>(this->size()));
    data_blocks_start_m.push_back(
        static_cast<char*>(::operator new(new_block_size)));
    block_size_m.push_back(new_block_size);
    open_block_m = data_blocks_start_m.size() - 1;
    free_location_m = data_blocks_start_m[open_block_m];
    open_block_end_m = free_location_m + new_block_size;
    #ifdef DEBUG_MEMORY
    std::cout << "creating new block " << static_cast<void*>(data_blocks_start_m.back()) << " of size " << new_block_size << " bytes or 0x" << std::hex << new_block_size << std::dec   << std::endl;
    #endif
  }

  SmallVector<char*, 8> data_blocks_start_m;
  SmallVector<std::size_t, 8> block_size_m;
  std::size_t open_block_m;
  char* free_location_m;
  char* open_block_end_m;
};

namespace half_edge_polytope {
namespace default_sizes {
static constexpr UnsignedIndex_t complete_block_vertices = 64;
static constexpr UnsignedIndex_t complete_block_faces = 32;
static constexpr UnsignedIndex_t complete_block_half_edges = 128;
}  // namespace default_sizes
}  // namespace half_edge_polytope

template <class PtType, class VertexType = Vertex<PtType>,
          class HalfEdgeType = HalfEdge<VertexType>,
          class FaceType = Face<HalfEdgeType>,
          UnsignedIndex_t kMaxHalfEdges =
              half_edge_polytope::default_sizes::complete_block_half_edges,
          UnsignedIndex_t kMaxVertices =
              half_edge_polytope::default_sizes::complete_block_vertices,
          UnsignedIndex_t kMaxFaces =
              half_edge_polytope::default_sizes::complete_block_faces>
class HalfEdgePolytope {
 public:
  using pt_type = PtType;
  using vertex_type = VertexType;
  using half_edge_type = HalfEdgeType;
  using face_type = FaceType;

  static constexpr UnsignedIndex_t maxHalfEdges = kMaxHalfEdges;
  static constexpr UnsignedIndex_t maxVertices = kMaxVertices;
  static constexpr UnsignedIndex_t maxFaces = kMaxFaces;

  HalfEdgePolytope(void)
      : initial_half_edge_storage_size_m(0),
        initial_vertex_storage_size_m(0),
        initial_face_storage_size_m(0),
        storage_m(),
        vertex_storage_m(),
        half_edge_storage_m(),
        face_storage_m() {}

  static HalfEdgePolytope fromKnownSizes(
      const UnsignedIndex_t a_number_of_half_edges,
      const UnsignedIndex_t a_number_of_vertices,
      const UnsignedIndex_t a_number_of_faces);

  void reset(void);

  void resize(const UnsignedIndex_t a_number_of_half_edges,
              const UnsignedIndex_t a_number_of_vertices,
              const UnsignedIndex_t a_number_of_faces);

  void resizeFor(const UnsignedIndex_t a_number_of_half_edges,
                 const UnsignedIndex_t a_number_of_vertices,
                 const UnsignedIndex_t a_number_of_faces);

  HalfEdgeType& getHalfEdge(const UnsignedIndex_t a_index);

  const HalfEdgeType& getHalfEdge(const UnsignedIndex_t a_index) const;

  VertexType& getVertex(const UnsignedIndex_t a_index);
  const VertexType& getVertex(const UnsignedIndex_t a_index) const;
  FaceType& getFace(const UnsignedIndex_t a_index);
  const FaceType& getFace(const UnsignedIndex_t a_index) const;

  HalfEdgeType* getNewHalfEdge(void);
  HalfEdgeType* getNewHalfEdge(const HalfEdgeType& a_half_edge);
  HalfEdgeType* getNewHalfEdge(HalfEdgeType&& a_half_edge);

  VertexType* getNewVertex(void);
  VertexType* getNewVertex(VertexType&& a_vertex);

  FaceType* getNewFace(void);
  FaceType* getNewFace(FaceType&& a_face);

  template <class GeometryType>
  void setVertexLocations(const GeometryType& a_geometry);

  std::size_t rawCapacity(void) const { return storage_m.capacity(); }
  std::size_t firstBlockCapacity(void) const { return storage_m.blockSize(0); }
  void resizeRawCapacity(const std::size_t a_raw_size) {
    storage_m.resize(a_raw_size);
  }

 protected:
  UnsignedIndex_t getNumberOfInitialFaces(void) const;
  UnsignedIndex_t getNumberOfInitialVertices(void) const;

  HalfEdgePolytope(const UnsignedIndex_t a_number_of_half_edges,
                   const UnsignedIndex_t a_number_of_vertices,
                   const UnsignedIndex_t a_number_of_faces);
  std::size_t initial_half_edge_storage_size_m;
  std::size_t initial_vertex_storage_size_m;
  std::size_t initial_face_storage_size_m;
  HalfEdgeStorage storage_m;
  std::vector<VertexType*> vertex_storage_m;
  std::vector<HalfEdgeType*> half_edge_storage_m;
  std::vector<FaceType*> face_storage_m;
};

// IO
// template <class PtType, class VertexType, class HalfEdgeType, class
// FaceType,
//           UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
//           UnsignedIndex_t kMaxFaces>
// inline std::ostream &
// operator<<(std::ostream &out,
//            const HalfEdgePolytope<PtType, VertexType, HalfEdgeType,
//            FaceType,
//                                   kMaxHalfEdges, kMaxVertices, kMaxFaces>
//               &a_polyhedron);
}  // namespace IRL

#include "irl/geometry/half_edge_structures/half_edge_polytope.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_H_
