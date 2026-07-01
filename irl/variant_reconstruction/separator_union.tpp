// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SEPARATOR_UNION_TPP_
#define IRL_SEPARATOR_UNION_TPP_

namespace IRL {

inline SeparatorUnion::SeparatorUnion(void)
    : type_m(SeparatorType::OnePlane), planes_m({Plane(), Plane()}) {}

inline SeparatorUnion::SeparatorUnion(const Plane& a_plane) {
  type_m = SeparatorType::OnePlane;
  planes_m[0] = a_plane;
};

inline SeparatorUnion::SeparatorUnion(const Plane& a_plane_0,
                                      const Plane& a_plane_1) {
  type_m = SeparatorType::TwoPlanes;
  planes_m[0] = a_plane_0;
  planes_m[1] = a_plane_1;
};

inline SeparatorUnion::SeparatorUnion(const Paraboloid& a_paraboloid) {
  type_m = SeparatorType::Paraboloid;
  paraboloid_m = a_paraboloid;
};

inline SeparatorUnion::SeparatorUnion(const Cylinder& a_cylinder) {
  type_m = SeparatorType::Cylinder;
  cylinder_m = a_cylinder;
};

inline SeparatorUnion& SeparatorUnion::operator=(const SeparatorUnion& other) {
  type_m = other.type();
  switch (type_m) {
    case SeparatorType::OnePlane:
      planes_m[0] = other.getPlanes()[0];
      break;
    case SeparatorType::TwoPlanes:
      planes_m = other.getPlanes();
      break;
    case SeparatorType::Paraboloid:
      paraboloid_m = other.getParaboloid();
      break;
    case SeparatorType::Cylinder:
      cylinder_m = other.getCylinder();
      break;
    default:
      std::runtime_error("Unrecognized reconstruction type in SeparatorUnion");
  }
  return *this;
};

// This just replaces the current separator union
// Only defined because it is required by AMReX
inline SeparatorUnion& SeparatorUnion::operator+=(const SeparatorUnion& other) {
  type_m = other.type();
  switch (type_m) {
    case SeparatorType::OnePlane:
      planes_m[0] = other.getPlanes()[0];
      break;
    case SeparatorType::TwoPlanes:
      planes_m = other.getPlanes();
      break;
    case SeparatorType::Paraboloid:
      paraboloid_m = other.getParaboloid();
      break;
    case SeparatorType::Cylinder:
      cylinder_m = other.getCylinder();
      break;
    default:
      std::runtime_error("Unrecognized reconstruction type in SeparatorUnion");
  }
  return *this;
};

inline SeparatorUnion& SeparatorUnion::operator=(const Plane& other) {
  type_m = SeparatorType::OnePlane;
  planes_m[0] = other;
  return *this;
};

inline SeparatorUnion& SeparatorUnion::operator=(const PlanarSeparator& other) {
  switch (other.getNumberOfPlanes()) {
    case 0:
      type_m = SeparatorType::OnePlane;
      planes_m[0] = Plane();
      break;
    case 1:
      type_m = SeparatorType::OnePlane;
      planes_m[0] = other[0];
      break;
    case 2:
      type_m = SeparatorType::TwoPlanes;
      planes_m[0] = other[0];
      planes_m[1] = other[0];
      break;
    default:
      std::runtime_error("SeparatorUnion cannot contain more than 2 planes");
  }
  return *this;
};

inline SeparatorUnion& SeparatorUnion::operator=(const Paraboloid& other) {
  type_m = SeparatorType::Paraboloid;
  paraboloid_m = other;
  return *this;
};

inline SeparatorUnion& SeparatorUnion::operator=(const Cylinder& other) {
  type_m = SeparatorType::Cylinder;
  cylinder_m = other;
  return *this;
};

inline const SeparatorUnion::SeparatorType& SeparatorUnion::type(void) const {
  return type_m;
};
inline SeparatorUnion::SeparatorType& SeparatorUnion::type(void) {
  return type_m;
};
inline const Plane& SeparatorUnion::getPlane(void) const {
  return planes_m[0];
};
inline Plane& SeparatorUnion::getPlane(void) { return planes_m[0]; };
inline const Plane& SeparatorUnion::getPlane(
    const UnsignedIndex_t a_index) const {
  return planes_m[a_index];
};
inline Plane& SeparatorUnion::getPlane(const UnsignedIndex_t a_index) {
  return planes_m[a_index];
};
inline const std::array<Plane, 2>& SeparatorUnion::getPlanes(void) const {
  return planes_m;
};
inline std::array<Plane, 2>& SeparatorUnion::getPlanes(void) {
  return planes_m;
};
inline const Paraboloid& SeparatorUnion::getParaboloid(void) const {
  return paraboloid_m;
};
inline Paraboloid& SeparatorUnion::getParaboloid(void) { return paraboloid_m; };
inline const Cylinder& SeparatorUnion::getCylinder(void) const {
  return cylinder_m;
};
inline Cylinder& SeparatorUnion::getCylinder(void) { return cylinder_m; };

inline void SeparatorUnion::setPlane(const int a_index, const Plane& a_plane) {
  if ((type_m == SeparatorType::OnePlane && a_index == 0) ||
      (type_m == SeparatorType::TwoPlanes && a_index >= 0 && a_index <= 1)) {
    planes_m[a_index] = a_plane;
  }
}

inline void SeparatorUnion::setToPlane(void) {
  type_m = SeparatorType::OnePlane;
  planes_m[0] = Plane();
}
inline void SeparatorUnion::setToPlane(const Plane& a_plane) {
  type_m = SeparatorType::OnePlane;
  planes_m[0] = a_plane;
}
inline void SeparatorUnion::setToPlanes(const Plane& a_plane0,
                                        const Plane& a_plane1) {
  type_m = SeparatorType::TwoPlanes;
  planes_m[0] = a_plane0;
  planes_m[1] = a_plane1;
}
inline void SeparatorUnion::setToParaboloid(void) {
  type_m = SeparatorType::Paraboloid;
}
inline void SeparatorUnion::setToParaboloid(const Paraboloid& a_paraboloid) {
  type_m = SeparatorType::Paraboloid;
  paraboloid_m = a_paraboloid;
}
inline void SeparatorUnion::setToCylinder(void) {
  type_m = SeparatorType::Cylinder;
}
inline void SeparatorUnion::setToCylinder(const Cylinder& a_cylinder) {
  type_m = SeparatorType::Cylinder;
  cylinder_m = a_cylinder;
}

inline void SeparatorUnion::setToFull(void) {
  type_m = SeparatorType::Paraboloid;
  paraboloid_m.markAsAlwaysAbove();
}

inline void SeparatorUnion::setToEmpty(void) {
  type_m = SeparatorType::Paraboloid;
  paraboloid_m.markAsAlwaysBelow();
}

inline const bool SeparatorUnion::isFull(void) const {
  switch (this->type()) {
    case SeparatorType::OnePlane:
      return this->getPlane().distance() > 1.0e14;
    case SeparatorType::TwoPlanes:
      return false;
    case SeparatorType::Paraboloid:
      return this->getParaboloid().isAlwaysAbove();
    case SeparatorType::Cylinder:
      return this->getCylinder().isAlwaysAbove();
    default:
      std::runtime_error("SeparatorUnion type not recognized");
      return false;
  }
}

inline const bool SeparatorUnion::isEmpty(void) const {
  switch (this->type()) {
    case SeparatorType::OnePlane:
      return this->getPlane().distance() < -1.0e14;
    case SeparatorType::TwoPlanes:
      return false;
    case SeparatorType::Paraboloid:
      return this->getParaboloid().isAlwaysBelow();
    case SeparatorType::Cylinder:
      return this->getCylinder().isAlwaysBelow();
    default:
      std::runtime_error("SeparatorUnion type not recognized");
      return false;
  }
}

inline void SeparatorUnion::serialize(ByteBuffer* a_buffer) const {
  a_buffer->pack(&this->type(), 1);
  switch (this->type()) {
    case SeparatorType::OnePlane:
      this->getPlane().serialize(a_buffer);
      break;
    case SeparatorType::TwoPlanes:
      this->getPlanes()[0].serialize(a_buffer);
      this->getPlanes()[1].serialize(a_buffer);
      break;
    case SeparatorType::Paraboloid:
      this->getParaboloid().serialize(a_buffer);
      break;
    case SeparatorType::Cylinder:
      this->getCylinder().serialize(a_buffer);
      break;
    default:
      std::runtime_error("SeparatorUnion type cannot be serialized");
  }
};

inline void SeparatorUnion::unpackSerialized(ByteBuffer* a_buffer) {
  a_buffer->unpack(&this->type(), 1);
  switch (this->type()) {
    case SeparatorType::OnePlane:
      this->getPlane().unpackSerialized(a_buffer);
      break;
    case SeparatorType::TwoPlanes:
      this->getPlanes()[0].unpackSerialized(a_buffer);
      this->getPlanes()[1].unpackSerialized(a_buffer);
      break;
    case SeparatorType::Paraboloid:
      this->getParaboloid().unpackSerialized(a_buffer);
      break;
    case SeparatorType::Cylinder:
      this->getCylinder().unpackSerialized(a_buffer);
      break;
    default:
      std::runtime_error("SeparatorUnion type cannot be unpacked");
  }
};

inline void SeparatorUnion::shift(const Pt a_shift) {
  switch (type_m) {
    case SeparatorType::OnePlane:
      planes_m[0].distance() += planes_m[0].normal() * a_shift;
      break;
    case SeparatorType::TwoPlanes:
      planes_m[0].distance() += planes_m[0].normal() * a_shift;
      planes_m[1].distance() += planes_m[1].normal() * a_shift;
      break;
    case SeparatorType::Paraboloid:
      paraboloid_m.setDatum(paraboloid_m.getDatum() + a_shift);
      break;
    case SeparatorType::Cylinder:
      cylinder_m.setDatum(cylinder_m.getDatum() + a_shift);
      break;
    default:
      std::runtime_error("SeparatorUnion type cannot shift datum");
  }
};

inline void SeparatorUnion::reflect(const SeparatorUnion& a_ref,
                                    const int a_dir, const double a_loc) {
  type_m = a_ref.type();
  switch (type_m) {
    case SeparatorType::OnePlane: {
      planes_m[0].normal() = a_ref.getPlane().normal();
      planes_m[0].distance() = a_ref.getPlane().distance();
      planes_m[0].normal()[a_dir] *= -1.0;
      planes_m[0].distance() -= 2.0 * a_ref.getPlane().normal()[a_dir] * a_loc;
      break;
    }
    case SeparatorType::TwoPlanes: {
      for (UnsignedIndex_t i = 0; i <= 1; ++i) {
        planes_m[i].normal() = a_ref.getPlane(i).normal();
        planes_m[i].distance() = a_ref.getPlane(i).distance();
        planes_m[i].normal()[a_dir] *= -1.0;
        planes_m[i].distance() -=
            2.0 * a_ref.getPlane(i).normal()[a_dir] * a_loc;
      }
      break;
    }
    case SeparatorType::Paraboloid: {
      Pt p_datum = a_ref.getParaboloid().getDatum();
      ReferenceFrame p_frame = a_ref.getParaboloid().getReferenceFrame();
      p_datum[a_dir] = 2.0 * a_loc - p_datum[a_dir];
      for (UnsignedIndex_t i = 0; i < 3; ++i) {
        p_frame[i][a_dir] *= -1.0;
      }
      paraboloid_m.setDatum(p_datum);
      paraboloid_m.setReferenceFrame(p_frame);
      break;
    }
    case SeparatorType::Cylinder: {
      Pt c_datum = a_ref.getCylinder().getDatum();
      ReferenceFrame c_frame = a_ref.getCylinder().getReferenceFrame();
      c_datum[a_dir] = 2.0 * a_loc - c_datum[a_dir];
      for (UnsignedIndex_t i = 0; i < 3; ++i) {
        c_frame[i][a_dir] *= -1.0;
      }
      cylinder_m.setDatum(c_datum);
      cylinder_m.setReferenceFrame(c_frame);
      break;
    }
    default: {
      std::runtime_error("SeparatorUnion type cannot reflect");
    }
  }
};

inline std::ostream& operator<<(std::ostream& out,
                                const SeparatorUnion& a_reconstruction) {
  if (a_reconstruction.type() == SeparatorUnion::SeparatorType::OnePlane) {
    out << "Plane 0: " << a_reconstruction.getPlane();
  } else if (a_reconstruction.type() ==
             SeparatorUnion::SeparatorType::TwoPlanes) {
    out << "Plane 0: " << a_reconstruction.getPlane(0)
        << "\nPlane 1: " << a_reconstruction.getPlane(1);
  } else if (a_reconstruction.type() ==
             SeparatorUnion::SeparatorType::Paraboloid) {
    out << "Paraboloid: " << a_reconstruction.getParaboloid();
  } else if (a_reconstruction.type() ==
             SeparatorUnion::SeparatorType::Cylinder) {
    out << "Cylinder: " << a_reconstruction.getCylinder();
  } else {
    throw std::runtime_error("Union type cannot be printed");
  }
  return out;
}

inline SeparatorUnion operator*(const SeparatorUnion& a_sep1,
                                const SeparatorUnion& a_sep2) {
  return a_sep1;
}

inline SeparatorUnion operator*(const double a_double,
                                const SeparatorUnion& a_sep) {
  return a_sep;
}

inline SeparatorUnion operator*(const SeparatorUnion& a_sep,
                                const double a_double) {
  return a_double * a_sep;
}

}  // namespace IRL

#endif  // IRL_SEPARATOR_UNION_TPP_
