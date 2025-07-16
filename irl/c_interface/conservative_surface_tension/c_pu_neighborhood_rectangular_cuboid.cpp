#include "irl/c_interface/conservative_surface_tension/c_pu_neighborhood_rectangular_cuboid.h"

#include <cassert>

extern "C" {

    void c_PUSTNeigh_RectCub_new(c_PUSTNeigh_RectCub* a_self) {
        assert(a_self->obj_ptr == nullptr);
        a_self->obj_ptr = new IRL::PUSTNeighborhood<IRL::RectangularCuboid>;
    }

    void c_PUSTNeigh_RectCub_delete(c_PUSTNeigh_RectCub* a_self) {
        delete a_self->obj_ptr;
        a_self->obj_ptr = nullptr;
    }

    void c_PUSTNeigh_RectCub_setSize(c_PUSTNeigh_RectCub* a_self, const int* a_size){
        assert(a_self != nullptr);
        assert(a_self-> != nullptr);
        a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
    }

    void c_PUSTNeigh_RectCub_setMember(c_PUSTNeigh_RectCub* a_self,
                                        const int* a_index,
                                        const c_RectCub* a_rectangular_cuboid,
                                        const PlanarSeparator* a_planar_separator) {
        assert(a_self != nullptr);
        assert(a_self->obj_ptr != nullptr);
        assert(a_rectangular_cuboid != nullptr);
        assert(a_rectangular_cuboid->obj_ptr != nullptr);
        assert(*a_index >=0);
        assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
        a_self->obj_ptr->setMember(static_cast<IRL::UnsignedIndex_t>(*a_index),
                                    a_rectangular_cuboid->obj_ptr,
                                    a_planar_separator);
    }

    void c_PUSTNeigh_RectCub_addMember(c_PUSTNeigh_RectCub* a_self,
                                        const  c_RectCub* a_rectangular_cuboid,
                                        const PlanarSeparator* a_planar_separator) {
        assert(a_self != nullptr);
        assert(a_self->obj_ptr != nullptr);
        assert(a_rectangular_cuboid != nullptr);
        assert(a_rectangular_cuboid->obj_ptr != nullptr);
        assert(a_planar_separator != nullptr);
        a_self->obj_ptr->addMember(a_rectangular_cuboid->obj_ptr, a_planar_separator);
    }
    
    void c_PUSTNeigh_RectCub_emptyNeighborhood(c_PUSTNeigh_RectCub* a_self) {
        assert(a_self != nullptr);
        assert(a_self->obj_ptr != nullptr);
        a_self->obj_ptr->emptyNeighborhood();
    }

    void c_PUSTNeigh_RectCub_setCenterOfStencil(c_PUSTNeigh_RectCub* a_self,
                                                const int* a_center_cell_index) {
        assert(a_self != nullptr);
        assert(a_self->obj_ptr != nullptr);
        assert(*a_center_cell_index >= 0);
        assert(*a_center_cell_index < static_cast<int>(a_self->obj_ptr->size()));
        a_self->obj_ptr->setCenterOfStencil(
            static_cast<IRL::UnsignedIndex_t>(*a_center_cell_index));
    }
}