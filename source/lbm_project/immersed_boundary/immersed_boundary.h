//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary.h
* @author Zhengliang Liu
* @brief define classes to manage immersed boundary method.
* @date  2024-2-03
*/
#ifndef SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_H_
#define SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_H_
#include <memory>
#include <vector>
#include "./lbm_interface.h"
namespace rootproject {
namespace lbmproject {
class GeometryVertexImmersedBoundary : public amrproject::GeometryVertex {
 public:
    std::array<DefReal, 3> velocity_{};
    std::array<DefReal, 3> force_{};
    DefReal area_;
};
class FsiImmersedBoundary {
 public:
    void CalculateBodyForce(const GridInfoLbmInteface& grid_info,
        std::vector<std::unique_ptr<GeometryVertexImmersedBoundary>> *ptr_geo_vertices,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const;

 protected:
    DefInt stencil_dis_ = 2;
    DefReal ib_stiffness_ = 2.;
    DefInt k0IndexForceX_ = 0, k0IndexForceY_ = k0IndexForceX_ + 1, k0IndexForceZ_ = k0IndexForceY_ + 1;

    void DirectForcingScheme2D(const DefReal dt_lbm, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const amrproject::DomainInfo& domain_info, DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
        DefMap<DefInt>* const ptr_map_node_for_ib, GeometryVertexImmersedBoundary* const cptr_geo_vertex) const;
    void DirectForcingScheme3D(const DefReal dt_lbm, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const amrproject::DomainInfo& domain_info, DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
        DefMap<DefInt>* const ptr_map_node_for_ib, GeometryVertexImmersedBoundary* const cptr_geo_vertex) const;

 private:
    DefReal StencilDisOne(DefReal dis) const;
    DefReal StencilDisTwo(DefReal dis) const;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_IMMERSED_BOUNDARY_H_
