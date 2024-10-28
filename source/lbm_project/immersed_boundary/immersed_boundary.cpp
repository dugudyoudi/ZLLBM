//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary.cpp
* @author Zhengliang Liu
* @brief functions used to implement immersed boundary method.
* @date  2023-11-6
*/
#include "grid/grid_manager.h"
#include "./lbm_interface.h"
#include "immersed_boundary/immersed_boundary.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
void FsiImmersedBoundary::CalculateBodyForce(const GridInfoLbmInteface& grid_info,
    std::vector<std::unique_ptr<GeometryVertexImmersedBoundary>> *ptr_geo_vertices,
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_map_grid_nodes) const {
    const amrproject::GridManagerInterface* ptr_grid_manager = grid_info.GetPtrToParentGridManager();

    const DefInt i_level = grid_info.GetGridLevel();
    const DefInt dim = ptr_grid_manager->k0GridDims_;
    std::vector<DefSFBitset> domain_min_n_level(dim), domain_max_n_level(dim);
    const amrproject::SFBitsetAuxInterface& sfbitset_aux = *grid_info.GetPtrSFBitsetAux();
    DefMap<DefInt> map_node_for_ib;
    amrproject::DomainInfo domain_info = grid_info.GetDomainInfo();
    std::function<void(const amrproject::SFBitsetAuxInterface&,
        const amrproject::DomainInfo&, DefMap<std::unique_ptr<GridNodeLbm>>* const,
        DefMap<DefInt>* const, GeometryVertexImmersedBoundary* const)> func_lagrangian_force;
    const LbmCollisionOptInterface& collision_opt =
        dynamic_cast<SolverLbmInterface*>(grid_info.GetPtrToSolver())->GetCollisionOperator(grid_info.GetGridLevel());
    const DefReal dt_lbm =  collision_opt.GetDtLbm();
    DefReal area_factor = 1. / dt_lbm / dt_lbm;
    if (dim == 2) {
        func_lagrangian_force = [this, area_factor](const amrproject::SFBitsetAuxInterface& sfbitset_aux,
            const amrproject::DomainInfo& domain_info, DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
            DefMap<DefInt>* const ptr_map_node_for_ib, GeometryVertexImmersedBoundary* const ptr_geo_vertex) {
            DirectForcingScheme2D(area_factor, sfbitset_aux, domain_info, ptr_map_grid_nodes, ptr_map_node_for_ib,
                ptr_geo_vertex);
        };
    } else if (dim == 3) {
        area_factor/=dt_lbm;
        func_lagrangian_force = [this, area_factor] (const amrproject::SFBitsetAuxInterface& sfbitset_aux,
            const amrproject::DomainInfo& domain_info,  DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
            DefMap<DefInt>* const ptr_map_node_for_ib, GeometryVertexImmersedBoundary* const ptr_geo_vertex) {
            DirectForcingScheme3D(area_factor, sfbitset_aux, domain_info, ptr_map_grid_nodes, ptr_map_node_for_ib,
                ptr_geo_vertex);
        };
    }
    for (auto& geo_vertex : *ptr_geo_vertices) {
        func_lagrangian_force(sfbitset_aux, domain_info, ptr_map_grid_nodes, &map_node_for_ib, geo_vertex.get());
    }
}

void FsiImmersedBoundary::DirectForcingScheme2D(const DefReal area_factor,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const amrproject::DomainInfo& domain_info, DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
    DefMap<DefInt>* const ptr_map_node_for_ib, GeometryVertexImmersedBoundary* const ptr_geo_vertex) const {
    std::array<DefReal, 2> index_ref = {ptr_geo_vertex->coordinate_.at(kXIndex) / domain_info.grid_space_[kXIndex],
        ptr_geo_vertex->coordinate_.at(kYIndex) /  domain_info.grid_space_[kYIndex]};
    std::array<DefAmrLUint, 2> indices = {static_cast<DefAmrLUint>(index_ref.at(kXIndex)+kEps),
        static_cast<DefAmrLUint>(index_ref.at(kYIndex)+kEps)};
    std::array<DefReal, 2> dis_ref = {index_ref.at(kXIndex) - indices.at(kXIndex),
        index_ref.at(kYIndex) - indices.at(kYIndex)};
    DefSFBitset sfbitset_ref = static_cast<const amrproject::SFBitsetAux2D&>(sfbitset_aux).SFBitsetEncoding(indices);
    std::vector<DefSFBitset> nodes_in_region;
    DefInt valid_length = sfbitset_aux.FindNodesInPeriodicRegionCorner(sfbitset_ref, stencil_dis_,
        domain_info.periodic_min_, domain_info.periodic_max_,
        domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region);
    std::function<DefReal(const DefReal)> func_stencil = [this](const DefReal dis){
        return StencilDisOne(dis);
    };
    if (valid_length == 2) {
        func_stencil = [this](const DefReal dis){
            return StencilDisTwo(dis);
        };
    }
    DefInt vec_index_x, vec_index_y;
    DefInt total_length = stencil_dis_*2;
    DefReal u_ref = 0, v_ref = 0;
    std::vector<DefReal> vec_delta_x(valid_length*2), vec_delta_y(valid_length*2);
    DefReal delta_total;
    // calculate Lagrangian force
    for (DefInt iy = -valid_length; iy < valid_length; iy++) {
        vec_index_y = (stencil_dis_ + iy) * total_length + stencil_dis_;
        vec_delta_y.at(iy+valid_length) = func_stencil(std::fabs(dis_ref.at(kYIndex) - iy - 1));
        for (DefInt ix = -valid_length; ix < valid_length; ix++) {
            vec_index_x = vec_index_y + ix;
            vec_delta_x.at(ix+valid_length) = func_stencil(std::fabs(dis_ref.at(kXIndex) - ix - 1));
            const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
            if (sfbitset_tmp!= amrproject::SFBitsetAuxInterface::kInvalidSFbitset) {
                delta_total = vec_delta_x.at(ix+valid_length)*vec_delta_y.at(iy+valid_length);
                u_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kXIndex)*delta_total;
                v_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kYIndex)*delta_total;
                if (ptr_map_node_for_ib->find(sfbitset_tmp) == ptr_map_node_for_ib->end()) {
                    if (ptr_map_grid_nodes->at(sfbitset_tmp)->force_.size() <= k0IndexForceY_) {
                        amrproject::LogManager::LogError(
                            "Size of forces stored in LBM node is less than indices for immersed boundary forces");
                    } else {
                        // reset immersed boundary force on grid nodes as zeros
                        ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceX_) = 0.;
                        ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceY_) = 0.;
                    }
                }
            } else {
                amrproject::LogManager::LogError("Space filling code of node is invalid");
            }
        }
    }
    ptr_geo_vertex->force_.at(kXIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kXIndex) - u_ref);
    ptr_geo_vertex->force_.at(kYIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kYIndex) - v_ref);

    // distribute Lagrangian force to Eulerian grid
    const DefReal area = ptr_geo_vertex->area_ * area_factor;
    for (DefInt iy = -valid_length; iy < valid_length; iy++) {
        vec_index_y = (stencil_dis_ + iy) * total_length + stencil_dis_;
        for (DefInt ix = -valid_length; ix < valid_length; ix++) {
            vec_index_x = vec_index_y + ix;
            const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
            delta_total = vec_delta_x.at(ix+valid_length)*vec_delta_y.at(iy+valid_length)*area;
            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceX_)+=
                ptr_geo_vertex->force_.at(kXIndex)*delta_total;
            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceY_)+=
                ptr_geo_vertex->force_.at(kYIndex)*delta_total;
        }
    }
}
void FsiImmersedBoundary::DirectForcingScheme3D(const DefReal area_factor,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const amrproject::DomainInfo& domain_info, DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
    DefMap<DefInt>* const ptr_map_node_for_ib, GeometryVertexImmersedBoundary* const ptr_geo_vertex) const {
    std::array<DefReal, 3> index_ref = {ptr_geo_vertex->coordinate_.at(kXIndex) / domain_info.grid_space_[kXIndex],
        ptr_geo_vertex->coordinate_.at(kYIndex) / domain_info.grid_space_[kYIndex],
        ptr_geo_vertex->coordinate_.at(kZIndex) / domain_info.grid_space_[kZIndex]};
    std::array<DefAmrLUint, 3> indices = {static_cast<DefAmrLUint>(index_ref.at(kXIndex)+kEps),
        static_cast<DefAmrLUint>(index_ref.at(kYIndex)+kEps), static_cast<DefAmrLUint>(index_ref.at(kZIndex)+kEps)};
    std::array<DefReal, 3> dis_ref = {index_ref.at(kXIndex) - indices.at(kXIndex),
        index_ref.at(kYIndex) - indices.at(kYIndex), index_ref.at(kZIndex) - indices.at(kZIndex) };
    DefSFBitset sfbitset_ref = static_cast<const amrproject::SFBitsetAux3D&>(sfbitset_aux).SFBitsetEncoding(indices);
    std::vector<DefSFBitset> nodes_in_region;
    DefInt valid_length = sfbitset_aux.FindNodesInPeriodicRegionCorner(sfbitset_ref, stencil_dis_,
        domain_info.periodic_min_, domain_info.periodic_max_,
        domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region);
    std::function<DefReal(const DefReal)> func_stencil = [this](const DefReal dis){
        return StencilDisOne(dis);
    };
    if (valid_length == 2) {
        func_stencil = [this](const DefReal dis){
            return StencilDisTwo(dis);
        };
    }
    DefInt vec_index_x, vec_index_y, vec_index_z;
    DefInt total_length = stencil_dis_*2;
    DefReal u_ref = 0, v_ref = 0, w_ref = 0;
    std::vector<DefReal> vec_delta_x(valid_length*2), vec_delta_y(valid_length*2), vec_delta_z(valid_length*2);
    DefReal delta_total;
    for (DefInt iz = -valid_length; iz < valid_length; iz++) {
        vec_index_z = (stencil_dis_ + iz) * total_length + stencil_dis_;
        vec_delta_z.at(iz+valid_length) = func_stencil(std::fabs(dis_ref.at(kZIndex) - iz - 1));
        for (DefInt iy = -valid_length; iy < valid_length; iy++) {
            vec_index_y = (vec_index_z + iy) * total_length + stencil_dis_;
            vec_delta_y.at(iy+valid_length) = func_stencil(std::fabs(dis_ref.at(kYIndex) - iy - 1));
            for (DefInt ix = -valid_length; ix < valid_length; ix++) {
                vec_index_x = vec_index_y + ix;
                vec_delta_x.at(ix+valid_length) = func_stencil(std::fabs(dis_ref.at(kXIndex) - ix - 1));
                const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
                if (sfbitset_tmp!= amrproject::SFBitsetAuxInterface::kInvalidSFbitset) {
                    delta_total = vec_delta_x.at(ix+valid_length)
                        *vec_delta_y.at(iy+valid_length)*vec_delta_z.at(iz+valid_length);
                    u_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kXIndex)*delta_total;
                    v_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kYIndex)*delta_total;
                    w_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kZIndex)*delta_total;
                    if (ptr_map_node_for_ib->find(sfbitset_tmp) == ptr_map_node_for_ib->end()) {
                        if (ptr_map_grid_nodes->at(sfbitset_tmp)->force_.size() <= k0IndexForceZ_) {
                            amrproject::LogManager::LogError(
                                "Size of forces stored in LBM node is less than indices for immersed boundary forces");
                        } else {
                            // reset immersed boundary force on grid nodes as zeros
                            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceX_) = 0.;
                            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceY_) = 0.;
                            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceZ_) = 0.;
                        }
                    }
                } else {
                    amrproject::LogManager::LogError(
                        "Space filling code of node is invalid");
                }
            }
        }
    }
    ptr_geo_vertex->force_.at(kXIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kXIndex) - u_ref);
    ptr_geo_vertex->force_.at(kYIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kYIndex) - v_ref);
    ptr_geo_vertex->force_.at(kZIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kZIndex) - w_ref);

    const DefReal area = ptr_geo_vertex->area_ * area_factor;
    for (DefInt iz = -valid_length; iz < valid_length; iz++) {
        vec_index_z = (stencil_dis_ + iz) * total_length + stencil_dis_;
        vec_delta_z.at(iz+valid_length) = func_stencil(std::fabs(dis_ref.at(kZIndex) - iz - 1));
        for (DefInt iy = -valid_length; iy < valid_length; iy++) {
            vec_index_y = (vec_index_z + iy) * total_length + stencil_dis_;
            for (DefInt ix = -valid_length; ix < valid_length; ix++) {
                vec_index_x = vec_index_y + ix;
                const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
                delta_total = vec_delta_x.at(ix+valid_length)*vec_delta_y.at(iy+valid_length)
                    *vec_delta_z.at(iz+valid_length)*area;
                ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceX_)+=
                    ptr_geo_vertex->force_.at(kXIndex)*delta_total;
                ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceY_)+=
                    ptr_geo_vertex->force_.at(kYIndex)*delta_total;
                ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceZ_)+=
                    ptr_geo_vertex->force_.at(kZIndex)*delta_total;
            }
        }
    }
}

DefReal FsiImmersedBoundary::StencilDisOne(DefReal dist) const {
    if (std::fabs(dist) < 1.) {
        return 1. - std::fabs(dist);
    } else {
        return 0.;
    }
}
DefReal FsiImmersedBoundary::StencilDisTwo(DefReal dist) const {
    if (dist < 2.) {
        DefReal dist_abs = std::fabs(dist);
        DefReal dist_sq = dist*dist;
        if (dist_abs < 1.) {
            return 0.125 * (3. - 2. * dist_abs + sqrt(1. + 4. * dist_abs - 4. * dist_sq));
        } else {
            return 0.125 * (5. - 2. * dist_abs - sqrt(-7. + 12. * dist_abs - 4. * dist_sq));
        }
    } else {
        return 0.;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
