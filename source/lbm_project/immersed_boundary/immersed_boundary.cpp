//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file immersed_boundary.cpp
* @author Zhengliang Liu
* @brief functions used to implement immersed boundary method.
* @date  2023-11-6
*/
#include <unordered_map>
#include <utility>
#include <cstdio>
#include <string>
#include "mpi/mpi_manager.h"
#include "grid/grid_manager.h"
#include "./lbm_interface.h"
#include "immersed_boundary/immersed_boundary.h"
#include "io/log_write.h"
#include "io/input_parser.h"
namespace rootproject {
namespace lbmproject {
DefMap<DefInt> FsiImmersedBoundary::map_ib_node_for_reset_ = {};
void FsiImmersedBoundary::ClearNodesRecordForIB() {
    map_ib_node_for_reset_.clear();
}
/**
 * @brief function to compute Eulerian force distributed by the immersed boundary.
 * @param[in] grid_info information of LBM grid.
 * @param[in, out] ptr_geo_vertices pointer to vertices of a geometry coupled with immersed boundary method.
 * @param[out] ptr_map_grid_nodes point to LBM grid nodes.
 */
void FsiImmersedBoundary::CalculateBodyForce(const GridInfoLbmInteface& grid_info,
    std::unordered_map<DefSizet, std::unique_ptr<amrproject::GeometryVertex>>* const ptr_geo_vertices,
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_map_grid_nodes) {
    const amrproject::GridManagerInterface* ptr_grid_manager = grid_info.GetPtrToParentGridManager();

    const DefInt i_level = grid_info.GetGridLevel();
    const DefInt dim = ptr_grid_manager->k0GridDims_;
    std::vector<DefSFBitset> domain_min_n_level(dim), domain_max_n_level(dim);
    const amrproject::SFBitsetAuxInterface& sfbitset_aux = *grid_info.GetPtrSFBitsetAux();
    DefMap<DefInt> map_node_for_ib;
    amrproject::DomainInfo domain_info = grid_info.GetDomainInfo();
    std::function<void(const amrproject::SFBitsetAuxInterface&,
        const amrproject::DomainInfo&, DefMap<std::unique_ptr<GridNodeLbm>>* const,
        GeometryVertexImmersedBoundary* const)> func_ib_force;
    SolverLbmInterface* ptr_solver_lbm = nullptr;
    if (auto ptr = grid_info.GetPtrToSolver().lock()) {
        ptr_solver_lbm = dynamic_cast<SolverLbmInterface*>(ptr.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal dt_lbm = ptr_solver_lbm->GetCollisionOperator(grid_info.GetGridLevel()).GetDtLbm();
    const std::function<void(const DefReal, const std::vector<DefReal>&, const std::vector<DefReal>&,
        DefReal* const, std::vector<DefReal>* const)> func_macro_with_force = ptr_solver_lbm->func_macro_with_force_;

    if (dim == 2) {
        std::function<void(const DefReal, GridNodeLbm* const,
            DefReal* const, std::vector<DefReal>* const)> func_macro_without_ib_force =
            [this, func_macro_with_force, ptr_solver_lbm](const DefReal dt_lbm, GridNodeLbm* const ptr_node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            ptr_node->force_.at(k0IndexForceX_) = 0.;
            ptr_node->force_.at(k0IndexForceY_) = 0.;
            const std::vector<DefReal> force(ptr_solver_lbm->GetAllForcesForANode(*ptr_node));
            func_macro_with_force(dt_lbm, ptr_node->f_, force, ptr_rho, ptr_velocity);
        };
        func_ib_force = [this, dt_lbm, func_macro_without_ib_force](
            const amrproject::SFBitsetAuxInterface& sfbitset_aux,
            const amrproject::DomainInfo& domain_info, DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
            GeometryVertexImmersedBoundary* const ptr_geo_vertex) {
            DirectForcingScheme2D(dt_lbm, sfbitset_aux, domain_info,
                func_macro_without_ib_force, ptr_map_grid_nodes, ptr_geo_vertex);
        };
    } else if (dim == 3) {
        std::function<void(const DefReal, GridNodeLbm* const,
            DefReal* const, std::vector<DefReal>* const)> func_macro_without_ib_force =
            [this, func_macro_with_force, ptr_solver_lbm](const DefReal dt_lbm, GridNodeLbm* const ptr_node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            ptr_node->force_.at(k0IndexForceX_) = 0.;
            ptr_node->force_.at(k0IndexForceY_) = 0.;
            ptr_node->force_.at(k0IndexForceZ_) = 0.;
            const std::vector<DefReal> force(ptr_solver_lbm->GetAllForcesForANode(*ptr_node));
            func_macro_with_force(dt_lbm, ptr_node->f_, force, ptr_rho, ptr_velocity);
        };
        func_ib_force = [this, dt_lbm, func_macro_without_ib_force] (
            const amrproject::SFBitsetAuxInterface& sfbitset_aux,
            const amrproject::DomainInfo& domain_info,  DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
            GeometryVertexImmersedBoundary* const ptr_geo_vertex) {
            DirectForcingScheme3D(dt_lbm, sfbitset_aux, domain_info,
                func_macro_without_ib_force, ptr_map_grid_nodes, ptr_geo_vertex);
        };
    }
    for (auto& geo_vertex : *ptr_geo_vertices) {
        func_ib_force(sfbitset_aux, domain_info, ptr_map_grid_nodes,
            dynamic_cast<GeometryVertexImmersedBoundary*>(geo_vertex.second.get()));
    }
}
/**
 * @brief function to implement 2D direct forcing scheme.
 * @param[in] dt_lbm  time spacing of current level.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] domain_info information of computational domain.
 * @param[in] func_compute_macro function to compute macroscopic variables.
 * @param[out] ptr_map_grid_nodes point to LBM grid nodes.
 * @param[out] ptr_geo_vertex point to an immersed boundary vertex.
 */
void FsiImmersedBoundary::DirectForcingScheme2D(const DefReal dt_lbm,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux, const amrproject::DomainInfo& domain_info,
    const std::function<void(const DefReal, GridNodeLbm* const,
    DefReal* const, std::vector<DefReal>* const)>& func_compute_macro,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
    GeometryVertexImmersedBoundary* const ptr_geo_vertex) const {
    std::array<DefReal, 2> index_ref = {ptr_geo_vertex->coordinate_.at(kXIndex) / domain_info.grid_space_[kXIndex],
        ptr_geo_vertex->coordinate_.at(kYIndex) / domain_info.grid_space_[kYIndex]};
    std::array<DefAmrLUint, 2> indices = {static_cast<DefAmrLUint>(index_ref.at(kXIndex)+kEps),
        static_cast<DefAmrLUint>(index_ref.at(kYIndex)+kEps)};
    std::array<DefReal, 2> dis_ref = {index_ref.at(kXIndex) - indices.at(kXIndex),
        index_ref.at(kYIndex) - indices.at(kYIndex)};
    const amrproject::SFBitsetAux2D& sfbitset_aux2d = static_cast<const amrproject::SFBitsetAux2D&>(sfbitset_aux);
    DefSFBitset sfbitset_ref = sfbitset_aux2d.SFBitsetEncoding(indices);
    std::vector<DefSFBitset> nodes_in_region;
    std::vector<std::pair<DefAmrLUint, DefSFBitset>> nodes_overlap;
    DefInt valid_length = 0;
    if (domain_info.bool_periodic_domain_) {
        valid_length = sfbitset_aux2d.FindNodesInPeriodicRegionCornerOverlap(sfbitset_ref, stencil_dis_,
            domain_info.periodic_min_, domain_info.periodic_max_,
            domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region, &nodes_overlap);
    } else {
        valid_length = sfbitset_aux2d.FindNodesInPeriodicRegionCorner(sfbitset_ref, stencil_dis_,
            domain_info.periodic_min_, domain_info.periodic_max_,
            domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region);
    }
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
    DefReal delta_x, delta_y;
    std::vector<DefReal> vec_delta_total(4*valid_length*valid_length);
    // calculate Lagrangian force
    for (DefInt iy = -valid_length; iy < valid_length; iy++) {
        vec_index_y = (stencil_dis_ + iy) * total_length + stencil_dis_;
        delta_y = func_stencil(std::fabs(dis_ref.at(kYIndex) - iy - 1));
        for (DefInt ix = -valid_length; ix < valid_length; ix++) {
            vec_index_x = vec_index_y + ix;
            delta_x = func_stencil(std::fabs(dis_ref.at(kXIndex) - ix - 1));
            const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
            if (sfbitset_tmp!= amrproject::SFBitsetAuxInterface::kInvalidSFbitset) {
                if (ptr_map_grid_nodes->find(sfbitset_tmp) !=  ptr_map_grid_nodes->end()) {
                    if (map_ib_node_for_reset_.find(sfbitset_tmp) == map_ib_node_for_reset_.end()) {
                        if (static_cast<DefInt>(ptr_map_grid_nodes->at(sfbitset_tmp)->force_.size()) > k0IndexForceY_) {
                            // compute macro variables and reset immersed boundary force on grid nodes as zeros
                            func_compute_macro(dt_lbm, ptr_map_grid_nodes->at(sfbitset_tmp).get(),
                                &ptr_map_grid_nodes->at(sfbitset_tmp)->rho_,
                                &ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_);
                        } else {
                            amrproject::LogManager::LogError("size of forces stored in LBM node"
                                " is less than indices for immersed boundary forces");
                        }
                        map_ib_node_for_reset_.insert({sfbitset_tmp, 0});
                    }
                    vec_delta_total.at(vec_index_x) = delta_x*delta_y;
                    u_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kXIndex)*vec_delta_total.at(vec_index_x);
                    v_ref+=ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kYIndex)*vec_delta_total.at(vec_index_x);
                } else {
                    std::vector<DefReal> coordinates;
                    sfbitset_aux.SFBitsetComputeCoordinateVir(sfbitset_tmp, domain_info.grid_space_, &coordinates);
                    amrproject::LogManager::LogError("node (" + std::to_string(coordinates.at(kXIndex)) + ", "
                        + std::to_string(coordinates.at(kYIndex)) + ") does not exist.");
                }
            } else {
                amrproject::LogManager::LogError("Space filling code of node is invalid");
            }
        }
    }
    ptr_geo_vertex->force_.at(kXIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kXIndex) - u_ref);
    ptr_geo_vertex->force_.at(kYIndex) = ib_stiffness_ * (ptr_geo_vertex->velocity_.at(kYIndex) - v_ref);

    // distribute Lagrangian force to Eulerian grid
    const DefReal area = ptr_geo_vertex->area_ / dt_lbm / domain_info.grid_space_[kXIndex];
    for (DefInt iy = -valid_length; iy < valid_length; iy++) {
        vec_index_y = (stencil_dis_ + iy) * total_length + stencil_dis_;
        for (DefInt ix = -valid_length; ix < valid_length; ix++) {
            vec_index_x = vec_index_y + ix;
            const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
            const DefReal& delta_total = vec_delta_total.at(vec_index_x)*area;
            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceX_)+=
                ptr_geo_vertex->force_.at(kXIndex)*delta_total;
            ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceY_)+=
                ptr_geo_vertex->force_.at(kYIndex)*delta_total;
        }
    }
    if (domain_info.bool_periodic_domain_) {
        for (auto iter : nodes_overlap) {
            if (map_ib_node_for_reset_.find(iter.second) == map_ib_node_for_reset_.end()) {
                if (static_cast<DefInt>(ptr_map_grid_nodes->at(iter.second)->force_.size()) > k0IndexForceY_) {
                    // compute macro variables and reset immersed boundary force on grid nodes as zeros
                    func_compute_macro(dt_lbm, ptr_map_grid_nodes->at(iter.second).get(),
                        &ptr_map_grid_nodes->at(iter.second)->rho_,
                        &ptr_map_grid_nodes->at(iter.second)->velocity_);
                } else {
                    amrproject::LogManager::LogError("size of forces stored in LBM node"
                        " is less than indices for immersed boundary forces");
                }
                map_ib_node_for_reset_.insert({iter.second, 0});
            }
            const DefReal& delta_total = vec_delta_total.at(iter.first)*area;
            ptr_map_grid_nodes->at(iter.second)->force_.at(k0IndexForceX_) +=
                ptr_geo_vertex->force_.at(kXIndex)*delta_total;
            ptr_map_grid_nodes->at(iter.second)->force_.at(k0IndexForceY_) +=
                ptr_geo_vertex->force_.at(kYIndex)*delta_total;
        }
    }
}
/**
 * @brief function to implement 3D direct forcing scheme.
 * @param[in] dt_lbm  time spacing of current level.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] domain_info information of computational domain.
 * @param[in] func_compute_macro function to compute macroscopic variables.
 * @param[out] ptr_map_grid_nodes point to LBM grid nodes.
 * @param[out] ptr_map_grid_nodes point to LBM grid nodes.
 * @param[in, out] ptr_geo_vertex point to an immersed boundary vertex.
 */
void FsiImmersedBoundary::DirectForcingScheme3D(const DefReal dt_lbm,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux, const amrproject::DomainInfo& domain_info,
    const std::function<void(const DefReal, GridNodeLbm* const,
    DefReal* const, std::vector<DefReal>* const)>& func_compute_macro,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes,
    GeometryVertexImmersedBoundary* const ptr_geo_vertex) const {
    std::array<DefReal, 3> index_ref = {ptr_geo_vertex->coordinate_.at(kXIndex) / domain_info.grid_space_[kXIndex],
        ptr_geo_vertex->coordinate_.at(kYIndex) / domain_info.grid_space_[kYIndex],
        ptr_geo_vertex->coordinate_.at(kZIndex) / domain_info.grid_space_[kZIndex]};
    std::array<DefAmrLUint, 3> indices = {static_cast<DefAmrLUint>(index_ref.at(kXIndex)+kEps),
        static_cast<DefAmrLUint>(index_ref.at(kYIndex)+kEps), static_cast<DefAmrLUint>(index_ref.at(kZIndex)+kEps)};
    std::array<DefReal, 3> dis_ref = {index_ref.at(kXIndex) - indices.at(kXIndex),
        index_ref.at(kYIndex) - indices.at(kYIndex), index_ref.at(kZIndex) - indices.at(kZIndex) };
    const amrproject::SFBitsetAux3D& sfbitset_aux3d = static_cast<const amrproject::SFBitsetAux3D&>(sfbitset_aux);
    DefSFBitset sfbitset_ref = sfbitset_aux3d.SFBitsetEncoding(indices);
    std::vector<DefSFBitset> nodes_in_region;
    std::vector<std::pair<DefAmrLUint, DefSFBitset>> nodes_overlap;
    DefInt valid_length = 0;
    if (domain_info.bool_periodic_domain_) {
        valid_length = sfbitset_aux3d.FindNodesInPeriodicRegionCornerOverlap(sfbitset_ref, stencil_dis_,
            domain_info.periodic_min_, domain_info.periodic_max_,
            domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region, &nodes_overlap);
    } else {
        valid_length = sfbitset_aux3d.FindNodesInPeriodicRegionCorner(sfbitset_ref, stencil_dis_,
            domain_info.periodic_min_, domain_info.periodic_max_,
            domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region);
    }
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
    DefReal delta_x, delta_y, delta_z;
    std::vector<DefReal> vec_delta_total(8*valid_length*valid_length*valid_length);
    for (DefInt iz = -valid_length; iz < valid_length; iz++) {
        vec_index_z = (stencil_dis_ + iz) * total_length + stencil_dis_;
        delta_z = func_stencil(std::fabs(dis_ref.at(kZIndex) - iz - 1));
        for (DefInt iy = -valid_length; iy < valid_length; iy++) {
            vec_index_y = (vec_index_z + iy) * total_length + stencil_dis_;
            delta_y = func_stencil(std::fabs(dis_ref.at(kYIndex) - iy - 1));
            for (DefInt ix = -valid_length; ix < valid_length; ix++) {
                vec_index_x = vec_index_y + ix;
                delta_x = func_stencil(std::fabs(dis_ref.at(kXIndex) - ix - 1));
                const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
                if (sfbitset_tmp!= amrproject::SFBitsetAuxInterface::kInvalidSFbitset) {
                    if (ptr_map_grid_nodes->find(sfbitset_tmp) !=  ptr_map_grid_nodes->end()) {
                        if (map_ib_node_for_reset_.find(sfbitset_tmp) == map_ib_node_for_reset_.end()) {
                            if (static_cast<DefInt>(ptr_map_grid_nodes->at(sfbitset_tmp)->force_.size())
                                <= k0IndexForceZ_) {
                                amrproject::LogManager::LogError("size of forces stored in LBM node"
                                    " is less than indices for immersed boundary forces");
                            } else {
                                // compute macro variables and reset immersed boundary force on grid nodes as zeros
                                func_compute_macro(dt_lbm, ptr_map_grid_nodes->at(sfbitset_tmp).get(),
                                    &ptr_map_grid_nodes->at(sfbitset_tmp)->rho_,
                                    &ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_);
                            }
                            map_ib_node_for_reset_.insert({sfbitset_tmp, 0});
                        }
                        vec_delta_total.at(vec_index_x) = delta_x*delta_y*delta_z;
                        u_ref+=
                            ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kXIndex)*vec_delta_total.at(vec_index_x);
                        v_ref+=
                            ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kYIndex)*vec_delta_total.at(vec_index_x);
                        w_ref+=
                            ptr_map_grid_nodes->at(sfbitset_tmp)->velocity_.at(kZIndex)*vec_delta_total.at(vec_index_x);
                    } else {
                        std::vector<DefReal> coordinates;
                        sfbitset_aux.SFBitsetComputeCoordinateVir(sfbitset_tmp, domain_info.grid_space_, &coordinates);
                        amrproject::LogManager::LogError("node (" + std::to_string(coordinates.at(kXIndex)) + ", "
                            + std::to_string(coordinates.at(kYIndex))  + ", "
                            + std::to_string(coordinates.at(kZIndex)) + ") does not exist.");
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

    const DefReal area = ptr_geo_vertex->area_/dt_lbm/domain_info.grid_space_[kXIndex]/domain_info.grid_space_[kXIndex];
    for (DefInt iz = -valid_length; iz < valid_length; iz++) {
        vec_index_z = (stencil_dis_ + iz) * total_length + stencil_dis_;
        for (DefInt iy = -valid_length; iy < valid_length; iy++) {
            vec_index_y = (vec_index_z + iy) * total_length + stencil_dis_;
            for (DefInt ix = -valid_length; ix < valid_length; ix++) {
                vec_index_x = vec_index_y + ix;
                const DefSFBitset& sfbitset_tmp = nodes_in_region.at(vec_index_x);
                const DefReal& delta_total = vec_delta_total.at(vec_index_x)*area;
                ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceX_)+=
                    ptr_geo_vertex->force_.at(kXIndex)*delta_total;
                ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceY_)+=
                    ptr_geo_vertex->force_.at(kYIndex)*delta_total;
                ptr_map_grid_nodes->at(sfbitset_tmp)->force_.at(k0IndexForceZ_)+=
                    ptr_geo_vertex->force_.at(kZIndex)*delta_total;
            }
        }
    }

    if (domain_info.bool_periodic_domain_) {
        for (auto iter : nodes_overlap) {
            if (map_ib_node_for_reset_.find(iter.second) == map_ib_node_for_reset_.end()) {
                if (static_cast<DefInt>(ptr_map_grid_nodes->at(iter.second)->force_.size()) > k0IndexForceY_) {
                    // compute macro variables and reset immersed boundary force on grid nodes as zeros
                    func_compute_macro(dt_lbm, ptr_map_grid_nodes->at(iter.second).get(),
                        &ptr_map_grid_nodes->at(iter.second)->rho_,
                        &ptr_map_grid_nodes->at(iter.second)->velocity_);
                } else {
                    amrproject::LogManager::LogError("size of forces stored in LBM node"
                        " is less than indices for immersed boundary forces");
                }
                map_ib_node_for_reset_.insert({iter.second, 0});
            }
            const DefReal& delta_total = vec_delta_total.at(iter.first)*area;
            ptr_map_grid_nodes->at(iter.second)->force_.at(k0IndexForceX_) +=
                ptr_geo_vertex->force_.at(kXIndex)*delta_total;
            ptr_map_grid_nodes->at(iter.second)->force_.at(k0IndexForceY_) +=
                ptr_geo_vertex->force_.at(kYIndex)*delta_total;
            ptr_map_grid_nodes->at(iter.second)->force_.at(k0IndexForceZ_) +=
                ptr_geo_vertex->force_.at(kZIndex)*delta_total;
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
/**
 * @brief function to copy IB forces to the buffer.
 * @param[in] node_lbm  reference of a LBM node.
 * @param[out] ptr_node_buffer pointer to the buffer for a node.
 */
void FsiImmersedBoundary::CopyIBNodeToBuffer2D(const GridNodeLbm& node_lbm, char* const ptr_node_buffer) const {
    if (static_cast<DefInt>(node_lbm.force_.size()) > k0IndexForceY_) {
        constexpr DefSizet force_size = sizeof(DefReal);
        std::memcpy(ptr_node_buffer, &node_lbm.force_.at(k0IndexForceX_), force_size);
        std::memcpy(ptr_node_buffer + force_size, &node_lbm.force_.at(k0IndexForceY_), force_size);
    } else {
        amrproject::LogManager::LogError("force size is less than required for immersed boundary.");
    }
}
/**
 * @brief function to copy IB forces to the buffer.
 * @param[in] node_lbm  reference of a LBM node.
 * @param[out] ptr_node_buffer pointer to the buffer for a node.
 */
void FsiImmersedBoundary::CopyIBNodeToBuffer3D(const GridNodeLbm& node_lbm, char* const ptr_node_buffer) const {
    if (DefInt(node_lbm.force_.size()) > k0IndexForceZ_) {
        constexpr DefSizet force_size = sizeof(DefReal);
        std::memcpy(ptr_node_buffer, &node_lbm.force_.at(k0IndexForceX_), force_size);
        std::memcpy(ptr_node_buffer + force_size, &node_lbm.force_.at(k0IndexForceY_), force_size);
        std::memcpy(ptr_node_buffer + 2*force_size, &node_lbm.force_.at(k0IndexForceZ_), force_size);
    } else {
        amrproject::LogManager::LogError("force size is less than required for immersed boundary.");
    }
}
/**
 * @brief function to read IB forces from buffer and added to those stored in node.
 * @param[in] ptr_node_buffer  pointer to the buffer of a given node.
 * @param[out] ptr_node pointer to LBM node.
 */
void FsiImmersedBoundary::ReadIBNodeFromBuffer2D(
    const char* ptr_node_buffer, GridNodeLbm* const ptr_node_lbm) const {
    DefReal force_tmp = 0.;
    constexpr DefSizet force_size = sizeof(DefReal);
    std::memcpy(&force_tmp, ptr_node_buffer, force_size);
    ptr_node_lbm->force_.at(k0IndexForceX_) += force_tmp;
    std::memcpy(&force_tmp, ptr_node_buffer + force_size, force_size);
    ptr_node_lbm->force_.at(k0IndexForceY_) += force_tmp;
}
/**
 * @brief function to read IB forces from buffer and added to those stored in node.
 * @param[in] ptr_node_buffer  pointer to the buffer of a given node.
 * @param[out] ptr_node pointer to LBM node.
 */
void FsiImmersedBoundary::ReadIBNodeFromBuffer3D(
    const char* ptr_node_buffer, GridNodeLbm* const ptr_node_lbm) const {
    DefReal force_tmp = 0.;
    constexpr DefSizet force_size = sizeof(DefReal);
    std::memcpy(&force_tmp, ptr_node_buffer, force_size);
    ptr_node_lbm->force_.at(k0IndexForceX_) += force_tmp;
    std::memcpy(&force_tmp, ptr_node_buffer + force_size, force_size);
    ptr_node_lbm->force_.at(k0IndexForceY_) += force_tmp;
    std::memcpy(&force_tmp, ptr_node_buffer + 2 * force_size, force_size);
    ptr_node_lbm->force_.at(k0IndexForceZ_) += force_tmp;
}
/**
 * @brief function to send and receive node information used in immersed boundary method.
 * @param[in] dim  dimension of the immersed boundary forces.
 * @param[in] i_level refinement level.
 * @param[in] mpi_manager  class to manage mpi communication.
 * @param[in, out] ptr_map_grid_nodes pointer to container storing grid node infomation.
 */
void FsiImmersedBoundary::SendNReceiveNodesForImmersedBoundary(const DefInt dim,
    const DefInt i_level, const amrproject::MpiManager& mpi_manager,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) {
#ifdef ENABLE_MPI
    int ib_force_size = static_cast<int>(dim*sizeof(DefReal));
    std::vector<amrproject::MpiManager::BufferSizeInfo> send_buffer_info, receive_buffer_info;
    mpi_manager.SendNReceiveGridNodeBufferSize(ib_force_size,
        i_level, map_ib_node_for_mpi_send_, &send_buffer_info, &receive_buffer_info);

    std::function<void(const GridNodeLbm&, char* const)> func_copy_a_node_to_buffer;
    std::function<void(const char* , GridNodeLbm* const)> func_read_a_node_from_buffer;
    if (dim == 2) {
        func_copy_a_node_to_buffer = [this](const GridNodeLbm& node_ref, char* const ptr_node_buffer) {
            CopyIBNodeToBuffer2D(node_ref, ptr_node_buffer);
        };
        func_read_a_node_from_buffer = [this](const char* ptr_node_buffer, GridNodeLbm* const ptr_node) {
            ReadIBNodeFromBuffer2D(ptr_node_buffer, ptr_node);
        };
    } else if (dim == 3) {
        func_copy_a_node_to_buffer = [this](const GridNodeLbm& node_ref, char* const ptr_node_buffer) {
            CopyIBNodeToBuffer3D(node_ref, ptr_node_buffer);
        };
        func_read_a_node_from_buffer = [this](const char* ptr_node_buffer, GridNodeLbm* const ptr_node) {
            ReadIBNodeFromBuffer3D(ptr_node_buffer, ptr_node);
        };
    } else {
        amrproject::LogManager::LogError("dimension is not 2 or 3.");
    }

    int num_ranks = mpi_manager.GetNumOfRanks(), rank_id = mpi_manager.GetRankId();
    std::vector<std::vector<MPI_Request>> vec_vec_reqs_send(num_ranks);
    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send;
    for (int i = 1; i < num_ranks; ++i) {
        int i_rank_send = (rank_id + i) % num_ranks;
        int i_rank_receive = (rank_id - i + num_ranks)% num_ranks;
        std::unique_ptr<char[]> ptr_buffer_receive =
            mpi_manager.BlockingSendNReceiveGridNode<GridNodeLbm>(i_rank_send, i_rank_receive, ib_force_size,
            map_ib_node_for_mpi_send_, send_buffer_info, receive_buffer_info, func_copy_a_node_to_buffer,
            *ptr_map_grid_nodes, &vec_vec_reqs_send.at(i_rank_send), &vec_ptr_buffer_send);
        // decode buffer
        if (ptr_buffer_receive != nullptr) {
            DefSizet buffer_size = receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(1)
                + (receive_buffer_info.at(i_rank_receive).num_chunks_ - 1)
                *receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(0);
            const DefSizet key_size = sizeof(DefSFBitset);
            const DefSizet num_nodes = buffer_size/(sizeof(DefSFBitset) + ib_force_size);
            DefSizet position = 0;
            DefSFBitset key_code;
            const char* ptr_buffer = ptr_buffer_receive.get();
            for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
                std::memcpy(&key_code, ptr_buffer + position, key_size);
                position += key_size;
                if (ptr_map_grid_nodes->find(key_code) != ptr_map_grid_nodes->end()) {
                    if (map_ib_node_for_reset_.find(key_code) == map_ib_node_for_reset_.end()) {
                        ptr_map_grid_nodes->at(key_code)->force_.at(k0IndexForceX_) = 0.;
                        ptr_map_grid_nodes->at(key_code)->force_.at(k0IndexForceY_) = 0.;
                        if (dim == 3) {
                            ptr_map_grid_nodes->at(key_code)->force_.at(k0IndexForceZ_) = 0.;
                        }
                        map_ib_node_for_reset_.insert({key_code, 0});
                    }
                    func_read_a_node_from_buffer(ptr_buffer + position, ptr_map_grid_nodes->at(key_code).get());
                }
                position += ib_force_size;
                if (position > buffer_size) {
                    amrproject::LogManager::LogError("Buffer to store node information overflows, please check"
                        " size of node info for mpi communication");
                }
            }
        }
    }

    int i_send = 0;
    for (int i_rank = 0; i_rank < num_ranks; ++i_rank) {
        int i_rank_send = (rank_id + i_rank) % num_ranks;
        if (send_buffer_info.at(i_rank_send).bool_exist_) {
            MPI_Waitall(static_cast<int>(vec_vec_reqs_send.at(i_send).size()),
                vec_vec_reqs_send.at(i_send).data(), MPI_STATUSES_IGNORE);
            ++i_send;
        }
    }
#endif  //  ENABLE_MPI
}
/**
 * @brief function to setup geometry for mpi communication.
 * @param[in] mpi_manager  class managing MPI communication.
 * @param[in] grid_info class storting grid node information on current rank.
 */
void GeometryInfoImmersedBoundary::SetupGeometryInfo(const DefReal time,
    const amrproject::MpiManager& mpi_manager,
    const amrproject::GridInfoInterface& grid_info) {
    this->amrproject::GeometryInfoInterface::SetupGeometryInfo(time, mpi_manager, grid_info);
#ifdef ENABLE_MPI
    // find nodes near immersed boundary and on mpi outer layers
    SolverLbmInterface* ptr_solver_lbm = nullptr;
    if (auto ptr = grid_info.GetPtrToSolver().lock()) {
        ptr_solver_lbm = dynamic_cast<SolverLbmInterface*>(ptr.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal dt_lbm = ptr_solver_lbm->GetCollisionOperator(grid_info.GetGridLevel()).GetDtLbm();
    amrproject::DomainInfo domain_info = grid_info.GetDomainInfo();
    const amrproject::SFBitsetAuxInterface& sfbitset_aux = *grid_info.GetPtrSFBitsetAux();
    std::vector<DefReal> grid_space_background = sfbitset_aux.GetBackgroundGridSpacing();
    const DefInt i_level = grid_info.GetGridLevel();

    std::vector<DefSFBitset> nodes_in_region;
    DefSFBitset sfbitset_tmp;
    DefSFCodeToUint code_background, code_min = mpi_manager.GetSFBitsetMinCurrentRank().to_ullong(),
        code_max = mpi_manager.GetSFBitsetMaxCurrentRank().to_ullong();
    std::vector<DefSFCodeToUint>::iterator iter_index;
    std::vector<DefSFCodeToUint> ull_max = mpi_manager.GetSFCodeMaxAllRanks();
    int i_rank;
    for (const auto& iter_vertex : map_vertices_info_) {
        sfbitset_tmp = sfbitset_aux.SFBitsetEncodingCoordi(domain_info.grid_space_,
            {iter_vertex.second->coordinate_.at(kXIndex),
            iter_vertex.second->coordinate_.at(kYIndex), iter_vertex.second->coordinate_.at(kZIndex)});
        sfbitset_aux.FindNodesInPeriodicRegionCorner(sfbitset_tmp, stencil_dis_,
            domain_info.periodic_min_, domain_info.periodic_max_,
            domain_info.domain_min_n_level_, domain_info.domain_max_n_level_, &nodes_in_region);
        for (const auto& iter_sfbitset : nodes_in_region) {
            code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_sfbitset).to_ullong();
            if (code_background < code_min || code_background > code_max) {
                iter_index = std::lower_bound(ull_max.begin(), ull_max.end(), code_background);
                i_rank = static_cast<int>(iter_index - ull_max.begin());
                if (map_ib_node_for_mpi_send_.find(i_rank) == map_ib_node_for_mpi_send_.end()) {
                    map_ib_node_for_mpi_send_.insert(std::make_pair(i_rank, DefMap<DefInt>()));
                } else {
                    map_ib_node_for_mpi_send_.at(i_rank).insert({iter_sfbitset, 0});
                }
            }
        }
    }
#endif  //  ENABLE_MPI
}
/**
 * @brief function to read and set geometry parameters.
 * @param[in] level default geometry level.
 * @param[in, out] ptr_geo_parameters map storing geometry parameters.
 */
void GeometryInfoImmersedBoundary::ReadAndSetGeoParameters(const DefInt level,
    std::map<std::string, amrproject::ParserData>* const ptr_geo_parameters) {
    amrproject::InputParser input_parser;
    input_parser.GetValue<bool>("ib.output_ib_force", ptr_geo_parameters, &write_ib_force_);
    amrproject::GeometryInfoInterface::ReadAndSetGeoParameters(level, ptr_geo_parameters);
}
/**
 * @brief function to write time history of Lagrangian force acting on the geometry.
 * @param[in] time  current time.
 */
void GeometryInfoImmersedBoundary::WriteTimeHisLagrangianForce(const DefReal time, const DefReal dx_background) const {
    if (!write_ib_force_) {
        return;
    }
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
    FILE* fp = nullptr;

    std::string filename;
    if (GetName().empty()) {
        filename = "ib_force_for_geo_" + std::to_string(i_geo_) + "_rank" +std::to_string(rank_id) + ".txt";
    } else {
        filename = "ib_force_for_geo_" + GetName() + "_rank" +std::to_string(rank_id) + ".txt";
    }

    if (time < 1 + kEps) {
        std::remove(filename.c_str());
        fopen_s(&fp, filename.c_str(), "w");
        if (fp) {
            fprintf_s(fp, "%s %s %s %s\n", "time", "force_x", "force_y", "force_z");
        } else {
            amrproject::LogManager::LogError("Failed to open the file for writing IB forces.\n");
        }
    } else {
        fopen_s(&fp, filename.c_str(), "a");
    }
    if (fp) {
        DefReal fx = 0., fy = 0., fz = 0.;
        DefReal scale = 1. / dx_background;
        if (k0GeoDim_ == 3) {
            scale/= dx_background;
        }
        for (const auto& iter_vertex : map_vertices_info_) {
            GeometryVertexImmersedBoundary* ptr_vertex_info =
                dynamic_cast<GeometryVertexImmersedBoundary*>(iter_vertex.second.get());
            fx += ptr_vertex_info->force_.at(kXIndex) * ptr_vertex_info->area_ * scale;
            fy += ptr_vertex_info->force_.at(kYIndex) * ptr_vertex_info->area_ * scale;
            fz += ptr_vertex_info->force_.at(kZIndex) * ptr_vertex_info->area_ * scale;
        }
        fprintf_s(fp, "%f %f %f %f\n", time - 1., fx, fy, fz);
        fclose(fp);
    } else {
        amrproject::LogManager::LogError("Failed to open the file for writing IB forces.\n");
    }
}
/**
 * @brief function to instantiate class of geometry information based on geometry type.
 * @param[in] dims dimension of the geometry.
 * @param[in] geo_type type of geometry.
 * @return share pointer of geometry information
 */
std::shared_ptr<amrproject::GeometryInfoInterface> GeoIBTypeReader::ReadGeoType(
    const DefInt dims, const std::string& geo_type) const {

    if (geo_type == "origin_ib") {
        return std::make_shared<GeometryInfoImmersedBoundary>(dims);
    } else {
        return amrproject::GeoTypeReader::ReadGeoType(dims, geo_type);
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject


