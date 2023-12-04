//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_periodic.cpp
* @author Zhengliang Liu
* @brief functions used to implement bounce back boundary condition.
* @date  2023-11-6
*/
#include "boundary_conditions/lbm_boundary_conditions.h"
#include "lbm_interface.h"
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to calculate periodic boundary condition for 2 dimensional models.
 * @param[in] boundary_type type of the boundary.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryPeriodic2D::CalBoundaryCondition(const ELbmBoundaryType boundary_type,
    const DefMap<DefAmrIndexUint>& boundary_nodes,
    GridInfoLbmInteface* const ptr_grid_info) const {
    const SolverLbmInterface& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmInterface>(ptr_grid_info->ptr_solver_));
    DefSFBitset boundary_counterpart, set_coordinate;
    const amrproject::GridManager2D& grid_manager2d =
        *dynamic_cast<amrproject::GridManager2D*>(ptr_grid_info->ptr_solver_->ptr_grid_manager_);
    switch (boundary_type) {
    case ELbmBoundaryType::kBoundaryXNeg:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kXIndex);
        set_coordinate = grid_manager2d.k0SFBitsetTakeXRef_.at(grid_manager2d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryXPos:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kXIndex);
        set_coordinate = grid_manager2d.k0SFBitsetTakeXRef_.at(grid_manager2d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYNeg:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kYIndex);
        set_coordinate = grid_manager2d.k0SFBitsetTakeYRef_.at(grid_manager2d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYPos:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kYIndex);
        set_coordinate = grid_manager2d.k0SFBitsetTakeYRef_.at(grid_manager2d.kRefOthers_);
        break;
    default:
        break;
    }
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->ptr_lbm_grid_;
    std::vector<DefAmrIndexUint> indices, inverse_indices;
    GetBoundaryNInverseIndices(boundary_type, lbm_solver, &indices, &inverse_indices);
    DefSFBitset sfbitset_counterpart;
    DefAmrIndexUint i, num_q_one_direction = static_cast<DefAmrIndexUint>(indices.size());
    for (const auto& iter_node : boundary_nodes) {
        sfbitset_counterpart = (iter_node.first & set_coordinate) | boundary_counterpart;
        if (grid_node.find(sfbitset_counterpart) != grid_node.end()) {
            GridNodeLbm& node = *grid_node.at(iter_node.first),
                node_counterpart =  *grid_node.at(sfbitset_counterpart);
            for (i = 0; i < num_q_one_direction; ++i) {
                node.f_.at(indices[i]) = node_counterpart.f_.at(indices[i]);
                node.f_collide_.at(indices[i]) = node_counterpart.f_collide_.at(indices[i]);
            }
        } else {
            std::array<DefReal, 2> coordinates, coordinates2,
                grid_space = {ptr_grid_info->grid_space_[kXIndex], ptr_grid_info->grid_space_[kYIndex]};
            grid_manager2d.SFBitsetComputeCoordinate(iter_node.first, grid_space, &coordinates);
            grid_manager2d.SFBitsetComputeCoordinate(sfbitset_counterpart, grid_space, &coordinates);
            amrproject::LogManager::LogError("Node (" + std::to_string(coordinates2[kXIndex]) + ", "
                + std::to_string(coordinates2[kYIndex]) + ") at periodic boundary is not found for node ("
                + std::to_string(coordinates[kXIndex]) + ", "+ std::to_string(coordinates[kYIndex]) + ") in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
    }
}
/**
 * @brief function to calculate periodic boundary condition for 3 dimensional models.
 * @param[in] boundary_type type of the boundary.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryPeriodic3D::CalBoundaryCondition(const ELbmBoundaryType boundary_type,
    const DefMap<DefAmrIndexUint>& boundary_nodes,
    GridInfoLbmInteface* const ptr_grid_info) const {
    const SolverLbmInterface& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmInterface>(ptr_grid_info->ptr_solver_));
    DefSFBitset boundary_counterpart, set_coordinate;
    const amrproject::GridManager3D& grid_manager3d =
        *dynamic_cast<amrproject::GridManager3D*>(ptr_grid_info->ptr_solver_->ptr_grid_manager_);
    switch (boundary_type) {
    case ELbmBoundaryType::kBoundaryXNeg:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kXIndex);
        set_coordinate = grid_manager3d.k0SFBitsetTakeXRef_.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryXPos:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kXIndex);
        set_coordinate = grid_manager3d.k0SFBitsetTakeXRef_.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYNeg:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kYIndex);
        set_coordinate = grid_manager3d.k0SFBitsetTakeYRef_.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYPos:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kYIndex);
        set_coordinate = grid_manager3d.k0SFBitsetTakeYRef_.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryZNeg:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kZIndex);
        set_coordinate = grid_manager3d.k0SFBitsetTakeZRef_.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryZPos:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kZIndex);
        set_coordinate = grid_manager3d.k0SFBitsetTakeZRef_.at(grid_manager3d.kRefOthers_);
        break;
    default:
        break;
    }
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->ptr_lbm_grid_;
    std::vector<DefAmrIndexUint> indices, inverse_indices;
    GetBoundaryNInverseIndices(boundary_type, lbm_solver, &indices, &inverse_indices);
    DefSFBitset sfbitset_counterpart;
    DefAmrIndexUint i, num_q_one_direction = static_cast<DefAmrIndexUint>(indices.size());
    for (const auto& iter_node : boundary_nodes) {
        sfbitset_counterpart = (iter_node.first & set_coordinate) | boundary_counterpart;
        if (grid_node.find(sfbitset_counterpart) != grid_node.end()) {
            GridNodeLbm& node = *grid_node.at(iter_node.first),
                node_counterpart =  *grid_node.at(sfbitset_counterpart);
            for (i = 0; i < num_q_one_direction; ++i) {
                node.f_.at(indices[i]) = node_counterpart.f_.at(indices[i]);
                node.f_collide_.at(indices[i]) = node_counterpart.f_collide_.at(indices[i]);
            }
        } else {
            std::array<DefReal, 3> coordinates, coordinates2,
                grid_space = {ptr_grid_info->grid_space_[kXIndex], ptr_grid_info->grid_space_[kYIndex],
                ptr_grid_info->grid_space_[kZIndex]};
            grid_manager3d.SFBitsetComputeCoordinate(iter_node.first, grid_space, &coordinates);
            grid_manager3d.SFBitsetComputeCoordinate(sfbitset_counterpart, grid_space, &coordinates);
            amrproject::LogManager::LogError("Node (" + std::to_string(coordinates2[kXIndex]) + ", "
                + std::to_string(coordinates2[kYIndex])  + ", " + std::to_string(coordinates2[kZIndex])
                + ") at periodic boundary is not found for node ("
                + std::to_string(coordinates[kXIndex]) + ", "+ std::to_string(coordinates[kYIndex])
                + ", "+ std::to_string(coordinates[kZIndex]) + ") in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject