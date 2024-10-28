//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_periodic.cpp
* @author Zhengliang Liu
* @brief functions used to implement bounce back boundary condition.
* @date  2023-11-6
*/
#include "boundary_conditions/lbm_boundary_conditions.h"
#include "./lbm_interface.h"
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
    const DefMap<DefInt>& boundary_nodes, GridInfoLbmInteface* const ptr_grid_info) const {
    const SolverLbmInterface& lbm_solver = *(dynamic_cast<SolverLbmInterface*>(ptr_grid_info->GetPtrToSolver()));
    DefSFBitset boundary_counterpart, set_coordinate;
    const amrproject::GridManager2D& grid_manager2d =
        *dynamic_cast<amrproject::GridManager2D*>(ptr_grid_info->GetPtrToParentGridManager());
    const std::array<DefSFBitset, 2>& take_xref = grid_manager2d.GetTakeXRef(),
        take_yref =  grid_manager2d.GetTakeYRef();
    switch (boundary_type) {
    case ELbmBoundaryType::kBoundaryXMin:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kXIndex);
        set_coordinate = take_xref.at(grid_manager2d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryXMax:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kXIndex);
        set_coordinate = take_xref.at(grid_manager2d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYMin:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kYIndex);
        set_coordinate = take_yref.at(grid_manager2d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYMax:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kYIndex);
        set_coordinate = take_yref.at(grid_manager2d.kRefOthers_);
        break;
    default:
        break;
    }
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();
    std::vector<DefInt> indices, inverse_indices;
    GetBoundaryNInverseIndices(boundary_type, lbm_solver, &indices, &inverse_indices);
    DefSFBitset sfbitset_counterpart;
    DefInt i, num_q_one_direction = static_cast<DefInt>(indices.size());
    for (const auto& iter_node : boundary_nodes) {
        sfbitset_counterpart = (iter_node.first & set_coordinate) | boundary_counterpart;
        if (grid_node.find(sfbitset_counterpart) != grid_node.end()) {
            GridNodeLbm& node = *grid_node.at(iter_node.first),
                node_counterpart =  *grid_node.at(sfbitset_counterpart);
            for (i = 0; i < num_q_one_direction; ++i) {
                node.f_.at(indices[i]) = node_counterpart.f_.at(indices[i]);
            }
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
    const DefMap<DefInt>& boundary_nodes,
    GridInfoLbmInteface* const ptr_grid_info) const {
    const SolverLbmInterface& lbm_solver = *(dynamic_cast<SolverLbmInterface*>(ptr_grid_info->GetPtrToSolver()));
    DefSFBitset boundary_counterpart, set_coordinate;
    const amrproject::GridManager3D& grid_manager3d =
        *dynamic_cast<amrproject::GridManager3D*>(lbm_solver.GetPtrToParentGridManager());
    const std::array<DefSFBitset, 2>& take_xref = grid_manager3d.GetTakeXRef(),
        take_yref =  grid_manager3d.GetTakeYRef(), take_zref = grid_manager3d.GetTakeZRef();
    switch (boundary_type) {
    case ELbmBoundaryType::kBoundaryXMin:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kXIndex);
        set_coordinate = take_xref.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryXMax:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kXIndex);
        set_coordinate = take_xref.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYMin:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kYIndex);
        set_coordinate = take_yref.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryYMax:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kYIndex);
        set_coordinate = take_yref.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryZMin:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMax_.at(kZIndex);
        set_coordinate = take_zref.at(grid_manager3d.kRefOthers_);
        break;
    case ELbmBoundaryType::kBoundaryZMax:
        boundary_counterpart = ptr_grid_info->k0VecBitsetDomainMin_.at(kZIndex);
        set_coordinate = take_zref.at(grid_manager3d.kRefOthers_);
        break;
    default:
        break;
    }
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();
    std::vector<DefInt> indices, inverse_indices;
    GetBoundaryNInverseIndices(boundary_type, lbm_solver, &indices, &inverse_indices);
    DefSFBitset sfbitset_counterpart;
    DefInt i, num_q_one_direction = static_cast<DefInt>(indices.size());
    for (const auto& iter_node : boundary_nodes) {
        sfbitset_counterpart = (iter_node.first & set_coordinate) | boundary_counterpart;
        if (grid_node.find(sfbitset_counterpart) != grid_node.end()) {
            GridNodeLbm& node = *grid_node.at(iter_node.first),
                node_counterpart =  *grid_node.at(sfbitset_counterpart);
            for (i = 0; i < num_q_one_direction; ++i) {
                node.f_.at(indices[i]) = node_counterpart.f_.at(indices[i]);
            }
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject