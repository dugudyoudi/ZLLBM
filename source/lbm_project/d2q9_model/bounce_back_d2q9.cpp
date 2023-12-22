//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file bounce_back_d2q9.cpp
* @author Zhengliang Liu
* @brief functions to manage boundary condition for D2Q9 model.
* @date  2023-9-30
*/
#include <string>
#include <array>
#include "d2q9_model/lbm_d2q9.h"
#include "io/log_write.h"
#include "grid/grid_manager.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to calculate bounce back boundary condition for D2Q9 model.
 * @param[in] boundary_type type of the boundary.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryBounceBackD2Q9::CalBoundaryCondition(const ELbmBoundaryType boundary_type,
    const DefMap<DefAmrIndexUint>& boundary_nodes, GridInfoLbmInteface* const ptr_grid_info) const {
    const SolverLbmD2Q9& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmD2Q9>(ptr_grid_info->ptr_solver_));
    const DefReal rho0 = lbm_solver.k0Rho_;
    const std::array<DefReal, 2>& velocity = boundary_velocity_;
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->ptr_lbm_grid_;
    switch (boundary_type) {
    case ELbmBoundaryType::kBoundaryXNeg:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFXpY0Z0] = node.f_collide_[lbm_solver.kFXnY0Z0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnY0Z0)
                * lbm_solver.k0Cx_.at(lbm_solver.kFXnY0Z0) * velocity[kXIndex];
            node.f_[lbm_solver.kFXpYpZ0] = node.f_collide_[lbm_solver.kFXnYnZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYnZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYnZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYnZ0) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXpYnZ0] = node.f_collide_[lbm_solver.kFXnYpZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYpZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYpZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYpZ0) * velocity[kYIndex]);
        }
        break;
    case ELbmBoundaryType::kBoundaryXPos:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFXnY0Z0] = node.f_collide_[lbm_solver.kFXpY0Z0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpY0Z0)
                * lbm_solver.k0Cx_.at(lbm_solver.kFXpY0Z0) * velocity[kXIndex];
            node.f_[lbm_solver.kFXnYpZ0] = node.f_collide_[lbm_solver.kFXpYnZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYnZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYnZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYnZ0) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXnYnZ0] = node.f_collide_[lbm_solver.kFXpYpZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYpZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYpZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYpZ0) * velocity[kYIndex]);
        }
        break;
    case ELbmBoundaryType::kBoundaryYNeg:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFX0YpZ0] = node.f_collide_[lbm_solver.kFX0YnZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFX0YnZ0)
                * lbm_solver.k0Cx_.at(lbm_solver.kFX0YnZ0) * velocity[kXIndex];
            node.f_[lbm_solver.kFXpYpZ0] = node.f_collide_[lbm_solver.kFXnYnZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYnZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYnZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYnZ0) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXnYpZ0] = node.f_collide_[lbm_solver.kFXpYnZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYnZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYnZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYnZ0) * velocity[kYIndex]);
        }
        break;
    case ELbmBoundaryType::kBoundaryYPos:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFX0YnZ0] = node.f_collide_[lbm_solver.kFX0YpZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFX0YpZ0)
                * lbm_solver.k0Cx_.at(lbm_solver.kFX0YpZ0) * velocity[kXIndex];
            node.f_[lbm_solver.kFXpYnZ0] = node.f_collide_[lbm_solver.kFXnYpZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYpZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYpZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYpZ0) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXnYnZ0] = node.f_collide_[lbm_solver.kFXpYpZ0]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYpZ0)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYpZ0) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYpZ0) * velocity[kYIndex]);
        }
        break;
    default:
        amrproject::LogManager::LogWarning("type of the boundary is not supported.");
        break;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject