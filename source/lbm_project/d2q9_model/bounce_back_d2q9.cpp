//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
    const DefMap<DefInt>& boundary_nodes, GridInfoLbmInteface* const ptr_grid_info) const {
    const SolverLbmD2Q9& lbm_solver = *(dynamic_cast<SolverLbmD2Q9*>(ptr_grid_info->GetPtrSolver()));
    const DefReal rho0 = lbm_solver.GetDefaultDensity();
    const std::array<DefReal, 2>& velocity = boundary_velocity_;
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->ptr_lbm_grid_nodes_;
    switch (boundary_type) {
    case ELbmBoundaryType::kBoundaryXMin:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFXpY0Z0_] = node.f_collide_[lbm_solver.kFXnY0Z0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnY0Z0_)
                * lbm_solver.k0Cx_.at(lbm_solver.kFXnY0Z0_) * velocity[kXIndex];
            node.f_[lbm_solver.kFXpYpZ0_] = node.f_collide_[lbm_solver.kFXnYnZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYnZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYnZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYnZ0_) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXpYnZ0_] = node.f_collide_[lbm_solver.kFXnYpZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYpZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYpZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYpZ0_) * velocity[kYIndex]);
        }
        break;
    case ELbmBoundaryType::kBoundaryXMax:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFXnY0Z0_] = node.f_collide_[lbm_solver.kFXpY0Z0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpY0Z0_)
                * lbm_solver.k0Cx_.at(lbm_solver.kFXpY0Z0_) * velocity[kXIndex];
            node.f_[lbm_solver.kFXnYpZ0_] = node.f_collide_[lbm_solver.kFXpYnZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYnZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYnZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYnZ0_) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXnYnZ0_] = node.f_collide_[lbm_solver.kFXpYpZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYpZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYpZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYpZ0_) * velocity[kYIndex]);
        }
        break;
    case ELbmBoundaryType::kBoundaryYMin:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFX0YpZ0_] = node.f_collide_[lbm_solver.kFX0YnZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFX0YnZ0_)
                * lbm_solver.k0Cx_.at(lbm_solver.kFX0YnZ0_) * velocity[kXIndex];
            node.f_[lbm_solver.kFXpYpZ0_] = node.f_collide_[lbm_solver.kFXnYnZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYnZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYnZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYnZ0_) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXnYpZ0_] = node.f_collide_[lbm_solver.kFXpYnZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYnZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYnZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYnZ0_) * velocity[kYIndex]);
        }
        break;
    case ELbmBoundaryType::kBoundaryYMax:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[lbm_solver.kFX0YnZ0_] = node.f_collide_[lbm_solver.kFX0YpZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFX0YpZ0_)
                * lbm_solver.k0Cx_.at(lbm_solver.kFX0YpZ0_) * velocity[kXIndex];
            node.f_[lbm_solver.kFXpYnZ0_] = node.f_collide_[lbm_solver.kFXnYpZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXnYpZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXnYpZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXnYpZ0_) * velocity[kYIndex]);
            node.f_[lbm_solver.kFXnYnZ0_] = node.f_collide_[lbm_solver.kFXpYpZ0_]
                - 6 * rho0 * lbm_solver.k0Weights_.at(lbm_solver.kFXpYpZ0_)
                * (lbm_solver.k0Cx_.at(lbm_solver.kFXpYpZ0_) * velocity[kXIndex]
                + lbm_solver.k0Cy_.at(lbm_solver.kFXpYpZ0_) * velocity[kYIndex]);
        }
        break;
    default:
        amrproject::LogManager::LogWarning("type of the boundary is not supported.");
        break;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject