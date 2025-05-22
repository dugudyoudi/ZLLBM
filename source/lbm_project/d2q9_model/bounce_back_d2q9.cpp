//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
 * @param[in] boundary_dir  boundary direction.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryBounceBackD2Q9::CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
    const DefMap<DefInt>& boundary_nodes, GridInfoLbmInteface* const ptr_grid_info) const {
    SolverLbmD2Q9* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_grid_info->GetPtrToSolver().lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmD2Q9*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal rho0 = ptr_lbm_solver->GetDefaultDensity();
    const std::array<DefReal, 2>& velocity = boundary_velocity_;
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();
    switch (boundary_dir) {
    case amrproject::EDomainBoundaryDirection::kBoundaryXMin:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[ptr_lbm_solver->kFXpY0Z0_] = node.f_collide_[ptr_lbm_solver->kFXnY0Z0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXnY0Z0_)
                * ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXnY0Z0_) * velocity[kXIndex];
            node.f_[ptr_lbm_solver->kFXpYpZ0_] = node.f_collide_[ptr_lbm_solver->kFXnYnZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXnYnZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXnYnZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXnYnZ0_) * velocity[kYIndex]);
            node.f_[ptr_lbm_solver->kFXpYnZ0_] = node.f_collide_[ptr_lbm_solver->kFXnYpZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXnYpZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXnYpZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXnYpZ0_) * velocity[kYIndex]);
        }
        break;
    case amrproject::EDomainBoundaryDirection::kBoundaryXMax:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[ptr_lbm_solver->kFXnY0Z0_] = node.f_collide_[ptr_lbm_solver->kFXpY0Z0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXpY0Z0_)
                * ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXpY0Z0_) * velocity[kXIndex];
            node.f_[ptr_lbm_solver->kFXnYpZ0_] = node.f_collide_[ptr_lbm_solver->kFXpYnZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXpYnZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXpYnZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXpYnZ0_) * velocity[kYIndex]);
            node.f_[ptr_lbm_solver->kFXnYnZ0_] = node.f_collide_[ptr_lbm_solver->kFXpYpZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXpYpZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXpYpZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXpYpZ0_) * velocity[kYIndex]);
        }
        break;
    case amrproject::EDomainBoundaryDirection::kBoundaryYMin:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[ptr_lbm_solver->kFX0YpZ0_] = node.f_collide_[ptr_lbm_solver->kFX0YnZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFX0YnZ0_)
                * ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFX0YnZ0_) * velocity[kXIndex];
            node.f_[ptr_lbm_solver->kFXpYpZ0_] = node.f_collide_[ptr_lbm_solver->kFXnYnZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXnYnZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXnYnZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXnYnZ0_) * velocity[kYIndex]);
            node.f_[ptr_lbm_solver->kFXnYpZ0_] = node.f_collide_[ptr_lbm_solver->kFXpYnZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXpYnZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXpYnZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXpYnZ0_) * velocity[kYIndex]);
        }
        break;
    case amrproject::EDomainBoundaryDirection::kBoundaryYMax:
        for (const auto& iter_node : boundary_nodes) {
            GridNodeLbm& node = *grid_node.at(iter_node.first);
            node.f_[ptr_lbm_solver->kFX0YnZ0_] = node.f_collide_[ptr_lbm_solver->kFX0YpZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFX0YpZ0_)
                * ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFX0YpZ0_) * velocity[kXIndex];
            node.f_[ptr_lbm_solver->kFXpYnZ0_] = node.f_collide_[ptr_lbm_solver->kFXnYpZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXnYpZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXnYpZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXnYpZ0_) * velocity[kYIndex]);
            node.f_[ptr_lbm_solver->kFXnYnZ0_] = node.f_collide_[ptr_lbm_solver->kFXpYpZ0_]
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(ptr_lbm_solver->kFXpYpZ0_)
                * (ptr_lbm_solver->k0Cx_.at(ptr_lbm_solver->kFXpYpZ0_) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(ptr_lbm_solver->kFXpYpZ0_) * velocity[kYIndex]);
        }
        break;
    default:
        amrproject::LogManager::LogWarning("type of the boundary is not supported.");
        break;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
