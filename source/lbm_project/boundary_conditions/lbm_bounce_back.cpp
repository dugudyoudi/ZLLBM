//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_bounce_back.cpp
* @author Zhengliang Liu
* @brief functions used to implement bounce back boundary condition.
* @date  2023-11-6
*/
#include "boundary_conditions/lbm_boundary_conditions.h"
#include "./lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
void BoundaryBounceBack2D::SetValues(const std::vector<DefReal> values) {
    if (values.size() != 2) {
        amrproject::LogManager::LogWarning("Size of given values should be 2 in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    boundary_velocity_ = {values.at(0), values.at(1)};
}
/**
 * @brief function to calculate bounce back boundary condition for 2 dimensional models.
 * @param[in] boundary_dir boundary direction.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryBounceBack2D::CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
    const DefMap<DefInt>& boundary_nodes,
    GridInfoLbmInteface* const ptr_grid_info) const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_grid_info->GetPtrToSolver().lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal rho0 = ptr_lbm_solver->GetDefaultDensity();
    const std::array<DefReal, 2>& velocity = boundary_velocity_;
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();
    std::vector<DefInt> indices, inverse_indices;
    GetBoundaryNInverseIndices(boundary_dir, *ptr_lbm_solver, &indices, &inverse_indices);
    DefInt i, num_q_one_direction = static_cast<DefInt>(indices.size());
    for (const auto& iter_node : boundary_nodes) {
        GridNodeLbm& node = *grid_node.at(iter_node.first);
        for (i = 0; i < num_q_one_direction; ++i) {
            node.f_.at(indices[i]) = node.f_collide_.at(inverse_indices[i])
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(inverse_indices[i])
                * (ptr_lbm_solver->k0Cx_.at(inverse_indices[i]) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(inverse_indices[i]) * velocity[kYIndex]);
        }
    }
}
void BoundaryBounceBack3D::SetValues(const std::vector<DefReal> values) {
    if (values.size() != 3) {
        amrproject::LogManager::LogWarning("Size of given values should be 3 in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    boundary_velocity_ = {values.at(0), values.at(1), values.at(2)};
}
/**
 * @brief function to calculate bounce back boundary condition for 3 dimensional models.
 * @param[in] boundary_dir boundary direction.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryBounceBack3D::CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
    const DefMap<DefInt>& boundary_nodes,
    GridInfoLbmInteface* const ptr_grid_info) const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_grid_info->GetPtrToSolver().lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal rho0 = ptr_lbm_solver->GetDefaultDensity();
    const std::array<DefReal, 3>& velocity = boundary_velocity_;
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();
    std::vector<DefInt> indices, inverse_indices;
    GetBoundaryNInverseIndices(boundary_dir, *ptr_lbm_solver, &indices, &inverse_indices);
    DefInt i, num_q_one_direction = static_cast<DefInt>(indices.size());
    for (const auto& iter_node : boundary_nodes) {
        GridNodeLbm& node = *grid_node.at(iter_node.first);
        for (i = 0; i < num_q_one_direction; ++i) {
            node.f_.at(indices[i]) = node.f_collide_.at(inverse_indices[i])
                - 6 * rho0 * ptr_lbm_solver->k0Weights_.at(inverse_indices[i])
                * (ptr_lbm_solver->k0Cx_.at(inverse_indices[i]) * velocity[kXIndex]
                + ptr_lbm_solver->k0Cy_.at(inverse_indices[i]) * velocity[kYIndex]
                + ptr_lbm_solver->k0Cz_.at(inverse_indices[i]) * velocity[kZIndex]);
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
