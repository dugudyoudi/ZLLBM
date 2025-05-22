//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_noneq_extrapolation.cpp
* @author Zhengliang Liu
* @brief functions used to implement non-equilibrium extrapolation boundary condition.
*/
#include "boundary_conditions/lbm_boundary_conditions.h"
#include "./lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
void BoundaryNonEqExtrapolation2D::SetValues(const std::vector<DefReal> values) {
    if (values.size() != 2) {
        amrproject::LogManager::LogWarning("Size of given values should be 2 in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    boundary_velocity_ = {values.at(0), values.at(1)};
}
/**
 * @brief function to calculate non-equilibrium extrapolation condition for 2 dimensional models.
 * @param[in] boundary_dir boundary direction.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryNonEqExtrapolation2D::CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
    const DefMap<DefInt>& boundary_nodes, GridInfoLbmInteface* const ptr_grid_info) const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_grid_info->GetPtrToSolver().lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal rho0 = ptr_lbm_solver->GetDefaultDensity();
    std::vector<DefReal> velocity = {boundary_velocity_[0], boundary_velocity_[1]};
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();

    std::function<DefSFBitset(const DefSFBitset&)> func_find_inner = nullptr;
    amrproject::SFBitsetAux2D aux2d;
    switch (boundary_dir) {
        case amrproject::EDomainBoundaryDirection::kBoundaryXMin:
            func_find_inner = [&aux2d](const DefSFBitset& bitset) { return aux2d.FindXPos(bitset); };
            break;
        case amrproject::EDomainBoundaryDirection::kBoundaryXMax:
            func_find_inner = [&aux2d](const DefSFBitset& bitset) { return aux2d.FindXNeg(bitset); };
                break;
        case amrproject::EDomainBoundaryDirection::kBoundaryYMin:
            func_find_inner = [&aux2d](const DefSFBitset& bitset) { return aux2d.FindYPos(bitset); };
                break;
        case amrproject::EDomainBoundaryDirection::kBoundaryYMax:
            func_find_inner = [&aux2d](const DefSFBitset& bitset) { return aux2d.FindYNeg(bitset); };
                break;
        default:
            amrproject::LogManager::LogWarning("Unknown boundary type");
            break;
    }

    DefReal dt = ptr_lbm_solver->GetCollisionOperator(ptr_grid_info->GetGridLevel()).GetDtLbm();
    DefReal rho_inner = 1.;
    const DefInt nq = ptr_lbm_solver->k0NumQ_;
    std::vector<DefReal> velocity_inner(2);
    std::vector<DefReal> feq_inner(nq), feq(nq);
    DefInt iq;
    for (const auto& iter_node : boundary_nodes) {
        GridNodeLbm& node = *grid_node.at(iter_node.first);
        const GridNodeLbm& node_inner = *grid_node.at(func_find_inner(iter_node.first));
        ptr_lbm_solver->func_macro_(dt, node_inner, &rho_inner, &velocity_inner);
        ptr_lbm_solver->func_cal_feq_(rho_inner, velocity_inner, &feq_inner);
        ptr_lbm_solver->func_cal_feq_(rho0, velocity, &feq);
        for (iq = 0; iq < nq; ++iq) {
            node.f_.at(iq) = feq.at(iq) + node_inner.f_.at(iq) - feq_inner.at(iq);
        }
    }
}
void BoundaryNonEqExtrapolation3D::SetValues(const std::vector<DefReal> values) {
    if (values.size() != 3) {
        amrproject::LogManager::LogWarning("Size of given values should be 3 in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    boundary_velocity_ = {values.at(0), values.at(1), values.at(2)};
}
/**
 * @brief function to calculate  non-equilibrium extrapolation boundary condition for 3 dimensional models.
 * @param[in] boundary_dir boundary direction.
 * @param[in] boundary_nodes space filling codes of nodes on the boundary.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void BoundaryNonEqExtrapolation3D::CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
    const DefMap<DefInt>& boundary_nodes,
    GridInfoLbmInteface* const ptr_grid_info) const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_grid_info->GetPtrToSolver().lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    const DefReal rho0 = ptr_lbm_solver->GetDefaultDensity();
    std::vector<DefReal> velocity(boundary_velocity_.begin(), boundary_velocity_.end());
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_node = *ptr_grid_info->GetPtrToLbmGrid();

    std::function<DefSFBitset(const DefSFBitset&)> func_find_inner = nullptr;
    amrproject::SFBitsetAux3D aux3d;
    switch (boundary_dir) {
        case amrproject::EDomainBoundaryDirection::kBoundaryXMin:
            func_find_inner = [&aux3d](const DefSFBitset& bitset) { return aux3d.FindXPos(bitset); };
            break;
        case amrproject::EDomainBoundaryDirection::kBoundaryXMax:
            func_find_inner = [&aux3d](const DefSFBitset& bitset) { return aux3d.FindXNeg(bitset); };
                break;
        case amrproject::EDomainBoundaryDirection::kBoundaryYMin:
            func_find_inner = [&aux3d](const DefSFBitset& bitset) { return aux3d.FindYPos(bitset); };
                break;
        case amrproject::EDomainBoundaryDirection::kBoundaryYMax:
            func_find_inner = [&aux3d](const DefSFBitset& bitset) { return aux3d.FindYNeg(bitset); };
                break;
        case amrproject::EDomainBoundaryDirection::kBoundaryZMin:
            func_find_inner = [&aux3d](const DefSFBitset& bitset) { return aux3d.FindZPos(bitset); };
                break;
        case amrproject::EDomainBoundaryDirection::kBoundaryZMax:
            func_find_inner = [&aux3d](const DefSFBitset& bitset) { return aux3d.FindZNeg(bitset); };
                break;
        default:
            amrproject::LogManager::LogWarning("Unknown boundary type");
            break;
    }

    DefReal dt = ptr_lbm_solver->GetCollisionOperator(ptr_grid_info->GetGridLevel()).GetDtLbm();
    DefReal rho_inner = 1.;
    const DefInt nq = ptr_lbm_solver->k0NumQ_;
    std::vector<DefReal> velocity_inner(3);
    std::vector<DefReal> feq_inner(nq), feq(nq);
    DefInt iq;
    for (const auto& iter_node : boundary_nodes) {
        GridNodeLbm& node = *grid_node.at(iter_node.first);
        const GridNodeLbm& node_inner = *grid_node.at(func_find_inner(iter_node.first));
        ptr_lbm_solver->func_macro_(dt, node_inner, &rho_inner, &velocity_inner);
        ptr_lbm_solver->func_cal_feq_(rho_inner, velocity_inner, &feq_inner);
        ptr_lbm_solver->func_cal_feq_(rho0, velocity, &feq);
        for (iq = 0; iq < nq; ++iq) {
            node.f_.at(iq) = feq.at(iq) + node_inner.f_.at(iq) - feq_inner.at(iq);
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
