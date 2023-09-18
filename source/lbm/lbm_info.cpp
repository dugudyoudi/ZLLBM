#include "lbm/lbm_info.h"
//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include "lbm/lbm_info.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbm {
void GridInfoLBM::SetNumberOfVecElements() {
    k0NumRealForEachNode_ = static_cast<DefAmrIndexUint>(
        std::dynamic_pointer_cast<SolverLbmDnQn>(ptr_solver_)->k0FzIndex_) + 1;
}
void GridInfoLBM::initial_grid_node(const DefSFBitset& bit_set_in) {
    std::shared_ptr<SolverLbmDnQn> ptr_solver =
        std::dynamic_pointer_cast<SolverLbmDnQn>(ptr_solver_);
    map_grid_node_.at(bit_set_in).vec_real[ptr_solver->k0RhoIndex_] = 0.;
    map_grid_node_.at(bit_set_in).vec_real[ptr_solver->k0UxIndex_] = 0.2;
    map_grid_node_.at(bit_set_in).vec_real[ptr_solver->k0UyIndex_] = 0.;
}
void SolverLbmDnQn::SolverInitial() {
    io::LogManager::LogError("Need to specify LBM DnQn Model"
     + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
}
void SolverLbmD2Q9::SolverInitial() {
    if (k0SolverDims_ != 2) {
        io::LogManager::LogError("k0SolverDims_ for LBM D2Q9 Model should be"
            "2 rather than: " + std::to_string(k0SolverDims_)
             + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    k0NumQ_ = 9;
    k0Cx_ = { 0., 1., 0., -1., 0., 1., -1., -1., 1. };
    k0Cy_ = { 0., 0., 1., 0., -1., 1., 1., -1., -1. };
    k0RhoIndex_ = k0NumQ_ * num_distribution_func_set_;
    k0UxIndex_ = k0RhoIndex_ + 1;
    k0UyIndex_ = k0UxIndex_ + 1;
    k0UzIndex_ = k0UyIndex_;
    DefAmrIndexUint k0FxIndex, k0FyIndex, k0FzIndex;
    if (bool_vec_forces_) {
        k0FxIndex = k0UzIndex_ + 1;
        k0FyIndex = k0FxIndex_ + 1;
        k0FzIndex = k0FyIndex_;
    } else {
        k0FxIndex = k0UzIndex_;
        k0FyIndex = k0UzIndex_;
        k0FzIndex = k0UzIndex_;
    }
}
}  // end namespace lbm
}  // end namespace rootproject