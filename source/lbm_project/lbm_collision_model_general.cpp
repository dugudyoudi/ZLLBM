//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_collision_model_general.cpp
* @author Zhengliang Liu
* @brief functions used for managing general collision models for LBM.
* @date  2023-10-30
*/
#include "lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
void GridInfoLbmInteface::SetCollisionOperator() {
    const SolverLbmInterface& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_));
    switch (k0CollisionOperatorType_) {
    case ELbmCollisionOperatorType::kLbmSrt: {
        if (bool_forces_) {
            ptr_collision_operator_ = LbmSrtForceCollisionOptCreator(lbm_solver);
        } else {
            ptr_collision_operator_ = LbmSrtCollisionOptCreator(lbm_solver);
        }
        break;
    }
    default:
        amrproject::LogManager::LogError("Collision operator is not defined in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        break;
    }
}
LbmCollisionOptInterface::LbmCollisionOptInterface(const SolverLbmInterface& lbm_solver) {
    viscosity_lbm_ = lbm_solver.k0LbmViscosity_;
}
void LbmStrCollisionOpt::CalRelaxationTime() {
    tau_srt_ = viscosity_lbm_ * SolverLbmInterface::kCsSqReciprocal / dt_lbm_ + 0.5;
}
void LbmStrCollisionOpt::CalRelaxationTimeNode(const GridNodeLbm& node) {
    tau_srt_ = viscosity_lbm_ * SolverLbmInterface::kCsSqReciprocal / dt_lbm_ + 0.5;
}
void LbmStrCollisionOpt::CollisionOperator(
    const SolverLbmInterface& lbm_solver, GridNodeLbm* const ptr_node) const {
    DefAmrIndexUint num_q = lbm_solver.k0NumQ_;
    std::vector<DefReal> feq(num_q, 0.);
    lbm_solver.func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
    for (DefAmrIndexUint iq = 0; iq < num_q; ++iq) {
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - (ptr_node->f_[iq] - feq[iq]) / tau_srt_;
    }
}
void LbmStrForceCollisionOpt::CollisionOperator(
    const SolverLbmInterface& lbm_solver, GridNodeLbm* const ptr_node) const {
    DefAmrIndexUint num_q = lbm_solver.k0NumQ_;
    std::vector<DefReal> feq(num_q, 0.);
    lbm_solver.func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
    DefReal forcing_term;
    for (DefAmrIndexUint iq = 0; iq < num_q; ++iq) {
        forcing_term = (lbm_solver.*(lbm_solver.ptr_func_cal_force_iq_))(iq, *ptr_node);
        ptr_node->f_collide_[iq] = (1. - 1./ tau_srt_) * ptr_node->f_[iq] + feq[iq] / tau_srt_
             + (1. - 0.5 / tau_srt_) * forcing_term * dt_lbm_;

    }
}
}  // end namespace lbmproject
}  // end namespace rootproject