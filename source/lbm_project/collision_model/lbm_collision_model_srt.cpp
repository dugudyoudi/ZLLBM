//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_collision_model_general.cpp
* @author Zhengliang Liu
* @brief functions used for managing general collision models for LBM.
* @date  2023-10-30
*/
#include "./lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to conduct SRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmStrCollisionOpt::CollisionOperator(const SolverLbmInterface& lbm_solver,
    const std::vector<DefReal>& /*force*/, GridNodeLbm* const ptr_node) const {
    DefInt num_q = lbm_solver.k0NumQ_;
    std::vector<DefReal> feq(num_q, 0.);
    lbm_solver.func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
    for (DefInt iq = 0; iq < num_q; ++iq) {
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - (ptr_node->f_[iq] - feq[iq]) / tau_;
    }
}
/**
 * @brief function to conduct SRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[in] force body force considered for the current node.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmStrForceCollisionOpt::CollisionOperator(const SolverLbmInterface& lbm_solver,
    const std::vector<DefReal>& force, GridNodeLbm* const ptr_node) const {
    DefInt num_q = lbm_solver.k0NumQ_;
    std::vector<DefReal> feq(num_q, 0.);
    lbm_solver.func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
    DefReal forcing_term;
    for (DefInt iq = 0; iq < num_q; ++iq) {
        forcing_term = (lbm_solver.*(lbm_solver.ptr_func_cal_force_iq_))(iq, *ptr_node, force);
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - (ptr_node->f_[iq] - feq[iq]) / tau_
             + (1. - 0.5 / tau_) * forcing_term * dt_lbm_;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
