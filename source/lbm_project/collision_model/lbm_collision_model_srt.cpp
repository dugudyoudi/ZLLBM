//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
 * @param[in] feq equilibrium distribution function.
 * @param[out] ptr_node pointer to LBM node information.
 */
void LbmStrCollisionOpt::CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
    GridNodeLbm* const ptr_node) const {
    for (DefInt iq = 0; iq < num_q; ++iq) {
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - (ptr_node->f_[iq] - feq[iq]) / tau_eff_;
    }
}
/**
 * @brief function to conduct SRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[in] feq equilibrium distribution function.
 * @param[in] forcing_term foring term distribution function.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmStrCollisionOpt::CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
    const std::vector<DefReal>& forcing_term, GridNodeLbm* const ptr_node) const {
    for (DefInt iq = 0; iq < num_q; ++iq) {
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - (ptr_node->f_[iq] - feq[iq]) / tau_eff_
            + (1. - 0.5 / tau_eff_) * forcing_term[iq] * dt_lbm_;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
