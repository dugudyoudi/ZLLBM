//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_collision_model_mrt.cpp
* @author Zhengliang Liu
* @brief functions used for managing MRT collision models.
*/
#include "./lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to calculate relaxation time and matrix.
 * @param[in] vis_lbm viscosity scaling for LBM solver.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::CalRelaxationTime(const DefReal vis_lbm, const SolverLbmInterface& lbm_solver) {
    LbmCollisionOptInterface::CalRelaxationTime(vis_lbm, lbm_solver);
    const DefReal relax_tau = 1./tau0_;
    SetImSMMatrix(relax_tau, lbm_solver);
    SetImDMMatrix(relax_tau, lbm_solver);
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] tau_eff effective relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::SetEffectiveRelaxationTime(
    const DefReal tau_eff, const SolverLbmInterface& lbm_solver) {
    tau_eff_ = tau_eff;
    SetImSMMatrix(1./tau_eff, lbm_solver);
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] tau_eff effective relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::SetEffectiveRelaxationTimeForForcing(
    const DefReal tau_eff, const SolverLbmInterface& lbm_solver) {
    tau_eff_ = tau_eff;
    const DefReal relax_tau = 1./tau_eff;
    SetImSMMatrix(relax_tau, lbm_solver);
    SetImDMMatrix(relax_tau, lbm_solver);
}
LbmMrtCollisionOpt::LbmMrtCollisionOpt(const DefInt i_level, const SolverLbmInterface& lbm_solver)
    : LbmCollisionOptInterface(i_level, lbm_solver) {
    LbmCollisionOptInterface::CalRelaxationTime(viscosity_lbm_, lbm_solver);
    LbmCollisionOptInterface::CalRelaxationTimeRatio(viscosity_lbm_, lbm_solver);
    matrix_m_ = lbm_solver.GetMrtMMatrix();
    matrix_im_ = lbm_solver.GetMrtImMatrix();
    const DefReal relax_tau = 1./tau0_;
    diag_s_ = lbm_solver.InitialMrtSMatrix(relax_tau);
    SetImSMMatrix(relax_tau, lbm_solver);
    diag_d_ = lbm_solver.InitialMrtDMatrix(relax_tau);
    SetImDMMatrix(relax_tau, lbm_solver);
}
/**
 * @brief function to calculate matrix M^{-1}SM for MRT collision.
 * @param[in] relax_tau reciprocal of relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::SetImSMMatrix(const DefReal relax_tau, const SolverLbmInterface& lbm_solver) {
    lbm_solver.UpdateMrtSMatrix(relax_tau, &diag_s_);
    const DefInt num_q = lbm_solver.k0NumQ_;
    std::vector<std::vector<DefReal>> diag_s_m(num_q, std::vector<DefReal>(num_q, 0.));
    for (DefInt j = 0; j < num_q; ++j) {
        for (DefInt i = 0; i < num_q; ++i) {
            for (DefInt k = 0; k < num_q; ++k) {
                diag_s_m[i][j] = diag_s_m[i][j] + diag_s_[i][k] * matrix_m_[k][j];
            }
        }
    }
    matrix_im_s_m_.assign(num_q, std::vector<DefReal>(num_q, 0.));
    for (DefInt j = 0; j < num_q; ++j) {
        for (DefInt i = 0; i < num_q; ++i) {
            for (DefInt k = 0; k < num_q; ++k) {
                matrix_im_s_m_[i][j] = matrix_im_s_m_[i][j] + matrix_im_[i][k] * diag_s_m[k][j];
            }
        }
    }
}
/**
 * @brief function to calculate matrix M^{-1}DM for MRT collision.
 * @param[in] relax_tau reciprocal of relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::SetImDMMatrix(const DefReal relax_tau, const SolverLbmInterface& lbm_solver) {
    lbm_solver.UpdateMrtDMatrix(relax_tau, &diag_d_);
    const DefInt num_q = lbm_solver.k0NumQ_;
    std::vector<std::vector<DefReal>> diag_d_m(num_q, std::vector<DefReal>(num_q, 0.));
    for (DefInt j = 0; j < num_q; ++j) {
        for (DefInt i = 0; i < num_q; ++i) {
            for (DefInt k = 0; k < num_q; ++k) {
                diag_d_m[i][j] = diag_d_m[i][j] + diag_d_[i][k] * matrix_m_[k][j];
            }
        }
    }
    matrix_im_d_m_.assign(num_q, std::vector<DefReal>(num_q, 0.));
    for (DefInt j = 0; j < num_q; ++j) {
        for (DefInt i = 0; i < num_q; ++i) {
            for (DefInt k = 0; k < num_q; ++k) {
                matrix_im_d_m_[i][j] = matrix_im_d_m_[i][j] + matrix_im_[i][k] * diag_d_m[k][j];
            }
        }
    }
}
/**
 * @brief function to conduct SRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[in] feq equilibrium distribution function.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmMrtCollisionOpt::CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
    GridNodeLbm* const ptr_node) const {
    DefReal f_mrt = 0.;
    for (DefInt iq = 0; iq < num_q; ++iq) {
        f_mrt = 0.;
        for (DefInt is = 0; is < num_q; ++is) {
            f_mrt += matrix_im_s_m_.at(iq).at(is) * (ptr_node->f_[is] - feq[is]);
        }
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - f_mrt;
    }
}
/**
 * @brief function to conduct SRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[in] feq equilibrium distribution function.
 * @param[in] forcing_term foring term distribution function.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmMrtCollisionOpt::CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
    const std::vector<DefReal>& forcing_term, GridNodeLbm* const ptr_node) const {
    DefReal f_mrt = 0.;
    for (DefInt iq = 0; iq < num_q; ++iq) {
        f_mrt = 0.;
        for (DefInt is = 0; is < num_q; ++is) {
            f_mrt += matrix_im_s_m_.at(iq).at(is) * (ptr_node->f_[is] - feq[is])
                - matrix_im_d_m_.at(iq).at(is) * forcing_term[is] * dt_lbm_;
        }
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - f_mrt;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
