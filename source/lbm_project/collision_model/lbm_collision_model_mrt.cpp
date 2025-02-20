//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_collision_model_mrt.cpp
* @author Zhengliang Liu
* @brief functions used for managing MRT collision models.
* @date  2025-2-19
*/
#include "./lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to calculate relaxation time and matrix.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::CalRelaxationTime(const SolverLbmInterface& lbm_solver) {
    LbmCollisionOptInterface::CalRelaxationTime(lbm_solver);
    SetImSMMatrix(1./tau_, lbm_solver);
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] node reference to LBM node information
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtCollisionOpt::CalRelaxationTimeNode(const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) {
    LbmCollisionOptInterface::CalRelaxationTimeNode(node, lbm_solver);
    SetImSMMatrix(1./tau_, lbm_solver);
}
LbmMrtCollisionOpt::LbmMrtCollisionOpt(const DefInt i_level, const SolverLbmInterface& lbm_solver)
    : LbmCollisionOptInterface(i_level, lbm_solver) {
    LbmCollisionOptInterface::CalRelaxationTime(lbm_solver);
    LbmCollisionOptInterface::CalRelaxationTimeRatio(lbm_solver);
    matrix_m_ = lbm_solver.GetMrtMMatrix();
    matrix_im_ = lbm_solver.GetMrtImMatrix();
    const DefReal relax_tau = 1./tau_;
    diag_s_ = lbm_solver.InitialMrtSMatrix(relax_tau);
    SetImSMMatrix(relax_tau, lbm_solver);
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
 * @brief function to conduct MRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmMrtCollisionOpt::CollisionOperator(const SolverLbmInterface& lbm_solver,
    const std::vector<DefReal>& /*force*/, GridNodeLbm* const ptr_node) const {
    const DefInt num_q = lbm_solver.k0NumQ_;
    std::vector<DefReal> feq(num_q, 0.);
    DefReal f_mrt = 0.;
    lbm_solver.func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
    for (DefInt iq = 0; iq < num_q; ++iq) {
        f_mrt = 0.;
        for (DefInt is = 0; is < num_q; ++is) {
            f_mrt += matrix_im_s_m_.at(iq).at(is) * (ptr_node->f_[is] - feq[is]);
        }
        ptr_node->f_collide_[iq] = ptr_node->f_[iq] - f_mrt;
    }
}
/**
 * @brief function to calculate relaxation time and matrix.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtForceCollisionOpt::CalRelaxationTime(const SolverLbmInterface& lbm_solver) {
    LbmMrtCollisionOpt::CalRelaxationTime(lbm_solver);
    SetImDMMatrix(1./tau_, lbm_solver);
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] node reference to LBM node information
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtForceCollisionOpt::CalRelaxationTimeNode(const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) {
    LbmMrtCollisionOpt::CalRelaxationTimeNode(node, lbm_solver);
    SetImDMMatrix(1./tau_, lbm_solver);
}
LbmMrtForceCollisionOpt::LbmMrtForceCollisionOpt(const DefInt i_level, const SolverLbmInterface& lbm_solver)
    : LbmMrtCollisionOpt(i_level, lbm_solver) {
    const DefReal relax_tau = 1./tau_;
    diag_d_ = lbm_solver.InitialMrtDMatrix(relax_tau);
    SetImDMMatrix(relax_tau, lbm_solver);
}
/**
 * @brief function to calculate matrix M^{-1}DM for MRT collision.
 * @param[in] relax_tau reciprocal of relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmMrtForceCollisionOpt::SetImDMMatrix(const DefReal relax_tau, const SolverLbmInterface& lbm_solver) {
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
 * @brief function to conduct MRT collision procedure.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[in] force body force considered for the current node.
 * @param[out] lbm_solver pointer to LBM node information.
 */
void LbmMrtForceCollisionOpt::CollisionOperator(const SolverLbmInterface& lbm_solver,
    const std::vector<DefReal>& force, GridNodeLbm* const ptr_node) const {
    const DefInt num_q = lbm_solver.k0NumQ_;
    std::vector<DefReal> feq(num_q, 0.), forcing_term(num_q, 0.);
    lbm_solver.func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
    for (DefInt iq = 0; iq < num_q; ++iq) {
        forcing_term[iq] = (lbm_solver.*(lbm_solver.ptr_func_cal_force_iq_))(iq, *ptr_node, force);
    }
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
