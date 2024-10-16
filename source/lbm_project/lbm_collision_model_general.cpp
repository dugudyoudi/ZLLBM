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
void GridInfoLbmInteface::SetCollisionOperator() {
    const SolverLbmInterface& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_));
    switch (k0CollisionOperatorType_) {
    case ELbmCollisionOperatorType::kLbmSrt: {
        if (lbm_solver.bool_forces_) {
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
    viscosity_lbm_ = lbm_solver.GetDefaultViscosity();
}
/**
 * @brief function to calculate relaxation time ratio between coarse and fine grid.
 */
void LbmCollisionOptInterface::CalRelaxationTimeRatio() {
    DefReal relax_tau_f = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ * 2. + 0.5;
    DefReal relax_tau_c = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
    tau_collision_c2f_ = 0.5 * (relax_tau_f - 1)/(relax_tau_c - 1);
    tau_stream_c2f_ = 0.5 * relax_tau_f/relax_tau_c;
    relax_tau_f = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
    relax_tau_c = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ / 2 + 0.5;
    tau_collision_f2c_ = 2. * (relax_tau_c - 1.)/(relax_tau_f - 1.);
    tau_stream_f2c_ = 2.* relax_tau_c/relax_tau_f;
}
/**
 * @brief function to convert post collision distribution functions on coarse grid to fine grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (coarse grid).
 * @param[in] lbm_solver reference of lbm solver.
 * @param[in] f_collide_coarse post collision distribution functions of a coarse grid.
 * @param[out] ptr_f_collide_fine pointer to post collision distribution functions of a fine grid.
 */
void LbmCollisionOptInterface::PostCollisionCoarse2Fine(const DefReal dt_lbm, const std::vector<DefReal>& feq,
    const std::vector<DefReal>& f_collide_coarse, std::vector<DefReal>* const ptr_f_collide_fine) {
    DefSizet num_q = feq.size();
    ptr_f_collide_fine->resize(num_q);
    //  DefReal relax_tau_f = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm * 2. + 0.5;
    //  DefReal relax_tau_c = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm + 0.5;
    //  DefReal tau_ratio = 0.5 * (relax_tau_f - 1)/ (relax_tau_c - 1);

    for (int iq = 0; iq < num_q; ++iq) {
        ptr_f_collide_fine->at(iq) = feq.at(iq) + tau_collision_c2f_ * (f_collide_coarse[iq] - feq.at(iq));
    }
}
/**
 * @brief function to convert post collision distribution functions on fine grid to coarse grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (fine grid).
 * @param[in] feq equilibrium distribution functions.
 * @param[in] node_fine  a node on fine grid.
 * @param[out] ptr_node_coarse pointer to a node on coarse grid.
 */
void LbmCollisionOptInterface::PostCollisionFine2Coarse(const DefReal dt_lbm, const std::vector<DefReal>& feq,
    const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) {
    ptr_node_coarse->rho_ = node_fine.rho_;
    ptr_node_coarse->velocity_ = node_fine.velocity_;
    DefSizet num_q = feq.size();
    ptr_node_coarse->f_collide_.resize(num_q);
    // DefReal relax_tau_f = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm + 0.5;
    // DefReal relax_tau_c = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm / 2 + 0.5;
    //  DefReal tau_ratio = 2 * (relax_tau_c - 1)/ (relax_tau_f - 1);

    for (int iq = 0; iq < num_q; ++iq) {
        ptr_node_coarse->f_collide_[iq] = feq.at(iq) + tau_collision_f2c_ * (node_fine.f_collide_[iq] - feq.at(iq));
    }
}
/**
 * @brief function to convert post collision distribution functions on coarse grid to fine grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (coarse grid).
 * @param[in] feq reference of lbm solver.
 * @param[in] lbm_solver reference of lbm solver.
 * @param[in] ptr_func_cal_force_iq pointer to function calculating lbm forcing term.
 * @param[in] node_coarse a node on coarse grid.
 * @param[out] ptr_f_collide_fine pointer to post collision distribution functions of a fine grid.
 */
void LbmCollisionOptInterface::PostCollisionCoarse2FineForce(const DefReal dt_lbm, const std::vector<DefReal>& feq,
    const SolverLbmInterface& lbm_solver,
    DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
    const GridNodeLbm& node_coarse, std::vector<DefReal>* const ptr_f_collide_fine) {
    DefReal force_tmp;
    DefSizet num_q = feq.size();
    ptr_f_collide_fine->resize(num_q);
    for (int iq = 0; iq <num_q; ++iq) {
        ptr_f_collide_fine->at(iq) = feq.at(iq) + tau_collision_c2f_ * (node_coarse.f_collide_[iq] - feq.at(iq));
        force_tmp = (lbm_solver.*ptr_func_cal_force_iq)(iq, node_coarse);
        ptr_f_collide_fine->at(iq) +=  dt_lbm * force_tmp / 4. - tau_collision_c2f_ *force_tmp * dt_lbm / 2.;
    }
}
/**
 * @brief function to convert post collision distribution functions on fine grid to coarse grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (fine grid).
 * @param[in] feq equilibrium distribution functions.
 * @param[in] lbm_solver reference of lbm solver.
 * @param[in] ptr_func_cal_force_iq pointer to function calculating lbm forcing term.
 * @param[in] node_fine  a node on fine grid.
 * @param[out] ptr_node_coarse pointer to a node on coarse grid.
 */
void LbmCollisionOptInterface::PostCollisionFine2CoarseForce(const DefReal dt_lbm, const std::vector<DefReal>& feq,
    const SolverLbmInterface& lbm_solver,
    DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
    const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) {
    DefReal force_tmp;
    DefSizet num_q = feq.size();
    ptr_node_coarse->rho_ = node_fine.rho_;
    ptr_node_coarse->velocity_ = node_fine.velocity_;
    bool bool_forces_coarse = (node_fine.force_.size() != 0);
    for (int iq = 0; iq < num_q; ++iq) {
        ptr_node_coarse->f_collide_[iq] = feq.at(iq) + tau_collision_f2c_ * (node_fine.f_collide_[iq] - feq.at(iq));
        force_tmp = (lbm_solver.*ptr_func_cal_force_iq)(iq, node_fine);
        ptr_node_coarse->f_collide_[iq] +=  dt_lbm * force_tmp - tau_collision_f2c_ * force_tmp * dt_lbm / 2.;
    }
}
/**
 * @brief function to convert post stream distribution functions on coarse grid to fine grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (coarse grid).
 * @param[in] lbm_solver reference of lbm solver.
 * @param[in] f_coarse post stream distribution functions of a coarse grid.
 * @param[out] ptr_f_fine pointer to post stream distribution functions of a fine grid.
 */
void LbmCollisionOptInterface::PostStreamCoarse2Fine(const DefReal dt_lbm, const std::vector<DefReal>& feq,
    const std::vector<DefReal>& f_coarse, std::vector<DefReal>* const ptr_f_fine) {
    DefSizet num_q = feq.size();
    ptr_f_fine->resize(num_q);
    for (int iq = 0; iq < num_q; ++iq) {
        ptr_f_fine->at(iq) = feq.at(iq) + tau_stream_c2f_ * (f_coarse[iq] - feq.at(iq));
    }
}
/**
 * @brief function to convert post stream distribution functions on fine grid to coarse grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (fine grid).
 * @param[in] feq equilibrium distribution functions.
 * @param[in] node_fine  a node on fine grid.
 * @param[out] ptr_node_coarse pointer to a node on coarse grid.
 */
void LbmCollisionOptInterface::PostStreamFine2Coarse(const DefReal dt_lbm, const std::vector<DefReal>& feq,
    const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) {
    ptr_node_coarse->rho_ = node_fine.rho_;
    ptr_node_coarse->velocity_ = node_fine.velocity_;
    DefSizet num_q = feq.size();
    for (int iq = 0; iq < num_q; ++iq) {
        ptr_node_coarse->f_[iq] = feq.at(iq) + tau_stream_f2c_ * (node_fine.f_[iq] - feq.at(iq));
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
