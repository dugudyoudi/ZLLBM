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
 * @brief function to instantiate collision operator based on input type.
 * @param[in] i_level input level.
 * @param[in] collision_operator_type predefine type of collision operator.
 */
void SolverLbmInterface::SetCollisionOperator(const DefInt i_level,
    const ELbmCollisionOperatorType collision_operator_type) {
    if (collision_operators_.find(i_level) != collision_operators_.end()) {
        collision_operators_.erase(i_level);
    }
    switch (collision_operator_type) {
        case ELbmCollisionOperatorType::kLbmSrt: {
            if (bool_forces_) {
                collision_operators_.insert({i_level, std::make_unique<LbmStrForceCollisionOpt>(i_level, *this)});
            } else {
                collision_operators_.insert({i_level, std::make_unique<LbmStrCollisionOpt>(i_level, *this)});
            }
        }
        break;
        case ELbmCollisionOperatorType::kLbmMrt: {
            if (bool_forces_) {
                collision_operators_.insert({i_level, std::make_unique<LbmMrtForceCollisionOpt>(i_level, *this)});
            } else {
                collision_operators_.insert({i_level, std::make_unique<LbmMrtCollisionOpt>(i_level, *this)});
            }
        }
        break;
        default:
        amrproject::LogManager::LogError("Collision operator at level " + std::to_string(i_level) + " is not defined");
            break;
    }
}
/**
 * @brief function to get collision operator at a given level.
 * @param[in] i_level input level.
 */
const LbmCollisionOptInterface& SolverLbmInterface::GetCollisionOperator(DefInt i_level) const {
    if (collision_operators_.find(i_level) != collision_operators_.end()) {
        return *collision_operators_.at(i_level).get();
    } else {
        amrproject::LogManager::LogError("Collision operator at level "
            + std::to_string(i_level) + " is not instantiated");
    }
    return *collision_operators_.at(0).get();
}
LbmCollisionOptInterface::LbmCollisionOptInterface(const DefInt i_level, const SolverLbmInterface& lbm_solver) {
    viscosity_lbm_ = lbm_solver.GetDefaultViscosity();
    dt_lbm_ = 1./ static_cast<DefReal>(TwoPowerN(i_level));
}
/**
 * @brief function to calculate relaxation time and matrix.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmCollisionOptInterface::CalRelaxationTime(const SolverLbmInterface& lbm_solver) {
    tau_ = viscosity_lbm_ * lbm_solver.kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] node reference to LBM node information
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmCollisionOptInterface::CalRelaxationTimeNode(
    const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) {
    tau_ = viscosity_lbm_ * lbm_solver.kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
}
/**
 * @brief function to calculate relaxation time ratio between coarse and fine grid.
 * @param[in] lbm_solver reference to LBM solver
 */
void LbmCollisionOptInterface::CalRelaxationTimeRatio(const SolverLbmInterface& lbm_solver) {
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
void LbmCollisionOptInterface::PostCollisionCoarse2Fine(const std::vector<DefReal>& feq,
    const std::vector<DefReal>& f_collide_coarse, std::vector<DefReal>* const ptr_f_collide_fine) const {
    DefInt num_q = static_cast<DefInt>(feq.size());
    ptr_f_collide_fine->resize(num_q);
    //  DefReal relax_tau_f = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm * 2. + 0.5;
    //  DefReal relax_tau_c = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm + 0.5;
    //  DefReal tau_ratio = 0.5 * (relax_tau_f - 1)/ (relax_tau_c - 1);

    for (DefInt iq = 0; iq < num_q; ++iq) {
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
void LbmCollisionOptInterface::PostCollisionFine2Coarse(const std::vector<DefReal>& feq,
    const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) const {
    ptr_node_coarse->rho_ = node_fine.rho_;
    ptr_node_coarse->velocity_ = node_fine.velocity_;
    DefInt num_q = static_cast<DefInt>(feq.size());
    ptr_node_coarse->f_collide_.resize(num_q);
    // DefReal relax_tau_f = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm + 0.5;
    // DefReal relax_tau_c = viscosity_lbm_ * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm / 2 + 0.5;
    //  DefReal tau_ratio = 2 * (relax_tau_c - 1)/ (relax_tau_f - 1);

    for (DefInt iq = 0; iq < num_q; ++iq) {
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
void LbmCollisionOptInterface::PostCollisionCoarse2FineForce(const std::vector<DefReal>& feq,
    const SolverLbmInterface& lbm_solver,
    DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
    const GridNodeLbm& node_coarse, std::vector<DefReal>* const ptr_f_collide_fine) const {
    DefReal force_tmp;
    DefInt num_q = static_cast<DefInt>(feq.size());
    ptr_f_collide_fine->resize(num_q);
    for (DefInt iq = 0; iq <num_q; ++iq) {
        ptr_f_collide_fine->at(iq) = feq.at(iq) + tau_collision_c2f_ * (node_coarse.f_collide_[iq] - feq.at(iq));
        force_tmp = (lbm_solver.*ptr_func_cal_force_iq)(iq, node_coarse);
        ptr_f_collide_fine->at(iq) +=  dt_lbm_ * force_tmp / 4. - tau_collision_c2f_ *force_tmp * dt_lbm_ / 2.;
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
void LbmCollisionOptInterface::PostCollisionFine2CoarseForce(const std::vector<DefReal>& feq,
    const SolverLbmInterface& lbm_solver,
    DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
    const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) const {
    DefReal force_tmp;
    DefInt num_q = static_cast<DefInt>(feq.size());
    ptr_node_coarse->rho_ = node_fine.rho_;
    ptr_node_coarse->velocity_ = node_fine.velocity_;
    bool bool_forces_coarse = (node_fine.force_.size() != 0);
    for (DefInt iq = 0; iq < num_q; ++iq) {
        ptr_node_coarse->f_collide_[iq] = feq.at(iq) + tau_collision_f2c_ * (node_fine.f_collide_[iq] - feq.at(iq));
        force_tmp = (lbm_solver.*ptr_func_cal_force_iq)(iq, node_fine);
        ptr_node_coarse->f_collide_[iq] +=  dt_lbm_ * force_tmp - tau_collision_f2c_ * force_tmp * dt_lbm_ / 2.;
    }
}
/**
 * @brief function to convert post stream distribution functions on coarse grid to fine grid.
 * @param[in] dt_lbm  LBM time step at current refinement level (coarse grid).
 * @param[in] lbm_solver reference of lbm solver.
 * @param[in] f_coarse post stream distribution functions of a coarse grid.
 * @param[out] ptr_f_fine pointer to post stream distribution functions of a fine grid.
 */
void LbmCollisionOptInterface::PostStreamCoarse2Fine(const std::vector<DefReal>& feq,
    const std::vector<DefReal>& f_coarse, std::vector<DefReal>* const ptr_f_fine) const {
    DefInt num_q = static_cast<DefInt>(feq.size());
    ptr_f_fine->resize(num_q);
    for (DefInt iq = 0; iq < num_q; ++iq) {
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
void LbmCollisionOptInterface::PostStreamFine2Coarse(const std::vector<DefReal>& feq,
    const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) const {
    ptr_node_coarse->rho_ = node_fine.rho_;
    ptr_node_coarse->velocity_ = node_fine.velocity_;
    DefInt num_q = static_cast<DefInt>(feq.size());
    for (DefInt iq = 0; iq < num_q; ++iq) {
        ptr_node_coarse->f_[iq] = feq.at(iq) + tau_stream_f2c_ * (node_fine.f_[iq] - feq.at(iq));
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
