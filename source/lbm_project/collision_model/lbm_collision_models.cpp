//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_collision_models.cpp
* @author Zhengliang Liu
* @brief functions used for managing general collision models for LBM.
* @date  2023-10-30
*/
#include "./lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to push back type of collision operator based on input key word to container member for recording.
 * @param[in] operator_type identifier of a collision operator.
 */
void SolverLbmInterface::ReadCollisionOperator(const std::string& operator_type) {
    if (operator_type == "srt") {
        collision_type_.push_back(ELbmCollisionOperatorType::kLbmSrt);
    } else if (operator_type == "mrt") {
        collision_type_.push_back(ELbmCollisionOperatorType::kLbmMrt);
    } else {
        amrproject::LogManager::LogError("Collision model " + operator_type + "  is not defined, "
            "please choosing from: srt, and mrt");
    }
}
/**
 * @brief function to assign types of collision operator to container member for recording.
 * @param[in] operator_types identifiers of collision operators.
 */
void SolverLbmInterface::ReadCollisionOperator(const std::vector<std::string>& operator_types) {
    if (!operator_types.empty()) {
        collision_type_.clear();
        for (const auto& iter : operator_types) {
            ReadCollisionOperator(iter);
        }
    }
}
/**
 * @brief function to instantiate collision operator based on input type.
 * @param[in] i_level input level.
 * @param[in] collision_operator_type predefine type of collision operator.
 */
void SolverLbmInterface::InstantiateCollisionOperator(const DefInt i_level,
    const ELbmCollisionOperatorType collision_operator_type) {
    if (collision_operators_.find(i_level) != collision_operators_.end()) {
        collision_operators_.erase(i_level);
    }
    switch (collision_operator_type) {
        case ELbmCollisionOperatorType::kLbmSrt:
            collision_operators_.insert({i_level, std::make_unique<LbmStrCollisionOpt>(i_level, *this)});
            break;
        case ELbmCollisionOperatorType::kLbmMrt:
            collision_operators_.insert({i_level, std::make_unique<LbmMrtCollisionOpt>(i_level, *this)});
            break;
        default:
            amrproject::LogManager::LogError(
                "Collision operator at level " + std::to_string(i_level) + " is not defined");
            break;
    }
}
/**
 * @brief function to instantiate collision operator based on input type.
 * @param[in] i_level input level.
 */
void SolverLbmInterface::InstantiateCollisionOperator(const DefInt i_level) {
    if (collision_type_.empty()) {  // default is SRT collision model
        InstantiateCollisionOperator(i_level, ELbmCollisionOperatorType::kLbmSrt);
    } else {
        if (static_cast<DefInt>(collision_type_.size()) > i_level) {
            InstantiateCollisionOperator(i_level, collision_type_.at(i_level));
        } else {
            InstantiateCollisionOperator(i_level, collision_type_.back());
        }
    }
}
/**
 * @brief function to set collision operator based on input type.
 * @param[in] i_level input level.
 * @param[in] collision_operator_type predefine type of collision operator.
 */
void SolverLbmInterface::SetCollisionOperator(const DefInt i_level,
    const ELbmCollisionOperatorType collision_operator_type) {
    if (static_cast<DefInt>(collision_type_.size()) > i_level) {
        collision_type_.at(i_level) = collision_operator_type;
    } else {
        collision_type_.resize(i_level + 1);
        collision_type_.at(i_level) = collision_operator_type;
    }
}
/**
 * @brief function to get collision operator at a given level.
 * @param[in] i_level input level.
 */
LbmCollisionOptInterface& SolverLbmInterface::GetCollisionOperator(DefInt i_level) const {
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
 * @param[in] vis_lbm viscosity scaling for LBM solver.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmCollisionOptInterface::CalRelaxationTime(const DefReal vis_lbm,
    const SolverLbmInterface& lbm_solver) {
    tau0_ = vis_lbm * lbm_solver.kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
    tau_eff_ = tau0_;
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] tau_eff effective relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmCollisionOptInterface::SetEffectiveRelaxationTime(
    const DefReal tau_eff, const SolverLbmInterface& lbm_solver) {
    tau_eff_ = tau_eff;
}
/**
 * @brief function to calculate relaxation time and matrix based on node information.
 * @param[in] tau_eff effective relaxation time.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmCollisionOptInterface::SetEffectiveRelaxationTimeForForcing(
    const DefReal tau_eff, const SolverLbmInterface& lbm_solver) {
    tau_eff_ = tau_eff;
}
/**
 * @brief function to calculate relaxation time ratio between coarse and fine grid.
* @param[in] vis_lbm viscosity scaling for LBM solver.
 * @param[in] lbm_solver reference to LBM solver.
 */
void LbmCollisionOptInterface::CalRelaxationTimeRatio(const DefReal vis_lbm, const SolverLbmInterface& lbm_solver) {
    DefReal relax_tau_f = vis_lbm * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ * 2. + 0.5;
    DefReal relax_tau_c = vis_lbm * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
    tau_collision_c2f_ = 0.5 * (relax_tau_f - 1)/(relax_tau_c - 1);
    tau_stream_c2f_ = 0.5 * relax_tau_f/relax_tau_c;
    relax_tau_f = vis_lbm * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ + 0.5;
    relax_tau_c = vis_lbm * SolverLbmInterface::kCs_Sq_Reciprocal_ / dt_lbm_ / 2 + 0.5;
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
