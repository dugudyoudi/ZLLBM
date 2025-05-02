//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_solver_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM solver interface.
* @date  2023-9-30
*/
#include "./lbm_interface.h"
#include "grid/grid_manager.h"
#include "grid/grid_info_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to initialize LBM solver.
 */
void SolverLbmInterface::SolverInitial() {
    if (k0LbmViscosity_ < kEps) {
        amrproject::LogManager::LogWarning("LBM viscosity is less than predefined kEps("
        + std::to_string(kEps) + ") in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (k0SolverDims_ == 2) {
        SetDefault2DFunctions();
    } else {
        SetDefault3DFunctions();
    }
    if (this->bool_forces_) {
        this->func_macro_ = [this](const DefReal dt_lbm, const GridNodeLbm& lbm_node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_v) {
            std::vector<DefReal> force(this->GetAllForcesForANode(lbm_node));
            this->func_macro_with_force_(dt_lbm, lbm_node.f_, force, ptr_rho, ptr_v);
        };
    } else {
        this->func_macro_ = [this](const DefReal dt_lbm, const GridNodeLbm& lbm_node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_v) {
            this->func_macro_without_force_(lbm_node.f_, ptr_rho, ptr_v);
        };
    }
    SetFunctionForNodeCollision();
    InitialModelDependencies();
}
/**
* @brief function to read abd set input parameters for LBM solver.
* @param[in, out] ptr_input_parser pointer to class for input parsing. 
*/
void SolverLbmInterface::ReadAndSetupSolverParameters(amrproject::InputParser* const ptr_input_parser) {
    std::string key_for_this_func =  "solver.lbm";
    auto& input_map = ptr_input_parser->GetNestedMapInput();

    if (input_map.find(key_for_this_func) != input_map.end()) {
        if (input_map.at(key_for_this_func).find(name_) == input_map.at(key_for_this_func).end()) {
            amrproject::LogManager::LogError("key (" + name_
                + ") does not match any solver name in input file");
        } else {
            std::string values_str;
            auto ptr_solver_parameters = &input_map.at(key_for_this_func).at(solver_type_);
            if (ptr_input_parser->GetValue<DefReal>("lbm.viscosity", ptr_solver_parameters, &k0LbmViscosity_)) {
                if (ptr_input_parser->print_values_when_read_) {
                    amrproject::LogManager::LogInfo("Read and set LBM viscosity: " + std::to_string(k0LbmViscosity_));
                }
            }
            if (ptr_input_parser->GetValue<DefReal>("lbm.rho_ini", ptr_solver_parameters, &k0Rho_)) {
                if (ptr_input_parser->print_values_when_read_) {
                    amrproject::LogManager::LogInfo("Read and set initial density: " + std::to_string(k0Rho_));
                }
            }
            if (ptr_input_parser->GetValue<DefReal>("lbm.velocity_ini", ptr_solver_parameters, &k0Velocity_)) {
                if (ptr_input_parser->print_values_when_read_) {
                    values_str = ptr_input_parser->ValuesToOutputStr<DefReal>(k0Velocity_);
                    amrproject::LogManager::LogInfo("Read and set initial velocity: " + values_str);
                }
            }
            std::vector<DefReal> force_in;
            if (ptr_input_parser->GetValue<DefReal>("lbm.force_ini", ptr_solver_parameters, &force_in)) {
                SetDefaultForce(force_in);
                if (ptr_input_parser->print_values_when_read_) {
                    values_str = ptr_input_parser->ValuesToOutputStr<DefReal>(force_in);
                    amrproject::LogManager::LogInfo("Read and set initial force: " + values_str);
                }
            }
            if (ptr_input_parser->GetValue<DefReal>("lbm.const_force", ptr_solver_parameters, &force_in)) {
                SetConstantForce(force_in);
                if (ptr_input_parser->print_values_when_read_) {
                    values_str = ptr_input_parser->ValuesToOutputStr<DefReal>(force_in);
                    amrproject::LogManager::LogInfo("Read and set constant force: " + values_str);
                }
            }
            DefInt num_force = 0;
            if (ptr_input_parser->GetValue<DefInt>("lbm.force_size", ptr_solver_parameters, &num_force)) {
                SetNumForces(num_force);
                if (ptr_input_parser->print_values_when_read_) {
                    amrproject::LogManager::LogInfo("Read and set number of forces : " +  std::to_string(num_force));
                }
            }
            std::string str_tmp;
            if (ptr_input_parser->GetValue<std::string>("lbm.les_model", ptr_solver_parameters, &str_tmp)) {
                SetLesModel(str_tmp);
            }
            std::vector<std::string> str_vec_tmp;
            if (ptr_input_parser->GetValue<std::string>("lbm.collision_models", ptr_solver_parameters, &str_vec_tmp)) {
                ReadCollisionOperator(str_vec_tmp);
                if (ptr_input_parser->print_values_when_read_) {
                    values_str.clear();
                    for (const auto& iter_str : str_vec_tmp) {
                        values_str += " " + iter_str;
                    }
                    amrproject::LogManager::LogInfo("Read and set collision model: " + values_str);
                }
            }
        }
    } else {
        amrproject::LogManager::LogError("key (solver.lbm) is not found for the LBM solver in input file.");
    }
}
/**
 * @brief function to set default force for initializing node.
 */
void SolverLbmInterface::SetDefaultForce(const std::vector<DefReal>& force_in) {
    bool_forces_ = true;
    k0Force_ = force_in;
}
/**
 * @brief function to set constant forces.
 */
void SolverLbmInterface::SetConstantForce(const std::vector<DefReal>& force_in) {
    bool_forces_ = true;
    k0ConstForce_ = force_in;
}
std::vector<DefReal> SolverLbmInterface::GetAllForcesForANode(
    const GridNodeLbm& node) const {
    std::vector<DefReal> force(k0ConstForce_);
    force.resize(k0SolverDims_);
    for (DefInt i = 0; i< static_cast<DefInt>(node.force_.size())/k0SolverDims_; ++i) {
        force.at(kXIndex) += node.force_.at(i*k0SolverDims_);
        force.at(kYIndex) += node.force_.at(i*k0SolverDims_ + 1);
        if (k0SolverDims_ == 3) {
            force.at(kZIndex) += node.force_.at(i*k0SolverDims_ + 2);
        }
    }
    return force;
}
/**
 * @brief function to set pointers to the default 2D member functions.
 */
void SolverLbmInterface::SetDefault2DFunctions() {
    if (k0BoolCompressible_) {
        this->func_macro_without_force_ = [this](const std::vector<DefReal>& f,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro2DCompressible(f, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const std::vector<DefReal>& f,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce2DCompressible(dt_lbm, f, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
            this->CalFeq2DCompressible(rho, velocity, ptr_feq);
        };
    } else {
        this->func_macro_without_force_ = [this](const std::vector<DefReal>& f,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro2DIncompressible(f, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const std::vector<DefReal>& f,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce2DIncompressible(dt_lbm, f, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
            this->CalFeq2DIncompressible(rho, velocity, ptr_feq);
        };
    }
    this->ptr_func_cal_force_iq_ = &SolverLbmInterface::CalForcingTerming2D;
}
/**
 * @brief function to set pointers to the default 3D member functions.
 */
void SolverLbmInterface::SetDefault3DFunctions() {
    if (k0BoolCompressible_) {
        this->func_macro_without_force_ = [this](const std::vector<DefReal>& f,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro3DCompressible(f, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const std::vector<DefReal>& f,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce3DCompressible(dt_lbm, f, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
                this->CalFeq3DCompressible(rho, velocity, ptr_feq);
        };
    } else {
        this->func_macro_without_force_ = [this](const std::vector<DefReal>& f,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro3DIncompressible(f, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const std::vector<DefReal>& f,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce3DIncompressible(dt_lbm, f, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
                this->CalFeq3DIncompressible(rho, velocity, ptr_feq);
        };
    }
    this->ptr_func_cal_force_iq_ = &SolverLbmInterface::CalForcingTerming3D;
}
/**
 * @brief function to resize model related vectors for better performance.
 */
void SolverLbmInterface::ResizeModelRelatedVectors() {
    this->k0Velocity_.resize(k0SolverDims_);
    this->k0Velocity_.shrink_to_fit();
}
/**
 * @brief function for marching LBM time step at grid of a given refinement level.
 * @param[in] time_step_level time step at current grid refinement level in one background step.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[out] ptr_grid_info pointer to class storing LBM grid information.
 */
void SolverLbmInterface::InformationFromGridOfDifferentLevel(
    const DefInt time_step_level, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    amrproject::GridInfoInterface* const ptr_grid_info) {
    const DefInt i_level = ptr_grid_info->GetGridLevel();
    if (i_level > 0 && (time_step_level%2 == 0)) {
        GridInfoLbmInteface* ptr_lbm_grid_info_coarse = dynamic_cast<GridInfoLbmInteface*>(
            ptr_grid_manager_->vec_ptr_grid_info_.at(i_level - 1).get());

        ptr_grid_info->TransferInfoFromCoarseGrid(*ptr_grid_manager_->GetPtrToSFBitsetAux(),
            amrproject::NodeBitStatus::kNodeStatusCoarse2FineGhost_, *ptr_lbm_grid_info_coarse);

        ptr_grid_info->TransferInfoToCoarseGrid(*ptr_grid_manager_->GetPtrToSFBitsetAux(),
            amrproject::NodeBitStatus::kNodeStatusCoarse2FineGhost_, ptr_lbm_grid_info_coarse);
    }
}
/**
 * @brief function for marching LBM time step at grid of a given refinement level.
 * @param[in] time_scheme enum class to identify time stepping scheme used in computation.
 * @param[in] time_step_current time step at current grid refinement level in one background step.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void SolverLbmInterface::RunSolverOnGivenGrid(const amrproject::ETimeSteppingScheme time_scheme,
    const DefInt time_step_level, const DefReal time_step_current,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux, amrproject::GridInfoInterface* const ptr_grid_info) {
    GridInfoLbmInteface* ptr_lbm_grid_nodes_info = dynamic_cast<GridInfoLbmInteface*>(ptr_grid_info);
    if (ptr_lbm_grid_nodes_info->GetPtrToLbmGrid() != nullptr) {
        DefMap<std::unique_ptr<GridNodeLbm>>& grid_nodes = *ptr_lbm_grid_nodes_info->GetPtrToLbmGrid();

        Collision(ptr_lbm_grid_nodes_info->GetNodeFlagNotCollision(), ptr_lbm_grid_nodes_info);

        Stream(ptr_lbm_grid_nodes_info->GetNodeFlagNotStream(), sfbitset_aux, &grid_nodes);

        // this function is used to reset flags that have change in InformationFromGridOfDifferentLevel
        // since some nodes should not be calculated after transferring information between different levels
        // otherwise the calculated ones will overwrite correct ones
        ptr_lbm_grid_nodes_info->InitialNotComputeNodeFlag();
    }
}
/**
 * @brief function to set collision step according to settings.
 */
void SolverLbmInterface::SetFunctionForNodeCollision() {
    if (bool_forces_) {
        this->func_collision_node_ = [this](const DefReal dt_lbm,
            LbmCollisionOptInterface* const ptr_collision_opt, GridNodeLbm* const ptr_node) {
            std::vector<DefReal> force = GetAllForcesForANode(*ptr_node);
            func_macro_with_force_(dt_lbm, ptr_node->f_, force, &ptr_node->rho_, &ptr_node->velocity_);
            const DefInt& num_q = k0NumQ_;
            std::vector<DefReal> feq(num_q, 0.), forcing_term(num_q, 0.);
            func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
            (this->*(this->ptr_func_cal_force_iq_))(*ptr_node, force, &forcing_term);
            if (bool_les_model_) {
                DefReal tau_sgs = ptr_les_model_->CalSgsRelaxationTimeWithForce(
                    ptr_collision_opt->GetDtLbm(), ptr_collision_opt->GetRelaxationTime(),
                    feq, force, *ptr_node, *this);
                ptr_collision_opt->SetEffectiveRelaxationTimeForForcing(tau_sgs, *this);
            }
            ptr_collision_opt->CollisionOperator(num_q, feq, forcing_term, ptr_node);
        };
    } else {
        this->func_collision_node_ = [this](const DefReal dt_lbm,
            LbmCollisionOptInterface* const ptr_collision_opt, GridNodeLbm* const ptr_node) {
            func_macro_without_force_(ptr_node->f_, &ptr_node->rho_, &ptr_node->velocity_);
            const DefInt& num_q = k0NumQ_;
            std::vector<DefReal> feq(num_q, 0.);
            func_cal_feq_(ptr_node->rho_, ptr_node->velocity_, &feq);
            if (bool_les_model_) {
                DefReal tau_sgs = ptr_les_model_->CalSgsRelaxationTimeWithoutForce(
                    ptr_collision_opt->GetDtLbm(), ptr_collision_opt->GetRelaxationTime(),
                    feq, *ptr_node, *this);
                ptr_collision_opt->SetEffectiveRelaxationTime(tau_sgs, *this);
            }
            ptr_collision_opt->CollisionOperator(num_q, feq, ptr_node);
        };
    }
}
/**
 * @brief function to perform collision step in the LBM simulation.
 * @param[in] flag_not_compute flag indicating whether to compute or not.
 * @param[out] ptr_lbm_grid_nodes_info pointer to class storing LBM grid information.
 */
void SolverLbmInterface::Collision(
    const DefInt flag_not_compute, GridInfoLbmInteface* const ptr_lbm_grid_nodes_info) const {
    // choose function to compute macroscopic variables based on if the forces are considered
    std::function<void(const DefReal, const std::vector<DefReal>&, const std::vector<DefReal>&,
        DefReal* const, std::vector<DefReal>* const)> func_macro;
    LbmCollisionOptInterface& collision_opt = GetCollisionOperator(ptr_lbm_grid_nodes_info->GetGridLevel());
    const DefReal dt_lbm = collision_opt.GetDtLbm();
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_nodes = *ptr_lbm_grid_nodes_info->GetPtrToLbmGrid();

    if (ptr_lbm_grid_nodes_info->GetPtrToLbmGrid() != nullptr) {
        for (auto& iter_node : grid_nodes) {
            if (iter_node.second->flag_status_ & flag_not_compute) {
            } else {
                this->func_collision_node_(dt_lbm, &collision_opt, iter_node.second.get());
            }
        }
    }
}
/**
 * @brief function to perform steam step in the LBM simulation.
 * @param[in] flag_not_compute flag indicating whether to compute or not.
 * @param[out] ptr_lbm_grid_nodes_info pointer to class storing LBM grid information.
 */
void SolverLbmInterface::Stream(const DefInt flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    if (ptr_map_grid_nodes != nullptr) {
        for (auto& iter_node : *ptr_map_grid_nodes) {
            if (!(iter_node.second->flag_status_ & flag_not_compute)) {
                StreamInForAGivenNode(iter_node.first, sfbitset_aux, ptr_map_grid_nodes);
            }
        }
    }
}
/**
 * @brief function to set initial distribution functions based on the reference macroscopic variables.
 * @param[out] ptr_f pointer to a vector storing distribution function after streaming.
 * @param[out] ptr_f_collide pointer to a vector storing distribution function after collision.
 */
void SolverLbmInterface::SetInitialDisFuncBasedOnReferenceMacros(
    std::vector<DefReal>* const ptr_f, std::vector<DefReal>* const ptr_f_collide) const {;
    if (static_cast<DefInt>(k0Velocity_.size()) == k0SolverDims_) {
        CalInitialDisFuncAsFeq(k0Rho_, k0Velocity_, ptr_f, ptr_f_collide);
        ptr_f->resize(k0NumQ_);
        ptr_f_collide->resize(k0NumQ_);
    } else {
        amrproject::LogManager::LogError("size of initial velocity (k0Velocity_) is not equal to"
            " solver dimension (k0SolverDims_) in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
 * @brief function to calculate the initial distribution function as equilibrium distribution functions.
 * @param[in] rho density.
 * @param[in] velocity velocities.
 * @param[out] ptr_f  pointer to a vector storing distribution function after streaming.
 * @param[out] ptr_f_collide pointer to a vector storing distribution function after collision.
 */
void SolverLbmInterface::CalInitialDisFuncAsFeq(const DefReal rho, const std::vector<DefReal>& velocity,
    std::vector<DefReal>* const ptr_f, std::vector<DefReal>* const ptr_f_collide) const {
    func_cal_feq_(rho, velocity, ptr_f);
    *ptr_f_collide = *ptr_f;
}
/**
 * @brief function to calculate the equilibrium distribution functions for a 2D grid node.
 * @param[in] rho  density at equilibrium status.
 * @param[in] velocity  velocity at equilibrium status.
 * @param[out] ptr_feq pointer to the vector to store the calculated equilibrium distribution functions.
 */
void SolverLbmInterface::CalFeq2DCompressible(const DefReal rho,
    const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const {
    ptr_feq->resize(this->k0NumQ_);
    DefReal c_uv = 0.;
    DefReal uv_sq = Square(velocity.at(kXIndex)) + Square(velocity.at(kYIndex));
    for (DefInt iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = velocity.at(kXIndex) * k0Cx_.at(iq) + velocity.at(kYIndex) * k0Cy_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) * rho * (1. + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
/**
 * @brief function to calculate the equilibrium distribution functions for a 2D grid node (incompressible).
 * @param[in] rho  density at equilibrium status.
 * @param[in] velocity  velocity at equilibrium status.
 * @param[out] ptr_feq pointer to the vector to store the calculated equilibrium distribution functions.
 */
void SolverLbmInterface::CalFeq2DIncompressible(const DefReal rho,
    const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const {
    ptr_feq->resize(this->k0NumQ_);
    DefReal c_uv = 0.;
    DefReal uv_sq = Square(velocity.at(kXIndex)) + Square(velocity.at(kYIndex));
    for (DefInt iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = velocity.at(kXIndex) * k0Cx_.at(iq) + velocity.at(kYIndex) * k0Cy_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) * (rho + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
/**
 * @brief function to calculate the equilibrium distribution functions for a 32D grid node.
 * @param[in] rho  density at equilibrium status.
 * @param[in] velocity  velocity at equilibrium status.
 * @param[out] ptr_feq pointer to the vector to store the calculated equilibrium distribution functions.
 */
void SolverLbmInterface::CalFeq3DCompressible(const DefReal rho,
    const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const {
    ptr_feq->resize(this->k0NumQ_);
    DefReal c_uv = 0.;
    DefReal uv_sq = Square(velocity.at(kXIndex)) + Square(velocity.at(kYIndex))
        + Square(velocity.at(kZIndex));
    for (DefInt iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = velocity.at(kXIndex) * k0Cx_.at(iq) + velocity.at(kYIndex) * k0Cy_.at(iq)
            + velocity.at(kZIndex) * k0Cz_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) * rho * (1. + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
/**
 * @brief function to calculate the equilibrium distribution functions for a 2D grid node (incompressible).
 * @param[in] rho  density at equilibrium status.
 * @param[in] velocity  velocity at equilibrium status.
 * @param[out] ptr_feq pointer to the vector to store the calculated equilibrium distribution functions.
 */
void SolverLbmInterface::CalFeq3DIncompressible(const DefReal rho,
    const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const {
    ptr_feq->resize(this->k0NumQ_);
    DefReal c_uv = 0.;
    DefReal uv_sq = Square(velocity.at(kXIndex)) + Square(velocity.at(kYIndex))
        + Square(velocity.at(kZIndex));
    for (DefInt iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = velocity.at(kXIndex) * k0Cx_.at(iq) + velocity.at(kYIndex) * k0Cy_.at(iq)
            + velocity.at(kZIndex) * k0Cz_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) *  (rho + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
/**
 * @brief function to calculate the force term for a given lattice direction in 2D.
 * @param[in] iq the ith lattice direction. 
 * @param[in] node grid node containing LBM related information.
 * @param[in] force body force acting on the node.
 * @param[out] ptr_forcing_term pointer to forcing term disribution functions.
 */
void SolverLbmInterface::CalForcingTerming2D(const GridNodeLbm& node, const std::vector<DefReal>& force,
    std::vector<DefReal>* const ptr_forcing_term) const {
    ptr_forcing_term->resize(k0NumQ_);
    for (DefInt iq = 0; iq < k0NumQ_; ++iq) {
        ptr_forcing_term->at(iq) = k0Weights_[iq] * (3.*((k0Cx_[iq] - node.velocity_[kXIndex]) * force[kXIndex]
            + (k0Cy_[iq] - node.velocity_[kYIndex]) * force[kYIndex])
            + 9. * (node.velocity_.at(kXIndex) * k0Cx_.at(iq) + node.velocity_.at(kYIndex) * k0Cy_.at(iq))
            * (k0Cx_[iq] * force[kXIndex] + k0Cy_[iq] * force[kYIndex]));
    }
}
/**
 * @brief function to calculate the force term for a given lattice direction in 3D.
 * @param[in] node grid node containing LBM related information.
 * @param[in] force body force acting on the node.
 * @param[out] ptr_forcing_term pointer to forcing term disribution functions.
 */
void SolverLbmInterface::CalForcingTerming3D(const GridNodeLbm& node, const std::vector<DefReal>& force,
    std::vector<DefReal>* const ptr_forcing_term) const {
    ptr_forcing_term->resize(k0NumQ_);
    for (DefInt iq = 0; iq < k0NumQ_; ++iq) {
        ptr_forcing_term->at(iq) = k0Weights_[iq] * (3.*((k0Cx_[iq] - node.velocity_[kXIndex]) * force[kXIndex]
            + (k0Cy_[iq] - node.velocity_[kYIndex]) * force[kYIndex]
            + (k0Cz_[iq] - node.velocity_[kZIndex]) * force[kZIndex])
            + 9. * (node.velocity_.at(kXIndex) * k0Cx_.at(iq) + node.velocity_.at(kYIndex) * k0Cy_.at(iq)
            + node.velocity_.at(kZIndex) * k0Cz_.at(iq))
            * (k0Cx_[iq] * force[kXIndex] + k0Cy_[iq] * force[kYIndex] + k0Cz_[iq] * force[kZIndex]));
    }
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] f distribution function.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro2DCompressible(const std::vector<DefReal>& f,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacro2DIncompressible(f, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] f distribution function.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro2DIncompressible(const std::vector<DefReal>& f,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    rho = f.at(0);
    (*ptr_velocity) = {0, 0};
    for (DefInt iq = 1; iq < k0NumQ_; ++iq) {
        (*ptr_rho) += f.at(iq);
        ptr_velocity->at(kXIndex) += k0Cx_[iq] * f.at(iq);
        ptr_velocity->at(kYIndex) += k0Cy_[iq] * f.at(iq);
    }
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] f distribution function.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro3DCompressible(const std::vector<DefReal>& f,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacro3DIncompressible(f, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
    ptr_velocity->at(kZIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] f distribution function.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro3DIncompressible(const std::vector<DefReal>& f,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    rho = f.at(0);
    (*ptr_velocity) = {0, 0, 0};
    for (DefInt iq = 1; iq < k0NumQ_; ++iq) {
        (*ptr_rho) += f.at(iq);
        ptr_velocity->at(kXIndex) += k0Cx_[iq] * f.at(iq);
        ptr_velocity->at(kYIndex) += k0Cy_[iq] * f.at(iq);
        ptr_velocity->at(kZIndex) += k0Cz_[iq] * f.at(iq);
    }
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] f distribution function.
 * @param[in] force body force acting on the node.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce2DCompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacroForce2DIncompressible(dt_lbm, f, force, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] f distribution function.
 * @param[in] force body force acting on the node.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce2DIncompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    CalMacro2DIncompressible(f, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex) += 0.5 * force.at(kXIndex) * dt_lbm;
    ptr_velocity->at(kYIndex) += 0.5 * force.at(kYIndex) * dt_lbm;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] f distribution function.
 * @param[in] force body force acting on the node.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce3DCompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacroForce3DIncompressible(dt_lbm, f, force, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
    ptr_velocity->at(kZIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] f distribution function.
 * @param[in] force body force acting on the node.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce3DIncompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    CalMacro3DIncompressible(f, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex) += 0.5 * force.at(kXIndex) * dt_lbm;
    ptr_velocity->at(kYIndex) += 0.5 * force.at(kYIndex) * dt_lbm;
    ptr_velocity->at(kZIndex) += 0.5 * force.at(kZIndex) * dt_lbm;
}
}  // end namespace lbmproject
}  // end namespace rootproject
