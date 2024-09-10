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
    InitialModelDependencies();
}
std::vector<DefReal> SolverLbmInterface::GetAllForcesForANode(
    const DefInt dims, const GridNodeLbm& node) const {
    return node.force_;
}
/**
 * @brief function to set pointers to the default 2D member functions.
 */
void SolverLbmInterface::SetDefault2DFunctions() {
    if (k0BoolCompressible_) {
        this->func_macro_without_force_ = [this](const GridNodeLbm& node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro2DCompressible(node, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const GridNodeLbm& node,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce2DCompressible(dt_lbm, node, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
            this->CalFeq2DCompressible(rho, velocity, ptr_feq);
        };
    } else {
        this->func_macro_without_force_ = [this](const GridNodeLbm& node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro2DIncompressible(node, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const GridNodeLbm& node,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce2DIncompressible(dt_lbm, node, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
            this->CalFeq2DIncompressible(rho, velocity, ptr_feq);
        };
    }
    this->ptr_func_cal_force_iq_ = &SolverLbmInterface::CalForceIq2D;
}
/**
 * @brief function to set pointers to the default 3D member functions.
 */
void SolverLbmInterface::SetDefault3DFunctions() {
    if (k0BoolCompressible_) {
        this->func_macro_without_force_ = [this](const GridNodeLbm& node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro3DCompressible(node, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const GridNodeLbm& node,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce3DCompressible(dt_lbm, node, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
                this->CalFeq3DCompressible(rho, velocity, ptr_feq);
        };
    } else {
        this->func_macro_without_force_ = [this](const GridNodeLbm& node,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacro3DIncompressible(node, ptr_rho, ptr_velocity);
        };
        this->func_macro_with_force_ = [this](const DefReal dt_lbm, const GridNodeLbm& node,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroForce3DIncompressible(dt_lbm, node, force, ptr_rho, ptr_velocity);
        };
        this->func_cal_feq_ = [this](const DefReal rho, const std::vector<DefReal>& velocity,
            std::vector<DefReal>* const ptr_feq) {
                this->CalFeq3DIncompressible(rho, velocity, ptr_feq);
        };
    }
    this->ptr_func_cal_force_iq_ = &SolverLbmInterface::CalForceIq3D;
}
/**
 * @brief function to resize model related vectors for better performance.
 */
void SolverLbmInterface::ResizeModelRelatedVectors() {
    this->k0Velocity_.resize(k0SolverDims_);
    this->k0Velocity_.shrink_to_fit();
    this->k0Force_.resize(k0SolverDims_);
    this->k0Force_.shrink_to_fit();
    k0Cx_.resize(k0NumQ_);
    k0Cx_.shrink_to_fit();
    k0Cy_.resize(k0NumQ_);
    k0Cy_.shrink_to_fit();
    k0Weights_.resize(k0NumQ_);
    k0Weights_.shrink_to_fit();
    k0QIndicesNeg_.resize(k0SolverDims_);
    k0QIndicesNeg_.shrink_to_fit();
    for (auto i = 0; i < k0SolverDims_; ++i) {
        k0QIndicesNeg_.at(i).resize(k0NumQInOneDirection_);
        k0QIndicesNeg_.at(i).shrink_to_fit();
    }
    k0QIndicesPos_.resize(k0SolverDims_);
    k0QIndicesPos_.shrink_to_fit();
     for (auto i = 0; i < k0SolverDims_; ++i) {
        k0QIndicesPos_.at(i).resize(k0NumQInOneDirection_);
        k0QIndicesPos_.at(i).shrink_to_fit();
    }
}
/**
 * @brief function to call domain boundary conditions.
 */
void SolverLbmInterface::CallDomainBoundaryCondition(const amrproject::ETimeSteppingScheme time_scheme,
    const DefInt time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    amrproject::GridInfoInterface* const ptr_grid_info) {
    GridInfoLbmInteface* ptr_lbm_grid_nodes_info = dynamic_cast<GridInfoLbmInteface*>(ptr_grid_info);
    ptr_lbm_grid_nodes_info->ComputeDomainBoundaryCondition();
}
/**
 * @brief function for marching LBM time step at grid of a given refinement level.
 * @param[in] timing indicator for timing in one time step.
 * @param[in] time_scheme enum class to identify time stepping scheme used in computation.
 * @param[in] time_step_current time step at current grid refinement level in one background step.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[out] ptr_grid_info pointer to class storing LBM grid information.
 * @return 0 run successfully, otherwise some error occurred in this function.
 */
int SolverLbmInterface::InformationFromGridOfDifferentLevel(
    const amrproject::ETimingInOneStep timing, const amrproject::ETimeSteppingScheme time_scheme,
    const DefInt time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    amrproject::GridInfoInterface* const ptr_grid_info) {
    DefInt i_level = ptr_grid_info->i_level_;
    if (i_level > 0 && (time_step_current%2 == 0)) {
        GridInfoLbmInteface* ptr_lbm_grid_info = dynamic_cast<GridInfoLbmInteface*>(ptr_grid_info);
        GridInfoLbmInteface* ptr_lbm_grid_info_coarse = dynamic_cast<GridInfoLbmInteface*>(
            ptr_grid_manager_->vec_ptr_grid_info_.at(i_level - 1).get());

        ptr_grid_info->TransferInfoFromCoarseGrid(*ptr_grid_manager_->GetSFBitsetAuxPtr(),
            amrproject::NodeBitStatus::kNodeStatusCoarse2FineGhost_, *ptr_lbm_grid_info_coarse);

        ptr_grid_info->TransferInfoToCoarseGrid(*ptr_grid_manager_->GetSFBitsetAuxPtr(),
            amrproject::NodeBitStatus::kNodeStatusCoarse2FineGhost_, ptr_lbm_grid_info_coarse);
    }
    return 0;
}
/**
 * @brief function for marching LBM time step at grid of a given refinement level.
 * @param[in] time_scheme enum class to identify time stepping scheme used in computation.
 * @param[in] time_step_current time step at current grid refinement level in one background step.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[out] ptr_grid_info pointer to class storing grid information.
 */
void SolverLbmInterface::RunSolverForNodesOnNormalGrid(const amrproject::ETimeSteppingScheme time_scheme,
    const DefInt time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    amrproject::GridInfoInterface* const ptr_grid_info) {
    DefInt i_level = ptr_grid_info->i_level_;
    GridInfoLbmInteface* ptr_lbm_grid_nodes_info = dynamic_cast<GridInfoLbmInteface*>(ptr_grid_info);
    if (ptr_lbm_grid_nodes_info->GetPointerToLbmGrid() != nullptr) {
        DefMap<std::unique_ptr<GridNodeLbm>>& grid_nodes = *ptr_lbm_grid_nodes_info->GetPointerToLbmGrid();

        Collision(ptr_lbm_grid_nodes_info->NodeFlagNotCollision_, ptr_lbm_grid_nodes_info);

        Stream(ptr_lbm_grid_nodes_info->NodeFlagNotStream_, sfbitset_aux, &grid_nodes);

        // this function is used to reset flags that have change in InformationFromGridOfDifferentLevel
        // since some nodes should not be calculated after transferring information between different levels
        // otherwise the calculated ones will overwrite correct ones
        ptr_lbm_grid_nodes_info->InitialNotComputeNodeFlag();
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
    std::function<void(const DefReal, const GridNodeLbm&, const std::vector<DefReal>&,
        DefReal* const, std::vector<DefReal>* const)> func_macro;
    DefReal dt_lbm = ptr_lbm_grid_nodes_info->ptr_collision_operator_->dt_lbm_;
    DefMap<std::unique_ptr<GridNodeLbm>>& grid_nodes = *ptr_lbm_grid_nodes_info->GetPointerToLbmGrid();

    if (ptr_lbm_grid_nodes_info->ptr_lbm_grid_nodes_ != nullptr) {
        if (ptr_lbm_grid_nodes_info->bool_forces_) {
            if (grid_nodes.begin()->second->force_.size() != k0SolverDims_) {
                amrproject::LogManager::LogError("Size of forces should be " + std::to_string(k0SolverDims_)
                    + "rather than " + std::to_string(grid_nodes.begin()->second->force_.size()) + " in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            func_macro = func_macro_with_force_;
        } else {
            func_macro = [this](const DefReal dt_lbm, const GridNodeLbm& node, const std::vector<DefReal>& /*force*/,
                DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
                func_macro_without_force_(node, ptr_rho, ptr_velocity);};
        }
        std::vector<DefReal> force;
        for (auto& iter_node : grid_nodes) {
            if (iter_node.second->flag_status_ & flag_not_compute) {
            } else {
                force = GetAllForcesForANode(k0SolverDims_, *iter_node.second.get());
                func_macro(dt_lbm, *iter_node.second.get(), force,
                    &iter_node.second->rho_, &iter_node.second->velocity_);
                ptr_lbm_grid_nodes_info->ptr_collision_operator_->CollisionOperator(
                    *this, force, iter_node.second.get());
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
    if (k0Velocity_.size() == k0SolverDims_) {
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
    for (int iq = 0; iq < this->k0NumQ_; ++iq) {
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
    for (int iq = 0; iq < this->k0NumQ_; ++iq) {
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
    for (int iq = 0; iq < this->k0NumQ_; ++iq) {
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
    for (int iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = velocity.at(kXIndex) * k0Cx_.at(iq) + velocity.at(kYIndex) * k0Cy_.at(iq)
            + velocity.at(kZIndex) * k0Cz_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) *  (rho + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
/**
 * @brief function to calculate the force term for a given lattice direction in 2D.
 * @param[in] iq the ith lattice direction. 
 * @param[in] node grid node containing LBM related information.
 * @param[in] force force should be added.
 * @return the calculated LBM body force term
 */
DefReal SolverLbmInterface::CalForceIq2D(const int iq, const GridNodeLbm& node,
    const std::vector<DefReal>& force) const {
    return k0Weights_[iq] * (3.*((k0Cx_[iq] - node.velocity_[kXIndex]) * force[kXIndex]
        + (k0Cy_[iq] - node.velocity_[kYIndex]) * force[kYIndex])
        + 9. * (node.velocity_.at(kXIndex) * k0Cx_.at(iq) + node.velocity_.at(kYIndex) * k0Cy_.at(iq))
        * (k0Cx_[iq] * force[kXIndex] + k0Cy_[iq] * force[kYIndex]));
}
/**
 * @brief function to calculate the force term for a given lattice direction in 3D.
 * @param[in] iq the ith lattice direction. 
 * @param[in] node grid node containing LBM related information.
 * @param[in] force force should be added.
 * @return the calculated LBM body force term
 */
DefReal SolverLbmInterface::CalForceIq3D(const int iq, const GridNodeLbm& node,
    const std::vector<DefReal>& force) const {
    return k0Weights_[iq] * (3.*((k0Cx_[iq] - node.velocity_[kXIndex]) * force[kXIndex]
        + (k0Cy_[iq] - node.velocity_[kYIndex]) * force[kYIndex]
        + (k0Cz_[iq] - node.velocity_[kZIndex]) * force[kZIndex])
        + 9. * (node.velocity_.at(kXIndex) * k0Cx_.at(iq) + node.velocity_.at(kYIndex) * k0Cy_.at(iq)
            + node.velocity_.at(kZIndex) * k0Cz_.at(iq))
        * (k0Cx_[iq] * force[kXIndex] + k0Cy_[iq] * force[kYIndex] + k0Cz_[iq] * force[kZIndex]));
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] node LBM node information.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_rho pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro2DCompressible(const GridNodeLbm& node,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacro2DIncompressible(node, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] node LBM node information.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro2DIncompressible(const GridNodeLbm& node,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    rho = node.f_[0];
    (*ptr_velocity) = {0, 0};
    for (DefInt iq = 1; iq < k0NumQ_; ++iq) {
        (*ptr_rho) += node.f_[iq];
        ptr_velocity->at(kXIndex) += k0Cx_[iq] * node.f_[iq];
        ptr_velocity->at(kYIndex) += k0Cy_[iq] * node.f_[iq];
    }
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] node LBM node information.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro3DCompressible(const GridNodeLbm& node,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacro3DIncompressible(node, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
    ptr_velocity->at(kZIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (without forcing term).
 * @param[in] node LBM node information.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacro3DIncompressible(const GridNodeLbm& node,
    DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    rho = node.f_[0];
    (*ptr_velocity) = {0, 0, 0};
    for (DefInt iq = 1; iq < k0NumQ_; ++iq) {
        (*ptr_rho) += node.f_[iq];
        ptr_velocity->at(kXIndex) += k0Cx_[iq] *  node.f_[iq];
        ptr_velocity->at(kYIndex) += k0Cy_[iq] *  node.f_[iq];
        ptr_velocity->at(kZIndex) += k0Cz_[iq] *  node.f_[iq];
    }
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] node LBM node information.
 * @param[in] force forcing term.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce2DCompressible(const DefReal dt_lbm, const GridNodeLbm& node,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacroForce2DIncompressible(dt_lbm, node, force, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] node LBM node information.
 * @param[in] force forcing term.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce2DIncompressible(const DefReal dt_lbm, const GridNodeLbm& node,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacro2DIncompressible(node, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex) += 0.5 * node.force_[kXIndex] * dt_lbm;
    ptr_velocity->at(kYIndex) += 0.5 * node.force_[kYIndex] * dt_lbm;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] node LBM node information.
 * @param[in] force forcing term.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce3DCompressible(const DefReal dt_lbm, const GridNodeLbm& node,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacroForce3DIncompressible(dt_lbm, node, force, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex)/=rho;
    ptr_velocity->at(kYIndex)/=rho;
    ptr_velocity->at(kZIndex)/=rho;
}
/**
 * @brief function to calculate macroscopic variables based on distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param[in] node LBM node information.
 * @param[in] force forcing term.
 * @param[out] ptr_rho pointer to fluid density.
 * @param[out] ptr_velocity pointer to fluid velocity.
 */
void SolverLbmInterface::CalMacroForce3DIncompressible(const DefReal dt_lbm, const GridNodeLbm& node,
    const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const {
    DefReal& rho = (*ptr_rho);
    CalMacro3DIncompressible(node, ptr_rho, ptr_velocity);
    ptr_velocity->at(kXIndex) += 0.5 * node.force_[kXIndex] * dt_lbm;
    ptr_velocity->at(kYIndex) += 0.5 * node.force_[kYIndex] * dt_lbm;
    ptr_velocity->at(kZIndex) += 0.5 * node.force_[kZIndex] * dt_lbm;
}
}  // end namespace lbmproject
}  // end namespace rootproject
