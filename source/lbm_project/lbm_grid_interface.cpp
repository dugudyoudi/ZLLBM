//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_grid_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM grid interface.
* @date  2023-9-30
*/
#include <mpi.h>
#include "lbm_interface.h"
#include "io/log_write.h"
#include "io/vtk_writer.h"
namespace rootproject {
namespace lbmproject {
    /**
 * @brief function to set variables stored in the node as zeroes.
 * @param[out] ptr_node pointer to a node.
 */
void GridInfoLbmInteface::SetNodeVariablesAsZeros(amrproject::GridNode* const ptr_node) {
    SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    GridNodeLbm* ptr_lbm_node = dynamic_cast<GridNodeLbm*>(ptr_node);
    ptr_lbm_node->rho_ = 0.;
    ptr_lbm_node->velocity_.assign(ptr_lbm_solver->k0SolverDims_, 0.);
    ptr_lbm_node->f_collide_.assign(ptr_lbm_solver->k0NumQ_, 0.);
}
/**
 * @brief function to create grid node.
 */
std::unique_ptr<amrproject::GridNode> GridInfoLbmInteface::GridNodeCreator() {
    const SolverLbmInterface& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_));
    if (!bool_forces_) {
        return std::make_unique<GridNodeLbm>(
            lbm_solver.k0Rho_, lbm_solver.k0Velocity_, f_ini_, f_collide_ini_);
    } else {
        return std::make_unique<GridNodeLbm>(
            lbm_solver.k0Rho_, lbm_solver.k0Velocity_, lbm_solver.k0Force_, f_ini_, f_collide_ini_);
    }
}
/**
 * @brief function to initializes the grid information.
 */
void GridInfoLbmInteface::InitialGridInfo() {
    SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    ptr_lbm_solver->SetInitialDisFuncBasedOnReferenceMacros(&f_ini_, &f_collide_ini_);
    k0SizeOfAllDistributionFunctions_ = ptr_lbm_solver->k0NumQ_*sizeof(DefReal);
    SetCollisionOperator();
    ptr_collision_operator_->viscosity_lbm_ = ptr_lbm_solver->k0LbmViscosity_;
    ptr_collision_operator_->dt_lbm_ = 1./ static_cast<DefReal>(TwoPowerN(i_level_));
    ptr_collision_operator_->CalRelaxationTime();
    ptr_collision_operator_->CalRelaxationTimeRatio();
    switch (interp_method_) {
    case amrproject::EInterpolationMethod::kLinear:
        max_interp_length_ = 1;
        if (ptr_lbm_solver->k0SolverDims_ == 2) {
            func_node_interp_ = [this](const DefAmrIndexLUint interp_length, const DefAmrIndexLUint region_length,
                const DefAmrUint flag_not_for_interp_coarse,
                const DefSFBitset& sfbitset_in, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
                const std::vector<DefSFBitset>& sfbitset_region, const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_fine,
                const amrproject::GridInfoInterface& coarse_grid_info,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_coarse, GridNodeLbm* const ptr_node) {
                    return this->InterpolationLinear2D<GridNodeLbm>(
                        region_length, flag_not_for_interp_coarse, sfbitset_in, sfbitset_aux,
                        sfbitset_region, nodes_fine, coarse_grid_info, nodes_coarse, ptr_node);
            };
        } else {
            func_node_interp_ = [this](const DefAmrIndexLUint interp_length, const DefAmrIndexLUint region_length,
                const DefAmrUint flag_not_for_interp_coarse,
                const DefSFBitset& sfbitset_in, const  amrproject::SFBitsetAuxInterface& sfbitset_aux,
                const std::vector<DefSFBitset>& sfbitset_region,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_fine,
                const amrproject::GridInfoInterface& coarse_grid_info,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_coarse, GridNodeLbm* const ptr_node) {
                    return this->InterpolationLinear3D<GridNodeLbm>(
                        region_length, flag_not_for_interp_coarse, sfbitset_in, sfbitset_aux,
                        sfbitset_region, nodes_fine, coarse_grid_info, nodes_coarse, ptr_node);
            };
        }
        break;
    case amrproject::EInterpolationMethod::kLagrangian:
        if (ptr_lbm_solver->k0SolverDims_ == 2) {
            func_node_interp_ = [this](const DefAmrIndexLUint interp_length, const DefAmrIndexLUint region_length,
                const DefAmrUint flag_not_for_interp_coarse,
                const DefSFBitset& sfbitset_in, const  amrproject::SFBitsetAuxInterface& sfbitset_aux,
                const std::vector<DefSFBitset>& sfbitset_region,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_fine,
                const amrproject::GridInfoInterface& coarse_grid_info,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_coarse, GridNodeLbm* const ptr_node) {
                    return this->InterpolationLagrangian2D<GridNodeLbm>(interp_length,
                        region_length, flag_not_for_interp_coarse, sfbitset_in, sfbitset_aux,
                        sfbitset_region, nodes_fine, coarse_grid_info, nodes_coarse, ptr_node);
            };
        } else {
            func_node_interp_ = [this](const DefAmrIndexLUint interp_length, const DefAmrIndexLUint region_length,
                const DefAmrUint flag_not_for_interp_coarse,
                const DefSFBitset& sfbitset_in, const  amrproject::SFBitsetAuxInterface& sfbitset_aux,
                const std::vector<DefSFBitset>& sfbitset_region,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_fine,
                const amrproject::GridInfoInterface& coarse_grid_info,
                const DefMap<std::unique_ptr<GridNodeLbm>>& nodes_coarse, GridNodeLbm* const ptr_node) {
                    return this->InterpolationLagrangian3D<GridNodeLbm>(interp_length,
                        region_length, flag_not_for_interp_coarse, sfbitset_in, sfbitset_aux,
                        sfbitset_region, nodes_fine, coarse_grid_info, nodes_coarse, ptr_node);
            };
        }
    default:
        break;
    }
}
/**
* @brief  function to reinterpret type of grid nodes as LBM node type.
*/
void GridInfoLbmInteface::SetPointerToCurrentNodeType() {
    if (!map_grid_node_.empty()) {
        auto& first_element = map_grid_node_.begin()->second;
        if (dynamic_cast<GridNodeLbm*>(first_element.get())) {
            // The elements in map_nodes are of type GridNodeLbm,
            // assuming all nodes in map_grid_node_ are the same type.
            ptr_lbm_grid_nodes_ = reinterpret_cast<DefMap<std::unique_ptr<GridNodeLbm>>*>(&map_grid_node_);
        } else {
            std::string msg = "type of nodes stored in map_grid_node_ is not GridNodeLbm, "
            "please check if appropriate node creator is available in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__);
            amrproject::LogManager::LogError(msg);
        }
    } else {
        // noting that the map is empty
        ptr_lbm_grid_nodes_ = reinterpret_cast<DefMap<std::unique_ptr<GridNodeLbm>>*>(&map_grid_node_);
    }
}
/**
 * @brief function to get pointer to the map store LBM nodes.
 */
DefMap<std::unique_ptr<GridNodeLbm>>* GridInfoLbmInteface::GetPointerToLbmGrid() {
    if (ptr_lbm_grid_nodes_ == nullptr) {
        SetPointerToCurrentNodeType();
    }
    return ptr_lbm_grid_nodes_;
}
/**
 * @brief function to call assigned boundary conditions for each domain boundary.
 */
void GridInfoLbmInteface::ComputeDomainBoundaryCondition() {
    for (auto i =0; i < domain_boundary_min_.size(); ++i) {
        switch (i) {
        case kXIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXNeg)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXNeg)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryXNeg, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for x minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kYIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYNeg)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYNeg)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryYNeg, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for y minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kZIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZNeg)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZNeg)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryZNeg, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for z minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        default:
            amrproject::LogManager::LogError("Type of boundary not defined for the " + std::to_string(i)
                + "th minimum domain boundary in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            break;
        }
    }
    for (auto i =0; i < domain_boundary_max_.size(); ++i) {
        switch (i) {
        case kXIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXPos)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXPos)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryXPos, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for x maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kYIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYPos)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYPos)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryYPos, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for y maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kZIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZPos)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZPos)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryZPos, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for z maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        default:
            amrproject::LogManager::LogError("Type of boundary not defined for the " + std::to_string(i)
                + "th maximum domain boundary in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            break;
        }
    }
}
/**
 * @brief function to setup grid information at the beginning of each time step.
 * @param[in] time_step count of time steps.
 * @param[out] ptr_grid_manager point to class manages mesh containing grids at different refinement levels.
 */
void GridInfoLbmInteface::SetUpGridAtBeginningOfTimeStep(const DefAmrIndexUint time_step,
    amrproject::GridManagerInterface* const ptr_grid_manager) {
    // reinterpret pointer to grid information at one level lower and one level higher
    if (i_level_ > 0) {
        std::dynamic_pointer_cast<GridInfoLbmInteface>(
            ptr_grid_manager->vec_ptr_grid_info_.at(i_level_ - 1))->GetPointerToLbmGrid();
    }
    if (i_level_  < ptr_grid_manager->vec_ptr_grid_info_.size() - 1) {
        std::dynamic_pointer_cast<GridInfoLbmInteface>(
            ptr_grid_manager->vec_ptr_grid_info_.at(i_level_ + 1))->GetPointerToLbmGrid();
    }
}
/**
 * @brief function to initialize node that does not need to be computed.
 * @param[in] grid_manager class manages mesh containing flags to classify nodes.
 */
void GridInfoLbmInteface::InitialNotComputeNodeFlag() {
    NodeFlagNotCollision_ = amrproject::NodeBitStatus::kNodeStatusMpiPartitionOutside_;
    NodeFlagNotStream_ = amrproject::NodeBitStatus::kNodeStatusMpiPartitionOutside_;
}
/**
 * @brief function to transfer information on the interface from the coarse grid to fine grid.
 * @param sfbitset_aux class to manage functions for space filling code computation.
 * @param node_flag flag indicates coarse node status, which will be not used for interpolation.
 * @param grid_info_coarse grid information on the coarse grid.
 * @return 0 run successfully, 1 does not do anything since find or coarse grid is not available.
 * @note information on the layers from 1 to k0NumFine2CoarseLayer_ is transferred from the coarse grid.
 */
int GridInfoLbmInteface::TransferInfoFromCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const DefAmrUint node_flag_not_interp, const amrproject::GridInfoInterface& grid_info_coarse) {
    DefAmrIndexUint dims = ptr_solver_->k0SolverDims_;
    if (GetPointerToLbmGrid() == nullptr ||
        (dynamic_cast<const GridInfoLbmInteface&>(grid_info_coarse).ptr_lbm_grid_nodes_ == nullptr)) {
        return 1;
    }
    // identify periodic boundaries
    std::vector<bool> periodic_min(dims), periodic_max(dims);
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXNeg) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXNeg)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_min.at(kXIndex) = true;
    } else {
        periodic_min.at(kXIndex) = false;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYNeg) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYNeg)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_min.at(kYIndex) = true;
    } else {
        periodic_min.at(kYIndex) = false;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXPos) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXPos)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_max.at(kXIndex) = true;
    } else {
        periodic_max.at(kXIndex) = false;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYPos) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYPos)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_max.at(kYIndex) = true;
    } else {
        periodic_max.at(kYIndex) = false;
    }
    if (dims == 3) {
        if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZNeg) != domain_boundary_condition_.end()
            && (domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZNeg)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic)) {
            periodic_min.at(kZIndex) = true;
        } else {
            periodic_min.at(kZIndex) = false;
        }
        if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZPos) != domain_boundary_condition_.end()
            && (domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZPos)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic)) {
            periodic_max.at(kZIndex) = true;
        } else {
            periodic_max.at(kZIndex) = false;
        }
    }

    // calculate domain min and max at one level lower
    std::vector<DefSFBitset> domain_min(dims), domain_max(dims);
    for (DefAmrIndexUint i_dims = 0; i_dims < dims; ++i_dims) {
        domain_min.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(
            1, k0VecBitsetDomainMin_.at(i_dims));
        domain_max.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(
            1, k0VecBitsetDomainMax_.at(i_dims));
    }

    DefSFBitset sfbitset_coarse;
    std::vector<DefSFBitset> nodes_in_region;
    DefAmrIndexLUint valid_length = max_interp_length_;
    SolverLbmInterface* lbm_solver_manager = dynamic_cast<SolverLbmInterface*>(ptr_solver_.get());
    const GridInfoLbmInteface& lbm_grid_coarse =
        dynamic_cast<const GridInfoLbmInteface&>(grid_info_coarse);
    for (auto& iter_interface : map_ptr_interface_layer_info_) {
        for (DefAmrIndexUint i_layer = 1; i_layer < k0NumFine2CoarseLayer_; ++i_layer) {
            // layers on inner interfaces
            for (auto& iter_node : iter_interface.second->vec_inner_fine2coarse_.at(i_layer)) {
                sfbitset_coarse = sfbitset_aux.SFBitsetToNLowerLevelVir(1, iter_node.first);
                valid_length = sfbitset_aux.FindNodesInPeriodicReginOfGivenLength(sfbitset_coarse, max_interp_length_,
                    periodic_min, periodic_max, domain_min, domain_max, &nodes_in_region);
                SetNodeVariablesAsZeros(ptr_lbm_grid_nodes_->at(iter_node.first).get());
                func_node_interp_(valid_length, max_interp_length_, node_flag_not_interp, iter_node.first, sfbitset_aux,
                    nodes_in_region, *ptr_lbm_grid_nodes_, grid_info_coarse, *lbm_grid_coarse.ptr_lbm_grid_nodes_,
                    ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
            // layers on outer interfaces
            for (auto& iter_node : iter_interface.second->vec_outer_fine2coarse_.at(i_layer)) {
                sfbitset_coarse = sfbitset_aux.SFBitsetToNLowerLevelVir(1, iter_node.first);
                valid_length = sfbitset_aux.FindNodesInPeriodicReginOfGivenLength(sfbitset_coarse, max_interp_length_,
                    periodic_min, periodic_max, domain_min, domain_max, &nodes_in_region);
                SetNodeVariablesAsZeros(ptr_lbm_grid_nodes_->at(iter_node.first).get());
                func_node_interp_(valid_length, max_interp_length_, node_flag_not_interp, iter_node.first, sfbitset_aux,
                    nodes_in_region, *ptr_lbm_grid_nodes_, grid_info_coarse, *lbm_grid_coarse.ptr_lbm_grid_nodes_,
                    ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
        }
    }
    return 0;
}
/**
 * @brief function to transfer information on the interface from the fine grid to coarse grid.
 * @param sfbitset_aux class to manage functions for space filling code computation.
 * @param node_flag node indicator does not used in this implementation.
 * @param grid_info_fine grid information on the fine grid.
 * @return 0 run successfully, 1 does not do anything since find or coarse grid is not available.
 * @note information on the layers from 1 to k0NumCoarse2FineLayer_ is transferred from the fine grid.
 */
int GridInfoLbmInteface::TransferInfoFromFineGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const DefAmrUint node_flag, const amrproject::GridInfoInterface& grid_info_fine) {
    DefSFBitset sfbitset_fine;
    if (GetPointerToLbmGrid() == nullptr ||
        (dynamic_cast<const GridInfoLbmInteface&>(grid_info_fine).ptr_lbm_grid_nodes_ == nullptr)) {
        return 1;
    }
    for (auto& iter_interface : map_ptr_interface_layer_info_) {
        for (DefAmrIndexUint i_layer = 1; i_layer < k0NumCoarse2FineLayer_; ++i_layer) {
            // layers on inner interfaces
            for (auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(i_layer)) {
                sfbitset_fine = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_node.first);
                grid_info_fine.NodeInfoFine2Coarse(*(dynamic_cast<const GridInfoLbmInteface&>(grid_info_fine)
                    .ptr_lbm_grid_nodes_->at(sfbitset_fine).get()), ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
            // layers on outer interfaces
            for (auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(i_layer)) {
                sfbitset_fine = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_node.first);
                grid_info_fine.NodeInfoFine2Coarse(*(dynamic_cast<const GridInfoLbmInteface&>(grid_info_fine)
                    .ptr_lbm_grid_nodes_->at(sfbitset_fine).get()), ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
        }
    }
    return 0;
}
/**
 * @brief function to transfer information stored in coarse LBM node to fine node.
 * @param[in]  coarse_base_node   reference of coarse LBM node.
 * @param[in]  ptr_base_fine_node   pointer to fine LBM node.
 */
void GridInfoLbmInteface::NodeInfoCoarse2fine(const amrproject::GridNode& coarse_base_node,
    amrproject::GridNode* const ptr_base_fine_node) const {
    const GridNodeLbm& coarse_node = dynamic_cast<const GridNodeLbm&>(coarse_base_node);
    GridNodeLbm* ptr_fine_node = dynamic_cast<GridNodeLbm*>(ptr_base_fine_node);
    std::vector<DefReal> feq;
    const SolverLbmInterface& lbm_solver = *std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    lbm_solver.func_cal_feq_(coarse_node.rho_, coarse_node.velocity_, &feq);
    if (bool_forces_) {
        ptr_collision_operator_->Coarse2FineForce(ptr_collision_operator_->dt_lbm_,
            feq, lbm_solver, lbm_solver.ptr_func_cal_force_iq_, coarse_node,  ptr_fine_node);
    } else {
        ptr_collision_operator_->Coarse2Fine(ptr_collision_operator_->dt_lbm_,
            feq, coarse_node,  ptr_fine_node);
    }
}
/**
 * @brief function to transfer information stored in fine LBM node to coarse node.
 * @param[in]  fine_base_node   reference of fine LBM node.
 * @param[in]  ptr_base_coarse_node   pointer to coarse LBM node.
 */
void GridInfoLbmInteface::NodeInfoFine2Coarse(const amrproject::GridNode& fine_base_node,
    amrproject::GridNode* const ptr_base_coarse_node) const {
    const GridNodeLbm& fine_node = dynamic_cast<const GridNodeLbm&>(fine_base_node);
    GridNodeLbm* ptr_coarse_node = dynamic_cast<GridNodeLbm*>(ptr_base_coarse_node);
    std::vector<DefReal> feq;
    const SolverLbmInterface& lbm_solver = *std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    lbm_solver.func_cal_feq_(fine_node.rho_, fine_node.velocity_, &feq);
    if (bool_forces_) {
        ptr_collision_operator_->Fine2CoarseForce(ptr_collision_operator_->dt_lbm_,
            feq, lbm_solver, lbm_solver.ptr_func_cal_force_iq_, fine_node,  ptr_coarse_node);
    } else {
        ptr_collision_operator_->Fine2Coarse(ptr_collision_operator_->dt_lbm_,
            feq, fine_node,  ptr_coarse_node);
    }
}
/**
 * @brief function to setup output infomation on variables of LBM node.
 */
void GridInfoLbmInteface::SetupOutputVariables() {
    for (const auto& iter_str : output_variable_name_) {
        if (iter_str == "rho") {
            output_variables_.push_back(std::make_unique<OutputLBMNodeVariableInfo>());
            OutputLBMNodeVariableInfo& output_info_temp =
                *dynamic_cast<OutputLBMNodeVariableInfo*>(output_variables_.back().get());
            output_info_temp.variable_dims_ = 1;
            output_info_temp.func_get_var_ = [] (const GridNodeLbm& node)->std::vector<DefReal> {
                return {node.rho_};
            };
            output_info_temp.output_name_ = "rho";
        } else if (iter_str == "velocity") {
            output_variables_.push_back(std::make_unique<OutputLBMNodeVariableInfo>());
            OutputLBMNodeVariableInfo& output_info_temp =
                *dynamic_cast<OutputLBMNodeVariableInfo*>(output_variables_.back().get());
            output_info_temp.variable_dims_ = ptr_solver_->k0SolverDims_;
            output_info_temp.func_get_var_ = [] (const GridNodeLbm& node)->std::vector<DefReal> {
                return node.velocity_;
            };
            output_info_temp.output_name_ = "velocity";
        }
    }
    GetPointerToLbmGrid();
}
/**
* @brief   function to write scalar and vector variables.
* @param[in]  fp   pointer to output file.
* @param[in]  base64_instance reference to class to convert data to uint8 type.
* @param[in]  bool_binary write data in binary or ascii format.
* @param[in]  output_data_format output data (real or integer) format.
* @param[in]  map_node_index indices of nodes for unstructured grid.
*/
void GridInfoLbmInteface::WriteOutputScalarAndVectors(
    FILE* const fp, const bool bool_binary,
    const amrproject::Base64Utility& base64_instance,
    const amrproject::OutputDataFormat& output_data_format,
    const DefMap<DefSizet>& map_node_index) const {
    for (auto& iter_var : output_variables_) {
        OutputOneVariable(fp, bool_binary, base64_instance,
            *dynamic_cast<OutputLBMNodeVariableInfo*>(iter_var.get()), output_data_format, map_node_index);
    }
}
/**
* @brief function to write a scalar of each node.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in binary or ascii format.
* @param[in]  output_info output infomation of a node variable.
* @param[in]  base64_instance reference to class to convert data to uint8 type.
* @param[in]  output_data_format output data (real or integer) format.
* @param[in]  map_node_index indices of nodes for unstructured grid.
*/
int GridInfoLbmInteface::OutputOneVariable(
    FILE* const fp, const bool bool_binary,
    const amrproject::Base64Utility& base64_instance,
    const OutputLBMNodeVariableInfo& output_info,
    const amrproject::OutputDataFormat& output_data_format,
    const DefMap<DefSizet>& map_node_index) const {
    const DefMap<std::unique_ptr<GridNodeLbm>>& map_grid_node = *ptr_lbm_grid_nodes_;
    std::string str_format, str_temp;
    DefAmrIndexUint dims = output_info.variable_dims_;
    if (bool_binary) {
        str_format =  "binary";
    } else {
        str_format =  "ascii";
    }
    str_temp.assign("      <DataArray NumberOfComponents=\"" + std::to_string(dims)
        + "\" type=\"" + output_data_format.output_real_.format_name_
        + "\" Name=\"" + output_info.output_name_ + "\" format=\"" + str_format + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            const std::vector<DefReal>& vec_var =  output_info.func_get_var_(*map_grid_node.at(iter->first).get());
            for (DefAmrIndexUint i_dims = 0;  i_dims < dims; ++i_dims) {
                base64_instance.AddToVecChar(output_data_format.output_real_.CastType(vec_var.at(i_dims)), &vec_uint8);
            }
        }
        base64_instance.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = "     "
            + output_data_format.output_real_.printf_format_;
        for (auto iter = map_node_index.begin(); iter != map_node_index.end(); ++iter) {
            const std::vector<DefReal>& vec_var =  output_info.func_get_var_(*map_grid_node.at(iter->first).get());
            for (DefAmrIndexUint i_dims = 0; i_dims < dims; ++i_dims) {
                fprintf_s(fp, "  ");
                fprintf_s(fp, str_format.c_str(), vec_var.at(i_dims));
            }
            fprintf_s(fp, "\n");
        }
    }
    fprintf_s(fp, "      </DataArray>\n");
    return 0;
}
}  // end namespace lbmproject
}  // end namespace rootproject