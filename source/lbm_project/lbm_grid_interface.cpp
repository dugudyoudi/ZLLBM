//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_grid_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM grid interface.
* @date  2023-9-30
*/
#include <iomanip>
#include "./lbm_interface.h"
#include "io/log_write.h"
#include "io/vtk_writer.h"
#include "criterion/criterion_manager.h"
#include "criterion/geometry_default_shape.h"
#include "grid/grid_info_interface.h"
#include "mpi/mpi_manager.h"
#include "immersed_boundary/immersed_boundary.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to addition assign coefficient for interpolation.
 * @param[in] node_in input grid node.
 * @param[in] coefficient coefficient for interpolation.
 */
void GridNodeLbm::InterpolationAdditionAssignCoefficient(
    const amrproject::GridNode& node_in, const DefReal coefficient) {
    GridNodeLbm& this_node = dynamic_cast<GridNodeLbm&>(*this);
    const GridNodeLbm& lbm_node_in = dynamic_cast<const GridNodeLbm&>(node_in);
    this_node.rho_ += lbm_node_in.rho_ * coefficient;
    this_node.velocity_.resize(lbm_node_in.velocity_.size());
    std::transform(lbm_node_in.velocity_.begin(), lbm_node_in.velocity_.end(),
        this_node.velocity_.begin(), this_node.velocity_.begin(),
        [coefficient](DefReal num1, DefReal num2) { return num2 + coefficient*num1;});

    this_node.f_.resize(lbm_node_in.f_.size());
    std::transform(lbm_node_in.f_.begin(), lbm_node_in.f_.end(),
        this_node.f_.begin(), this_node.f_.begin(),
        [coefficient](DefReal num1, DefReal num2) { return num2 + coefficient*num1;});
}
/**
 * @brief function to set variables stored in the node as zeroes.
 * @param[out] ptr_node pointer to a node.
 */
void GridInfoLbmInteface::SetNodeVariablesAsZeros(amrproject::GridNode* const ptr_node) {
    if (auto ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_.lock())) {
        GridNodeLbm* ptr_lbm_node = dynamic_cast<GridNodeLbm*>(ptr_node);
        ptr_lbm_node->rho_ = 0.;
        ptr_lbm_node->velocity_.assign(ptr_lbm_solver->GetSolverDim(), 0.);
        ptr_lbm_node->f_.assign(ptr_lbm_solver->k0NumQ_, 0.);
    } else {
        amrproject::LogManager::LogError("Point to solver has been expired");
    }
}
/**
 * @brief function to advancing simulation at current time step.
 * @param[in] time_step_current total time step at current grid refinement level.
 * @param[in] mpi_manager class managing MPI communication.
 * @param[out] ptr_criterion_manager pointer to managing criterion for refinement.
 */
void GridInfoLbmInteface::UpdateCriterion(const DefReal time_step_current,
    const amrproject::MpiManager& mpi_manager, amrproject::CriterionManager* const ptr_criterion_manager) {
    bool first_ib = true;
    std::vector<GeometryInfoImmersedBoundary*> vec_ib_geo;
    for (auto& iter_geo : ptr_criterion_manager->vec_ptr_geometries_) {
        // update when level of geometry and grid are the same
        if (iter_geo->GetLevel() == i_level_) {
            if (iter_geo->GetStatus() == amrproject::EGeometryStatus::kMoving) {
                iter_geo->UpdateGeometry(time_step_current);
            }
            if (auto derived_ptr = dynamic_cast<GeometryInfoImmersedBoundary*>(iter_geo.get())) {
                if (first_ib) {
                    derived_ptr->ClearNodesRecordForIB();
                    first_ib = false;
                }
                derived_ptr->CalculateBodyForce(*this, &(iter_geo->map_vertices_info_), GetPtrToLbmGrid());
                vec_ib_geo.emplace_back(dynamic_cast<GeometryInfoImmersedBoundary*>(iter_geo.get()));
            }
        }
        if (i_level_ == 0) {
            if (auto derived_ptr = dynamic_cast<GeometryInfoImmersedBoundary*>(iter_geo.get())) {
                derived_ptr->WriteTimeHisLagrangianForce(time_step_current, grid_space_[kXIndex]);
            }
        }
    }
#ifdef ENABLE_MPI
    if (!first_ib) {
        DefInt dim = ptr_parent_grid_manager_->k0GridDims_;
        for (auto& iter_geo : vec_ib_geo) {
            iter_geo->SendNReceiveNodesForImmersedBoundary(dim,
                i_level_, mpi_manager, ptr_lbm_grid_nodes_);
        }
    }
#endif  //  ENABLE_MPI
}
/**
 * @brief function to advancing simulation at current time step.
 * @param[in] time_scheme enum class to identify time stepping scheme used in computation.
 * @param[in] time_step_level time step at current grid refinement level in one background step.
 * @param[in] time_step_current total time step at current grid refinement level.
 * @param[in] ptr_mpi_manager pointer to class managing MPI communication.
 * @param[out] ptr_criterion_manager pointer to managing criterion for refinement.
 */
void GridInfoLbmInteface::AdvancingAtCurrentTime(const amrproject::ETimeSteppingScheme time_scheme,
    const DefInt time_step_level, const DefReal time_step_current,
    amrproject::MpiManager* const ptr_mpi_manager,
    amrproject::CriterionManager* const ptr_criterion_manager) {
    const DefInt i_level = i_level_;
    GetPtrToLbmGrid();
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }

    UpdateCriterion(time_step_current, *ptr_mpi_manager, ptr_criterion_manager);

#ifdef ENABLE_MPI
    std::vector<amrproject::MpiManager::BufferSizeInfo> send_buffer_info, receive_buffer_info;
    std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send, vec_ptr_buffer_receive;
    if (i_level > 0) {
        ComputeInfoInInterpMpiLayers(interp_nodes_inner_layer_);
    }
    ptr_mpi_manager->ProcessMpiLayersInfoAndCommunicate(&send_buffer_info, &receive_buffer_info,
        &vec_vec_reqs_send, &vec_vec_reqs_receive, &vec_ptr_buffer_send, &vec_ptr_buffer_receive, this);
#endif  //  ENABLE_MPI

    ptr_lbm_solver->RunSolverOnGivenGrid(time_scheme, time_step_level, time_step_current, *ptr_sfbitset_aux_, this);

#ifdef ENABLE_MPI
    std::function<void(const char*,  amrproject::GridNode* const)> func_read_a_node_from_buffer =
        [](const char* const ptr_node_buffer, amrproject::GridNode* const ptr_node) {
        ptr_node->ReadANodeFromBufferForMpi(ptr_node_buffer);
    };
    ptr_mpi_manager->WaitAndReadGridNodesFromBuffer(send_buffer_info,
        receive_buffer_info, vec_ptr_buffer_receive, func_read_a_node_from_buffer,
        &vec_vec_reqs_send, &vec_vec_reqs_receive, this);
#endif  //  ENABLE_MPI

    ComputeDomainBoundaryCondition();

    // use information in current time step
#ifdef ENABLE_MPI
    if (i_level > 0) {
        ptr_mpi_manager->MpiCommunicationForInterpolation(*ptr_sfbitset_aux_,
            *GetPtrToParentGridManager()->vec_ptr_grid_info_.at(i_level - 1).get(), this);
    }
#endif  //  ENABLE_MPI

    ptr_lbm_solver->InformationFromGridOfDifferentLevel(time_step_level, *ptr_sfbitset_aux_, this);
}
/**
 * @brief function to create grid node.
 */
std::unique_ptr<amrproject::GridNode> GridInfoLbmInteface::GridNodeCreator() const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    if (!ptr_lbm_solver->bool_forces_) {
        return std::make_unique<GridNodeLbm>(
            ptr_lbm_solver->GetDefaultDensity(), ptr_lbm_solver->GetDefaultVelocity(), f_ini_, f_collide_ini_);
    } else {
        return std::make_unique<GridNodeLbm>(
            ptr_lbm_solver->GetDefaultDensity(), ptr_lbm_solver->GetDefaultVelocity(),
            ptr_lbm_solver->GetDefaultForce(), f_ini_, f_collide_ini_);
    }
}
/**
 * @brief function to initializes the grid information.
 */
void GridInfoLbmInteface::InitialGridInfo(const DefInt dims) {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    ptr_lbm_solver->SetInitialDisFuncBasedOnReferenceMacros(&f_ini_, &f_collide_ini_);
    ptr_lbm_solver->InstantiateCollisionOperator(i_level_);
    ChooseInterpolationMethod(dims);
    if (ptr_lbm_solver->GetNumForces() != static_cast<DefInt>(ptr_lbm_solver->GetDefaultForce().size())) {
        amrproject::LogManager::LogWarning("Dimension of forces (k0NumForces_) is not equal to"
            " the dimension of default forces (k0Force_) in the solver for"
            " initialization in GridInfoLbmInteface::InitialGridInfo.");
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
        // if the map is empty
        ptr_lbm_grid_nodes_ = reinterpret_cast<DefMap<std::unique_ptr<GridNodeLbm>>*>(&map_grid_node_);
    }
}
/**
 * @brief function to get pointer to the map store LBM nodes.
 */
DefMap<std::unique_ptr<GridNodeLbm>>* GridInfoLbmInteface::GetPtrToLbmGrid() {
    if (ptr_lbm_grid_nodes_ == nullptr) {
        SetPointerToCurrentNodeType();
    }
    return ptr_lbm_grid_nodes_;
}
/**
 * @brief function to call assigned boundary conditions for each domain boundary.
 */
void GridInfoLbmInteface::ComputeDomainBoundaryCondition() {
    for (DefInt i = 0; i < static_cast<DefInt>(domain_boundary_min_.size()); ++i) {
        switch (i) {
        case kXIndex: {
                if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryXMin)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(
                        amrproject::EDomainBoundaryDirection::kBoundaryXMin)->CalBoundaryCondition(
                        amrproject::EDomainBoundaryDirection::kBoundaryXMin, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for x minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kYIndex: {
                if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryYMin)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(
                        amrproject::EDomainBoundaryDirection::kBoundaryYMin)->CalBoundaryCondition(
                        amrproject::EDomainBoundaryDirection::kBoundaryYMin, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for y minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kZIndex: {
                if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryZMin)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(
                        amrproject::EDomainBoundaryDirection::kBoundaryZMin)->CalBoundaryCondition(
                        amrproject::EDomainBoundaryDirection::kBoundaryZMin, domain_boundary_min_.at(i), this);
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
    for (DefInt i =0; i < static_cast<DefInt>(domain_boundary_max_.size()); ++i) {
        switch (i) {
        case kXIndex: {
                if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryXMax)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(
                        amrproject::EDomainBoundaryDirection::kBoundaryXMax)->CalBoundaryCondition(
                        amrproject::EDomainBoundaryDirection::kBoundaryXMax, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for x maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kYIndex: {
                if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryYMax)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(
                        amrproject::EDomainBoundaryDirection::kBoundaryYMax)->CalBoundaryCondition(
                        amrproject::EDomainBoundaryDirection::kBoundaryYMax, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for y maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kZIndex: {
                if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryZMax)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(
                        amrproject::EDomainBoundaryDirection::kBoundaryZMax)->CalBoundaryCondition(
                        amrproject::EDomainBoundaryDirection::kBoundaryZMax, domain_boundary_max_.at(i), this);
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
 * @brief function to call assigned boundary conditions for a given node.
 * @param[in] flag_node_boundary flag indicting node on which boundary.
 * @param[in] bitset_in space filling code of a given node.
 */
void GridInfoLbmInteface:: ComputeDomainBoundaryConditionForANode(
    int flag_node_boundary, const DefSFBitset bitset_in) {
    if (flag_node_boundary > kFlagInsideDomain_) {
        if ((flag_node_boundary & kFlagXMinBoundary_) == kFlagXMinBoundary_) {
            if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryXMin)
                != domain_boundary_condition_.end()) {
                domain_boundary_condition_.at(
                    amrproject::EDomainBoundaryDirection::kBoundaryXMin)->CalBoundaryCondition(
                    amrproject::EDomainBoundaryDirection::kBoundaryXMin, {std::make_pair(bitset_in, 0)}, this);
            }
        }
        if ((flag_node_boundary & kFlagXMaxBoundary_) == kFlagXMaxBoundary_) {
            if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryXMax)
                != domain_boundary_condition_.end()) {
                domain_boundary_condition_.at(
                    amrproject::EDomainBoundaryDirection::kBoundaryXMax)->CalBoundaryCondition(
                    amrproject::EDomainBoundaryDirection::kBoundaryXMax, {std::make_pair(bitset_in, 0)}, this);
            }
        }
        if ((flag_node_boundary & kFlagYMinBoundary_) == kFlagYMinBoundary_) {
            if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryYMin)
                != domain_boundary_condition_.end()) {
                domain_boundary_condition_.at(
                    amrproject::EDomainBoundaryDirection::kBoundaryYMin)->CalBoundaryCondition(
                    amrproject::EDomainBoundaryDirection::kBoundaryYMin, {std::make_pair(bitset_in, 0)}, this);
            }
        }
        if ((flag_node_boundary & kFlagYMaxBoundary_) == kFlagYMaxBoundary_) {
            if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryYMax)
                != domain_boundary_condition_.end()) {
                domain_boundary_condition_.at(
                    amrproject::EDomainBoundaryDirection::kBoundaryYMax)->CalBoundaryCondition(
                    amrproject::EDomainBoundaryDirection::kBoundaryYMax, {std::make_pair(bitset_in, 0)}, this);
            }
        }
        if ((flag_node_boundary & kFlagZMinBoundary_)== kFlagZMinBoundary_) {
            if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryZMin)
                != domain_boundary_condition_.end()) {
                domain_boundary_condition_.at(
                    amrproject::EDomainBoundaryDirection::kBoundaryZMin)->CalBoundaryCondition(
                    amrproject::EDomainBoundaryDirection::kBoundaryZMin, {std::make_pair(bitset_in, 0)}, this);
            }
        }
        if ((flag_node_boundary & kFlagZMaxBoundary_) == kFlagZMaxBoundary_) {
            if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryZMax)
                != domain_boundary_condition_.end()) {
                domain_boundary_condition_.at(
                    amrproject::EDomainBoundaryDirection::kBoundaryZMax)->CalBoundaryCondition(
                    amrproject::EDomainBoundaryDirection::kBoundaryZMax, {std::make_pair(bitset_in, 0)}, this);
            }
        }
    }
}
/**
 * @brief function to setup grid information at the beginning of each time step.
 * @param[in] time_step count of time steps.
 */
void GridInfoLbmInteface::SetUpGridAtBeginningOfTimeStep(const DefInt time_step) {
    amrproject::GridManagerInterface* ptr_grid_manager = GetPtrToParentGridManager();
    // reinterpret pointer to grid information at one level lower and one level higher
    if (i_level_ > 0) {
        std::dynamic_pointer_cast<GridInfoLbmInteface>(
            ptr_grid_manager->vec_ptr_grid_info_.at(i_level_ - 1))->GetPtrToLbmGrid();
    }
    if (i_level_  < static_cast<DefInt>(ptr_grid_manager->vec_ptr_grid_info_.size()) - 1) {
        std::dynamic_pointer_cast<GridInfoLbmInteface>(
            ptr_grid_manager->vec_ptr_grid_info_.at(i_level_ + 1))->GetPtrToLbmGrid();
    }
}
/**
 * @brief function to initialize node that does not need to be computed.
 * @param[in] grid_manager class manages mesh containing flags to classify nodes.
 */
void GridInfoLbmInteface::InitialNotComputeNodeFlag() {
    NodeFlagNotCollision_ = amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_
        |amrproject::NodeBitStatus::kNodeStatusMpiInterpInner_;
    NodeFlagNotStream_ = amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_
        |amrproject::NodeBitStatus::kNodeStatusMpiInterpInner_;
}
/**
 * @brief function to transfer information on the interface from the coarse grid to fine grid.
 * @param sfbitset_aux class to manage functions for space filling code computation.
 * @param node_flag flag indicates coarse node status, which will be not used for interpolation.
 * @param grid_info_coarse grid information on the coarse grid.
 * @return 0 run successfully, otherwise some error occurred in this function.
 */
int GridInfoLbmInteface::TransferInfoFromCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const DefInt node_flag_not_interp, const amrproject::GridInfoInterface& grid_info_coarse) {
    const DefInt dims = ptr_parent_grid_manager_->k0GridDims_;
    if (GetPtrToLbmGrid() == nullptr ||
        (dynamic_cast<const GridInfoLbmInteface&>(grid_info_coarse).ptr_lbm_grid_nodes_ == nullptr)) {
        amrproject::LogManager::LogError("pointer to LBM grid nodes is null");
    }
    // identify periodic boundaries
    std::vector<bool> periodic_min(dims), periodic_max(dims);
    // CheckIfPeriodicDomainRequired(dims, &periodic_min, &periodic_max);
    if (domain_boundary_condition_.find(
        amrproject::EDomainBoundaryDirection::kBoundaryXMin) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryXMin)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_min.at(kXIndex) = true;
    } else {
        periodic_min.at(kXIndex) = false;
    }
    if (domain_boundary_condition_.find(
        amrproject::EDomainBoundaryDirection::kBoundaryYMin) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryYMin)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_min.at(kYIndex) = true;
    } else {
        periodic_min.at(kYIndex) = false;
    }
    if (domain_boundary_condition_.find(
        amrproject::EDomainBoundaryDirection::kBoundaryXMax) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryXMax)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_max.at(kXIndex) = true;
    } else {
        periodic_max.at(kXIndex) = false;
    }
    if (domain_boundary_condition_.find(
        amrproject::EDomainBoundaryDirection::kBoundaryYMax) != domain_boundary_condition_.end()
        && (domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryYMax)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic)) {
        periodic_max.at(kYIndex) = true;
    } else {
        periodic_max.at(kYIndex) = false;
    }
    if (dims == 3) {
        if (domain_boundary_condition_.find(
            amrproject::EDomainBoundaryDirection::kBoundaryZMin) != domain_boundary_condition_.end()
            && (domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryZMin)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic)) {
            periodic_min.at(kZIndex) = true;
        } else {
            periodic_min.at(kZIndex) = false;
        }
        if (domain_boundary_condition_.find(
            amrproject::EDomainBoundaryDirection::kBoundaryZMax) != domain_boundary_condition_.end()
            && (domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryZMax)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic)) {
            periodic_max.at(kZIndex) = true;
        } else {
            periodic_max.at(kZIndex) = false;
        }
    }

    // calculate domain min and max at one level lower
    std::vector<DefSFBitset> domain_min(dims), domain_max(dims);
    for (DefInt i_dims = 0; i_dims < dims; ++i_dims) {
        domain_min.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(
            1, k0VecBitsetDomainMin_.at(i_dims));
        domain_max.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(
            1, k0VecBitsetDomainMax_.at(i_dims));
    }

    DefSFBitset sfbitset_coarse;
    std::vector<DefSFBitset> nodes_in_region;
    DefAmrLUint valid_length = max_interp_length_;
    DefInt start_layer = k0NumFine2CoarseLayer_ - k0NumFine2CoarseGhostLayer_;

    for (auto& iter_interface : map_ptr_interface_layer_info_) {
        for (DefInt i_layer = start_layer; i_layer < k0NumFine2CoarseLayer_; ++i_layer) {
            // layers on inner interfaces
            for (auto& iter_node : iter_interface.second->vec_inner_fine2coarse_.at(i_layer)) {
                sfbitset_coarse = sfbitset_aux.SFBitsetToNLowerLevelVir(1, iter_node.first);
                valid_length = sfbitset_aux.FindNodesInPeriodicRegionCorner(sfbitset_coarse, max_interp_length_,
                    periodic_min, periodic_max, domain_min, domain_max, &nodes_in_region);
                if (ptr_lbm_grid_nodes_->find(iter_node.first) == ptr_lbm_grid_nodes_->end()) {
                    amrproject::LogManager::LogError("node for interpolation is not found");
                }
                func_node_interp_(valid_length, max_interp_length_, node_flag_not_interp, iter_node.first, sfbitset_aux,
                    nodes_in_region, map_grid_node_, grid_info_coarse, grid_info_coarse.map_grid_node_,
                    ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
            // layers on outer interfaces
            for (auto& iter_node : iter_interface.second->vec_outer_fine2coarse_.at(i_layer)) {
                sfbitset_coarse = sfbitset_aux.SFBitsetToNLowerLevelVir(1, iter_node.first);
                valid_length = sfbitset_aux.FindNodesInPeriodicRegionCorner(sfbitset_coarse, max_interp_length_,
                    periodic_min, periodic_max, domain_min, domain_max, &nodes_in_region);
                if (ptr_lbm_grid_nodes_->find(iter_node.first) != ptr_lbm_grid_nodes_->end()) {
                    SetNodeVariablesAsZeros(ptr_lbm_grid_nodes_->at(iter_node.first).get());
                } else {
                    amrproject::LogManager::LogError("node for interpolation is not found");
                }
                func_node_interp_(valid_length, max_interp_length_, node_flag_not_interp, iter_node.first, sfbitset_aux,
                    nodes_in_region, map_grid_node_, grid_info_coarse, grid_info_coarse.map_grid_node_,
                    ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
        }
    }

    return 0;
}
/**
 * @brief function to transfer information on the interface from the fine grid to coarse grid.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] node_flag node indicator does not used in this implementation.
 * @param[out] ptr_grid_info_coarse grid information on the fine grid.
 * @return 0 run successfully,  otherwise some error occurred in this function.
 * @note information on the layers from k0NumCoarse2FineGhostLayer_ to k0NumCoarse2FineLayer_ is transferred from the fine grid.
 */
int GridInfoLbmInteface::TransferInfoToCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const DefInt node_flag, amrproject::GridInfoInterface* const ptr_grid_info_coarse) {
    DefSFBitset sfbitset_fine;
    if (GetPtrToLbmGrid() == nullptr ||
        (dynamic_cast<GridInfoLbmInteface*>(ptr_grid_info_coarse)->ptr_lbm_grid_nodes_ == nullptr)) {
        amrproject::LogManager::LogError("pointer to LBM grid is null");
    }
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_grid_coarse =
        dynamic_cast<const GridInfoLbmInteface*>(ptr_grid_info_coarse)->ptr_lbm_grid_nodes_;
    DefInt num_c2f_ghost_layer = ptr_grid_info_coarse->GetNumCoarse2FineGhostLayer(),
        num_c2f_layer = ptr_grid_info_coarse->GetNumCoarse2FineLayer();
    for (auto& iter_interface : ptr_grid_info_coarse->map_ptr_interface_layer_info_) {
        for (DefInt i_layer = num_c2f_ghost_layer; i_layer < num_c2f_layer; ++i_layer) {
            // layers on inner interfaces
            for (auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(i_layer)) {
                sfbitset_fine = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_node.first);
                if (ptr_lbm_grid_nodes_->find(sfbitset_fine) != ptr_lbm_grid_nodes_->end()) {
                    NodeInfoFine2Coarse(*(ptr_lbm_grid_nodes_->at(sfbitset_fine).get()),
                        ptr_lbm_grid_coarse->at(iter_node.first).get());
                } else if (interp_nodes_outer_layer_.find(sfbitset_fine) != interp_nodes_outer_layer_.end()
                    && ptr_lbm_grid_coarse->find(iter_node.first) != ptr_lbm_grid_coarse->end()) {
                    NodeInfoFine2Coarse(*(interp_nodes_outer_layer_.at(sfbitset_fine).get()),
                        ptr_lbm_grid_coarse->at(iter_node.first).get());
                }
            }
            // layers on outer interfaces
            for (auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(i_layer)) {
                sfbitset_fine = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_node.first);
                if (ptr_lbm_grid_nodes_->find(sfbitset_fine) != ptr_lbm_grid_nodes_->end()) {
                    NodeInfoFine2Coarse(*(ptr_lbm_grid_nodes_->at(sfbitset_fine).get()),
                        ptr_lbm_grid_coarse->at(iter_node.first).get());
                } else if (interp_nodes_outer_layer_.find(sfbitset_fine) != interp_nodes_outer_layer_.end()
                    && ptr_lbm_grid_coarse->find(iter_node.first) != ptr_lbm_grid_coarse->end()) {
                    NodeInfoFine2Coarse(*(interp_nodes_outer_layer_.at(sfbitset_fine).get()),
                        ptr_lbm_grid_coarse->at(iter_node.first).get());
                }
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
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    std::vector<DefReal> feq(ptr_lbm_solver->k0NumQ_);
    ptr_fine_node->velocity_.resize(ptr_lbm_solver->GetSolverDim());
    ptr_fine_node->rho_ = coarse_node.rho_;
    ptr_fine_node->velocity_.assign(coarse_node.velocity_.begin(), coarse_node.velocity_.end());

    ptr_lbm_solver->func_macro_without_force_(coarse_node.f_, &ptr_fine_node->rho_, &ptr_fine_node->velocity_);
    ptr_lbm_solver->func_cal_feq_(ptr_fine_node->rho_, ptr_fine_node->velocity_, &feq);
    ptr_lbm_solver->GetCollisionOperator(i_level_).PostStreamCoarse2Fine(feq, coarse_node.f_, &ptr_fine_node->f_);
}
/**
 * @brief function to transfer information stored in fine LBM node to coarse node.
 * @param[in]  fine_base_node   reference of fine LBM node.
 * @param[out]  ptr_base_coarse_node   pointer to coarse LBM node.
 */
void GridInfoLbmInteface::NodeInfoFine2Coarse(const amrproject::GridNode& fine_base_node,
    amrproject::GridNode* const ptr_base_coarse_node) const {
    const GridNodeLbm& fine_node = dynamic_cast<const GridNodeLbm&>(fine_base_node);
    GridNodeLbm* ptr_coarse_node = dynamic_cast<GridNodeLbm*>(ptr_base_coarse_node);
    std::vector<DefReal> feq;
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    ptr_coarse_node->rho_ = fine_node.rho_;
    ptr_coarse_node->velocity_.resize(ptr_lbm_solver->GetSolverDim());
    ptr_coarse_node->velocity_.assign(fine_node.velocity_.begin(), fine_node.velocity_.end());
    ptr_lbm_solver->func_cal_feq_(fine_node.rho_, fine_node.velocity_, &feq);
    ptr_lbm_solver->GetCollisionOperator(i_level_).PostStreamFine2Coarse(feq, fine_node, ptr_coarse_node);
}
/**
 * @brief function to setup output infomation on variables of LBM node.
 */
void GridInfoLbmInteface::SetupOutputVariables() {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    DefReal dt = ptr_lbm_solver->GetCollisionOperator(i_level_).GetDtLbm();
    std::function<void(const DefReal dt_lbm, const std::vector<DefReal>& f, const std::vector<DefReal>& force,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity)> func_macro;
    if (ptr_lbm_solver->bool_forces_) {
        func_macro = ptr_lbm_solver->func_macro_with_force_;
    } else {
        func_macro = [ptr_lbm_solver](const DefReal dt_lbm, const std::vector<DefReal>& f,
            const std::vector<DefReal>& /*force*/, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            ptr_lbm_solver->func_macro_without_force_(f, ptr_rho, ptr_velocity);
        };
    }
    for (const auto& iter_str : output_variable_name_) {
        if (iter_str == "rho") {
            output_variables_.emplace_back(std::make_unique<OutputLBMNodeVariableInfo>());
            OutputLBMNodeVariableInfo& output_info_tmp =
                *dynamic_cast<OutputLBMNodeVariableInfo*>(output_variables_.back().get());
            output_info_tmp.variable_dims_ = 1;
            output_info_tmp.func_get_var_ = [] (GridNodeLbm* const ptr_node)->std::vector<DefReal> {
                return {ptr_node->rho_};
            };
            output_info_tmp.output_name_ = "rho";
        } else if (iter_str == "velocity") {
            output_variables_.emplace_back(std::make_unique<OutputLBMNodeVariableInfo>());
            OutputLBMNodeVariableInfo& output_info_tmp =
                *dynamic_cast<OutputLBMNodeVariableInfo*>(output_variables_.back().get());
            output_info_tmp.variable_dims_ = ptr_lbm_solver->GetSolverDim();
            output_info_tmp.func_get_var_ = [ptr_lbm_solver, dt, func_macro]
                (GridNodeLbm* const ptr_node)->std::vector<DefReal> {
                std::vector<DefReal> force(ptr_lbm_solver->GetAllForcesForANode(*ptr_node));
                func_macro(dt, ptr_node->f_, force, &ptr_node->rho_, &ptr_node->velocity_);
                return ptr_node->velocity_;
            };
            output_info_tmp.output_name_ = "velocity";
        }
    }
    GetPtrToLbmGrid();
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
    std::string str_format, str_tmp;
    DefInt dims = output_info.variable_dims_;
    if (bool_binary) {
        str_format =  "binary";
    } else {
        str_format =  "ascii";
    }
    str_tmp.assign("      <DataArray NumberOfComponents=\"" + std::to_string(dims)
        + "\" type=\"" + output_data_format.output_real_.format_name_
        + "\" Name=\"" + output_info.output_name_ + "\" format=\"" + str_format + "\">\n");
    fprintf(fp, "%s", str_tmp.c_str());

    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            const std::vector<DefReal>& vec_var = output_info.func_get_var_(map_grid_node.at(iter->first).get());
            for (DefInt i_dims = 0;  i_dims < dims; ++i_dims) {
                base64_instance.AddToVecChar(output_data_format.output_real_.CastType(vec_var.at(i_dims)), &vec_uint8);
            }
        }
        base64_instance.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf(fp, "%c", iter);
        }
        fprintf(fp, "\n");
    } else {
        std::string str_format = "     "
            + output_data_format.output_real_.printf_format_;
        for (auto iter = map_node_index.begin(); iter != map_node_index.end(); ++iter) {
            const std::vector<DefReal>& vec_var =  output_info.func_get_var_(map_grid_node.at(iter->first).get());
            for (DefInt i_dims = 0; i_dims < dims; ++i_dims) {
                fprintf(fp, "  ");
                fprintf(fp, str_format.c_str(), vec_var.at(i_dims));
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "      </DataArray>\n");
    return 0;
}
/**
 * @brief function to calculate the size of the grid node information needed for checkpoint IO.
 * @return size of the grid node information for checkpoint IO.
 */
int GridInfoLbmInteface::GetSizeOfGridNodeInfoForCheckPoint() const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    int size_info = sizeof(DefInt) + static_cast<int>(ptr_lbm_solver->k0NumQ_*sizeof(DefReal));
    if (ptr_lbm_solver->bool_forces_) {
        size_info += static_cast<int>(ptr_lbm_solver->GetNumForces()*sizeof(DefReal));
    }
    return size_info;
}
}  // end namespace lbmproject
}  // end namespace rootproject
