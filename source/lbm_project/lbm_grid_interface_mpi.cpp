//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_grid_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM grid interface.
* @date  2023-9-30
*/
#include <mpi.h>
#include "./lbm_interface.h"
#include "io/log_write.h"
#include "io/vtk_writer.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to check domain boundary is periodic.
 * @param[in] dims dimensions of grid.
 * @param[out] ptr_periodic_min pointer to vector indicating periodic boundaries at minimum.
 * @param[out] ptr_periodic_max pointer to vector indicating periodic boundaries at maximum.
 * @return if true, at least one of the domain boundaries is periodic.
 * @note initially designed for finding periodic boundary in order to setup communication layers for them.
 */
bool GridInfoLbmInteface::CheckIfPeriodicDomainRequired(const DefAmrIndexUint dims,
    std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const {
    ptr_periodic_min->assign(dims, false);
    ptr_periodic_max->assign(dims, false);
    bool bool_has_periodic = false;
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXNeg)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXNeg)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_min->at(kXIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYNeg)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYNeg)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_min->at(kYIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXPos)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXPos)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_max->at(kXIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYPos)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYPos)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_max->at(kYIndex) = true;
        bool_has_periodic = true;
    }
    if (dims == 3) {
        if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZNeg)
            != domain_boundary_condition_.end()
            && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZNeg)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic) {
            ptr_periodic_min->at(kZIndex) = true;
            bool_has_periodic = true;
        }
        if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZPos)
            != domain_boundary_condition_.end()
            && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZPos)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic) {
            ptr_periodic_max->at(kZIndex) = true;
            bool_has_periodic = true;
        }
    }
    return bool_has_periodic;
}
/**
 * @brief function to calculate the size of the grid node information needed for MPI communication.
 * @return size of the grid node information for MPI communication.
 */
int GridInfoLbmInteface::SizeOfGridNodeInfoForMpiCommunication() const {
    int size_info = k0SizeOfAllDistributionFunctions_;
    if (bool_forces_) {
        size_info += ptr_solver_->k0SolverDims_*sizeof(DefReal);
    }
    return size_info;
}
/**
 * @brief function to compute information of node in inner mpi layer before communication.
 * @param[in] sfbitset_in space filling code of the input node.
 * @param[in] lbm_solver class to manage LBM solver.
 */
void GridInfoLbmInteface::ComputeNodeInfoBeforeMpiCommunication(
    const DefSFBitset sfbitset_in, const SolverLbmInterface& lbm_solver) {
    lbm_solver.StreamInForAGivenNode(sfbitset_in, *ptr_sfbitset_aux_, ptr_lbm_grid_nodes_);
}
/**
 * @brief function to copy node information to a buffer.
 * @param map_nodes container storing space filling codes of the nodes need to be copied.
 * @param ptr_buffer pointer to the buffer storing node information.
 * @note this function is a non-constant function since some information on nodes will be calculated before sending.
 */
void GridInfoLbmInteface::CopyNodeInfoToBuffer(
    const DefMap<DefAmrIndexUint>& map_nodes, char* const ptr_buffer) {
    GetPointerToLbmGrid();
    if (ptr_lbm_grid_nodes_ == nullptr) {
        std::string msg = "pointer to lbm grid nodes is null in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__);
        amrproject::LogManager::LogError(msg);
    } else  {
        DefSizet position = 0;
        SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
        for (const auto& iter : map_nodes) {
            if (ptr_lbm_grid_nodes_->find(iter.first) != ptr_lbm_grid_nodes_->end()) {
                // calculate distribution functions before copy
                ComputeNodeInfoBeforeMpiCommunication(iter.first, *ptr_lbm_solver);

                std::memcpy(ptr_buffer + position, &(iter.first), sizeof(DefSFBitset));
                position+=sizeof(DefSFBitset);
                std::memcpy(ptr_buffer + position, ptr_lbm_grid_nodes_->at(iter.first)->f_.data(),
                    k0SizeOfAllDistributionFunctions_);
                position+=k0SizeOfAllDistributionFunctions_;
                if (bool_forces_) {
                    int force_size = ptr_solver_->k0SolverDims_*sizeof(DefReal);
                    std::memcpy(ptr_buffer + position, ptr_lbm_grid_nodes_->at(iter.first)->force_.data(), force_size);
                    position += force_size;
                }
            } else {
                std::string msg = "grid node does not exist for copying to a buffer in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__);
                amrproject::LogManager::LogError(msg);
            }
        }
    }
}
/**
 * @brief function to compute macroscopic variables based on distribution functions in the last time step.
 * @param map_inner_nodes container storing space filling codes of inner mpi communication layers will be sent to other ranks.
 * @param map_outer_nodes container storing space filling codes of outer mpi communication layer of the current rank.
 */
void GridInfoLbmInteface::ComputeLocalInfoOnMpiLayers(
    const std::map<int, DefMap<DefAmrIndexUint>>& map_inner_nodes,
    const DefMap<DefAmrIndexUint>& map_outer_nodes) {
    DefReal dt_lbm = ptr_collision_operator_->dt_lbm_;
    std::function<void(const DefReal, GridNodeLbm* const)> func_macro;
    SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    GetPointerToLbmGrid();
    if (bool_forces_) {
        if (ptr_lbm_grid_nodes_->begin()->second->force_.size() != ptr_solver_->k0SolverDims_) {
            amrproject::LogManager::LogError("Size of forces should be "
                + std::to_string(ptr_solver_->k0SolverDims_)
                + "rather than " + std::to_string(ptr_lbm_grid_nodes_->begin()->second->force_.size()) + " in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        func_macro = ptr_lbm_solver->func_macro_with_force_;
    } else {
        func_macro = ptr_lbm_solver->func_macro_without_force_;
    }
    DefAmrIndexUint flag_not_compute =
        amrproject::NodeBitStatus::kNodeStatusFine2Coarse0_|amrproject::NodeBitStatus::kNodeStatusFine2CoarseM1_;
    // collision for nodes in outer and inner MPI communication layers
    for (const auto& iter_node : map_outer_nodes) {
        if (!(ptr_lbm_grid_nodes_->at(iter_node.first)->flag_status_&flag_not_compute)) {
            func_macro(dt_lbm, ptr_lbm_grid_nodes_->at(iter_node.first).get());
            ptr_collision_operator_->CollisionOperator(*ptr_lbm_solver,
                ptr_lbm_grid_nodes_->at(iter_node.first).get());
        }
    }
    for (const auto& iter_layer : map_inner_nodes) {
        for (const auto& iter_node : iter_layer.second) {
            if (!(ptr_lbm_grid_nodes_->at(iter_node.first)->flag_status_&flag_not_compute)) {
                func_macro(dt_lbm, ptr_lbm_grid_nodes_->at(iter_node.first).get());
                ptr_collision_operator_->CollisionOperator(*ptr_lbm_solver,
                    ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
        }
    }
}
/**
 * @brief function to compute information of node in inner mpi layer after communication.
 * @param[in] sfbitset_in space filling code of the input node.
 * @param[in] lbm_solver class to manage LBM solver.
 */
void GridInfoLbmInteface::ComputeNodeInfoAfterMpiCommunication(
    const DefSFBitset sfbitset_in, const SolverLbmInterface& lbm_solver) {
    lbm_solver.StreamOutForAGivenNode(sfbitset_in, *ptr_sfbitset_aux_, ptr_lbm_grid_nodes_);
}
/**
 * @brief function to read node information from a buffer consisting all chunks.
 * @param buffer_size total size of the buffer.
 * @param buffer  pointer to the buffer storing node information.
 */
void GridInfoLbmInteface::ReadNodeInfoFromBuffer(
    const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    int node_info_size = SizeOfGridNodeInfoForMpiCommunication();
    DefSizet num_nodes = buffer_size/(key_size + node_info_size);
    int force_size = ptr_solver_->k0SolverDims_*sizeof(DefReal);
    // deserialize data stored in buffer
    DefSizet position = 0;
    DefSFBitset key_code;
    SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();

    for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        if (ptr_lbm_grid_nodes_->find(key_code) != ptr_lbm_grid_nodes_->end()) {
            std::memcpy(ptr_lbm_grid_nodes_->at(key_code)->f_.data(),
                ptr_buffer + position, k0SizeOfAllDistributionFunctions_);
            position += k0SizeOfAllDistributionFunctions_;
            if (bool_forces_) {
                std::memcpy(ptr_lbm_grid_nodes_->at(key_code)->force_.data(),
                    ptr_buffer + position, force_size);
                 position += force_size;
            }
             int i_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
        } else {
            std::string msg = "grid node does not exist for copying from a buffer in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__);
            amrproject::LogManager::LogError(msg);
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject