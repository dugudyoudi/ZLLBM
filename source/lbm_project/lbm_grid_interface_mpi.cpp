//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_grid_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM grid interface.
*/
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
bool GridInfoLbmInteface::CheckIfPeriodicDomainRequired(const DefInt dims,
    std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const {
    ptr_periodic_min->assign(dims, false);
    ptr_periodic_max->assign(dims, false);
    bool bool_has_periodic = false;
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXMin)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXMin)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_min->at(kXIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYMin)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYMin)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_min->at(kYIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXMax)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXMax)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_max->at(kXIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYMax)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYMax)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_max->at(kYIndex) = true;
        bool_has_periodic = true;
    }
    if (dims == 3) {
        if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZMin)
            != domain_boundary_condition_.end()
            && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZMin)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic) {
            ptr_periodic_min->at(kZIndex) = true;
            bool_has_periodic = true;
        }
        if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZMax)
            != domain_boundary_condition_.end()
            && domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZMax)->boundary_scheme_
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
    SolverLbmInterface& lbm_solver = *std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    if (lbm_solver.bool_forces_) {
        size_info += lbm_solver.GetNumForces()*sizeof(DefReal);
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
 * @brief function to copy information of node needed for interpolation to a buffer.
 * @param[in] coarse_grid_info class storting grid information at lower level.
 * @param[in] map_nodes container storing space filling codes of the nodes need to be copied.
 * @param[out] ptr_buffer pointer to the buffer storing node information.
 * @note node not exist at current level will be find in the lower level.
 */
int GridInfoLbmInteface::CopyInterpolationNodeInfoToBuffer(const GridInfoInterface& coarse_grid_info,
    const DefMap<DefInt>& map_nodes, char* const ptr_buffer) {
    GetPointerToLbmGrid();
    SolverLbmInterface& lbm_solver = *std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_);
    if (ptr_lbm_grid_nodes_ == nullptr) {
        amrproject::LogManager::LogError("pointer to lbm grid nodes is null");
    }
    const GridInfoLbmInteface& coarse_grid_info_lbm = dynamic_cast<const GridInfoLbmInteface&>(coarse_grid_info);
    if (coarse_grid_info_lbm.ptr_lbm_grid_nodes_ == nullptr) {
        amrproject::LogManager::LogError("pointer to lbm grid nodes at one lower refinement level is null");
    }
    DefSizet position = 0;
    const DefMap<std::unique_ptr<GridNodeLbm>>& lbm_coarse_grid_nodes = *coarse_grid_info_lbm.ptr_lbm_grid_nodes_;
    for (const auto& iter : map_nodes) {
        if (ptr_lbm_grid_nodes_->find(iter.first) != ptr_lbm_grid_nodes_->end()) {
            std::memcpy(ptr_buffer + position, &(iter.first), sizeof(DefSFBitset));
            position+=sizeof(DefSFBitset);
            std::memcpy(ptr_buffer + position, ptr_lbm_grid_nodes_->at(iter.first)->f_.data(),
                k0SizeOfAllDistributionFunctions_);
            position+=k0SizeOfAllDistributionFunctions_;
            if (lbm_solver.bool_forces_) {
                int force_size = lbm_solver.GetNumForces()*sizeof(DefReal);
                std::memcpy(ptr_buffer + position, ptr_lbm_grid_nodes_->at(iter.first)->force_.data(), force_size);
                position += force_size;
            }
        } else {
            DefSFBitset sfbitset_lower = ptr_sfbitset_aux_->SFBitsetToNLowerLevelVir(1, iter.first);
            if (lbm_coarse_grid_nodes.find(sfbitset_lower) != lbm_coarse_grid_nodes.end()) {
                GridNodeLbm node_coarse2fine;
                coarse_grid_info.NodeInfoCoarse2fine(
                    *lbm_coarse_grid_nodes.at(sfbitset_lower).get(), &node_coarse2fine);
                std::memcpy(ptr_buffer + position, &(iter.first), sizeof(DefSFBitset));
                position+=sizeof(DefSFBitset);

                std::memcpy(ptr_buffer + position, node_coarse2fine.f_.data(),
                    k0SizeOfAllDistributionFunctions_);

                position+=k0SizeOfAllDistributionFunctions_;
                if (lbm_solver.bool_forces_) {
                    int force_size = lbm_solver.GetNumForces()*sizeof(DefReal);
                    std::memcpy(ptr_buffer + position,
                        lbm_coarse_grid_nodes.at(sfbitset_lower)->force_.data(), force_size);
                    position += force_size;
                }
            } else {
                std::vector<DefReal> indices;
                ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(iter.first, grid_space_, &indices);
                std::string msg;
                if (indices.size() == 2) {
                    msg = "grid node (" + std::to_string(indices[kXIndex]) + ", " + std::to_string(indices[kYIndex])
                        + ") at " + std::to_string(i_level_) + " at level not exist for copying to a buffer";
                } else {
                    msg = "grid node (" + std::to_string(indices[kXIndex]) + ", " + std::to_string(indices[kYIndex])
                        + std::to_string(indices[kZIndex]) +  + ") at " + std::to_string(i_level_)
                        + " level does not exist for copying to a buffer";
                }
                amrproject::LogManager::LogError(msg);
                return -1;
            }
        }
    }
    return 0;
}
/**
 * @brief function to compute macroscopic variables based on distribution functions in the last time step.
 * @param map_inner_nodes container storing space filling codes of inner mpi communication layers will be sent to other ranks.
 * @param map_outer_nodes container storing space filling codes of outer mpi communication layer of the current rank.
 */
void GridInfoLbmInteface::ComputeInfoInMpiLayers(
    const std::map<int, DefMap<DefInt>>& map_inner_nodes,
    const DefMap<DefInt>& map_outer_nodes) {
    DefReal dt_lbm = ptr_collision_operator_->dt_lbm_;
    DefInt dims = ptr_solver_->GetSolverDim();
    std::function<void(const DefReal, const GridNodeLbm&, const std::vector<DefReal>&,
        DefReal* const, std::vector<DefReal>* const)> func_macro;
    SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_grid_nodes = GetPointerToLbmGrid();
    if (ptr_lbm_solver->bool_forces_) {
        func_macro = ptr_lbm_solver->func_macro_with_force_;
    } else {
        func_macro = [ptr_lbm_solver](const DefReal dt_lbm, const GridNodeLbm& node,
            const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            ptr_lbm_solver->func_macro_without_force_(node, ptr_rho, ptr_velocity);};
    }

    // collision for nodes in outer and inner MPI communication layers
    DefInt flag_not_collide = NodeFlagNotCollision_&(~(amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_));
    ptr_lbm_solver->CollisionForGivenNodes<DefInt>(flag_not_collide, map_outer_nodes, this);
    DefMap<DefInt> map_one_layer_near_inner;
    std::vector<DefSFBitset> vec_neighbor;
    for (const auto& iter_layer : map_inner_nodes) {
        for (const auto& iter_node : iter_layer.second) {
            if (!(ptr_lbm_grid_nodes_->at(iter_node.first)->flag_status_&flag_not_collide)) {
                const std::vector<DefReal>& force
                    = ptr_lbm_solver->GetAllForcesForANode(*ptr_lbm_grid_nodes_->at(iter_node.first).get());
                func_macro(dt_lbm, *ptr_lbm_grid_nodes_->at(iter_node.first).get(), force,
                    &ptr_lbm_grid_nodes_->at(iter_node.first)->rho_,
                    &ptr_lbm_grid_nodes_->at(iter_node.first)->velocity_);
                ptr_collision_operator_->CollisionOperator(*ptr_lbm_solver, force,
                    ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
            ptr_sfbitset_aux_->SFBitsetFindAllNeighborsVir(iter_node.first, &vec_neighbor);
            for (const auto& iter_neighbor : vec_neighbor) {
                if (map_one_layer_near_inner.find(iter_neighbor) == map_one_layer_near_inner.end()) {
                    map_one_layer_near_inner.insert({iter_neighbor, 0});
                }
            }
        }
    }
    // Since stream step will be performed before mpi communication, post-collision distribution
    // functions of neighboring nodes are need
    ptr_lbm_solver->CollisionForGivenNodes<DefInt>(flag_not_collide, map_one_layer_near_inner, this);

    DefInt flag_not_stream = NodeFlagNotStream_&(~(amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_));
    for (const auto& iter_layer : map_inner_nodes) {
        ptr_lbm_solver->StreamInForGivenNodes<DefInt>(
            flag_not_stream, *ptr_sfbitset_aux_, iter_layer.second, this);
    }
}
/**
 * @brief function to compute information of node in inner mpi layer after communication.
 * @param[in] sfbitset_in space filling code of the input node.
 * @param[in] lbm_solver class to manage LBM solver.
 */
void GridInfoLbmInteface::ComputeNodeInfoAfterMpiCommunication(
    const DefSFBitset sfbitset_in, const SolverLbmInterface& lbm_solver) {
}
/**
 * @brief function to read node information from a buffer consisting all chunks.
 * @param buffer_size total size of the buffer.
 * @param buffer  pointer to the buffer storing node information.
 * @return 0 run successfully, -1 failure since at least one node is not found. 
 */
int GridInfoLbmInteface::ReadInterpolationNodeInfoFromBuffer(
    const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    int node_info_size = SizeOfGridNodeInfoForMpiCommunication();
    DefSizet num_nodes = buffer_size/(key_size + node_info_size);
    SolverLbmInterface& lbm_solver = *std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_);
    int force_size = lbm_solver.GetNumForces()*sizeof(DefReal);
    // deserialize data stored in buffer
    DefSizet position = 0;
    DefSFBitset key_code;
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_interp_nodes =
        reinterpret_cast<DefMap<std::unique_ptr<GridNodeLbm>>*>(&interp_nodes_outer_layer_);
    for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        if (ptr_lbm_interp_nodes->find(key_code) != ptr_lbm_interp_nodes->end()) {
            std::memcpy(ptr_lbm_interp_nodes->at(key_code)->f_.data(),
                ptr_buffer + position, k0SizeOfAllDistributionFunctions_);
            position += k0SizeOfAllDistributionFunctions_;
            if (lbm_solver.bool_forces_) {
                std::memcpy(ptr_lbm_grid_nodes_->at(key_code)->force_.data(),
                    ptr_buffer + position, force_size);
                 position += force_size;
            }
        } else {
            std::vector<DefReal> coordinates;
            std::vector<DefReal> grid_spacing_background = ptr_sfbitset_aux_->GetBackgroundGridSpacing();
            ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(key_code, grid_space_, &coordinates);
            std::string msg;
            if (coordinates.size() == 2) {
                msg = "grid node (" + std::to_string(coordinates[kXIndex]) + ", " + std::to_string(coordinates[kYIndex])
                    + ") at " + std::to_string(i_level_) + " level does not exist for copying from a buffer";
            } else {
                msg = "grid node (" + std::to_string(coordinates[kXIndex]) + ", " + std::to_string(coordinates[kYIndex])
                    + ", " + std::to_string(coordinates[kZIndex]) + ") at " + std::to_string(i_level_)
                    + " level does not exist for copying from a buffer";
            }
            amrproject::LogManager::LogError(msg);
            return -1;
        }
    }
    return 0;
}
}  // end namespace lbmproject
}  // end namespace rootproject
