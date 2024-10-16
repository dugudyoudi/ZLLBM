//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_info_mpi.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid informaiton for mpi communication.
* @date  2022-6-7
*/
#include <mpi.h>

#include <string>
#include "grid/grid_info_interface.h"
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to set mpi communication region for periodic boundary conditions.
 * @param[in] dims  dimensions of the grid.
 * @param[in] sfbitset_aux class to manage space filling codes.
 * @param[in] bool_periodic_min booleans indicating if the boundary is periodic at minimum domain boundaries.
 * @param[in] bool_periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[out] ptr_partition_interface pointer to container storing nodes on partition information.
 */
void GridInfoInterface::SetPeriodicBoundaryAsPartitionInterface(DefInt dims,
    const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
    DefMap<DefInt>* const ptr_partition_interface) {
    DefInt flag0 = ~0;
    DefInt i_level = i_level_;
    if (bool_periodic_min.size() != dims) {
        LogManager::LogError("size of bool_periodic_min " + std::to_string(bool_periodic_min.size())
            + " is not equal to the grid dimension " + std::to_string(dims) + "in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (bool_periodic_max.size() != dims) {
        LogManager::LogError("size of bool_periodic_max " + std::to_string(bool_periodic_min.size())
            + " is not equal to the grid dimension " + std::to_string(dims) + "in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    for (DefInt i_dims = 0; i_dims < dims; ++i_dims) {
        if (bool_periodic_min.at(i_dims)) {
            ptr_partition_interface->insert(
                domain_boundary_min_.at(i_dims).begin(), domain_boundary_min_.at(i_dims).end());
        }
        if (bool_periodic_max.at(i_dims)) {
            ptr_partition_interface->insert(
                domain_boundary_max_.at(i_dims).begin(), domain_boundary_max_.at(i_dims).end());
        }
    }
}
/**
 * @brief function to add nodes not exist on current rank but used in interpolation.
 * @param[in] bool_periodic_min booleans indicating if the boundary is periodic at minimum domain boundaries.
 * @param[in] bool_periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] sfbitset_aux class to manage space filling codes.
 * @param[in] refinement_interface nodes on the refinement interface.
 * @param[in] map_nodes_lower grid nodes at one lower refinement level.
 */
int GridInfoInterface::AddGhostNodesForInterpolation(const std::vector<bool>& bool_periodic_min,
    const std::vector<bool>& bool_periodic_max, const SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<DefInt>& refinement_interface, const DefMap<std::unique_ptr<GridNode>>& map_nodes_lower) {
    DefSFBitset sfbitset_lower, sfbitset_current;
    std::vector<DefSFBitset> vec_in_region;
    size_t dims = k0VecBitsetDomainMin_.size();
    std::vector<DefSFBitset> domain_min(dims), domain_max(dims);
    for (DefInt i_dims = 0; i_dims < dims; ++i_dims) {
        domain_min.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(1, k0VecBitsetDomainMin_.at(i_dims));
        domain_max.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(1, k0VecBitsetDomainMax_.at(i_dims));
    }
    for (auto& iter_node : refinement_interface) {
        sfbitset_lower = sfbitset_aux.SFBitsetToNLowerLevelVir(1, iter_node.first);
        sfbitset_aux.FindNodesInPeriodicRegionCorner(sfbitset_lower,
            max_interp_length_, bool_periodic_min, bool_periodic_max,
            domain_min, domain_max, &vec_in_region);
        for (const auto& iter_node_region : vec_in_region) {
            if (iter_node_region != SFBitsetAuxInterface::kInvalidSFbitset
                && (map_nodes_lower.find(iter_node_region) == map_nodes_lower.end()
                || ((map_nodes_lower.at(iter_node_region)->flag_status_ & NodeBitStatus::kNodeStatusCoarse2FineGhost_))
                == NodeBitStatus::kNodeStatusCoarse2FineGhost_)) {
                sfbitset_current = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_node_region);
                if (map_grid_node_.find(sfbitset_current) == map_grid_node_.end()
                    ||((map_grid_node_.at(sfbitset_current)->flag_status_ & NodeBitStatus::kNodeStatusFine2CoarseGhost_)
                    == NodeBitStatus::kNodeStatusFine2CoarseGhost_)) {
                    interp_nodes_outer_layer_.insert({sfbitset_current, GridNodeCreator()});
                }
            }
        }
    }
    return 0;
}

void GridInfoInterface::RemoveUnnecessaryF2CNodesOnMpiOuterLayer(const DefSFCodeToUint code_min,
    const DefSFCodeToUint code_max, const DefInt num_outer_layer, DefMap<DefInt>* const ptr_mpi_outer_layer) {
    std::vector<DefSFBitset> nodes_in_region;
    DefSFCodeToUint code_current;
    SFBitsetAuxInterface& sfbitset_aux = *ptr_sfbitset_aux_;
    DefInt dims = GetPtrToParentGridManager()->k0GridDims_;
    std::vector<bool> periodic_min(dims, false), periodic_max(dims, false);
    CheckIfPeriodicDomainRequired(dims, &periodic_min, &periodic_max);
    std::vector<DefAmrLUint> indices_min = GetPtrToParentGridManager()->GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = GetPtrToParentGridManager()->GetMaxIndexOfBackgroundNodeArrAsVec();
    std::vector<DefSFBitset> domain_min_n_level(dims), domain_max_n_level(dims);
    sfbitset_aux.GetMinAtGivenLevel(i_level_, indices_min, &domain_min_n_level);
    sfbitset_aux.GetMaxAtGivenLevel(i_level_, indices_max, &domain_max_n_level);
    for (const auto& iter_interface : map_ptr_interface_layer_info_) {
        for (auto& iter_layer : iter_interface.second->vec_inner_fine2coarse_) {
            std::vector<DefSFBitset> vec_node_to_remove;
            for (const auto& iter_node : iter_layer) {
                code_current = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level_, iter_node.first).to_ullong();
                if (code_current > code_max || code_current < code_min) {
                    bool bool_remove = true;
                    sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_node.first, num_outer_layer,
                        periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &nodes_in_region);
                    for (auto& iter_node_region : nodes_in_region) {
                        if (map_grid_node_.find(iter_node_region) != map_grid_node_.end()
                            && (map_grid_node_.at(iter_node_region)->flag_status_
                            & NodeBitStatus::kNodeStatusMpiPartitionInner_)) {
                            bool_remove = false;
                            break;
                        }
                    }
                    if (bool_remove) {
                        vec_node_to_remove.push_back(iter_node.first);
                    }
                }
            }
            for (const auto& iter_node : vec_node_to_remove) {
                iter_layer.erase(iter_node);
                if (ptr_mpi_outer_layer->find(iter_node)!= ptr_mpi_outer_layer->end()) {
                    ptr_mpi_outer_layer->erase(iter_node);
                }
                if (map_grid_node_.find(iter_node) != map_grid_node_.end()) {
                    map_grid_node_.erase(iter_node);
                }
            }
        }
        for (auto& iter_layer : iter_interface.second->vec_outer_fine2coarse_) {
            for (const auto& iter_node : iter_layer) {
                std::vector<DefSFBitset> vec_node_to_remove;
                code_current = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level_, iter_node.first).to_ullong();
                if (code_current > code_max || code_current < code_min) {
                    bool bool_remove = true;
                    sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_node.first, num_outer_layer,
                        periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &nodes_in_region);
                    for (auto& iter_node_region : nodes_in_region) {
                        if (map_grid_node_.find(iter_node_region) != map_grid_node_.end()
                            && (map_grid_node_.at(iter_node_region)->flag_status_
                            & NodeBitStatus::kNodeStatusMpiPartitionInner_)) {
                            bool_remove = false;
                            break;
                        }
                    }
                    if (bool_remove) {
                        vec_node_to_remove.push_back(iter_node.first);
                    }
                }
                for (const auto& iter_node : vec_node_to_remove) {
                    iter_layer.erase(iter_node);
                    if (ptr_mpi_outer_layer->find(iter_node)!= ptr_mpi_outer_layer->end()) {
                        ptr_mpi_outer_layer->erase(iter_node);
                    }
                    if (map_grid_node_.find(iter_node) != map_grid_node_.end()) {
                        map_grid_node_.erase(iter_node);
                    }
                }
            }
        }
    }
}
/**
 * @brief function to copy node information to a buffer.
 * @param[in] map_nodes container storing space filling codes of the nodes need to be copied.
 * @param[out] ptr_buffer pointer to the buffer storing node information.
 * @note this function is a non-constant function since some information on nodes will be calculated before sending.
 */
int GridInfoInterface::CopyNodeInfoToBuffer(
    const DefMap<DefInt>& map_nodes, char* const ptr_buffer) const {
    DefSizet position = 0;
    int node_info_size = SizeOfGridNodeInfoForMpiCommunication();
    DefSizet buffer_size = (node_info_size + sizeof(DefSFBitset)) * map_nodes.size();
    for (const auto& iter : map_nodes) {
        if (map_grid_node_.find(iter.first) != map_grid_node_.end()) {
            std::memcpy(ptr_buffer + position, &(iter.first), sizeof(DefSFBitset));
            position+=sizeof(DefSFBitset);
            map_grid_node_.at(iter.first)->CopyNodInfoToBuffer(ptr_buffer + position);
            position+=node_info_size;
            if (position > buffer_size) {
                LogManager::LogError("Buffer to store node information overflows (buffer size is "
                    + std::to_string(buffer_size) + " ), please check"
                    " size of node info for mpi communication");
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
            LogManager::LogError(msg);
        }
    }
    return 0;
}
/**
 * @brief function to read node information from a buffer consisting all chunks.
 * @param buffer_size total size of the buffer.
 * @param buffer  pointer to the buffer storing node information.
 */
void GridInfoInterface::ReadNodeInfoFromBuffer(
    const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    int node_info_size = SizeOfGridNodeInfoForMpiCommunication();
    DefSizet num_nodes = buffer_size/(key_size + node_info_size);
    // deserialize data stored in buffer
    DefSizet position = 0;
    DefSFBitset key_code;
    for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        if (map_grid_node_.find(key_code) != map_grid_node_.end()) {
            map_grid_node_.at(key_code)->ReadNodInfoFromBuffer(ptr_buffer + position);
            position+=node_info_size;
            if (position > buffer_size) {
                LogManager::LogError("Buffer to store node information overflows, please check"
                    " size of node info for mpi communication");
            }
        }   // may receive nodes that do not exist in current rank on c2f interface
    }
}
/**
 * @brief function to copy information of node needed for interpolation to a buffer.
 * @param[in] coarse_grid_info class storting grid information at lower level.
 * @param[in] map_nodes container storing space filling codes of the nodes need to be copied.
 * @param[out] ptr_buffer pointer to the buffer storing node information.
 */
int GridInfoInterface::CopyNodeVariablesToBuffer(
    const std::function<void(const GridNode& node_ref, char* const)>& func_copy_buffer,
    const DefMap<DefInt>& map_nodes, char* const ptr_buffer) {
    DefSizet position = 0;
    DefSizet var_size = GetVariablesSizeForCopy();
    for (const auto& iter : map_nodes) {
        if (map_grid_node_.find(iter.first) != map_grid_node_.end()) {
            std::memcpy(ptr_buffer + position, &(iter.first), sizeof(DefSFBitset));
            position+=sizeof(DefSFBitset);
            func_copy_buffer(*map_grid_node_.at(iter.first).get(), ptr_buffer+position);
            position+=var_size;
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
    return 0;
}
int GridInfoInterface::ReadNodeVariablesFromBuffer(const DefSizet buffer_size,
        const std::function<void(const char* ptr_node_buffer, const GridNode* const ptr_node)> func_read_buffer,
        const std::unique_ptr<char[]>& buffer,  DefMap<std::unique_ptr<GridNode>>* const ptr_grid_nodes) {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    DefSizet var_size = GetVariablesSizeForCopy();
    if (buffer_size%(key_size + var_size) != 0) {
        LogManager::LogError("Buffer size is not divisible by size of variables,"
            " i.e. number of nodes if not an integer");
    }
    DefSizet num_nodes = buffer_size/(key_size + var_size);
    // deserialize data stored in buffer
    DefSizet position = 0;
    DefSFBitset key_code;
    for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        if (map_grid_node_.find(key_code) != map_grid_node_.end()) {
            func_read_buffer(ptr_buffer + position, map_grid_node_.at(key_code).get());
            position += var_size;
        } else {
            std::vector<DefReal> indices;
            ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(key_code, grid_space_, &indices);
            std::string msg;
            if (indices.size() == 2) {
                msg = "grid node (" + std::to_string(indices[kXIndex]) + ", " + std::to_string(indices[kYIndex])
                    + ") at " + std::to_string(i_level_) + " level does not exist for copying from a buffer";
            } else {
                msg = "grid node (" + std::to_string(indices[kXIndex]) + ", " + std::to_string(indices[kYIndex])
                    + ", " + std::to_string(indices[kZIndex]) + ") at " + std::to_string(i_level_)
                    + " level does not exist for copying from a buffer";
            }
            amrproject::LogManager::LogError(msg);
            return -1;
        }
    }
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject

