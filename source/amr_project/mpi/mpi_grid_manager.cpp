//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is enabled.
* @date  2023-5-2
*/
#include <limits>
#include "mpi/mpi_manager.h"
#include "grid/grid_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to serialize data (a DefSFBitset and a DefAmrUint) and save it to buffer
 * @param[in] map_nodes node information
 * @param[out] ptr_buffer_size pointer to size of the buffer in bytes.
 * @return sunique pointer to a char array to store the serialized data.
 */
std::unique_ptr<char[]> MpiManager::SerializeNodeStoreUint(const DefMap<DefAmrUint>& map_nodes,
    int* const ptr_buffer_size) const {
    int key_size = sizeof(DefSFBitset), node_size = sizeof(DefAmrUint);
    int num_nodes = 1;
    if  (sizeof(int) + map_nodes.size()*(key_size + node_size) > (std::numeric_limits<int>::max)()) {
        LogManager::LogError("size of the buffer is greater"
         " than the maximum of int in MpiManager::SerializeData(DefMap<DefAmrUint>)");
    } else {
        num_nodes = static_cast<int>(map_nodes.size());
    }
    int& buffer_size = *ptr_buffer_size;
    buffer_size = sizeof(int) + (key_size + node_size) * num_nodes;
    // allocation buffer to store the serialized data
    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_nodes, sizeof(int));
    position += sizeof(int);
    // serialize data stored in nodes
    DefSFCodeToUint key_host;
    for (auto& iter : map_nodes) {
        // convert bitset into unsigned long long and save it to buffer
        key_host = iter.first.to_ullong();
        std::memcpy(ptr_buffer + position, &key_host, key_size);
        position += key_size;
        // save node data to buffer
        std::memcpy(ptr_buffer + position, &iter.second, node_size);
        position += node_size;
    }
    return buffer;
}
/**
 * @brief function to deserialize data (a DefSFBitset and a DefAmrUint) from a buffer.
 * @param[in] buffer buffer has received data.
 * @param[out] map_nodes node information.
 */
void MpiManager::DeserializeNodeStoreUint(
    const std::unique_ptr<char[]>& buffer, DefMap<DefAmrUint>* const map_nodes) const {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset), node_size = sizeof(DefAmrUint);
    // number of nodes
    int num_nodes;
    int position = 0;
    std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
    DefAmrUint node_data;
    // deserialize data stored in buffer
    position += sizeof(int);
    DefSFCodeToUint key_code;
    DefSFBitset key_host;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        key_host = static_cast<DefSFBitset>(key_code);
        position += key_size;
        std::memcpy(&node_data, ptr_buffer + position, node_size);
        position += node_size;
        map_nodes->insert({ key_host, node_data });
    }
}
/**
 * @brief function to serialize space filling code (DefSFBitset) and save it to buffer
 * @param[in] map_nodes node information
 * @param[out] ptr_buffer_size pointer to size of the buffer in bytes.
 * @return unique pointer to a char array to store the serialized data.
 */
std::unique_ptr<char[]> MpiManager::SerializeNodeSFBitset(
    const DefMap<DefAmrIndexUint>& map_nodes, int* const ptr_buffer_size) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes = 1;
    if  (sizeof(int) + map_nodes.size() *(key_size) > (std::numeric_limits<int>::max)()) {
        LogManager::LogError("size of the buffer is greater than the"
         " maximum of int in MpiManager::SerializeData(DefMap<DefAmrUint>)");
    } else {
        num_nodes = static_cast<int>(map_nodes.size());
    }
    int& buffer_size = *ptr_buffer_size;
    buffer_size = sizeof(int) + (key_size) * num_nodes;
    // allocation buffer to store the serialized data
    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_nodes, sizeof(int));
    position += sizeof(int);
    // serialize data stored in nodes
    DefSFCodeToUint key_host;
    for (auto& iter : map_nodes) {
        // convert bitset into unsigned long long and save it to buffer
        key_host = iter.first.to_ullong();
        std::memcpy(ptr_buffer + position, &key_host, key_size);
        position += key_size;
    }
    return buffer;
}
/**
 * @brief function to deserialize space filling code (DefSFBitset) from a buffer
 * @param[in] flag_node a given value assign to the value corresponding to space filling code
 * @param[in] buffer buffer has received data
 * @param[out] map_nodes node information
 */
void MpiManager::DeserializeNodeSFBitset(const DefAmrIndexUint flag_node,
    const std::unique_ptr<char[]>& buffer, DefMap<DefAmrIndexUint>* const map_nodes) const {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    // number of nodes
    int num_nodes;
    int position = 0;
    std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
    // deserialize data stored in buffer
    position += sizeof(int);
    DefSFCodeToUint key_code;
    DefSFBitset key_host;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        key_host = static_cast<DefSFBitset>(key_code);
        position += key_size;
        map_nodes->insert({ key_host, flag_node });
    }
}
/**
 * @brief function to send and receive partitioned grid for each rank.
 * @param[in] flag_size0 flag of existing node.
 * @param[in] flag_fine2coarse_outmost flag indicting node is on the outmost fine to coarse refinement layer.
 * @param[in] bitset_min minimum space filling code for each partition.
 * @param[in] bitset_max maximum space filling code for each partition.
 * @param[in] indices_min minimum indices of the computational domain.
 * @param[in] indices_max maximum indices of the computational domain.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[in] vec_sfbitset_rank0  space-filling codes of nodes at all refinement levels on rank 0.
 * @param[out] ptr_sfbitset_each nodes on current rank at all levels.
 * @param[out] ptr_sfbitset_ghost_each nodes near refinement interface in a given number of extended mpi communication layers.
 * @param[out] ptr_vec_grid_info grid information at all ranks.
 * @note sent nodes are whose space filling code is >= bitset_min and <= bitset_max,
 *       as well as those in outer mpi communication layers.
 *       Nodes coincide with those in ptr_sfbitset_each are not included in ptr_sfbitset_ghost_each
 */
void MpiManager::IniSendNReceivePartitionedGrid(const DefAmrIndexUint dims,
    const DefAmrIndexUint flag_size0, const DefAmrUint flag_fine2coarse_outmost,
    const std::vector<DefSFBitset>& bitset_min, const std::vector<DefSFBitset>& bitset_max,
    const std::vector<DefAmrIndexLUint>& indices_min, const std::vector<DefAmrIndexLUint>& indices_max,
    const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset_rank0,
    std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_each,
    std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_ghost_each,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const {
    if (vec_sfbitset_rank0.size() - 1 > INT_MAX) {
        LogManager::LogError("size of vector vec_sfbitset_rank0 is larger than INT_MAX in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    int max_level = static_cast<int>(vec_sfbitset_rank0.size() - 1);
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_min(bitset_min.size()), ull_max(bitset_max.size());
    std::array<DefAmrIndexLUint, 2> indices_min_2d, indices_max_2d;
    if (dims == 2) {
        if (indices_min.size() != 2) {
            LogManager::LogError("size of vector indices_min shoule be 2 in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        if (indices_max.size() != 2) {
            LogManager::LogError("size of vector indices_max shoule be 2 in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        indices_min_2d = { indices_min.at(kXIndex), indices_min.at(kYIndex) };
        indices_max_2d = { indices_max.at(kXIndex), indices_max.at(kYIndex) };
    }
    std::array<DefAmrIndexLUint, 3> indices_min_3d, indices_max_3d;
    if (dims == 3) {
        if (indices_min.size() != 3) {
            LogManager::LogError("size of vector indices_min shoule be 3 in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        if (indices_max.size() != 3) {
            LogManager::LogError("size of vector indices_max shoule be 3 in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        indices_min_3d = { indices_min.at(kXIndex), indices_min.at(kYIndex), indices_min.at(kZIndex)};
        indices_max_3d = { indices_max.at(kXIndex), indices_max.at(kYIndex), indices_max.at(kZIndex)};
    }

    if (rank_id == 0) {
        if (vec_sfbitset_rank0.size() != ptr_sfbitset_each->size()) {
            LogManager::LogError("size of the input vector (vec_sfbitset_rank0) is different from the"
                " size of the output vector (ptr_sfbitset_each) in MpiManager::IniSendNReceivePartitionedGrid in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        for (DefSizet i = 0; i < bitset_max.size(); ++i) {
            ull_min.at(i) = bitset_min.at(i).to_ullong();
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefSFBitset) - 1;
    DefSFCodeToUint background_code;
    int int_interface_status;
    DefAmrIndexUint num_ghost_lower = k0NumPartitionOuterLayers_/2,
     num_ghost_upper = (k0NumPartitionOuterLayers_ + 1)/2;
    DefSizet num_max = ull_max.size();
    std::vector<DefMap<DefAmrIndexUint>> partition_interface_background(num_ranks);
    if (rank_id == 0) {  // partition nodes on rank 0
        for (auto i_rank = 0; i_rank < num_ranks; ++i_rank) {
            if (dims == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
                IniFindInterfaceForPartitionFromMinNMax(
                     bitset_min.at(i_rank), bitset_max.at(i_rank), indices_min_2d, indices_max_2d,
                    dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux), &partition_interface_background.at(i_rank));
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
            } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
                IniFindInterfaceForPartitionFromMinNMax(
                    bitset_min.at(i_rank), bitset_max.at(i_rank), indices_min_3d, indices_max_3d,
                    dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux), &partition_interface_background.at(i_rank));
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
            }
        }
    }
    for (int i_level = max_level; i_level > 0; --i_level) {
        DefMap<std::set<int>> outmost_for_all_ranks;
        int i_level_lower = i_level - 1;
        if (rank_id == 0) {  // partition nodes on rank 0
            std::vector<std::vector<DefMap<DefAmrIndexUint>>> vec_nodes_ranks(num_ranks);
            int i_rank = 0;
            // min and max code at i_level_lower refinement level
            std::vector<DefSFCodeToUint> code_min_current(num_ranks), code_max_current(num_ranks);
            std::vector<DefSFBitset> bitset_cell_lower_ghost, bitset_all_ghost;
            std::vector<DefSFBitset> corresponding_ones(dims),
                domain_min_m1_n_level(dims), domain_max_p1_n_level(dims);
            std::vector<DefMap<DefAmrIndexUint>> map_ghost_lower_tmp_ranks(num_ranks),
                map_ghost_upper_tmp_ranks(num_ranks);
            if (dims == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
                const SFBitsetAux2D sfbitset_aux_2d = dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux);
                for (auto i = 0; i < num_ranks; ++i) {
                    code_min_current.at(i) = sfbitset_aux_2d.SFBitsetToNHigherLevel(
                        i_level_lower, ull_min.at(i)).to_ullong();
                    code_max_current.at(i) = sfbitset_aux_2d.SFBitsetToNHigherLevel(
                        i_level_lower, ull_max.at(i)).to_ullong();
                }
                GetNLevelCorrespondingOnes2D(i_level_lower, sfbitset_aux_2d, &corresponding_ones);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
            } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
                const SFBitsetAux3D sfbitset_aux_3d = dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux);
                for (auto i = 0; i < num_ranks; ++i) {
                    code_min_current.at(i) = sfbitset_aux_3d.SFBitsetToNHigherLevel(
                        i_level_lower, ull_min.at(i)).to_ullong();
                    code_max_current.at(i) = sfbitset_aux_3d.SFBitsetToNHigherLevel(
                        i_level_lower, ull_max.at(i)).to_ullong();
                }
                GetNLevelCorrespondingOnes3D(i_level_lower, sfbitset_aux_3d, &corresponding_ones);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
            }
            sfbitset_aux.GetMinM1AtGivenLevel(i_level_lower, indices_min, &domain_min_m1_n_level);
            sfbitset_aux.GetMaxP1AtGivenLevel(i_level_lower, indices_max, &domain_max_p1_n_level);
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            std::vector<DefSFCodeToUint>::iterator iter_index;
            std::vector<DefSFBitset> nodes_in_region;
            DefAmrIndexUint num_extend_coarse2fine = ptr_vec_grid_info->at(i_level_lower)->k0NumCoarse2FineLayer_ - 1;
            DefMap<DefAmrIndexUint> map_extend_coarse2fine;
            DefAmrIndexLUint num_layer_near_interface = 1;
            if (k0NumPartitionOuterLayers_ > num_extend_coarse2fine) {
                num_layer_near_interface = k0NumPartitionOuterLayers_ - num_extend_coarse2fine;
            }
            std::vector<DefMap<DefAmrIndexUint>> map_interface_near_coarse2fine(num_ranks);
            bool bool_near_coarse2fine;
            for (auto& iter_node : vec_sfbitset_rank0.at(i_level)) {
                background_code = (sfbitset_aux.SFBitsetToNLowerLevelVir(i_level_lower, iter_node.first)).to_ullong();
                iter_index = std::lower_bound(ull_max.begin(), ull_max.end(), background_code);
                i_rank = static_cast<int>(iter_index - ull_max.begin());
#ifdef DEBUG_CHECK_GRID
                if (i_rank == num_max) {
                    // lower_bound returns the next element of last in ull_max if not found the desired one,
                    // which means that the space filling code of the node exceeds the maximum given by ull_max
                    LogManager::LogError("nodes is out of computational domain"
                        " in MpiManager::IniSendNReceivePartitionedGrid in "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
#endif  // DEBUG_CHECK_GRID
                if (dims == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
                    // search potential ghost nodes near partition interface
                    int_interface_status = CheckNodeOnOuterBoundaryOfBackgroundCell2D(i_level_lower,
                        code_min_current.at(i_rank), code_max_current.at(i_rank), iter_node.first,
                        dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux), domain_min_m1_n_level,
                        domain_max_p1_n_level, corresponding_ones, partition_interface_background.at(i_rank));
                    if ((int_interface_status&1) == 1) {  // lower interface
                        if (num_ghost_lower == 0) {  // only one ghost layer at given level
                            map_ghost_lower_tmp_ranks.at(i_rank).insert({iter_node.first, flag_size0});
                        } else {
                            SearchForGhostLayerForMinNMax2D(iter_node.first, num_ghost_lower,
                                code_min_current.at(i_rank), &MpiManager::CodeSmallerThanMin,
                                flag_size0, dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux),
                                domain_min_m1_n_level, domain_max_p1_n_level, &map_ghost_lower_tmp_ranks.at(i_rank));
                        }
                    }
                    if ((int_interface_status&2) == 2) {  // upper interface
                        SearchForGhostLayerForMinNMax2D(iter_node.first, num_ghost_upper,
                            code_max_current.at(i_rank), &MpiManager::CodeLargerThanMax,
                            flag_size0, dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux),
                            domain_min_m1_n_level, domain_max_p1_n_level, &map_ghost_upper_tmp_ranks.at(i_rank));
                    }
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
                } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
                    // search potential ghost nodes near partition interface
                    int_interface_status = CheckNodeOnOuterBoundaryOfBackgroundCell3D(i_level_lower,
                        code_min_current.at(i_rank), code_max_current.at(i_rank), iter_node.first,
                        dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux), domain_min_m1_n_level,
                        domain_max_p1_n_level, corresponding_ones, partition_interface_background.at(i_rank));
                    if ((int_interface_status&1) == 1) {  // lower interface
                        SearchForGhostLayerForMinNMax3D(iter_node.first, num_ghost_lower,
                            code_min_current.at(i_rank), &MpiManager::CodeSmallerThanMin,
                            flag_size0, dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux),
                            domain_min_m1_n_level, domain_max_p1_n_level, &map_ghost_lower_tmp_ranks.at(i_rank));
                    }
                    if ((int_interface_status&2) == 2) {  // upper interface
                        SearchForGhostLayerForMinNMax3D(iter_node.first, num_ghost_upper,
                            code_max_current.at(i_rank), &MpiManager::CodeLargerThanMax,
                            flag_size0, dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux),
                            domain_min_m1_n_level, domain_max_p1_n_level, &map_ghost_upper_tmp_ranks.at(i_rank));
                    }
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
                }
                if (int_interface_status) {
                    sfbitset_aux.FindNodesInReginOfGivenLength(iter_node.first,
                        num_layer_near_interface, domain_min_m1_n_level, domain_max_p1_n_level,
                        &nodes_in_region);
                    bool_near_coarse2fine = false;
                    for (const auto& iter_node_in_region : nodes_in_region) {
                        if ((vec_sfbitset_rank0.at(i_level).find(iter_node_in_region)
                            != vec_sfbitset_rank0.at(i_level).end())
                            && (vec_sfbitset_rank0.at(i_level).at(iter_node_in_region)&flag_fine2coarse_outmost)) {
                            // at least one node in the region is on the outmost coarse to fine layer
                            bool_near_coarse2fine = true;
                            break;
                        }
                    }
                    if (bool_near_coarse2fine) {
                        map_interface_near_coarse2fine.at(i_rank).insert({iter_node.first, 0});
                    }
                }
                if (i_rank == 0) {  // nodes on rank 0 do not need send to other ranks
                    ptr_sfbitset_each->at(i_level).insert(iter_node);
                } else {   // nodes on rank 0 need to send to other ranks
                    if (i_counts.at(i_rank) == 0) {
                        vec_nodes_ranks.at(i_rank).push_back({iter_node});
                        i_chunk_each_rank.at(i_rank) +=1;
                        ++i_counts.at(i_rank);
                    } else {
                        vec_nodes_ranks.at(i_rank).at(i_chunk_each_rank.at(i_rank)).insert(iter_node);
                        ++i_counts.at(i_rank);
                        if (i_counts.at(i_rank) == max_buffer) {
                            // check if size of send buffer exceeds limits of int
                            i_counts.at(i_rank) = 0;
                        }
                    }
                }
                // node on outmost fine to coarse refinement interface
                if (iter_node.second & flag_fine2coarse_outmost) {
                    if (outmost_for_all_ranks.find(iter_node.first) == outmost_for_all_ranks.end()) {
                        outmost_for_all_ranks.insert({iter_node.first, {i_rank}});
                    } else {
                        outmost_for_all_ranks.at(iter_node.first).insert(i_rank);
                    }
                    // find nodes in all coarse to fine layers
                    sfbitset_aux.FindNodesInReginOfGivenLength(iter_node.first,
                        num_extend_coarse2fine, domain_min_m1_n_level, domain_max_p1_n_level,
                        &nodes_in_region);
                    for (const auto& iter_region : nodes_in_region) {
                        if (map_extend_coarse2fine.find(iter_region) == map_extend_coarse2fine.end()
                            && vec_sfbitset_rank0.at(i_level).find(iter_region)
                            != vec_sfbitset_rank0.at(i_level).end()) {
                            map_extend_coarse2fine.insert({iter_region, 0});
                        }
                    }
                }
            }

            // ghost nodes on rank 0
            for (const auto& iter_ghost : map_ghost_upper_tmp_ranks.at(0)) {
                if (vec_sfbitset_rank0.at(i_level).find(iter_ghost.first)
                    !=vec_sfbitset_rank0.at(i_level).end()) {
                    if (vec_sfbitset_rank0.at(i_level).at(iter_ghost.first) & flag_fine2coarse_outmost) {
                        if (outmost_for_all_ranks.find(iter_ghost.first) == outmost_for_all_ranks.end()) {
                            outmost_for_all_ranks.insert({iter_ghost.first, {0}});
                        } else {
                            outmost_for_all_ranks.at(iter_ghost.first).insert(0);
                        }
                    }
                    if (ptr_sfbitset_each->at(i_level).find(iter_ghost.first) == ptr_sfbitset_each->at(i_level).end()) {
                        ptr_sfbitset_each->at(i_level).insert(iter_ghost);
                    }
                }
            }
            for (const auto& iter_ghost : map_ghost_lower_tmp_ranks.at(0)) {
                if (vec_sfbitset_rank0.at(i_level).find(iter_ghost.first)
                    !=vec_sfbitset_rank0.at(i_level).end()) {
                    if (vec_sfbitset_rank0.at(i_level).at(iter_ghost.first) & flag_fine2coarse_outmost) {
                        if (outmost_for_all_ranks.find(iter_ghost.first) == outmost_for_all_ranks.end()) {
                            outmost_for_all_ranks.insert({iter_ghost.first, {0}});
                        } else {
                            outmost_for_all_ranks.at(iter_ghost.first).insert(0);
                        }
                    }
                    if (ptr_sfbitset_each->at(i_level).find(iter_ghost.first) == ptr_sfbitset_each->at(i_level).end()) {
                        ptr_sfbitset_each->at(i_level).insert(iter_ghost);
                    }
                }
            }
            // nodes in both refinement and mpi communication layers
            DefMap<DefAmrIndexUint> map_ghost_n_refinement;
            DefSFCodeToUint code_tmp;
            for (const auto& iter_ghost : map_interface_near_coarse2fine.at(0)) {
                sfbitset_aux.FindNodesInReginOfGivenLength(iter_ghost.first,
                    k0NumPartitionOuterLayers_, domain_min_m1_n_level, domain_max_p1_n_level,
                    &nodes_in_region);
                for (const auto& iter_region : nodes_in_region) {
                    code_tmp = iter_region.to_ullong();
                    if ((code_tmp < code_min_current.at(0)
                        || code_tmp > code_max_current.at(0))
                        && (map_extend_coarse2fine.find(iter_region)
                        != map_extend_coarse2fine.end()
                        || vec_sfbitset_rank0.at(i_level).find(iter_region)
                        == vec_sfbitset_rank0.at(i_level).end())) {
                        ptr_sfbitset_ghost_each->at(i_level).insert({iter_region, 0});
                    }
                }
            }
            // ghost nodes on other ranks
            std::vector<std::vector<DefMap<DefAmrIndexUint>>> vec_ghost_nodes_ranks(num_ranks);
            std::vector<int> i_ghost_chunk_each_rank(num_ranks, -1), i_ghost_counts(num_ranks, 0);
            for (auto i_rank = 1; i_rank < num_ranks; ++i_rank) {
                std::vector<std::unique_ptr<char[]>> vec_ptr_buffer;
                std::vector<MPI_Request> reqs_send;
                // find nodes on extened ghost layers for mpi communication
                DefMap<DefAmrIndexUint> partition_interface_rank0_lower_level;
                for (const auto& iter_ghost : map_ghost_upper_tmp_ranks.at(i_rank)) {
                    if (partition_interface_rank0_lower_level.find(iter_ghost.first)
                        == partition_interface_rank0_lower_level.end()
                        &&vec_sfbitset_rank0.at(i_level).find(iter_ghost.first)
                        !=vec_sfbitset_rank0.at(i_level).end()) {
                        if (vec_sfbitset_rank0.at(i_level).at(iter_ghost.first) & flag_fine2coarse_outmost) {
                            if (outmost_for_all_ranks.find(iter_ghost.first) == outmost_for_all_ranks.end()) {
                                outmost_for_all_ranks.insert({iter_ghost.first, {i_rank}});
                            } else {
                                outmost_for_all_ranks.at(iter_ghost.first).insert(i_rank);
                            }
                        }
                        partition_interface_rank0_lower_level.insert(iter_ghost);
                        if (i_counts.at(i_rank) == 0) {
                            vec_nodes_ranks.at(i_rank).push_back({iter_ghost});
                            i_chunk_each_rank.at(i_rank) +=1;
                            ++i_counts.at(i_rank);
                        } else {
                            vec_nodes_ranks.at(i_rank).at(i_chunk_each_rank.at(i_rank)).insert(iter_ghost);
                            ++i_counts.at(i_rank);
                            if (i_counts.at(i_rank) == max_buffer) {
                                i_counts.at(i_rank) = 0;
                            }
                        }
                    }
                }
                for (const auto& iter_ghost : map_ghost_lower_tmp_ranks.at(i_rank)) {
                    if (partition_interface_rank0_lower_level.find(iter_ghost.first)
                        == partition_interface_rank0_lower_level.end()
                        &&vec_sfbitset_rank0.at(i_level).find(iter_ghost.first)
                        !=vec_sfbitset_rank0.at(i_level).end()) {
                        if (vec_sfbitset_rank0.at(i_level).at(iter_ghost.first) & flag_fine2coarse_outmost) {
                            if (outmost_for_all_ranks.find(iter_ghost.first) == outmost_for_all_ranks.end()) {
                                outmost_for_all_ranks.insert({iter_ghost.first, {i_rank}});
                            } else {
                                outmost_for_all_ranks.at(iter_ghost.first).insert(i_rank);
                            }
                        }
                        partition_interface_rank0_lower_level.insert(iter_ghost);
                        if (i_counts.at(i_rank) == 0) {
                            vec_nodes_ranks.at(i_rank).push_back({iter_ghost});
                            i_chunk_each_rank.at(i_rank) +=1;
                            ++i_counts.at(i_rank);
                        } else {
                            vec_nodes_ranks.at(i_rank).at(i_chunk_each_rank.at(i_rank)).insert(iter_ghost);
                            ++i_counts.at(i_rank);
                            if (i_counts.at(i_rank) == max_buffer) {
                                i_counts.at(i_rank) = 0;
                            }
                        }
                    }
                }
                // send grid nodes
                int buffer_size_send = 0;
                int num_chunks;
                num_chunks = static_cast<int>(vec_nodes_ranks.at(i_rank).size());
                vec_ptr_buffer.resize(num_chunks);
                reqs_send.resize(num_chunks);
                MPI_Send(&num_chunks, 1, MPI_INT, i_rank, i_level, MPI_COMM_WORLD);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    vec_ptr_buffer.at(i_chunk) = SerializeNodeSFBitset(vec_nodes_ranks.at(i_rank).at(i_chunk),
                        &buffer_size_send);
                    MPI_Send(&buffer_size_send, 1, MPI_INT, i_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, i_rank,
                    i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
                }
                MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
                // nodes used on both refinement layers and outer mpi communication layers
                map_ghost_n_refinement.clear();
                for (const auto& iter_ghost : map_interface_near_coarse2fine.at(i_rank)) {
                    sfbitset_aux.FindNodesInReginOfGivenLength(iter_ghost.first,
                        k0NumPartitionOuterLayers_, domain_min_m1_n_level, domain_max_p1_n_level,
                        &nodes_in_region);
                    for (const auto& iter_region : nodes_in_region) {
                        code_tmp = iter_region.to_ullong();
                        if ((code_tmp < code_min_current.at(i_rank)
                            || code_tmp > code_max_current.at(i_rank))
                            &&map_ghost_n_refinement.find(iter_region) == map_ghost_n_refinement.end()
                            && (map_extend_coarse2fine.find(iter_region)
                            != map_extend_coarse2fine.end()
                            || vec_sfbitset_rank0.at(i_level).find(iter_region)
                            == vec_sfbitset_rank0.at(i_level).end())) {
                            map_ghost_n_refinement.insert({iter_region, 0});
                        }
                    }
                }
                for (const auto& iter_ghost : map_ghost_n_refinement) {
                    if (i_ghost_counts.at(i_rank) == 0) {
                        vec_ghost_nodes_ranks.at(i_rank).push_back({iter_ghost});
                        i_ghost_chunk_each_rank.at(i_rank) +=1;
                        ++i_ghost_counts.at(i_rank);
                    } else {
                        vec_ghost_nodes_ranks.at(i_rank).at(i_ghost_chunk_each_rank.at(i_rank)).insert(iter_ghost);
                        ++i_ghost_counts.at(i_rank);
                        if (i_ghost_counts.at(i_rank) == max_buffer) {
                            i_ghost_counts.at(i_rank) = 0;
                        }
                    }
                }
                // send ghost nodes
                num_chunks = static_cast<int>(vec_ghost_nodes_ranks.at(i_rank).size());
                vec_ptr_buffer.resize(num_chunks);
                reqs_send.resize(num_chunks);
                MPI_Send(&num_chunks, 1, MPI_INT, i_rank, i_level, MPI_COMM_WORLD);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    vec_ptr_buffer.at(i_chunk) = SerializeNodeSFBitset(
                        vec_ghost_nodes_ranks.at(i_rank).at(i_chunk), &buffer_size_send);
                    MPI_Send(&buffer_size_send, 1, MPI_INT, i_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, i_rank,
                        i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
                }
                MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
            }
        } else {
            int num_chunks = 0;
            // receive grid nodes
            MPI_Recv(&num_chunks, 1, MPI_INT, 0, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int buffer_size_receive;
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
                MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                    i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                DeserializeNodeSFBitset(flag_size0, vec_ptr_buffer, &ptr_sfbitset_each->at(i_level));
            }

            // receive ghost nodes
            MPI_Recv(&num_chunks, 1, MPI_INT, 0, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                int buffer_size_receive;
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
                MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                    i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                DeserializeNodeSFBitset(flag_size0, vec_ptr_buffer, &ptr_sfbitset_ghost_each->at(i_level));
            }
        }
        // send and receive nodes on outmost coarse to fine refinement interfaces
        DefAmrUint flag0 = static_cast<DefAmrUint>(flag_size0);
        IniSendNReceiveCoarse2Fine0Interface(dims, i_level_lower,
            ptr_vec_grid_info->at(i_level_lower)->k0NumCoarse2FineLayer_, flag0,
            outmost_for_all_ranks, &ptr_vec_grid_info->at(i_level_lower)->map_ptr_interface_layer_info_);
    }
}
/**
 * @brief function to send, receive grid information for each rank at all levels.
 * @param[in] flag_size0 flag of existing node
 * @param[in] flag_fine2coarse_outmost flag indicting node is on the outmost coarse to fine refinement layer.
 * @param[in] dims dimensions of the grid.
 * @param[in] max_level the maximum level of the grid.
 * @param[in] bitset_domain_min the minimum space filling code of the computational domain
 * @param[in] bitset_domain_max the maximum space filling code of the computational domain
 * @param[in] vec_cost computational cost of the each node for all refinement levels.
 * @param[in] bitset_aux An interface for manipulating bitsets.
 * @param[in] vec_tracking_creator creator for tracking grid information.
 * @param[in] ini_sfbitset_one_lower_level_rank0 pointer to initial mesh generated on rank 0.
 * @param[out] ptr_sfbitset_bound_current pointer to lower and upper bound space filling code of current rank
 * @param[out] ptr_sfbitset_one_lower_level_current_rank pointer to nodes on the current rank.
 * @param[out] ptr_sfbitset_ghost_one_lower_level_current_rank  pointer to non-mpi communication ghost nodes 
 *             near coarse to fine refinement interface at one level lower on the current rank.
 * @param[out] ptr_vec_grid_info pointer to vector of grid information.
 */
void MpiManager::IniSendNReceiveGridInfoAtGivenLevels(const DefAmrIndexUint flag_size0,
    const DefAmrUint flag_fine2coarse_outmost, const DefAmrIndexUint dims, const DefAmrIndexUint max_level,
    const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
    const std::vector<DefAmrIndexLUint>& indices_min, const std::vector<DefAmrIndexLUint>& indices_max,
    const std::vector<DefAmrIndexLUint>& vec_cost, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_creator,
    const std::vector<DefMap<DefAmrIndexUint>> ini_sfbitset_one_lower_level_rank0,
    std::array<DefSFBitset, 2>* const ptr_sfbitset_bound_current,
    std::vector<DefMap<DefAmrIndexUint>>* const  ptr_sfbitset_one_lower_level_current_rank,
    std::vector<DefMap<DefAmrIndexUint>>* const  ptr_sfbitset_ghost_one_lower_level_current_rank,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) {
    std::vector<DefSFBitset> bitset_min(num_of_ranks_), bitset_max(num_of_ranks_);
    if (rank_id_ == 0) {
        DefSizet vec_size = vec_cost.size();
        if (ini_sfbitset_one_lower_level_rank0.size() != vec_size) {
            LogManager::LogError("size of the input vector (vec_cost) is different from the size of the input vector "
             "(ptr_ini_sfbitset_one_lower_level_rank0) in MpiManager::IniSendNReceiveGridInfoAtGivenLevels in "
             + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        if (dims == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTION
            IniTraverseBackgroundForPartitionRank0(bitset_domain_min, bitset_domain_max,
                vec_cost, ini_sfbitset_one_lower_level_rank0,
                dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux), &bitset_min, &bitset_max);
#endif  // DEBUG_DISABLE_2D_FUNCTION
        } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTION
            IniTraverseBackgroundForPartitionRank0(bitset_domain_min, bitset_domain_max,
                vec_cost, ini_sfbitset_one_lower_level_rank0,
                dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux), &bitset_min, &bitset_max);

#endif  // DEBUG_DISABLE_3D_FUNCTION
        }
    }

    // broadcast minimum and maximum space filling code of all ranks to each rank
    MPI_Bcast(&bitset_min[0], static_cast<int>(num_of_ranks_), MPI_CODE_UINT_TYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitset_max[0], static_cast<int>(num_of_ranks_), MPI_CODE_UINT_TYPE, 0, MPI_COMM_WORLD);
    vec_sfbitset_min_all_ranks_ = bitset_min;
    vec_sfbitset_max_all_ranks_ = bitset_max;
    sfbitset_min_current_rank_ = bitset_min.at(rank_id_);
    sfbitset_max_current_rank_ = bitset_max.at(rank_id_);

    ptr_sfbitset_bound_current->at(0) = bitset_min.at(rank_id_);
    ptr_sfbitset_bound_current->at(1) = bitset_max.at(rank_id_);

    // send and receive grid nodes and nodes on the second outmost coarse to fine refinement interface
    IniSendNReceivePartitionedGrid(dims, flag_size0, flag_fine2coarse_outmost,
        bitset_min, bitset_max, indices_min, indices_max,
        sfbitset_aux, ini_sfbitset_one_lower_level_rank0, ptr_sfbitset_one_lower_level_current_rank,
        ptr_sfbitset_ghost_one_lower_level_current_rank, ptr_vec_grid_info);

    // send and receive tracking nodes
    for (DefAmrIndexUint i_level = max_level; i_level > 0; --i_level) {
        GridInfoInterface& grid_info = *(ptr_vec_grid_info->at(i_level));
        IniSendNReceiveTracking(dims, i_level, bitset_max, sfbitset_aux,
           vec_tracking_creator, &grid_info.map_ptr_tracking_grid_info_);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
