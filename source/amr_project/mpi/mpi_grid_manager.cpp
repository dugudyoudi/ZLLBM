//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is enabled.
* @date  2023-5-2
*/
#include <limits>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to serialize data (a DefSFBitset and a DefAmrUint) and save it to buffer
 * @param[in] map_nodes node information
 * @param[out] buffer buffer to send data
 * @return size of the buffer 
 */
int MpiManager::SerializeNodeStoreUint(const DefMap<DefAmrUint>& map_nodes,
    std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset), node_size = sizeof(DefAmrUint);
    int num_nodes = 1;
    if  (sizeof(int) + map_nodes.size() *(key_size + node_size) > 0x7FFFFFFF) {
        LogManager::LogError("size of the buffer is greater"
         " than the maximum of int in MpiManager::SerializeData(DefMap<DefAmrUint>)");
    } else {
        num_nodes = static_cast<int>(map_nodes.size());
    }
    int buffer_size = sizeof(int) + (key_size + node_size) * num_nodes;
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
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
    return buffer_size;
}
/**
 * @brief function to deserialize data (a DefSFBitset and a DefAmrUint) from a buffer
 * @param[in] buffer buffer has received data
 * @param[out] map_nodes node information
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
    DefSFCodeToUint key_net;
    DefSFBitset key_host;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_net, ptr_buffer + position, key_size);
        key_host = static_cast<DefSFBitset>(key_net);
        position += key_size;
        std::memcpy(&node_data, ptr_buffer + position, node_size);
        position += node_size;
        map_nodes->insert({ key_host, node_data });
    }
}
/**
 * @brief function to serialize space filling code (DefSFBitset) and save it to buffer
 * @param[in] map_nodes node information
 * @param[out] buffer buffer to send data
 * @return size of the buffer 
 */
int MpiManager::SerializeNodeSFBitset(const DefMap<DefAmrIndexUint>& map_nodes,
    std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes = 1;
    if  (sizeof(int) + map_nodes.size() *(key_size) > 0x7FFFFFFF) {
        LogManager::LogError("size of the buffer is greater than the"
         " maximum of int in MpiManager::SerializeData(DefMap<DefAmrUint>)");
    } else {
        num_nodes = static_cast<int>(map_nodes.size());
    }
    int buffer_size = sizeof(int) + (key_size) * num_nodes;
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
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
    return buffer_size;
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
    DefSFCodeToUint key_net;
    DefSFBitset key_host;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_net, ptr_buffer + position, key_size);
        key_host = static_cast<DefSFBitset>(key_net);
        position += key_size;
        map_nodes->insert({ key_host, flag_node });
    }
}
/**
 * @brief function to send and receive partitioned grid for each rank
 * @param[in] flag_size0 flag of existing node
 * @param[in] bitset_max maximum space filling code for each partition.
 * @param[in] vec_sfbitset space-filling codes for nodes at given refinements levels.
 * @param[in] bitset_aux class manage space filling curves.
 * @param[out] ptr_sfbitset_each nodes on current rank at sizeof(vec_sfbitset) levels.
 */
void MpiManager::IniSendNReceivePartitionedGrid(const DefAmrIndexUint flag_size0,
    const std::vector<DefSFBitset>& bitset_max,
    const std::vector<DefMap<DefAmrIndexUint>>& vec_sfbitset,
    const SFBitsetAuxInterface& bitset_aux,
    std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_each) const {
    if (vec_sfbitset.size() - 1 > INT_MAX) {
        LogManager::LogError("size of vector vec_sfbitset is larger than INT_MAX");
    }
    int max_level = static_cast<int>(vec_sfbitset.size() - 1);
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());

    if (rank_id == 0) {
        if (vec_sfbitset.size() != ptr_sfbitset_each->size()) {
            LogManager::LogError("size of the input vector (vec_sfbitset) is different from the"
                " size of the output vector (ptr_sfbitset_each) in MpiManager::IniSendNReceivePartitionedGrid");
        }
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefAmrIndexUint) - 1;
    DefSFCodeToUint background_code;
    DefSizet num_max = ull_max.size();
    for (int i_level = max_level; i_level > 0; --i_level) {
        if (rank_id == 0) {  // partition nodes on rank 0
            std::vector<std::vector<DefMap<DefAmrIndexUint>>> vec_nodes_ranks(num_ranks);
            int index;
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            for (const auto& iter_node : vec_sfbitset.at(i_level)) {
                background_code = (bitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first)).to_ullong();
                auto iter_index = std::lower_bound(ull_max.begin(),
                    ull_max.end(), background_code);
                index = static_cast<int>(iter_index - ull_max.begin());
                if (index == 0) {
                    ptr_sfbitset_each->at(i_level).insert(iter_node);
                } else {
#ifdef DEBUG_CHECK_GRID
                    if (index == num_max) {
                        // lower_bound returns the next element of last in ull_max if not found the desired one,
                        // which means that the space filling code of the node exceeds the maximum given by ull_max
                        LogManager::LogError("nodes is out of computational domain"
                            " in MpiManager::IniSendNReceivePartitionedGrid");
                    }
#endif  // DEBUG_CHECK_GRID
                    if (i_counts.at(index) == 0) {
                        vec_nodes_ranks.at(index).push_back({iter_node});
                        i_chunk_each_rank.at(index) +=1;
                        ++i_counts.at(index);
                    } else {
                        vec_nodes_ranks.at(index).at(i_chunk_each_rank.at(index)).insert(iter_node);
                        ++i_counts.at(index);
                        if (i_counts.at(index) == max_buffer) {
                            // check if size of send buffer exceeds limits of int
                            i_counts.at(index) = 0;
                        }
                    }
                }
            }
            int buffer_size_send = 0;
            for (auto iter_rank = 1; iter_rank < num_ranks; ++iter_rank) {
                int num_chunks = static_cast<int>(vec_nodes_ranks.at(iter_rank).size());
                MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, i_level, MPI_COMM_WORLD);
                std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
                std::vector<MPI_Request> reqs_send(num_chunks);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    buffer_size_send = SerializeNodeSFBitset(vec_nodes_ranks.at(iter_rank).at(i_chunk),
                    vec_ptr_buffer.at(i_chunk));
                    MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                    i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
                }
                MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
            }
        } else {
            int num_chunks = 0;
            MPI_Recv(&num_chunks, 1, MPI_INT, 0, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                int buffer_size_receive;
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
                MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                    i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                DeserializeNodeSFBitset(flag_size0, vec_ptr_buffer, &ptr_sfbitset_each->at(i_level));
            }
        }
    }
}
/**
 * @brief function to send, receive grid information for each rank at all levels.
 * @param[in] flag_size0 flag of existing node
 * @param[in] dims dimensions of the grid.
 * @param[in] max_level the maximum level of the grid.
 * @param[in] bitset_domain_min the minimum space filling code of the computational domain
 * @param[in] bitset_domain_max the maximum space filling code of the computational domain
 * @param[in] vec_cost computational cost of the 
 * @param[in] bitset_aux An interface for manipulating bitsets.
 * @param[in] ini_sfbitset_one_lower_level_rank0 initial mesh generated on rank 0.
 * @param[in] vec_tracking_creator creator for tracking grid information.
 * @param[out] ptr_sfbitset_bound_current pointer to lower and upper bound space filling code of current rank
 * @param[out] ptr_sfbitset_one_lower_level_current_rank pointer to mesh on current rank.
 * @param[out] ptr_vec_grid_info vector of pointers to grid information.
 */
void MpiManager::SendNReceiveGridInfoAtGivenLevels(const DefAmrIndexUint flag_size0,
    const DefAmrIndexUint dims, const DefAmrIndexUint max_level,
    const DefSFBitset bitset_domain_min, const DefSFBitset bitset_domain_max,
    const std::vector<DefAmrIndexLUint>& vec_cost, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefMap<DefAmrIndexUint>>&  ini_sfbitset_one_lower_level_rank0,
    const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_creator,
    std::array<DefSFBitset, 2>* const ptr_sfbitset_bound_current,
    std::vector<DefMap<DefAmrIndexUint>>* const  ptr_sfbitset_one_lower_level_current_rank,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const {
    std::vector<DefSFBitset> bitset_min(num_of_ranks_), bitset_max(num_of_ranks_);
    if (rank_id_ == 0) {
        DefSizet vec_size = vec_cost.size();
        if (ini_sfbitset_one_lower_level_rank0.size() != vec_size) {
            LogManager::LogError("size of the input vector (vec_cost) is different from the size of the input vector "
                "(ini_sfbitset_one_lower_level_rank0) in MpiManager::SendNReceiveGridInfoAtGivenLevels");
        }
        if (dims == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTION
            TraverseBackgroundForPartitionRank0(bitset_domain_min, bitset_domain_max,
                vec_cost, ini_sfbitset_one_lower_level_rank0,
                dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux), &bitset_min, &bitset_max);
#endif  // DEBUG_DISABLE_2D_FUNCTION
        } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTION
            TraverseBackgroundForPartitionRank0(bitset_domain_min, bitset_domain_max,
                vec_cost, ini_sfbitset_one_lower_level_rank0,
                dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux), &bitset_min, &bitset_max);

#endif  // DEBUG_DISABLE_3D_FUNCTION
        }
    }
    MPI_Bcast(&bitset_min[0], static_cast<int>(dims), MPI_CODE_UINT_TYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitset_max[0], static_cast<int>(dims), MPI_CODE_UINT_TYPE, 0, MPI_COMM_WORLD);

    ptr_sfbitset_bound_current->at(0) = bitset_min.at(rank_id_);
    ptr_sfbitset_bound_current->at(1) = bitset_max.at(rank_id_);
    IniSendNReceivePartitionedGrid(flag_size0, bitset_max, ini_sfbitset_one_lower_level_rank0,
     sfbitset_aux, ptr_sfbitset_one_lower_level_current_rank);

    // send and receive tracking and interface nodes
    for (DefAmrIndexUint i_level = max_level; i_level > 0; --i_level) {
        GridInfoInterface& grid_info = *(ptr_vec_grid_info->at(i_level));
        IniSendNReceiveTracking(dims, i_level, bitset_max, sfbitset_aux,
           vec_tracking_creator, &grid_info.map_ptr_tracking_grid_info_);
        IniSendNReceiveRefinementInterface(dims, i_level, grid_info.k0NumCoarse2FineLayer_,
           bitset_max, sfbitset_aux, &grid_info.map_ptr_interface_layer_info_);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
