//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager_mpi.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is enabled.
* @date  2023-5-2
*/
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#include <limits>
#include "grid/grid_manager.h"
#include "io/log_write.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to partition grid using space-filling curves. 
 * @param[in] vec_sfbitset space-filling codes for nodes of entire mesh.
 * @param[in] mpi_manager An instance of MpiManager containing rank and number of
 *     ranks information.
 * @param[out] ptr_bitset_min  pointer to minimum space filling code for each partition.
 * @param[out] ptr_bitset_max  pointer to maximum space filling code for each partition.
 */
void GridManagerInterface::IniPartiteGridBySpaceFillingCurves(
    const std::vector<DefMap<DefUint>>& vec_sfbitset,
    const MpiManager& mpi_manager,
    std::vector<DefSFBitset>* const ptr_bitset_min,
    std::vector<DefSFBitset>* const ptr_bitset_max) {
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    ptr_bitset_min->resize(num_ranks);
    ptr_bitset_max->resize(num_ranks);
    DefMap<DefUint> background_occupied;
    DefSFBitset bitset_background;
    DefLUint num_bk_node = CalNumOfBackgroundNode();
    DefLUint bk_cost = vec_ptr_grid_info_.at(0)->computational_cost_;
    DefLUint sum_load = num_bk_node * bk_cost;
    // add computational cost of refined nodes
    for (DefSizet i_level = k0MaxLevel_; i_level > 0; --i_level) {
        DefLUint node_cost = vec_ptr_grid_info_
            .at(i_level)->computational_cost_;
        for (const auto& iter_low : vec_sfbitset.at(i_level)) {
            bitset_background = NodeAtNLowerLevel(i_level, iter_low.first);
            if (background_occupied.find(bitset_background) ==
                background_occupied.end()) {
                background_occupied.insert(
                    { bitset_background, node_cost });
                sum_load += (node_cost - bk_cost);
            } else {
                background_occupied.at(bitset_background) += node_cost;
                sum_load += node_cost;
            }
        }
    }
    // calculate loads on each rank
    DefLUint ave_load = static_cast<DefLUint>(sum_load / num_ranks) + 1;
    DefLUint load_rank0 = sum_load - (num_ranks - 1) * ave_load;
    std::vector<DefLUint> rank_load(num_ranks, ave_load);
    rank_load.at(0) = load_rank0;   // load at rank 0

    TraverseBackgroundForPartition(rank_load, background_occupied,
        ptr_bitset_min, ptr_bitset_max);
}
/**
 * @brief function to serialize data (a uint) and save it to buffer
 * @param[in] map_nodes node information
 * @param[out] buffer buffer to send data
 * @return size of the buffer 
 */
int GridManagerInterface::SerializeNode(const DefMap<DefUint>& map_nodes,
    std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset), node_size = sizeof(DefUint);
    int num_nodes;
    if  (sizeof(int) + map_nodes.size() *(key_size + node_size) > 0x7FFFFFFF) {
        LogError("size of the buffer is greater than the maximum of int in MpiManager::SerializeData(DefMap<DefUint>)");
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
 * @brief function to deserialize data (a uint) from a buffer
 * @param[in] buffer buffer has received data
 * @param[out] map_nodes node information
 */
void GridManagerInterface::DeserializeNode(const std::unique_ptr<char[]>& buffer,
    DefMap<DefUint>* const map_nodes) const {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset), node_size = sizeof(DefUint);
    // number of nodes
    int num_nodes;
    int position = 0;
    std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
    DefUint node_data;
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
 * @brief function to initialize partitioned grid for each rank
 * @param[in] bitset_max maximum space filling code for each partition.
 * @param[in] vec_sfbitset space-filling codes for nodes of entire mesh.
 * @param[in] mpi_manager reference to mpi manager.
 * @param[out] ptr_sfbitset_each nodes on current rank.
 */
void GridManagerInterface::IniSendNReceivePartitionedGrid(
    const std::vector<DefSFBitset>& bitset_max,
    const std::vector<DefMap<DefUint>>& vec_sfbitset,
    const MpiManager& mpi_manager,
    std::vector<DefMap<DefUint>>* const ptr_sfbitset_each) const {
    int max_level = static_cast<int>(vec_sfbitset.size()) - 1;
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefUint) - 1;
    DefSFCodeToUint background_code;
    DefSizet num_max = ull_max.size();
    for (int i_level = max_level; i_level > 0; --i_level) {
        if (rank_id == 0) {  // partition nodes on rank 0
            std::vector<std::vector<DefMap<DefUint>>> vec_nodes_ranks(num_ranks);
            int index;
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            for (const auto& iter_node : vec_sfbitset.at(i_level)) {
                background_code = (NodeAtNLowerLevel(i_level, iter_node.first)).to_ullong();
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
                        LogError("nodes is out of computational domain in MpiManager::IniSendNReceivePartitionedGrid");
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
            if (num_ranks > 1) {
                int buffer_size_send = 0;
                for (auto iter_rank = 1; iter_rank < num_ranks; ++iter_rank) {
                    int num_chunks = static_cast<int>(vec_nodes_ranks.at(iter_rank).size());
                    MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, i_level, MPI_COMM_WORLD);
                    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
                    std::vector<std::unique_ptr<MPI_Request>> reqs_send(num_chunks);
                    std::vector<std::unique_ptr<MPI_Status>> stats_send(num_chunks);
                    for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                        buffer_size_send = SerializeNode(vec_nodes_ranks.at(iter_rank).at(i_chunk),
                        vec_ptr_buffer.at(i_chunk));
                        reqs_send[i_chunk].reset(new MPI_Request);
                        stats_send[i_chunk].reset(new MPI_Status);
                        MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                        MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                        i_chunk, MPI_COMM_WORLD, reqs_send[i_chunk].get());
                    }
                    MPI_Waitall(num_chunks, reinterpret_cast<MPI_Request*>(reqs_send.data()),
                    reinterpret_cast<MPI_Status*>(stats_send.data()));
                }
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
                DeserializeNode(vec_ptr_buffer, &ptr_sfbitset_each->at(i_level));
            }
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
