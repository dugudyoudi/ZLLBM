//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "grid/grid_manager.h"
#include "io/log_write.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
void MpiManager::StartupMpi(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_ranks_);
}
void MpiManager::SetMpiParameters() {
}
/**
* @brief      function for partition mesh based on the space filling curve.
* @param[in]  sfbitset_one_lower_level  nodes at one lower level.
* @param[in]  grid_manager    reference to grid manager.
* @param[out]  ptr_bitset_min  minimum space filling code for each partition.
* @param[out]  ptr_bitset_max  maximum space filling code for each partition.
*/
void MpiManager::IniPartiteGridBySpaceFillingCurves(
    const std::vector<DefMap<DefUint>>& sfbitset_one_lower_level,
    GridManagerInterface const& grid_manager,
    std::vector<DefSFBitset>* const ptr_bitset_min,
    std::vector<DefSFBitset>* const ptr_bitset_max) {
    ptr_bitset_min->resize(num_of_ranks_);
    ptr_bitset_max->resize(num_of_ranks_);
    DefMap<DefUint> background_occupied;
    DefSFBitset bitset_background;
    DefLUint num_bk_node = grid_manager.CalNumOfBackgroundNode();
    DefLUint bk_cost =
        grid_manager.vec_ptr_grid_info_.at(0)->computational_cost_;
    DefLUint sum_load = num_bk_node * bk_cost;
    // add computational cost of refined nodes
    for (DefSizet i_level = grid_manager.k0MaxLevel_; i_level > 0; --i_level) {
        DefLUint node_cost = grid_manager.vec_ptr_grid_info_
            .at(i_level)->computational_cost_;
        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            bitset_background = grid_manager.NodeAtNLowerLevel(
                i_level, iter_low.first);
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
    DefLUint ave_load = static_cast<DefLUint>(sum_load / num_of_ranks_) + 1;
    DefLUint load_rank0 = sum_load - (num_of_ranks_ - 1) * ave_load;
    std::vector<DefLUint> rank_load(num_of_ranks_, ave_load);
    rank_load.at(0) = load_rank0;   // load at rank 0

    grid_manager.TraverseBackgroundForPartition(rank_load, background_occupied,
        ptr_bitset_min, ptr_bitset_max);
}
/**
 * @brief function to serialize data (a uint) and save it to buffer
 * @param[in] map_nodes node information
 * @param[out] buffer buffer to send data
 * @return size of the buffer 
 */
int MpiManager::SerializeData(const DefMap<DefUint>& map_nodes,
    std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset), node_size = sizeof(DefUint);
    int num_nodes;
    if  (sizeof(int) + map_nodes.size() *(key_size + node_size) > 0x7FFFFFFF) {
        LogError("size of the buffer is greater than the maximum of int in MpiManager::SerializeData");
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
    DefSFCodeToUint key_host, key_net;
    for (auto& iter : map_nodes) {
        // divide bitset into chunks and save to buffer
        key_host = iter.first.to_ullong();
        key_net = HtoNUint(key_host);
        std::memcpy(ptr_buffer + position, &key_net, key_size);
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
void MpiManager::DeserializeData(const std::unique_ptr<char[]>& buffer,
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
        key_host = NtoHUint(key_net);
        position += key_size;
        std::memcpy(&node_data, ptr_buffer + position, node_size);
        position += node_size;
        map_nodes->insert({ key_host, node_data });
    }
}
/**
 * @brief function to broadcast grid bounds of all ranks on rank 0 to other ranks.
 * @param[in] ptr_bitset_bounds pointer to bounds.
 */
void MpiManager::IniBroadcastBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds) {
    if (rank_id_ == 0) {  // bitset_bounds on rank has be calculated
        int bit_size = static_cast<int>(ptr_bitset_bounds->size()) * sizeof(DefSFBitset);
        MPI_Bcast(ptr_bitset_bounds->data(), bit_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
}
/**
 * @brief function to initialize partitioned grid for each rank
 * @param[in] bitset_max maximum space filling code for each partition.
 * @param[in] sfbitset_all_one_lower_level  nodes of entire mesh at one lower level.
 * @param[in] grid_manager reference to grid manager.
 * @param[out] ptr_sfbitset_each_one_lower_level nodes on current rank at one lower level.
 */
void MpiManager::IniSendNReceivePartitionedGrid(
    const std::vector<DefSFBitset>& bitset_max,
    const std::vector<DefMap<DefUint>>& sfbitset_all_one_lower_level,
    GridManagerInterface const& grid_manager,
    std::vector<DefMap<DefUint>>* const ptr_sfbitset_each_one_lower_level) const {
    int max_level = static_cast<int>(sfbitset_all_one_lower_level.size()) - 1;
    int num_ranks = num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id_ == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    DefSFCodeToUint background_code;
    DefSizet num_max = ull_max.size();
    for (int i_level = max_level; i_level > 0; --i_level) {
        if (rank_id_ == 0) {  // partition nodes on rank 0
            int buffer_size_send;
            std::vector<std::vector<DefMap<DefUint>>> vec_nodes_ranks(num_ranks);
            int index;
            std::vector<int> i_chunk_each_rank(num_ranks, 1), i_counts(num_ranks, 1);
            for (const auto& iter_node : sfbitset_all_one_lower_level.at(i_level)) {
                background_code = (grid_manager.NodeAtNLowerLevel(i_level, iter_node.first)).to_ullong();
                auto iter_index = std::lower_bound(ull_max.begin(),
                     ull_max.end(), background_code);
                index = static_cast<int>(iter_index - ull_max.begin());
                if (index == 0) {
                    ptr_sfbitset_each_one_lower_level->at(i_level).insert(iter_node);
                } else {
#ifdef DEBUG_CHECK_GRID
                    if (index == num_max) {
                        // lower_bound returns the next element of last in ull_max if not found the desired one,
                        // which means that the space filling code of the node exceeds the maximum given by ull_max
                        LogError("nodes is out of computational domain in MpiManager::IniSendNReceivePartitionedGrid");
                    }
#endif  // DEBUG_CHECK_GRID
                    if (i_counts.at(index) == 1) {
                        vec_nodes_ranks.at(index).push_back({iter_node});
                        i_chunk_each_rank.at(index) +=1;
                    } else {
                        vec_nodes_ranks.at(index).at(i_chunk_each_rank.at(index)).insert(iter_node);
                        ++i_counts.at(index);
                        if (i_counts.at(index) == (std::numeric_limits<int>::max)()) {
                            // check if size of send buffer exceeds limits of int
                            i_counts.at(index) = 1;
                        }
                    }
                }
            }
            if (num_ranks > 1) {
                for (auto iter_rank = 1; iter_rank < num_ranks; ++iter_rank) {
                    int num_chunks = static_cast<int>(vec_nodes_ranks.at(iter_rank).size());
                    MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, i_level, MPI_COMM_WORLD);
                    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
                    std::vector<std::unique_ptr<MPI_Request>> reqs_send(num_chunks);
                    std::vector<std::unique_ptr<MPI_Status>> stats_send(num_chunks);
                    for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                        buffer_size_send = SerializeData(vec_nodes_ranks.at(iter_rank).at(i_chunk),
                        vec_ptr_buffer.at(i_chunk));
                        reqs_send[i_chunk].reset(new MPI_Request);
                        stats_send[i_chunk].reset(new MPI_Status);
                        MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_level, MPI_COMM_WORLD);
                        MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_CHAR, iter_rank,
                        i_level, MPI_COMM_WORLD, reqs_send[i_chunk].get());
                    }
                    MPI_Waitall(num_chunks, reinterpret_cast<MPI_Request*>(reqs_send.data()),
                    reinterpret_cast<MPI_Status*>(stats_send.data()));
                }
            }
        } else {
            int num_chunks;
            MPI_Recv(&num_chunks, sizeof(int), MPI_INT, 0, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
            std::vector<std::unique_ptr<MPI_Request>> reqs_receive(num_chunks);
            std::vector<std::unique_ptr<MPI_Status>> stats_receive(num_chunks);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                int buffer_size_receive;
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Irecv(vec_ptr_buffer.at(i_chunk).get(), buffer_size_receive, MPI_CHAR, 0,
                    i_level, MPI_COMM_WORLD, reqs_receive[i_chunk].get());
                DeserializeData(vec_ptr_buffer.at(i_chunk), &ptr_sfbitset_each_one_lower_level->at(i_level));
            }
            MPI_Waitall(num_chunks, reinterpret_cast<MPI_Request*>(reqs_receive.data()),
                reinterpret_cast<MPI_Status*>(stats_receive.data()));
        }
    }
}
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
void MpiManager::IniSendNReceivePartitionedGeoCoordi(
    const GridManager2D& grid_manager2d,
    const std::vector<DefSFBitset>& bitset_max,
    Geometry2DInterface* const ptr_geo2d) const {
    std::vector<DefReal> background_space(2), coordi(2);
    background_space.at(kXIndex) = grid_manager2d.k0DomainDx_.at(kXIndex);
    background_space.at(kYIndex) = grid_manager2d.k0DomainDx_.at(kYIndex);
    DefSFBitset bitset_temp;
    if (rank_id_ == 0) {
        for (const auto& i_point : ptr_geo2d->coordinate_origin_) {
            grid_manager2d.SFBitsetEncodingCoordi()
        }
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
void MpiManager::IniSendNReceivePartitionedGeoCoordi(
    const GridManager3D& grid_manager3d,
    const std::vector<DefSFBitset>& bitset_max,
    Geometry3DInterface* const ptr_geo3d) const {

}
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
