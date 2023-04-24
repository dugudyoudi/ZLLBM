//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
void MpiManager::StartupMpi(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_ranks_);
}
void MpiManager::SetMpiParameters() {
}
void MpiManager::PartiteGridBySpaceFillingCurves(
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
DefSizet MpiManager::SerializeData(const DefMap<DefUint>& map_nodes,
    std::unique_ptr<char[]>& buffer) const {
    DefSizet num_nodes = map_nodes.size();
    DefSizet buffer_size = (sizeof(DefSFBitset) + sizeof(DefUint)) * num_nodes;
    buffer = std::make_unique<char[]>(buffer_size);
    // store
    DefSizet position = 0;
    std::memcpy(buffer.get() + position, &buffer_size, sizeof(DefSizet));
    position += sizeof(DefSizet);
    for (auto& kv_pair : data) {
        // serialize the key 
        int key_size = kv_pair.first.size();
        std::memcpy(buffer.get() + position, &key_size, sizeof(int));
        position += sizeof(int); std::memcpy(buffer.get() + position, kv_pair.first.c_str(), key_size); position += key_size; // serialize the value int value_size = sizeof(int); auto value_ptr = std::make_unique<char[]>(value_size); std::memcpy(value_ptr.get(), &(kv_pair.second), value_size); std::memcpy(buffer.get() + position, &value_size, sizeof(int)); position += sizeof(int); std::memcpy(buffer.get() + position, value_ptr.get(), value_size); position += value_size;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif ENABLE_MPI
