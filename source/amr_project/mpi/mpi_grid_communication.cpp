//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_grid_communication.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is enabled.
* @date  2024-2-16
*/
#include "mpi/mpi_manager.h"
#include "grid/grid_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
void MpiManager::SerializeNodeInformation(const GridInfoInterface& grid_info) const {
    int key_size = sizeof(DefSFBitset);
    int node_info_size = grid_info.SizeOfGridNodeForMpiCommunication();
    DefAmrIndexUint i_level = grid_info.i_level_;
    DefSizet num_node_all = 0;
    DefSizet num_trunks, last_num_nodes, other_num_nodes;
    DefSizet count = 0;
    for (const auto& iter_rank : mpi_communication_inner_layers_.at(i_level)) {
        num_node_all = mpi_communication_inner_layers_.size();
        num_trunks = (node_info_size + key_size)*num_node_all/(0x7FFFFFFF - sizeof(int));
        if (num_trunks < 1) {  // buffer is enough to store information of all nodes
            other_num_nodes = num_node_all;
            last_num_nodes = 0;
        } else {
            last_num_nodes = num_node_all%(num_trunks - 1);
            other_num_nodes = (num_node_all - last_num_nodes)/(num_trunks - 1);
        }
        for (auto& iter_node : iter_rank) {
            count = 0;
            std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size);
            char* ptr_buffer = buffer.get();
            int position = 0;
            std::memcpy(ptr_buffer + position, &num_nodes, sizeof(int));
        }

    }

    // DefSizet num_node_for_communication = grid_info.
    // if  (sizeof(int) + set_nodes.size() *key_size > 0x7FFFFFFF) {
    //     LogManager::LogError("size of the buffer is greater than"
    //      " the maximum of int in MpiManager::IniSerializeTrackingNode in "
    //      + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    // }
    // int& buffer_size = *ptr_buffer_size;
    // buffer_size = sizeof(int) + key_size * num_nodes;
}

}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
