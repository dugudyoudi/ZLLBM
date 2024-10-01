//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager_debug.cpp
* @author Zhengliang Liu
* @brief functions used to debug mpi communication.
*/
#include <vector>
#include "mpi/mpi_manager.h"
#ifdef DEBUG_UNIT_TEST
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
void MpiManager::DebugMpiForAllGrids(const GridManagerInterface& grid_manager) const {
    for (const auto& iter : grid_manager.vec_ptr_grid_info_) {
        CheckMpiNodesCorrespondence(*iter.get());
        CheckMpiPeriodicCorrespondence(*iter.get());
    }
}
/**
* @brief      function to check if mpi outer layer only receives information in inner layer from the unique corresponding rank
* @param[in]  grid_info   reference to grid information.
*/
void MpiManager::CheckMpiNodesCorrespondence(const GridInfoInterface& grid_info) const {
    DefInt i_level = grid_info.i_level_;
    std::vector<bool> vec_receive_ranks(IdentifyRanksReceivingGridNode(i_level));
    std::vector<BufferSizeInfo> send_buffer_info, receive_buffer_info;
    std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send, vec_ptr_buffer_receive(num_of_ranks_);
    int node_info_size = 0;
    int node_buffer_size = static_cast<int>(sizeof(DefSFBitset));
    SendNReceiveGridNodeBufferSize(i_level, node_info_size,
        vec_receive_ranks, &send_buffer_info, &receive_buffer_info);
    // receive node in mpi outer layer via mpi communication
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (receive_buffer_info.at(i_rank).bool_exist_) {
            vec_vec_reqs_receive.push_back({});
            const int& num_chunks = receive_buffer_info.at(i_rank).num_chunks_;
            vec_vec_reqs_receive.back().resize(num_chunks);
            const int& buffer_size_rest = receive_buffer_info.at(i_rank).array_buffer_size_.at(0);
            DefSizet buffer_size_total = sizeof(int) + (num_chunks - 1)*buffer_size_rest
                + receive_buffer_info.at(i_rank).array_buffer_size_.at(1);
            vec_ptr_buffer_receive.at(i_rank) = std::make_unique<char[]>(buffer_size_total);
            std::unique_ptr<char[]>& buffer = vec_ptr_buffer_receive.at(i_rank);
            int position = 0;
            for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
                MPI_Irecv(buffer.get()+i_chunk*buffer_size_rest, buffer_size_rest, MPI_BYTE, i_rank,
                    i_chunk, MPI_COMM_WORLD, &vec_vec_reqs_receive.back().at(i_chunk));
            }
            int i_chunk_last = num_chunks - 1;
            MPI_Irecv(buffer.get()+ (num_chunks - 1)*buffer_size_rest,
                receive_buffer_info.at(i_rank).array_buffer_size_.at(1), MPI_BYTE, i_rank,
                i_chunk_last, MPI_COMM_WORLD, &vec_vec_reqs_receive.back().at(i_chunk_last));
        }
    }
    // send node in mpi inner layer via mpi communication
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (send_buffer_info.at(i_rank).bool_exist_) {
            vec_vec_reqs_send.push_back({});
            if (mpi_communication_inner_layers_.at(i_level).find(i_rank)
                !=mpi_communication_inner_layers_.at(i_level).end()) {
                const int& num_chunks = send_buffer_info.at(i_rank).num_chunks_;
                vec_vec_reqs_send.back().resize(num_chunks);
                int size_tmp;
                vec_ptr_buffer_send.push_back(SerializeNodeSFBitset(
                    mpi_communication_inner_layers_.at(i_level).at(i_rank), &size_tmp));
                const int& buffer_size_send = send_buffer_info.at(i_rank).array_buffer_size_.at(0);
                for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
                    MPI_Isend(vec_ptr_buffer_send.back().get()+i_chunk*buffer_size_send,
                        send_buffer_info.at(i_rank).array_buffer_size_.at(0),
                        MPI_BYTE, i_rank, i_chunk, MPI_COMM_WORLD, &vec_vec_reqs_send.back().at(i_chunk));
                }
                int i_chunk_last = num_chunks - 1;
                MPI_Isend(vec_ptr_buffer_send.back().get()+i_chunk_last*buffer_size_send,
                    send_buffer_info.at(i_rank).array_buffer_size_.at(1), MPI_BYTE, i_rank, i_chunk_last,
                    MPI_COMM_WORLD, &vec_vec_reqs_send.back().at(i_chunk_last));
            } else {
                LogManager::LogError("inner nodes will be sent to other ranks do not exist on"
                    + std::to_string(i_rank) + " in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
        }
    }

    int i_rev = 0;
    DefInt flag_receive_once = 1, flag_receive_more = 2;
    DefMap<DefInt> map_received;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (receive_buffer_info.at(i_rank).bool_exist_) {
            DefMap<DefInt> map_received_tmp;
            MPI_Waitall(static_cast<int>(vec_vec_reqs_receive.at(i_rev).size()),
                vec_vec_reqs_receive.at(i_rev).data(), MPI_STATUSES_IGNORE);
            DefSizet buffer_size = sizeof(int) + receive_buffer_info.at(i_rank).array_buffer_size_.at(1)
                + (receive_buffer_info.at(i_rank).num_chunks_ - 1)
                *receive_buffer_info.at(i_rank).array_buffer_size_.at(0);
            DeserializeNodeSFBitset(flag_receive_once, vec_ptr_buffer_receive.at(i_rank), &map_received_tmp);
            for (const auto& iter : map_received_tmp) {
                if (map_received.find(iter.first) == map_received.end()) {
                    map_received.insert({iter.first, flag_receive_once});
                } else {
                    map_received.at(iter.first) = flag_receive_more;
                    LogManager::LogError("receive more nodes from rank " + std::to_string(i_rank) + " in "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            ++i_rev;
        }
    }

//     if (rank_id_ == 2) {
//         std::ofstream file1("test0.txt");
//         for (auto iter : map_received) {
//             amrproject::SFBitsetAux2D aux2d;
//             std::array<DefAmrLUint, 2> indices;
//             aux2d.SFBitsetComputeIndices(iter.first, &indices);
//             file1  << indices[0] << " " << indices[1] << std::endl;

//         }
//         file1.close();

//        file1.open("test1.txt");
//         for (auto iter : mpi_communication_outer_layers_.at(1)) {
//             amrproject::SFBitsetAux2D aux2d;
//             std::array<DefAmrLUint, 2> indices;
//             aux2d.SFBitsetComputeIndices(iter.first, &indices);
//             file1 << indices[0] << " " << indices[1] << std::endl;
//         }
//         file1.close();
// }


    if (mpi_communication_outer_layers_.at(i_level).size() != map_received.size()) {
        LogManager::LogError("the number of received nodes is not the same as the number of nodes"
            " in mpi outer layer on rank " + std::to_string(rank_id_) + " at level " + std::to_string(i_level));
    }
    for (const auto& iter : mpi_communication_outer_layers_.at(i_level)) {
        if (map_received.find(iter.first) == map_received.end()) {
            LogManager::LogError("at least one of node in mpi outer layer is not received in on rank "
            + std::to_string(rank_id_) + " at level " + std::to_string(i_level));
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
/**
* @brief      function to check if ghost node are instantiated for mpi communication at periodic boundaries.
* @param[in]  grid_info   reference to grid information.
*/
void MpiManager::CheckMpiPeriodicCorrespondence(const GridInfoInterface& grid_info) const {
    DefInt dims = grid_info.GetPtrToParentGridManager()->k0GridDims_;
    std::vector<bool> periodic_min(dims, false), periodic_max(dims, false);
    grid_info.CheckIfPeriodicDomainRequired(dims, &periodic_min, &periodic_max);
    std::vector<DefSFBitset> domain_min_n_level(dims), domain_max_n_level(dims);
    std::vector<DefAmrLUint> indices_min = grid_info.GetPtrToParentGridManager()->GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = grid_info.GetPtrToParentGridManager()->GetMaxIndexOfBackgroundNodeArrAsVec();
    grid_info.ptr_sfbitset_aux_->GetMinAtGivenLevel(grid_info.i_level_, indices_min, &domain_min_n_level);
    grid_info.ptr_sfbitset_aux_->GetMaxAtGivenLevel(grid_info.i_level_, indices_max, &domain_max_n_level);
    std::vector<DefSFBitset> vec_in_region;
        //  if (i_rank == 2) {
        //     std::ofstream file1("test0.txt");
        // for (auto& iter : map_ghost_n_refinement) {
        //     amrproject::SFBitsetAux2D aux2d;
        //     std::array<DefAmrLUint, 2> indices;
        //     aux2d.SFBitsetComputeIndices(iter.first, &indices);
        //     file1  << indices[0] << " " << indices[1] << " " << iter.second << std::endl;

        // }        file1.close();
        // }

    for (int i = 0; i< grid_info.domain_boundary_min_.size(); ++i) {
        if (periodic_min.at(i)) {
            for (const auto& iter_node : grid_info.domain_boundary_min_.at(i)) {
                if ((grid_info.map_grid_node_.at(iter_node.first)->flag_status_
                    &NodeBitStatus::kNodeStatusMpiPartitionOuter_)
                    != NodeBitStatus::kNodeStatusMpiPartitionOuter_) {
                    grid_info.ptr_sfbitset_aux_->FindNodesInPeriodicRegionCenter(iter_node.first,
                        k0NumPartitionOuterLayers_, periodic_min, periodic_max,
                        domain_min_n_level, domain_max_n_level, &vec_in_region);
                    for (const auto& iter_region : vec_in_region) {
                        if (iter_region != grid_info.ptr_sfbitset_aux_->kInvalidSFbitset
                            && grid_info.map_grid_node_.find(iter_region) == grid_info.map_grid_node_.end()) {
                            std::vector<DefReal> coordinates(dims), coordinates2(dims);
                            if (dims == 2) {
                                grid_info.ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(
                                    iter_region, grid_info.grid_space_, &coordinates);
                                grid_info.ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(
                                    iter_region, grid_info.grid_space_, &coordinates2);
                                const std::array<DefReal, 2>& domain_min = dynamic_cast<SFBitsetAux2D*>(
                                    grid_info.ptr_sfbitset_aux_)->k0RealMin_;
                                coordinates[kXIndex] -= domain_min[kXIndex];
                                coordinates[kYIndex] -= domain_min[kYIndex];
                                coordinates2[kXIndex] -= domain_min[kXIndex];
                                coordinates2[kYIndex] -= domain_min[kYIndex];
                                LogManager::LogError("Node (" + std::to_string(coordinates2[kXIndex]) + ", "
                                    + std::to_string(coordinates2[kYIndex])
                                    + ") at periodic boundary is not found for node ("
                                    + std::to_string(coordinates[kXIndex]) + ", "
                                    + std::to_string(coordinates[kYIndex]) + ") at level "
                                    + std::to_string(grid_info.i_level_));
                            } else {
                                grid_info.ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(
                                    iter_region, grid_info.grid_space_, &coordinates);
                                grid_info.ptr_sfbitset_aux_->SFBitsetComputeCoordinateVir(
                                    iter_region, grid_info.grid_space_, &coordinates2);
                                const std::array<DefReal, 3>& domain_min = dynamic_cast<SFBitsetAux3D*>(
                                    grid_info.ptr_sfbitset_aux_)->k0RealMin_;
                                coordinates[kXIndex] -= domain_min[kXIndex];
                                coordinates[kYIndex] -= domain_min[kYIndex];
                                coordinates[kZIndex] -= domain_min[kZIndex];
                                coordinates2[kXIndex] -= domain_min[kXIndex];
                                coordinates2[kYIndex] -= domain_min[kYIndex];
                                coordinates2[kZIndex] -= domain_min[kZIndex];
                                LogManager::LogError("Node (" + std::to_string(coordinates2[kXIndex]) + ", "
                                    + std::to_string(coordinates2[kYIndex])  + ", "
                                    + std::to_string(coordinates2[kZIndex])
                                    + ") at periodic boundary is not found for node ("
                                    + std::to_string(coordinates[kXIndex]) + ", "+ std::to_string(coordinates[kYIndex])
                                    + ", "+ std::to_string(coordinates[kZIndex]) + ") at level "
                                    + std::to_string(grid_info.i_level_));
                            }
                        }
                    }
                }
            }
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
#endif  // DEBUG_UNIT_TEST
