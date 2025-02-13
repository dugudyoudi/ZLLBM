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
    const DefInt i_level = grid_info.GetGridLevel();
    std::vector<bool> vec_receive_ranks(IdentifyRanksReceivingGridNode(i_level));
    std::vector<BufferSizeInfo> send_buffer_info, receive_buffer_info;
    std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send, vec_ptr_buffer_receive(num_of_ranks_);
    int node_info_size = 0;
    int node_buffer_size = static_cast<int>(sizeof(DefSFBitset));
    if (static_cast<DefInt>(mpi_communication_inner_layers_.size()) > i_level) {
        SendNReceiveGridNodeBufferSize(node_info_size, i_level, mpi_communication_inner_layers_.at(i_level),
            &send_buffer_info, &receive_buffer_info);
    } else {
        LogManager::LogError("number of mpi inner layers is less the given grid level on" + std::to_string(rank_id_));
    }

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
                vec_ptr_buffer_send.emplace_back(SerializeNodeSFBitset(
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

    // received nodes may be more than nodes in mpi outer layer in coarse to fine layers
    if (mpi_communication_outer_layers_.at(i_level).size() > map_received.size()) {
        LogManager::LogError("the number of received nodes is less than the number of nodes"
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
    const std::vector<DefAmrLUint>
        indices_min = grid_info.GetPtrToParentGridManager()->GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = grid_info.GetPtrToParentGridManager()->GetMaxIndexOfBackgroundNodeArrAsVec();
    const DefInt i_level = grid_info.GetGridLevel();
    grid_info.GetPtrSFBitsetAux()->GetMinAtGivenLevel(i_level, indices_min, &domain_min_n_level);
    grid_info.GetPtrSFBitsetAux()->GetMaxAtGivenLevel(i_level, indices_max, &domain_max_n_level);
    DefSFBitset sfbitset_max = ~0, sfbitset_min = ~0, sfbitset_tmp = ~0;
    std::vector<DefReal> coordinates(dims), coordinates2(dims);
    std::vector<DefSFBitset> vec_in_region;
    const std::array<DefSFBitset, 2>& take_xref = grid_info.GetPtrSFBitsetAux()->GetTakeXRef(),
        take_yref = grid_info.GetPtrSFBitsetAux()->GetTakeYRef(),
        take_zref = grid_info.GetPtrSFBitsetAux()->GetTakeZRef();
    for (int i = 0; i < dims; ++i) {
        if (periodic_min.at(i)) {
            std::vector<DefInt> neg_layers(dims, k0NumPartitionOuterLayers_),
                pos_layers(dims, k0NumPartitionOuterLayers_);
            pos_layers.at(i) = k0NumPartitionOuterLayers_ - 1;
            for (const auto& iter_node : grid_info.domain_boundary_min_.at(i)) {
                if (i == kXIndex) {
                    sfbitset_max = (take_xref[SFBitsetAuxInterface::kRefOthers_]&iter_node.first)|domain_max_n_level[i];
                } else if (i == kYIndex) {
                    sfbitset_max = (take_yref[SFBitsetAuxInterface::kRefOthers_]&iter_node.first)|domain_max_n_level[i];
                } else if (i == kZIndex) {
                    sfbitset_max = (take_zref[SFBitsetAuxInterface::kRefOthers_]&iter_node.first)|domain_max_n_level[i];
                }
                bool exist_inner_mpi_node = false;
                grid_info.GetPtrSFBitsetAux()->FindNodesInPeriodicRegionCenter(iter_node.first, neg_layers, pos_layers,
                    periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
                for (const auto iter_region : vec_in_region) {
                    if (iter_region != SFBitsetAuxInterface::kInvalidSFbitset
                        && grid_info.map_grid_node_.find(iter_region) != grid_info.map_grid_node_.end()
                        && ((grid_info.map_grid_node_.at(iter_region)->flag_status_
                        &NodeBitStatus::kNodeStatusMpiPartitionOuter_)
                        !=NodeBitStatus::kNodeStatusMpiPartitionOuter_)) {
                            exist_inner_mpi_node = true;
                            break;
                        }
                }
                if (exist_inner_mpi_node && grid_info.domain_boundary_max_.at(i).find(sfbitset_max)
                    == grid_info.domain_boundary_max_.at(i).end()) {
                    grid_info.GetPtrSFBitsetAux()->SFBitsetComputeCoordinateVir(
                        iter_node.first, grid_info.GetGridSpace(), &coordinates);
                    grid_info.GetPtrSFBitsetAux()->SFBitsetComputeCoordinateVir(
                        sfbitset_max, grid_info.GetGridSpace(), &coordinates2);
                    if (dims == 2) {
                        LogManager::LogError("Can not find counterpart node ("
                            + std::to_string(coordinates2[kXIndex]) + ", "
                            + std::to_string(coordinates2[kYIndex]) + ") in direction " + std::to_string(i)
                            + " on the maximum boundary for node ("
                            + std::to_string(coordinates[kXIndex]) + ", "
                            + std::to_string(coordinates[kYIndex]) + ") on the minimum boundary at level "
                            + std::to_string(i_level) + " on rank " + std::to_string(rank_id_));
                    } else if (dims == 3) {
                        LogManager::LogError("Can not find counterpart node ("
                            + std::to_string(coordinates2[kXIndex]) + ", "
                            + std::to_string(coordinates2[kYIndex]) + ", "
                            + std::to_string(coordinates2[kZIndex]) + ") in direction " + std::to_string(i)
                            + " on the maximum boundary for node ("
                            + std::to_string(coordinates[kXIndex]) + ", "+ std::to_string(coordinates[kYIndex])
                            + ", "+ std::to_string(coordinates[kZIndex]) + ") on the minimum boundary at level "
                            + std::to_string(i_level) + " on rank " + std::to_string(rank_id_));
                    }
                }
            }
        }
        if (periodic_max.at(i)) {
            std::vector<DefInt> neg_layers(dims, k0NumPartitionOuterLayers_),
                pos_layers(dims, k0NumPartitionOuterLayers_);
            neg_layers.at(i) = k0NumPartitionOuterLayers_ - 1;
            for (const auto& iter_node : grid_info.domain_boundary_max_.at(i)) {
                if (i == kXIndex) {
                    sfbitset_min = (take_xref[SFBitsetAuxInterface::kRefOthers_]&iter_node.first)|domain_min_n_level[i];
                } else if (i == kYIndex) {
                    sfbitset_min = (take_yref[SFBitsetAuxInterface::kRefOthers_]&iter_node.first)|domain_min_n_level[i];
                } else if (i == kZIndex) {
                    sfbitset_min = (take_zref[SFBitsetAuxInterface::kRefOthers_]&iter_node.first)|domain_min_n_level[i];
                }
                bool exist_inner_mpi_node = false;
                grid_info.GetPtrSFBitsetAux()->FindNodesInPeriodicRegionCenter(iter_node.first, neg_layers, pos_layers,
                    periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
                for (const auto iter_region : vec_in_region) {
                    if (iter_region != SFBitsetAuxInterface::kInvalidSFbitset
                        && grid_info.map_grid_node_.find(iter_region) != grid_info.map_grid_node_.end()
                        && ((grid_info.map_grid_node_.at(iter_region)->flag_status_
                        &NodeBitStatus::kNodeStatusMpiPartitionOuter_)
                        !=NodeBitStatus::kNodeStatusMpiPartitionOuter_)) {
                            exist_inner_mpi_node = true;
                            break;
                        }
                }
                if (exist_inner_mpi_node && grid_info.domain_boundary_min_.at(i).find(sfbitset_min)
                    == grid_info.domain_boundary_min_.at(i).end()) {
                    grid_info.GetPtrSFBitsetAux()->SFBitsetComputeCoordinateVir(
                        iter_node.first, grid_info.GetGridSpace(), &coordinates);
                    grid_info.GetPtrSFBitsetAux()->SFBitsetComputeCoordinateVir(
                        sfbitset_min, grid_info.GetGridSpace(), &coordinates2);
                    if (dims == 2) {
                        LogManager::LogError("Can not find counterpart node ("
                            + std::to_string(coordinates2[kXIndex]) + ", "
                            + std::to_string(coordinates2[kYIndex]) + ") in direction " + std::to_string(i)
                            + " on the minimum boundary for node ("
                            + std::to_string(coordinates[kXIndex]) + ", "
                            + std::to_string(coordinates[kYIndex]) + ") on the maximum boundary at level "
                            + std::to_string(i_level) + " on rank " + std::to_string(rank_id_));
                    } else if (dims == 3) {
                        LogManager::LogError("Can not find counterpart node ("
                            + std::to_string(coordinates2[kXIndex]) + ", "
                            + std::to_string(coordinates2[kYIndex]) + ", "
                            + std::to_string(coordinates2[kZIndex]) + ") in direction " + std::to_string(i)
                            + " on the minimum boundary for node ("
                            + std::to_string(coordinates[kXIndex]) + ", "+ std::to_string(coordinates[kYIndex])
                            + ", "+ std::to_string(coordinates[kZIndex]) + ") on the maximum boundary at level "
                            + std::to_string(i_level) + " on rank " + std::to_string(rank_id_));
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
