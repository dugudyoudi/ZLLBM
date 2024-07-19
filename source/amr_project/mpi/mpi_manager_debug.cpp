//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager_debug.cpp
* @author Zhengliang Liu
* @brief functions used to debug mpi communication.
*/
#include <vector>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
* @brief      function to check if mpi outer layer only receives information in inner layer from the unique corresponding rank
* @param[in]  grid_info   reference to grid information.
*/
void MpiManager::CheckMpiNodesCorrespondence(const GridInfoInterface& grid_info) const {
    DefAmrIndexUint i_level = grid_info.i_level_;
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
    int error_flag;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (send_buffer_info.at(i_rank).bool_exist_) {
            vec_vec_reqs_send.push_back({});
            if (mpi_communication_inner_layers_.at(i_level).find(i_rank)
                !=mpi_communication_inner_layers_.at(i_level).end()) {
                const int& num_chunks = send_buffer_info.at(i_rank).num_chunks_;
                vec_vec_reqs_send.back().resize(num_chunks);
                int size_tmp;
                vec_ptr_buffer_send.push_back(SerializeNodeSFBitset(
                    mpi_communication_inner_layers_.at(i_level).at(i_rank), &size_tmp, &error_flag));
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
    DefAmrIndexUint flag_receive_once = 1, flag_receive_more = 2;
    DefMap<DefAmrIndexUint> map_received;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (receive_buffer_info.at(i_rank).bool_exist_) {
            DefMap<DefAmrIndexUint> map_received_tmp;
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

    // if (rank_id_ == 0) {
    //     std::ofstream file1("test0.txt");
    //     for (auto iter : mpi_communication_inner_layers_.at(1).at(1)) {
    //         amrproject::SFBitsetAux2D aux2d;
    //         std::array<DefAmrIndexLUint, 2> indices;
    //         aux2d.SFBitsetComputeIndices(iter.first, &indices);
    //         file1  << indices[0] << " " << indices[1] << std::endl;

    //     }
    //     file1.close();
    // } else {
    //     std::ofstream file1("test1.txt");
    //     for (auto iter : mpi_communication_outer_layers_.at(1)) {
    //         amrproject::SFBitsetAux2D aux2d;
    //         std::array<DefAmrIndexLUint, 2> indices;
    //         aux2d.SFBitsetComputeIndices(iter.first, &indices);
    //         file1 << indices[0] << " " << indices[1] << std::endl;
    //     }
    //     file1.close();
    // }

    if (mpi_communication_outer_layers_.at(i_level).size() != map_received.size()) {
        LogManager::LogError("the number of received nodes is not the same as the number of nodes"
            " in mpi outer layer in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    for (const auto& iter : mpi_communication_outer_layers_.at(i_level)) {
        if (map_received.find(iter.first) == map_received.end()) {
            LogManager::LogError("at least one of node in mpi outer layer is not received in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
