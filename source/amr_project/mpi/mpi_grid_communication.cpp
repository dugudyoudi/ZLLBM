//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_grid_communication.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is enabled.
* @date  2024-2-26
*/
#include "mpi/mpi_manager.h"
#include "grid/grid_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to identify information of grid node from which ranks will be received.
 * @param[in] i_level refinement level corresponding to key of the container storing inner communication layers.
 */
std::vector<bool> MpiManager::IdentifyRanksReceivingGridNode(const DefAmrIndexUint i_level) const {
    std::vector<char> send_ranks(num_of_ranks_, 0), receive_ranks(num_of_ranks_, 0);
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (mpi_communication_inner_layers_.at(i_level).find(i_rank)
            != mpi_communication_inner_layers_.at(i_level).end()) {
            send_ranks.at(i_rank) = 1;
        }
    }
    MPI_Alltoall(send_ranks.data(), 1, MPI_CHAR,
        receive_ranks.data(), 1, MPI_CHAR, MPI_COMM_WORLD);
    std::vector<bool> recv_data(receive_ranks.size());
    std::copy(receive_ranks.begin(), receive_ranks.end(), recv_data.begin());
    return recv_data;
}
/**
 * @brief function to send grid node information to other ranks.
 * @param[in] i_level_receive refinement level of receiving information.
 * @param[in] rank_id_receive indicators for current rank from which ranks it will receive information.
 * @param[in] grid_info class storting grid node information on current rank.
 * @param[out] ptr_send_buffer_info pointer to information of buffer size will be sent.
 * @param[out] ptr_receive_buffer_info pointer to information of buffer size will be received.
 */
// Information sends in an periodic iteration pattern for all ranks.

void MpiManager::SendNReceiveGridNodeBufferSize(const DefAmrIndexUint i_level_receive,
    const std::vector<bool>& rank_id_receive, const GridInfoInterface& grid_info,
    std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
    std::vector<BufferSizeInfo>* const ptr_receive_buffer_info) const {
    int key_size = sizeof(DefSFBitset);
    int node_info_size = grid_info.SizeOfGridNodeInfoForMpiCommunication();
    DefAmrIndexUint i_level = grid_info.i_level_;
    if (i_level_receive != i_level) {
        LogManager::LogError("refinement level of receiving information is not equal to the "
            "refinement level current grid in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    int num_chunks;
    int node_buffer_size = key_size + node_info_size;
    int reqs_rank_count = 0;
    ptr_send_buffer_info->resize(num_of_ranks_);
    ptr_receive_buffer_info->resize(num_of_ranks_);
    for (int i = 1; i < num_of_ranks_; ++i) {
        DefSizet num_node_all = 0, last_num_nodes = 0, other_num_nodes = 0;
        int i_rank_send = (rank_id_ + i) % num_of_ranks_;
        int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
        bool bool_send_to_i_rank = false;
        if ((mpi_communication_inner_layers_.size() > i_level)
            && (mpi_communication_inner_layers_.at(i_level).find(i_rank_send)
            != mpi_communication_inner_layers_.at(i_level).end())) {
            bool_send_to_i_rank = true;
        }
        if (bool_send_to_i_rank) {
            ptr_send_buffer_info->at(i_rank_send).bool_exist_ = true;
            num_node_all = mpi_communication_inner_layers_.at(i_level).at(i_rank_send).size();
            num_chunks = static_cast<int>(node_buffer_size * num_node_all/
                ((std::numeric_limits<int>::max)() - sizeof(int)) + 1);
            if (num_chunks < 2) {  // buffer is enough to store information of all nodes
                other_num_nodes = 0;
                last_num_nodes = num_node_all;
            } else {
                last_num_nodes = num_node_all%(num_chunks - 1);
                other_num_nodes = (num_node_all - last_num_nodes)/(num_chunks - 1);
            }
        }

        // send and receive number of chunks if needed
        DefSizet buffer_size_send = 0;
        if (bool_send_to_i_rank) {
            ptr_send_buffer_info->at(i_rank_send).num_chunks_ = num_chunks;
            MPI_Send(&num_chunks, 1, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
        }
        if (rank_id_receive.at(i_rank_receive)) {
            ptr_receive_buffer_info->at(i_rank_receive).bool_exist_ = true;
            MPI_Recv(&ptr_receive_buffer_info->at(i_rank_receive).num_chunks_,
                1, MPI_INT, i_rank_receive, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // send and receive buffer size of each chunk if needed
        if (bool_send_to_i_rank) {
            buffer_size_send = node_buffer_size * other_num_nodes;
            if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.at(0) = static_cast<int>(buffer_size_send);
            } else {
                LogManager::LogError("size of the buffer is greater than the maximum value of int in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            buffer_size_send = node_buffer_size * last_num_nodes;
            if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.at(1) = static_cast<int>(buffer_size_send);
            } else {
                LogManager::LogError("size of the buffer is greater than the maximum value of int in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            MPI_Send(ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.data(),
                2, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
        }
        if (rank_id_receive.at(i_rank_receive)) {
            MPI_Recv(ptr_receive_buffer_info->at(i_rank_receive).array_buffer_size_.data(),
                2, MPI_INT, i_rank_receive, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}
/**
 * @brief function to calculate size of buffer storing grid nodes will be sent to a given rank.
 * @param[in] rank_send rank id will be sent to.
 * @param[in] grid_info class storting grid node information on current rank.
 * @return size of buffer storing grid nodes.
 */
DefSizet MpiManager::CalculateBufferSizeForGridNode(const int rank_send,
    const GridInfoInterface& grid_info) const {
    int key_size = sizeof(DefSFBitset);
    int node_info_size = grid_info.SizeOfGridNodeInfoForMpiCommunication();
    int node_buffer_size = key_size + node_info_size;
    return node_buffer_size*mpi_communication_inner_layers_.at(grid_info.i_level_).at(rank_send).size();
}
/**
 * @brief function to send grid node information to a given rank.
 * @param[in] rank_send rank id will be sent to.
 * @param[in] send_buffer_info information of buffer size will be sent.
 * @param[out] ptr_grid_info pointer to class storting grid node information on current rank.
 * @param[out] ptr_buffer_send pointer to buffer storing nodes need to be sent
 * @param[out] ptr_reqs_send pointer to status of mpi sending for each splitted chunks of the target rank.
 * @note information will only be sent to ranks existing in mpi_communication_inner_layers_.
 */
// Information sends in an periodic iteration pattern for all ranks.
void MpiManager::SendGridNodeInformation(const int rank_send, const BufferSizeInfo& send_buffer_info,
    GridInfoInterface* const ptr_grid_info, char* const ptr_buffer_send,
    std::vector<MPI_Request>* const ptr_reqs_send) {
    DefAmrIndexUint i_level = ptr_grid_info->i_level_;

    if (mpi_communication_inner_layers_.at(i_level).find(rank_send)
        !=mpi_communication_inner_layers_.at(i_level).end()) {
        if (!send_buffer_info.bool_exist_) {
            LogManager::LogError("sending buffer size information is not initialized properly in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        const int& num_chunks = send_buffer_info.num_chunks_;
        ptr_reqs_send->resize(num_chunks);
        // copy grid node information to buffer
        DefSizet position = 0;
        ptr_grid_info->CopyNodeInfoToBuffer(mpi_communication_inner_layers_.at(i_level).at(rank_send), ptr_buffer_send);

        // send grid node information chunk by chunk via non-blocking communication
        const int& buffer_size_send = send_buffer_info.array_buffer_size_.at(0);
        for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
            MPI_Isend(ptr_buffer_send+i_chunk*buffer_size_send, send_buffer_info.array_buffer_size_.at(0),
                MPI_BYTE, rank_send, i_chunk, MPI_COMM_WORLD, &ptr_reqs_send->at(i_chunk));
        }
        int i_chunk_last = num_chunks - 1;
        MPI_Isend(ptr_buffer_send+i_chunk_last*buffer_size_send,
            send_buffer_info.array_buffer_size_.at(1), MPI_BYTE, rank_send, i_chunk_last,
            MPI_COMM_WORLD, &ptr_reqs_send->at(i_chunk_last));
    }
}
/**
 * @brief function to receive grid node information from other ranks.
 * @param[in] rank_receive rank id wil be received from.
 * @param[in] receive_buffer_info  information of buffer will be received.
 * @param[out] ptr_grid_info pointer to class storting grid node information.
 * @param[out] ptr_reqs_receive pointer to status of mpi sending for each splitted pieces of target ranks.
 */
std::unique_ptr<char[]> MpiManager::ReceiveGridNodeInformation(const int rank_receive, const int node_info_size,
    const BufferSizeInfo& receive_buffer_info, std::vector<MPI_Request>* const ptr_reqs_receive) const {
    if (receive_buffer_info.bool_exist_) {
        int node_buffer_size = sizeof(DefSFBitset) + node_info_size;
        const int& num_chunks = receive_buffer_info.num_chunks_;
        ptr_reqs_receive->resize(num_chunks);
        const int& buffer_size_rest = receive_buffer_info.array_buffer_size_.at(0);
        DefSizet buffer_size_total = (num_chunks - 1)*buffer_size_rest
            + receive_buffer_info.array_buffer_size_.at(1);
        std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size_total);
        int position = 0;
        for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
            MPI_Irecv(buffer.get()+i_chunk*buffer_size_rest, buffer_size_rest, MPI_BYTE, rank_receive,
                i_chunk, MPI_COMM_WORLD, &ptr_reqs_receive->at(i_chunk));
        }
        int i_chunk_last = num_chunks - 1;
        MPI_Irecv(buffer.get()+ (num_chunks - 1)*buffer_size_rest,
            receive_buffer_info.array_buffer_size_.at(1), MPI_BYTE, rank_receive,
            i_chunk_last, MPI_COMM_WORLD, &ptr_reqs_receive->at(i_chunk_last));
        return buffer;
    } else {
        std::unique_ptr<char[]> buffer = std::make_unique<char[]>(0);
        return buffer;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
