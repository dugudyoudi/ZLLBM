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
std::vector<bool> MpiManager::IdentifyRanksReceivingGridNode(const DefInt i_level) const {
    std::vector<int> send_ranks(num_of_ranks_, 0), receive_ranks(num_of_ranks_, 0);
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (mpi_communication_inner_layers_.at(i_level).find(i_rank)
            != mpi_communication_inner_layers_.at(i_level).end()
            && mpi_communication_inner_layers_.at(i_level).at(i_rank).size() > 0) {
            send_ranks.at(i_rank) = 1;
        }
    }
    MPI_Alltoall(send_ranks.data(), 1, MPI_INT,
        receive_ranks.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<bool> recv_data(receive_ranks.size());
    std::copy(receive_ranks.begin(), receive_ranks.end(), recv_data.begin());
    return recv_data;
}
/**
 * @brief function to send and receive buffer size information.
 * @param[in] node_info_size size of each node information will be sent.
 * @param[in] i_level refinement level.
 * @param[in] map_send_nodes  spacing filling codes of nodes will be sent.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[out] grid_info class storting grid information.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_send_buffer_info pointer to information of buffer size will be sent.
 * @param[out] ptr_receive_buffer_info pointer to information of buffer size will be received.
 */
// Information sends in an periodic iteration pattern for all ranks.
void MpiManager::SendNReceiveGridNodeBufferSize(const int node_info_size, const DefInt i_level,
    const std::map<int, DefMap<DefInt>>& map_send_nodes,
    std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
    std::vector<BufferSizeInfo>* const ptr_receive_buffer_info) const {
    int key_size = sizeof(DefSFBitset);
    int node_buffer_size = key_size + node_info_size;
    std::vector<MPI_Request> reqs_receive(num_of_ranks_);
    ptr_send_buffer_info->resize(num_of_ranks_);
    ptr_receive_buffer_info->resize(num_of_ranks_);
    for (int i = 1; i < num_of_ranks_; ++i) {
        // send and receive number of chunks if needed
        int i_rank_send = (rank_id_ + i) % num_of_ranks_;
        int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
        bool bool_send_to_i_rank = (map_send_nodes.find(i_rank_send)
            != map_send_nodes.end());
        int num_node_all = 0, last_num_nodes = 0, other_num_nodes = 0;
        if (bool_send_to_i_rank) {
            num_node_all = static_cast<int>(map_send_nodes.at(i_rank_send).size());
            int num_chunks = static_cast<int>(node_buffer_size * num_node_all
                /((std::numeric_limits<int>::max)() - sizeof(int)))  + 1;
            ptr_send_buffer_info->at(i_rank_send).bool_exist_ = true;
            ptr_send_buffer_info->at(i_rank_send).num_chunks_ = num_chunks;
            if (num_chunks < 2) {  // buffer is enough to store information of all nodes
                other_num_nodes = 0;
                last_num_nodes = num_node_all;
            } else {
                last_num_nodes = num_node_all%(num_chunks - 1);
                other_num_nodes = (num_node_all - last_num_nodes)/(num_chunks - 1);
            }
            ptr_send_buffer_info->at(i_rank_send).num_chunks_ = num_chunks;
            MPI_Send(&num_chunks, 1, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
        } else {
            ptr_send_buffer_info->at(i_rank_send).bool_exist_ = false;
            ptr_send_buffer_info->at(i_rank_send).num_chunks_ = 0;
            MPI_Send(&ptr_send_buffer_info->at(i_rank_send).num_chunks_,
                1, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
        }

        MPI_Recv(&ptr_receive_buffer_info->at(i_rank_receive).num_chunks_,
            1, MPI_INT, i_rank_receive, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (ptr_receive_buffer_info->at(i_rank_receive).num_chunks_ == 0) {
            ptr_receive_buffer_info->at(i_rank_receive).bool_exist_ = false;
        } else {
            ptr_receive_buffer_info->at(i_rank_receive).bool_exist_ = true;
        }

        // send and receive buffer size of each chunk if needed
        if (bool_send_to_i_rank) {
            DefSizet buffer_size_send = node_buffer_size * other_num_nodes;
            if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.at(0) = static_cast<int>(buffer_size_send);
            } else {
                LogManager::LogError("size of the buffer is greater than the maximum value of int");
            }
            buffer_size_send = node_buffer_size * last_num_nodes;
            if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.at(1) = static_cast<int>(buffer_size_send);
            } else {
                LogManager::LogError("size of the buffer is greater than the maximum value of int");
            }
            MPI_Send(ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.data(),
                2, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
        }
        if (ptr_receive_buffer_info->at(i_rank_receive).bool_exist_) {
            MPI_Recv(ptr_receive_buffer_info->at(i_rank_receive).array_buffer_size_.data(),
                2, MPI_INT, i_rank_receive, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}
/**
 * @brief function to calculate size of buffer storing grid nodes will be sent to a given rank.
 * @param[in] num_nodes number of nodes will be sent.
 * @param[in] grid_info class storting grid node information on current rank.
 * @return size of buffer storing grid nodes.
 */
DefSizet MpiManager::CalculateBufferSizeForGridNode(
    const DefSizet num_nodes, const GridInfoInterface& grid_info) const {
    int key_size = sizeof(DefSFBitset);
    int node_info_size = grid_info.GetSizeOfGridNodeInfoForMpiCommunication();
    int node_buffer_size = key_size + node_info_size;
    return node_buffer_size*num_nodes;
}
/**
 * @brief function to send grid node information to a given rank.
 * @param[in] rank_send rank id will be sent to.
 * @param[in] send_buffer_info information of buffer size will be sent.
 * @param[in] nodes_to_send space filling codes of nodes need will be sent.
 * @param[out] ptr_grid_info pointer to class storting grid node information on current rank.
 * @param[out] ptr_buffer_send pointer to buffer storing nodes need to be sent
 * @param[out] ptr_reqs_send pointer to status of mpi sending for each splitted chunks of the target rank.
 * @return 0 run successfully, otherwise some error occurred in this function.
 * @note information will only be sent of nodes existing in nodes_to_send.
 */
// Information sends in a periodic iteration pattern for all ranks.
void MpiManager::SendGridNodeInformation(const int rank_send, const BufferSizeInfo& send_buffer_info,
    const std::map<int, DefMap<DefInt>>& nodes_to_send,
    const std::function<void(const DefMap<DefInt>& , char* const)> func_copy_node_to_buffer,
    GridInfoInterface* const ptr_grid_info, char* const ptr_buffer_send,
    std::vector<MPI_Request>* const ptr_reqs_send) const {
    if (nodes_to_send.find(rank_send) !=nodes_to_send.end()) {
        if (!send_buffer_info.bool_exist_) {
            LogManager::LogError("sending buffer size information is not initialized properly");
        }
        const int& num_chunks = send_buffer_info.num_chunks_;
        ptr_reqs_send->resize(num_chunks);
        // copy grid node information to buffer
        func_copy_node_to_buffer(nodes_to_send.at(rank_send), ptr_buffer_send);

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
/**
 * @brief function to send and receive grid nodes via blocking mpi communication.
 * @param[in] node_info_size size of each node information will be sent.
 * @param[in] map_send_nodes  spacing filling codes of nodes will be sent.
 * @param[in] send_buffer_info buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] func_write_buffer function to write grid node infomation to buffer.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive  pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_vec_ptr_buffer_send pointer to buffer storing grid node information for sending.
 * @param[out] ptr_vec_ptr_buffer_receive pointer to buffer storing received grid node information.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 */
void MpiManager::NonBlockingSendNReceiveGridNode(const int node_info_size,
    const std::map<int, DefMap<DefInt>>& map_send_nodes,
    const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::function<void(const GridNode& node_ref, char* const)>& func_write_buffer,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive,
    GridInfoInterface* const ptr_grid_info) const {
    // send node information
    std::function<void(const DefMap<DefInt>& , char* const)> func_copy_node_to_buffer =
        [func_write_buffer, &ptr_grid_info](const DefMap<DefInt>& map_send, char* const ptr_buffer) {
        ptr_grid_info->CopyNodeInfoToBuffer(func_write_buffer, map_send, ptr_buffer);
    };
    int node_buffer_size = sizeof(DefSFBitset) + node_info_size;

    for (int i = 1; i < num_of_ranks_; ++i) {
        int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
        if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
            ptr_vec_vec_reqs_receive->push_back({});
            ptr_vec_ptr_buffer_receive->emplace_back(ReceiveGridNodeInformation(
                i_rank_receive, node_info_size,
                receive_buffer_info.at(i_rank_receive), &ptr_vec_vec_reqs_receive->back()));
        }
    }

    for (int i = 1; i < num_of_ranks_; ++i) {
        int i_rank_send = (rank_id_ + i) % num_of_ranks_;
        if (send_buffer_info.at(i_rank_send).bool_exist_) {
            ptr_vec_vec_reqs_send->push_back({});
            if (map_send_nodes.find(i_rank_send) == map_send_nodes.end()) {
                LogManager::LogError("nodes need to be sent does not exist for sending.");
            }
            ptr_vec_ptr_buffer_send->emplace_back(std::make_unique<char[]>(
                node_buffer_size*map_send_nodes.at(i_rank_send).size()));

            SendGridNodeInformation(i_rank_send, send_buffer_info.at(i_rank_send),
                map_send_nodes, func_copy_node_to_buffer,
                ptr_grid_info, ptr_vec_ptr_buffer_send->back().get(), &ptr_vec_vec_reqs_send->back());
        }
    }
}

/**
 * @brief function to send and receive grid node information from other ranks.
 * @param[out] ptr_send_buffer_info pointer to buffer size information of sending.
 * @param[out] ptr_receive_buffer_info pointer to buffer size information of receiving.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive  pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_vec_ptr_buffer_send pointer to buffer storing grid node information for sending.
 * @param[out] ptr_vec_ptr_buffer_receive pointer to buffer storing received grid node information.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 */
void MpiManager::SendAndReceiveGridNodesOnAllMpiLayers(
    std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
    std::vector<BufferSizeInfo>* const ptr_receive_buffer_info,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive,
    GridInfoInterface* ptr_grid_info) const {
    ptr_vec_vec_reqs_send->clear();
    ptr_vec_ptr_buffer_send->clear();
    ptr_vec_vec_reqs_receive->clear();
    ptr_vec_ptr_buffer_receive->clear();
    const DefInt i_level = ptr_grid_info->GetGridLevel();
    ptr_grid_info->ComputeInfoInMpiLayers(mpi_communication_inner_layers_.at(i_level),
        mpi_communication_outer_layers_.at(i_level));
    std::vector<bool> vec_receive_ranks(IdentifyRanksReceivingGridNode(i_level));

    if (mpi_communication_inner_layers_.size() > i_level) {
        SendNReceiveGridNodeBufferSize(ptr_grid_info->GetSizeOfGridNodeInfoForMpiCommunication(),
            i_level, mpi_communication_inner_layers_.at(i_level), ptr_send_buffer_info, ptr_receive_buffer_info);
    } else {
        LogManager::LogError("number of mpi inner layers is less the given grid level on" + std::to_string(rank_id_));
    }

    std::function<void(const GridNode&, char* const)> func_copy_a_node_to_buffer =
        [](const GridNode& node_ref, char* const ptr_node_buffer) {
        node_ref.CopyANodeToBufferForMpi(ptr_node_buffer);
    };
    NonBlockingSendNReceiveGridNode(ptr_grid_info->GetSizeOfGridNodeInfoForMpiCommunication(),
        mpi_communication_inner_layers_.at(i_level),
        *ptr_send_buffer_info, *ptr_receive_buffer_info, func_copy_a_node_to_buffer,
        ptr_vec_vec_reqs_send, ptr_vec_vec_reqs_receive,
        ptr_vec_ptr_buffer_send, ptr_vec_ptr_buffer_receive, ptr_grid_info);
}
/**
 * @brief function to wait for mpi communication and read grid node information from received buffer.
 * @param[in] send_buffer_info  buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] vec_ptr_buffer_receive  buffer storing received grid node information.
 * @param[in] func_read_a_node_from_buffer function to read grid node infomation from buffer.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive  pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 */
void MpiManager::WaitAndReadGridNodesFromBuffer(const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
    const std::function<void(const char*,  GridNode* const)>& func_read_a_node_from_buffer,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    GridInfoInterface* const ptr_grid_info) const {
    int i_rev = 0;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        int i_rank_receive = (rank_id_ - i_rank + num_of_ranks_)% num_of_ranks_;
        if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
            MPI_Waitall(static_cast<int>(ptr_vec_vec_reqs_receive->at(i_rev).size()),
                ptr_vec_vec_reqs_receive->at(i_rev).data(), MPI_STATUSES_IGNORE);
            DefSizet buffer_size = receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(1)
                + (receive_buffer_info.at(i_rank_receive).num_chunks_ - 1)
                *receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(0);
            ptr_grid_info->ReadNodeInfoFromBuffer(func_read_a_node_from_buffer,
                buffer_size, vec_ptr_buffer_receive.at(i_rev));
            ++i_rev;
        }
    }
    int i_send = 0;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        int i_rank_send = (rank_id_ + i_rank) % num_of_ranks_;
        if (send_buffer_info.at(i_rank_send).bool_exist_) {
            MPI_Waitall(static_cast<int>(ptr_vec_vec_reqs_send->at(i_send).size()),
                ptr_vec_vec_reqs_send->at(i_send).data(), MPI_STATUSES_IGNORE);
            ++i_send;
        }
    }
}
/**
 * @brief function to send and receive space filling code requested by current rank.
 * @param[in] num_node_request_current number of nodes requested by current rank.
 * @param[in] num_node_request_others number of nodes requested by other ranks.
 * @param[in] requested_nodes space filling code of nodes requested by current rank.
 * @param[out] ptr_vec_map_nodes_need_send space filling code of nodes need to be sent to other ranks.
 * @return 0 run successfully, otherwise some error occurred in this function.
 */
int MpiManager::SendNReceiveRequestNodesSFbitset(const std::vector<int>& num_node_request_current,
    const std::vector<int>& num_node_request_others, const std::vector<DefMap<DefInt>> &requested_nodes,
    std::map<int, DefMap<DefInt>>* const ptr_map_nodes_need_send) const {
    std::vector<MPI_Request> reqs_send(num_of_ranks_);
    for (int i = 1; i < num_of_ranks_; ++i) {
        int i_rank_send = (rank_id_ + i) % num_of_ranks_;
        int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
        if (num_node_request_current.at(i_rank_send) > 0) {
            int buffer_size_send = 0;
            std::unique_ptr<char[]> buffer_request;
            buffer_request = SerializeNodeSFBitset(requested_nodes.at(i_rank_send), &buffer_size_send);
            MPI_Send(&buffer_size_send, 1, MPI_INT, i_rank_send, 0, MPI_COMM_WORLD);
            MPI_Isend(buffer_request.get(), buffer_size_send, MPI_BYTE, i_rank_send,
                1, MPI_COMM_WORLD, &reqs_send[i_rank_send]);
        }
        if (num_node_request_others.at(i_rank_receive) > 0) {
            int buffer_size_receive = 0;
            MPI_Recv(&buffer_size_receive, 1, MPI_INT, i_rank_receive, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::unique_ptr<char[]> ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
            MPI_Recv(ptr_buffer.get(), buffer_size_receive, MPI_BYTE, i_rank_receive,
                1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ptr_map_nodes_need_send->insert(std::make_pair(i_rank_receive, DefMap<DefInt>()));
            DeserializeNodeSFBitset(0, ptr_buffer, &ptr_map_nodes_need_send->at(i_rank_receive));
        }
    }
    return 0;
}
/**
 * @brief function to send and receive which need nodes are needed for interpolation
 * @param[in] num_node_request_current number of nodes requested by current rank for interpolation.
 * @param[in] num_node_request_others number of nodes requested from other ranks for interpolation.
 * @param[in] requested_nodes space filling code of nodes requested from other ranks for interpolation.
 * @param[out] ptr_vec_map_nodes_need_send space filling code of nodes need to be sent to other ranks.
 * @return 0 run successfully, otherwise some error occurred in this function.
 */
void MpiManager::SendNReceiveSFbitsetForInterpolation(const DefInt i_level,
    const SFBitsetAuxInterface& sfbitset_aux, const DefMap<std::unique_ptr<GridNode>>& map_nodes_outer_layer,
    std::vector<int>* const ptr_vec_num_recv, std::map<int, DefMap<DefInt>>* const ptr_node_inner_layers) const {
    DefSFCodeToUint code_background;
    std::vector<DefMap<DefInt>> requested_nodes(num_of_ranks_);
    int which_rank;
    for (const auto& iter_node : map_nodes_outer_layer) {
        code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first).to_ullong();
        which_rank = CheckNodeInWhichRank(code_background);
        if (which_rank != rank_id_) {
            if (which_rank >= num_of_ranks_) {
                LogManager::LogError("Requested rank id (" + std::to_string(which_rank)
                    + ") shoule be less than the number of ranks (" + std::to_string(num_of_ranks_) +")");
            } else {
                requested_nodes.at(which_rank).insert({ iter_node.first, 0 });
            }
        } else {
            LogManager::LogError("Node already exists in current rank");
        }
    }
    std::vector<int> num_node_request_others(num_of_ranks_, 0);
    ptr_vec_num_recv->resize(num_of_ranks_);
    ptr_vec_num_recv->assign(num_of_ranks_, 0);

    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        if (requested_nodes.at(i_rank).size() > (std::numeric_limits<int>::max)()) {
            LogManager::LogError("Number of requested nodes is lager max int value ("
                + std::to_string((std::numeric_limits<int>::max)()) + ")");
        }
        ptr_vec_num_recv->at(i_rank) = static_cast<int>(requested_nodes.at(i_rank).size());
    }
    MPI_Alltoall(ptr_vec_num_recv->data(), 1, MPI_INT,
        num_node_request_others.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // send and receive which need nodes are needed for interpolation
    ptr_node_inner_layers->clear();
    if (SendNReceiveRequestNodesSFbitset(*ptr_vec_num_recv, num_node_request_others,
        requested_nodes, ptr_node_inner_layers) != 0) {
        LogManager::LogError("send and receive nodes requested for interpolation failed");
    }
}
/**
 * @brief function to send and receive ghost nodes used for interpolation.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] grid_info_lower class storting grid information at lower level.
 * @param[out] ptr_send_buffer_info pointer to buffer size information of sending.
 * @param[out] ptr_receive_buffer_info pointer to buffer size information of receiving.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive  pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_vec_ptr_buffer_send pointer to buffer storing grid node information for sending.
 * @param[out] ptr_vec_ptr_buffer_receive pointer to buffer storing received grid node information.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 * @return 0 run successfully, otherwise some error occurred in this function.
 */
int MpiManager::SendGhostNodeForInterpolation(const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefInt>& vec_num_recv,
    const GridInfoInterface& grid_info_lower, GridInfoInterface* ptr_grid_info,
    std::vector<BufferSizeInfo>* const ptr_send_buffer_info,
    std::vector<BufferSizeInfo>* const ptr_receive_buffer_info,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive) const {
    ptr_send_buffer_info->resize(num_of_ranks_);
    ptr_receive_buffer_info->resize(num_of_ranks_);
    const DefInt i_level = ptr_grid_info->GetGridLevel();
    if (i_level > 0) {
        // send and receive which nodes are needed for interpolation on other ranks
        int key_size = sizeof(DefSFBitset);
        int node_size = ptr_grid_info->GetSizeOfGridNodeInfoForMpiCommunication();
        int node_buffer_size = key_size + node_size;
        std::vector<MPI_Request> reqs_receive(num_of_ranks_);
        for (int i = 1; i < num_of_ranks_; ++i) {
            // send and receive number of chunks if needed
            int i_rank_send = (rank_id_ + i) % num_of_ranks_;
            int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
            bool bool_send_to_i_rank = (ptr_grid_info->interp_nodes_inner_layer_.find(i_rank_send)
                != ptr_grid_info->interp_nodes_inner_layer_.end());
            int num_node_all = 0, last_num_nodes = 0, other_num_nodes = 0;
            if (bool_send_to_i_rank) {
                num_node_all = static_cast<int>(ptr_grid_info->interp_nodes_inner_layer_.at(i_rank_send).size());
                int num_chunks = static_cast<int>(node_buffer_size * num_node_all
                    /((std::numeric_limits<int>::max)() - sizeof(int)) + 1);
                ptr_send_buffer_info->at(i_rank_send).bool_exist_ = true;
                ptr_send_buffer_info->at(i_rank_send).num_chunks_ = num_chunks;
                if (num_chunks < 2) {  // buffer is enough to store information of all nodes
                    other_num_nodes = 0;
                    last_num_nodes = num_node_all;
                } else {
                    last_num_nodes = num_node_all%(num_chunks - 1);
                    other_num_nodes = (num_node_all - last_num_nodes)/(num_chunks - 1);
                }
                ptr_send_buffer_info->at(i_rank_send).num_chunks_ = num_chunks;
                MPI_Send(&num_chunks, 1, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
            }

            bool bool_receive_from_i_rank = (ptr_grid_info->interp_nodes_inner_layer_.find(i_rank_receive)
                != ptr_grid_info->interp_nodes_inner_layer_.end());
            if (vec_num_recv.at(i_rank_receive) > 0) {
                ptr_receive_buffer_info->at(i_rank_receive).bool_exist_ = true;
                MPI_Recv(&ptr_receive_buffer_info->at(i_rank_receive).num_chunks_,
                    1, MPI_INT, i_rank_receive, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // send and receive buffer size of each chunk if needed
            if (bool_send_to_i_rank) {
                DefSizet buffer_size_send = node_buffer_size * other_num_nodes;
                if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                    ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.at(0) = static_cast<int>(buffer_size_send);
                } else {
                    LogManager::LogError("size of the buffer is greater than the maximum value of int");
                    return -1;
                }
                buffer_size_send = node_buffer_size * last_num_nodes;
                if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                    ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.at(1) = static_cast<int>(buffer_size_send);
                } else {
                    LogManager::LogError("size of the buffer is greater than the maximum value of int");
                }
                MPI_Send(ptr_send_buffer_info->at(i_rank_send).array_buffer_size_.data(),
                    2, MPI_INT, i_rank_send, i_level, MPI_COMM_WORLD);
            }
            if (vec_num_recv.at(i_rank_receive) > 0) {
                MPI_Recv(ptr_receive_buffer_info->at(i_rank_receive).array_buffer_size_.data(),
                    2, MPI_INT, i_rank_receive, i_level, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // send node information
        std::function<void(const GridNode&, char* const)> func_copy_a_node_to_buffer =
            [](const GridNode& node_ref, char* const ptr_node_buffer) {
            node_ref.CopyAInterpNodeToBufferForMpi(ptr_node_buffer);
        };
        std::function<void(const DefMap<DefInt>& , char* const)> func_copy_node_to_buffer =
            [func_copy_a_node_to_buffer, &ptr_grid_info, &grid_info_lower](
            const DefMap<DefInt>& map_send, char* const ptr_buffer) {
            ptr_grid_info->CopyInterpolationNodeInfoToBuffer(
                func_copy_a_node_to_buffer, grid_info_lower, map_send, ptr_buffer);
        };
        for (int i = 1; i < num_of_ranks_; ++i) {
            int i_rank_send = (rank_id_ + i) % num_of_ranks_;
            int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;

            if (ptr_send_buffer_info->at(i_rank_send).bool_exist_) {
                ptr_vec_vec_reqs_send->push_back({});
                if (ptr_grid_info->interp_nodes_inner_layer_.find(i_rank_send)
                    == ptr_grid_info->interp_nodes_inner_layer_.end()) {
                    LogManager::LogError("nodes need to be sent does not exist");
                }
                ptr_vec_ptr_buffer_send->emplace_back(std::make_unique<char[]>(CalculateBufferSizeForGridNode(
                    ptr_grid_info->interp_nodes_inner_layer_.at(i_rank_send).size(), *ptr_grid_info)));
                SendGridNodeInformation(i_rank_send, ptr_send_buffer_info->at(i_rank_send),
                    ptr_grid_info->interp_nodes_inner_layer_, func_copy_node_to_buffer,
                    ptr_grid_info, ptr_vec_ptr_buffer_send->back().get(),
                    &ptr_vec_vec_reqs_send->back());
            }
            if (ptr_receive_buffer_info->at(i_rank_receive).bool_exist_) {
                ptr_vec_vec_reqs_receive->push_back({});
                ptr_vec_ptr_buffer_receive->emplace_back(ReceiveGridNodeInformation(i_rank_receive, node_size,
                    ptr_receive_buffer_info->at(i_rank_receive), &ptr_vec_vec_reqs_receive->back()));
            }
        }
    }
    return 0;
}
/**
 * @brief function to wait for mpi communication and read ghost nodes for interpolation.
 * @param[in] send_buffer_info  buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[out] vec_ptr_buffer_receive buffer storing received grid node information.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive  pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 * @return 0 run successfully, otherwise some error occurred in this function.
 */
int MpiManager::WaitAndReadGhostNodeForInterpolation(const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    GridInfoInterface* const ptr_grid_info) const {
    int i_rev = 0;
    std::function<void(const char*,  GridNode* const ptr_node)> func_read_a_node_from_buffer =
        [](const char* const ptr_node_buffer, GridNode* const ptr_node) {
        ptr_node->ReadAInterpNodeFromBufferForMpi(ptr_node_buffer);
    };
    for (int i_rank = 1; i_rank < num_of_ranks_; ++i_rank) {
        int i_rank_receive = (rank_id_ - i_rank + num_of_ranks_)% num_of_ranks_;;
        if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
            MPI_Waitall(static_cast<int>(ptr_vec_vec_reqs_receive->at(i_rev).size()),
                ptr_vec_vec_reqs_receive->at(i_rev).data(), MPI_STATUSES_IGNORE);
            DefSizet buffer_size = receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(1)
                + (receive_buffer_info.at(i_rank_receive).num_chunks_ - 1)
                *receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(0);
            ptr_grid_info->ReadInterpolationNodeInfoFromBuffer(
                func_read_a_node_from_buffer, buffer_size, vec_ptr_buffer_receive.at(i_rev));
            ++i_rev;
        }
    }
    int i_send = 0;
    for (int i_rank = 1; i_rank < num_of_ranks_; ++i_rank) {
        int i_rank_send = (rank_id_ + i_rank) % num_of_ranks_;
        if (send_buffer_info.at(i_rank_send).bool_exist_) {
            MPI_Waitall(static_cast<int>(ptr_vec_vec_reqs_send->at(i_send).size()),
                ptr_vec_vec_reqs_send->at(i_send).data(), MPI_STATUSES_IGNORE);
            ++i_send;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
