//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_info_mpi.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid information at each refinement level when mpi is enabled.
* @date  2023-5-2
*/
#include <stddef.h>
#include <utility>
#include <tuple>
#include <set>
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#include "grid/grid_info_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
struct CriterionIndexForMpi {
    ECriterionType enum_criterion;
    DefSizet index_criterion;
    DefUint index_creator;
};
/**
 * @brief function to create and commit MPI std::tuple<ECriterionType, unit64_t, DefUint> types.
 * @param[out] ptr_mpi_tracking_index_type pointer to the std::tuple<ECriterionType, unit64_t, DefUint> types.
 * @return if the MPI types were created and committed successfully.
 */
int GridInfoInterface::CreateAndCommitCriterionIndexType(MPI_Datatype* ptr_mpi_tracking_index_type) {
    MPI_Datatype mpi_tracking_enum_type;
    MPI_Type_contiguous(sizeof(ECriterionType), MPI_BYTE, &mpi_tracking_enum_type);
    MPI_Type_commit(&mpi_tracking_enum_type);
    std::tuple<ECriterionType, uint64_t, DefUint> tuple_tmp;
    // the last MPI_UNSIGNED store the index (k0IndexCreator) of the creator for tracking grid nodes
    MPI_Datatype mpi_create_type[3] = {mpi_tracking_enum_type, MPI_SIZET_DATA_TYPE, MPI_UINT_DATA_TYPE};
    int mpi_block_length[3] = {1, 1, 1};
    MPI_Aint mpi_disp[3];
    mpi_disp[0] = offsetof(CriterionIndexForMpi, enum_criterion);
    mpi_disp[1] = offsetof(CriterionIndexForMpi, index_criterion);
    mpi_disp[2] = offsetof(CriterionIndexForMpi, index_creator);
    MPI_Type_create_struct(3, mpi_block_length, mpi_disp, mpi_create_type, ptr_mpi_tracking_index_type);
    MPI_Type_commit(ptr_mpi_tracking_index_type);
    MPI_Type_free(&mpi_tracking_enum_type);
    return MPI_SUCCESS;
}
/**
 * @brief function to serializes a set of tracking nodes into a buffer.
 * @param[in] set_nodes a set of space filling code for tracking nodes to be serialized.
 * @param[out] buffer a pointer to a char array where the serialized data will be stored.
 * @return The size of the serialized data in bytes.
 */
int TrackingGridInfoInterface::SerializeTrackingNode(const std::set<DefSFCodeToUint>& set_nodes,
        std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes;
    if  (sizeof(int) + set_nodes.size() *key_size > 0x7FFFFFFF) {
        LogError("size of the buffer is greater than the maximum of int in GridInfoInterface::SerializeTrackingNode");
    } else {
        num_nodes = static_cast<int>(set_nodes.size());
    }
    int buffer_size = sizeof(int) + key_size * num_nodes;
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_nodes, sizeof(int));
    position += sizeof(int);
    // serialize data stored in nodes
    for (const auto& iter : set_nodes) {
        // convert bitset into unsigned long long and save it to buffer
        std::memcpy(ptr_buffer + position, &iter, key_size);
        position += key_size;
    }
    return buffer_size;
}
/**
 * @brief function to deserialize tracking node data and insert it into the associated map.
 * @param[in] criterion_type type of criterion.
 * @param[in]  index_criterion index of the criterion.
 * @param[in]  buffer buffer containing the tracking node data.
 */
void TrackingGridInfoInterface::DeserializeTrackingNode(const std::unique_ptr<char[]>& buffer) {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    // number of nodes
    int num_nodes;
    int position = 0;
    std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
    // deserialize data stored in buffer
    position += sizeof(int);
    DefSFCodeToUint key_net;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_net, ptr_buffer + position, key_size);
        position += key_size;
        map_tracking_node_.insert({static_cast<DefSFBitset>(key_net), k0TrackNodeInstance_});
    }
}
/**
 * @brief function to send/receive tracking grid information to/from other MPI processes.
 * @param[in] dims dimensions of the grid.
 * @param[in] bitset_max vector containing the maximum space filling code for each grid.
 * @param[in] mpi_manager manager object containing the mpi related information.
 * @param[in] bitset_aux manager object for manipulating space filling code.
 * @param[in] vec_tracking_info_creator vector of unique pointers to objects for creating instance of tracking grid.
 */
void GridInfoInterface::IniSendNReceiveTracking(const DefUint dims, const std::vector<DefSFBitset>& bitset_max,
    const MpiManager& mpi_manager, const SFBitsetAuxInterface& bitset_aux,
    const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator) {
    MPI_Datatype mpi_tracking_index_type;
    CreateAndCommitCriterionIndexType(&mpi_tracking_index_type);
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }

    DefSizet num_max = ull_max.size();
    DefSizet num_tracking = map_ptr_tracking_grid_info_.size();
    MPI_Bcast(&num_tracking, 1, MPI_SIZET_DATA_TYPE, 0, MPI_COMM_WORLD);
    CriterionIndexForMpi tracking_index_tmp;
    std::vector<std::pair<ECriterionType, DefSizet>> vec_tracking_indices;
    if (rank_id == 0) {
        // save tracking indices stored in unordered_map to vector
        for (const auto& iter_tracking : map_ptr_tracking_grid_info_) {
            vec_tracking_indices.push_back(iter_tracking.first);
        }
    }

    for (DefSizet i_tracking = 0; i_tracking < num_tracking; ++i_tracking) {
        // broadcast the current index of the tracking grid
        if (rank_id == 0) {
            tracking_index_tmp.enum_criterion = vec_tracking_indices.at(i_tracking).first;
            tracking_index_tmp.index_criterion = vec_tracking_indices.at(i_tracking).second;
            tracking_index_tmp.index_creator =
                map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0IndexCreator;
#ifdef DEBUG_CHECK_GRID
        if (map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterNeg.size() != dims) {
            LogError("size of k0ExtendOuterNeg " + std::to_string(
            map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterNeg.size())
            + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking");
        }
        if (map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterPos.size() != dims) {
            LogError("size of k0ExtendOuterPos " + std::to_string(
            map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterPos.size())
            + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking");
        }
        if (map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerNeg.size() != dims) {
            LogError("size of k0ExtendInnerNeg " + std::to_string(
            map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerNeg.size())
            + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking");
        }
        if (map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerPos.size() != dims) {
            LogError("size of k0ExtendInnerPos " + std::to_string(
            map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerPos.size())
            + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking");
        }
#endif  // DEBUG_CHECK_GRID
        }
        MPI_Bcast(&tracking_index_tmp, 1, mpi_tracking_index_type, 0, MPI_COMM_WORLD);
        std::pair<ECriterionType, DefSizet> pair_tracking = {
            tracking_index_tmp.enum_criterion, tracking_index_tmp.index_criterion};
        if (rank_id > 0) {
            DefUint index_creator = tracking_index_tmp.index_creator;
            if (map_ptr_tracking_grid_info_.find(pair_tracking) == map_ptr_tracking_grid_info_.end()) {
                map_ptr_tracking_grid_info_.insert({pair_tracking,
                    vec_tracking_info_creator.at(index_creator).get()->CreateTrackingGridInfo()});
                map_ptr_tracking_grid_info_.at(pair_tracking)->k0IndexCreator = index_creator;
            }
            map_ptr_tracking_grid_info_.at(pair_tracking)->k0ExtendInnerNeg.resize(dims);
            map_ptr_tracking_grid_info_.at(pair_tracking)->k0ExtendInnerPos.resize(dims);
            map_ptr_tracking_grid_info_.at(pair_tracking)->k0ExtendOuterNeg.resize(dims);
            map_ptr_tracking_grid_info_.at(pair_tracking)->k0ExtendOuterPos.resize(dims);
        }
        TrackingGridInfoInterface* ptr_tracking_info = map_ptr_tracking_grid_info_.at(pair_tracking).get();

         // broadcast number of extending layers
        MPI_Bcast(ptr_tracking_info->k0ExtendInnerNeg.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendInnerPos.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendOuterNeg.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendOuterPos.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);

        // broadcast type information
        DefTypeUint uint_type[2] = {0, 0};
        if (rank_id == 0) {
            uint_type[0] = ptr_tracking_info->computational_cost_;
            uint_type[1] = static_cast<DefTypeUint>(ptr_tracking_info->grid_extend_type_);
        }
        MPI_Bcast(uint_type, 2, MPI_TYPEUINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        if (rank_id > 0) {
            ptr_tracking_info->computational_cost_ = uint_type[0];
            ptr_tracking_info->grid_extend_type_ = static_cast<EGridExtendType>(uint_type[1]);
        }

        //  send and receive tracking grid information
        if (rank_id == 0 && num_ranks > 1) {
            DefSFCodeToUint background_code, key_ull;
            int max_buffer = (std::numeric_limits<int>::max)() / ptr_tracking_info->SizeOfEachNodeForMpi() - 1;
            std::vector<std::vector<std::set<DefSFCodeToUint>>> vec_nodes_ranks(num_ranks);
            vec_nodes_ranks.at(0).push_back({});
            int index;
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            for (const auto& iter_node : ptr_tracking_info->map_tracking_node_) {
                background_code = (bitset_aux.SFBitsetToNLowerLevelVir(i_level_, iter_node.first)).to_ullong();
                auto iter_index = std::lower_bound(ull_max.begin(),
                    ull_max.end(), background_code);
                index = static_cast<int>(iter_index - ull_max.begin());
                key_ull = iter_node.first.to_ullong();
                if (index != 0) {
                    vec_nodes_ranks.at(0).at(0).insert(key_ull);  // nodes need to be deleted on rank 0
#ifdef DEBUG_CHECK_GRID
                    if (index == num_max) {
                        // lower_bound returns the next element of last in ull_max if not found the desired one,
                        // which means that the space filling code of the node exceeds the maximum given by ull_max
                        LogError("nodes is out of computational domain in GridInfoInterface::IniSendNReceiveTracking");
                    }
#endif  // DEBUG_CHECK_GRID
                    if (i_counts.at(index) == 0) {
                        vec_nodes_ranks.at(index).push_back({key_ull});
                        i_chunk_each_rank.at(index) +=1;
                        ++i_counts.at(index);
                    } else {
                        vec_nodes_ranks.at(index).at(i_chunk_each_rank.at(index)).insert(key_ull);
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
                MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, 0, MPI_COMM_WORLD);
                std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
                std::vector<MPI_Request> reqs_send(num_chunks);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    buffer_size_send = ptr_tracking_info->SerializeTrackingNode(
                        vec_nodes_ranks.at(iter_rank).at(i_chunk),
                    vec_ptr_buffer.at(i_chunk));
                    MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                    i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
                }
                MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
            }
            for (const auto& iter_key : vec_nodes_ranks.at(0).at(0)) {
                ptr_tracking_info->map_tracking_node_.erase(static_cast<DefSFBitset>(iter_key));
            }
        } else if (rank_id > 0) {
            int num_chunks = 0;
            MPI_Recv(&num_chunks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                int buffer_size_receive;
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
                MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                    i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                ptr_tracking_info->DeserializeTrackingNode(vec_ptr_buffer);
            }
        }
    }
    MPI_Type_free(&mpi_tracking_index_type);
}
/**
 * @brief function to serializes a set of tracking nodes into a buffer.
 * @param[in] set_nodes a set of space filling code for tracking nodes to be serialized.
 * @param[out] buffer a pointer to a char array where the serialized data will be stored.
 * @return The size of the serialized data in bytes.
 */
int InterfaceLayerInfo::SerializeInterfaceNode(const std::set<DefSFCodeToUint>& set_nodes,
        std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes;
    if  (sizeof(int) + set_nodes.size() *key_size > 0x7FFFFFFF) {
        LogError("size of the buffer is greater than the maximum of int in GridInfoInterface::SerializeTrackingNode");
    } else {
        num_nodes = static_cast<int>(set_nodes.size());
    }
    int buffer_size = sizeof(int) + key_size * num_nodes;
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_nodes, sizeof(int));
    position += sizeof(int);
    // serialize data stored in nodes
    for (const auto& iter : set_nodes) {
        // convert bitset into unsigned long long and save it to buffer
        std::memcpy(ptr_buffer + position, &iter, key_size);
        position += key_size;
    }
    return buffer_size;
}
/**
 * @brief function to deserialize tracking node data and insert it into the associated map.
 * @param[in] criterion_type type of criterion.
 * @param[in]  index_criterion index of the criterion.
 * @param[in]  buffer buffer containing the tracking node data.
 */
void InterfaceLayerInfo::DeserializeInterfaceNode(
     const std::unique_ptr<char[]>& buffer, DefMap<DefUint>* ptr_map_interface_layer) {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    // number of nodes
    int num_nodes;
    int position = 0;
    std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
    // deserialize data stored in buffer
    position += sizeof(int);
    DefSFCodeToUint key_net;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_net, ptr_buffer + position, key_size);
        position += key_size;
        ptr_map_interface_layer->insert({key_net, 0});
    }
}
/**
 * @brief  function to send and receive nodes in a given level of the interface layer.
 * @param[in] i_level  level of the interface layer to send/receive.
 * @param[in] ull_max vector of the maximum space filling codes of partitioned grid.
 * @param[in] mpi_manager manager object containing the mpi related information.
 * @param[in] bitset_aux manager object for manipulating space filling code.
 * @param[out] ptr_map_interface_layer pointer to nodes in the interface layer.
 */
void InterfaceLayerInfo::SendNReiveOneLayer(const DefSizet i_level, const std::vector<DefSFCodeToUint>& ull_max,
    const MpiManager& mpi_manager, const SFBitsetAuxInterface& bitset_aux,
    DefMap<DefUint>* const ptr_map_interface_layer) {
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    DefSFCodeToUint background_code, key_ull;
    std::vector<std::vector<std::set<DefSFCodeToUint>>> vec_nodes_ranks(num_ranks);
    vec_nodes_ranks.at(0).push_back({});
    DefSizet num_max = ull_max.size();
    int index;
    std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
    int max_buffer = (std::numeric_limits<int>::max)() / SizeOfEachNodeForMpi() - 1;
    if (rank_id == 0 && num_ranks > 1) {
        for (const auto& iter_node : *ptr_map_interface_layer) {
            background_code = (bitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first)).to_ullong();
            auto iter_index = std::lower_bound(ull_max.begin(),
                ull_max.end(), background_code);
            index = static_cast<int>(iter_index - ull_max.begin());
            key_ull = iter_node.first.to_ullong();
            if (index != 0) {
                vec_nodes_ranks.at(0).at(0).insert(key_ull);  // nodes need to be deleted on rank 0
#ifdef DEBUG_CHECK_GRID
                if (index == num_max) {
                    // lower_bound returns the next element of last in ull_max if not found the desired one,
                    // which means that the space filling code of the node exceeds the maximum
                    LogError("node is out of computational domain in InterfaceLayerInfo::SendNReiveOneLayer");
                }
#endif  // DEBUG_CHECK_GRID
                if (i_counts.at(index) == 0) {
                    vec_nodes_ranks.at(index).push_back({key_ull});
                    i_chunk_each_rank.at(index) +=1;
                    ++i_counts.at(index);
                } else {
                    vec_nodes_ranks.at(index).at(i_chunk_each_rank.at(index)).insert(key_ull);
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
            MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, 0, MPI_COMM_WORLD);
            std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
            std::vector<MPI_Request> reqs_send(num_chunks);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                buffer_size_send = SerializeInterfaceNode(
                    vec_nodes_ranks.at(iter_rank).at(i_chunk), vec_ptr_buffer.at(i_chunk));
                MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
            }
            MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
        }
        for (const auto& iter_key : vec_nodes_ranks.at(0).at(0)) {
            ptr_map_interface_layer->erase(static_cast<DefSFBitset>(iter_key));
        }
    } else if (rank_id > 0) {
        int num_chunks = 0;
        MPI_Recv(&num_chunks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
            int buffer_size_receive;
            MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
            MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            DeserializeInterfaceNode(vec_ptr_buffer, ptr_map_interface_layer);
        }
    }
}
/**
 * @brief function to send/receive information of interface layers.
 * @param[in] dims dimension of the grid.
 * @param[in] bitset_max vector containing the maximum space filling code for each grid.
 * @param[in] mpi_manager manager object containing the mpi related information.
 * @param[in] bitset_aux manager object for manipulating space filling code.
 */
void GridInfoInterface::IniSendNReceiveInterface(const DefUint dims,
    const std::vector<DefSFBitset>& bitset_max, const MpiManager& mpi_manager,
    const SFBitsetAuxInterface& bitset_aux) {
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitCriterionIndexType(&mpi_interface_index_type);
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }
    DefSizet num_interface = map_ptr_interface_layer_info_.size();
    MPI_Bcast(&num_interface, 1, MPI_SIZET_DATA_TYPE, 0, MPI_COMM_WORLD);
    CriterionIndexForMpi interface_index_tmp;
    std::vector<std::pair<ECriterionType, DefSizet>> vec_interface_indices;
    if (rank_id == 0) {
        // save interface indices stored in unordered_map to vector
        for (const auto& iter_interface : map_ptr_interface_layer_info_) {
            vec_interface_indices.push_back(iter_interface.first);
        }
    }

    for (uint64_t i_interface = 0; i_interface < num_interface; ++i_interface) {
        // broadcast the current index of the interfaces
        if (rank_id == 0) {
            interface_index_tmp.enum_criterion = vec_interface_indices.at(i_interface).first;
            interface_index_tmp.index_criterion = static_cast<uint64_t>(vec_interface_indices.at(i_interface).second);
            interface_index_tmp.index_creator = 0;
        }
        MPI_Bcast(&interface_index_tmp, 1, mpi_interface_index_type, 0, MPI_COMM_WORLD);
        std::pair<ECriterionType, DefSizet> pair_interface = {
                interface_index_tmp.enum_criterion, interface_index_tmp.index_criterion};
        if (rank_id > 0) {  // create interface instance
            if (map_ptr_interface_layer_info_.find(pair_interface) == map_ptr_interface_layer_info_.end()) {
                map_ptr_interface_layer_info_.insert({pair_interface, std::make_shared<InterfaceLayerInfo>()});
            }
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendInnerNeg.resize(dims);
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendInnerPos.resize(dims);
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendOuterNeg.resize(dims);
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendOuterPos.resize(dims);
            map_ptr_interface_layer_info_.at(pair_interface)->vec_outer_coarse2fine_.resize(k0NumCoarse2FineLayer_);
            map_ptr_interface_layer_info_.at(pair_interface)->vec_inner_coarse2fine_.resize(k0NumCoarse2FineLayer_);
        } else {
#ifdef DEBUG_CHECK_GRID
        if (map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendOuterNeg.size() != dims) {
            LogError("size of k0ExtendOuterNeg " + std::to_string(
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendOuterNeg.size())
            + "is not equal to the dimension " + std::to_string(dims) + "in IniSendNReceiveInterface");
        }
        if (map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendOuterPos.size() != dims) {
            LogError("size of k0ExtendOuterPos " + std::to_string(
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendOuterPos.size())
            + "is not equal to the dimension " + std::to_string(dims) + "in IniSendNReceiveInterface");
        }
        if (map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendInnerNeg.size() != dims) {
            LogError("size of k0ExtendInnerNeg " + std::to_string(
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendInnerNeg.size())
            + "is not equal to the dimension " + std::to_string(dims) + "in IniSendNReceiveInterface");
        }
        if (map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendInnerPos.size() != dims) {
            LogError("size of k0ExtendInnerPos " + std::to_string(
            map_ptr_interface_layer_info_.at(pair_interface)->k0ExtendInnerPos.size())
            + "is not equal to the dimension " + std::to_string(dims) + "in IniSendNReceiveInterface");
        }
        if (map_ptr_interface_layer_info_.at(pair_interface)->vec_outer_coarse2fine_.size()
           != k0NumCoarse2FineLayer_) {
            LogError("size of vec_outer_coarse2fine_ " + std::to_string(
            map_ptr_interface_layer_info_.at(pair_interface)->vec_outer_coarse2fine_.size())
            + "is not equal to the number of layers " + std::to_string(k0NumCoarse2FineLayer_)
            + "in IniSendNReceiveInterface");
        }
        if (map_ptr_interface_layer_info_.at(pair_interface)->vec_inner_coarse2fine_.size()
           != k0NumCoarse2FineLayer_) {
            LogError("size of vec_outer_coarse2fine_ " + std::to_string(
            map_ptr_interface_layer_info_.at(pair_interface)->vec_inner_coarse2fine_.size())
            + "is not equal to the number of layers " + std::to_string(k0NumCoarse2FineLayer_)
            + "in IniSendNReceiveInterface");
        }
#endif  // DEBUG_CHECK_GRID
        }
        InterfaceLayerInfo* ptr_interface = map_ptr_interface_layer_info_.at(pair_interface).get();

        // broadcast number of extending layers
        MPI_Bcast(ptr_interface->k0ExtendInnerNeg.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_interface->k0ExtendInnerPos.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_interface->k0ExtendOuterNeg.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_interface->k0ExtendOuterPos.data(),
           static_cast<int>(dims), MPI_UINT_DATA_TYPE, 0, MPI_COMM_WORLD);

        //  send and receive outer layers
        for (DefSizet i_layer = 0; i_layer < k0NumCoarse2FineLayer_; ++i_layer) {
            ptr_interface->SendNReiveOneLayer(i_level_, ull_max,
                mpi_manager, bitset_aux, &ptr_interface->vec_outer_coarse2fine_.at(i_layer));
        }
        //  send and receive inner layers
        for (DefSizet i_layer = 0; i_layer < k0NumCoarse2FineLayer_; ++i_layer) {
            ptr_interface->SendNReiveOneLayer(i_level_, ull_max,
                mpi_manager, bitset_aux, &ptr_interface->vec_inner_coarse2fine_.at(i_layer));
        }
    }

    MPI_Type_free(&mpi_interface_index_type);
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
