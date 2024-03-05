//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_grid_info.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid information at each refinement level when mpi is enabled.
* @date  2023-5-2
*/
#include <stddef.h>
#include <utility>
#include <tuple>
#include <set>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
#include "io/debug_write.h"
namespace rootproject {
namespace amrproject {
struct CriterionIndexForMpi {
    ECriterionType enum_criterion;
    DefAmrIndexUint index_criterion;
    DefAmrIndexUint index_creator;
};
/**
 * @brief function to create and commit MPI <ECriterionType, DefAmrIndexUint, DefAmrIndexUint> types.
 * @param[out] ptr_mpi_tracking_index_type pointer to the <ECriterionType, DefAmrIndexUint, DefAmrIndexUint> types.
 * @return if the MPI types were created and committed successfully.
 */
int MpiManager::CreateAndCommitCriterionIndexType(MPI_Datatype* ptr_mpi_tracking_index_type) const {
    MPI_Datatype mpi_tracking_enum_type;
    MPI_Type_contiguous(sizeof(ECriterionType), MPI_BYTE, &mpi_tracking_enum_type);
    MPI_Type_commit(&mpi_tracking_enum_type);
    // the last MPI_UNSIGNED store the index (k0IndexCreator) of the creator for tracking grid nodes
    MPI_Datatype mpi_create_type[3] = {mpi_tracking_enum_type, MPI_AMR_INDEX_UINT_TYPE, MPI_AMR_INDEX_UINT_TYPE};
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
 * @param[out] ptr_buffer_size pointer to size of the buffer in bytes.
 * @return unique pointer to a char array to store the serialized data.
 */
std::unique_ptr<char[]> MpiManager::IniSerializeTrackingNode(
    const std::set<DefSFCodeToUint>& set_nodes, int* const ptr_buffer_size) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes = 0;
    if  (sizeof(int) + set_nodes.size() *key_size > 0x7FFFFFFF) {
        LogManager::LogError("size of the buffer is greater than"
         " the maximum of int in MpiManager::IniSerializeTrackingNode in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        num_nodes = static_cast<int>(set_nodes.size());
    }
    int& buffer_size = *ptr_buffer_size;
    buffer_size = sizeof(int) + key_size * num_nodes;
    // allocation buffer to store the serialized data
    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size);
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
    return buffer;
}
/**
 * @brief function to deserialize tracking node data and insert it into the associated map.
 * @param[in]  index_criterion index of the criterion.
 * @param[in]  buffer buffer containing the tracking node data.
 * @param[in]  tracking_node_instance instance of the tracking node for initialization.
 * @param[out] ptr_map_tracking pointer to a container store tracking nodes.
 */
void MpiManager::IniDeserializeTrackingNode(const std::unique_ptr<char[]>& buffer,
    const TrackingNode& tracking_node_instance, DefMap<TrackingNode>* const ptr_map_tracking) const {
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
        ptr_map_tracking->insert({static_cast<DefSFBitset>(key_net), tracking_node_instance});
    }
}
/**
 * @brief function to send/receive tracking grid information to/from other MPI processes.
 * @param[in] dims dimensions of the grid.
 * @param[in] i_level refinement level of the tracking grid.
 * @param[in] bitset_max vector containing the maximum space filling code for each grid.
 * @param[in] bitset_aux manager object for manipulating space filling code.
 * @param[in] tracking_node_instance instance of the tracking node for initialization.
 * @param[in] vec_tracking_info_creator vector of unique pointers to objects for creating instance of tracking grid.
 * @param[out] ptr_map_tracking_info map of unique pointers to tracking grid information.
 */
void MpiManager::IniSendNReceiveTracking(const DefAmrIndexUint dims, const DefAmrIndexUint i_level,
    const std::vector<DefSFBitset>& bitset_max, const SFBitsetAuxInterface& bitset_aux,
    const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator,
    std::map<std::pair<ECriterionType, DefAmrIndexUint>, std::shared_ptr<TrackingGridInfoInterface>>*
    const ptr_map_tracking_info) const {
    MPI_Datatype mpi_tracking_index_type;
    CreateAndCommitCriterionIndexType(&mpi_tracking_index_type);
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (DefSizet i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }

    DefSizet num_max = ull_max.size();
    DefAmrIndexUint num_tracking = DefAmrIndexUint(ptr_map_tracking_info->size());
    MPI_Bcast(&num_tracking, 1, MPI_AMR_INDEX_UINT_TYPE, 0, MPI_COMM_WORLD);
    CriterionIndexForMpi tracking_index_tmp;
    std::vector<std::pair<ECriterionType, DefAmrIndexUint>> vec_tracking_indices;
    if (rank_id == 0) {
        // save tracking indices stored in unordered_map to vector
        for (const auto& iter_tracking : *ptr_map_tracking_info) {
            vec_tracking_indices.push_back(iter_tracking.first);
        }
    }

    for (DefAmrIndexUint i_tracking = 0; i_tracking < num_tracking; ++i_tracking) {
        // broadcast the current index of the tracking grid
        if (rank_id == 0) {
            tracking_index_tmp.enum_criterion = vec_tracking_indices.at(i_tracking).first;
            tracking_index_tmp.index_criterion = vec_tracking_indices.at(i_tracking).second;
            tracking_index_tmp.index_creator =
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0IndexCreator;
#ifdef DEBUG_CHECK_GRID
            if (ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterNeg_.size() != dims) {
                LogManager::LogError("size of k0ExtendOuterNeg_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterNeg_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            if (ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterPos_.size() != dims) {
                LogManager::LogError("size of k0ExtendOuterPos_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterPos_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            if (ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerNeg_.size() != dims) {
                LogManager::LogError("size of k0ExtendInnerNeg_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerNeg_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            if (ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerPos_.size() != dims) {
                LogManager::LogError("size of k0ExtendInnerPos_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerPos_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
#endif  // DEBUG_CHECK_GRID
        }
        MPI_Bcast(&tracking_index_tmp, 1, mpi_tracking_index_type, 0, MPI_COMM_WORLD);
        std::pair<ECriterionType, DefAmrIndexUint> pair_tracking = {
            tracking_index_tmp.enum_criterion, tracking_index_tmp.index_criterion};
        if (rank_id > 0) {
            DefAmrIndexUint index_creator = tracking_index_tmp.index_creator;
            if (ptr_map_tracking_info->find(pair_tracking) == ptr_map_tracking_info->end()) {
                ptr_map_tracking_info->insert({pair_tracking,
                    vec_tracking_info_creator.at(index_creator).get()->CreateTrackingGridInfo()});
                ptr_map_tracking_info->at(pair_tracking)->k0IndexCreator = index_creator;
            }
            ptr_map_tracking_info->at(pair_tracking)->k0ExtendInnerNeg_.resize(dims);
            ptr_map_tracking_info->at(pair_tracking)->k0ExtendInnerPos_.resize(dims);
            ptr_map_tracking_info->at(pair_tracking)->k0ExtendOuterNeg_.resize(dims);
            ptr_map_tracking_info->at(pair_tracking)->k0ExtendOuterPos_.resize(dims);
        }
        TrackingGridInfoInterface* ptr_tracking_info = ptr_map_tracking_info->at(pair_tracking).get();

         // broadcast number of extending layers
        MPI_Bcast(ptr_tracking_info->k0ExtendInnerNeg_.data(),
           static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendInnerPos_.data(),
           static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendOuterNeg_.data(),
           static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendOuterPos_.data(),
           static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);

        // broadcast type information
        DefAmrTypeUint uint_type[2] = {0, 0};
        if (rank_id == 0) {
            uint_type[0] = ptr_tracking_info->computational_cost_;
            uint_type[1] = static_cast<DefAmrTypeUint>(ptr_tracking_info->grid_extend_type_);
        }
        MPI_Bcast(uint_type, 2, MPI_AMR_TYPE_UINT_TYPE, 0, MPI_COMM_WORLD);
        if (rank_id > 0) {
            ptr_tracking_info->computational_cost_ = uint_type[0];
            ptr_tracking_info->grid_extend_type_ = static_cast<EGridExtendType>(uint_type[1]);
        }

        //  send and receive tracking grid information
        if (rank_id == 0 && num_ranks > 1) {
            DefSFCodeToUint background_code, key_ull;
            int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefSFBitset) - 1;
            std::vector<std::vector<std::set<DefSFCodeToUint>>> vec_nodes_ranks(num_ranks);
            vec_nodes_ranks.at(0).push_back({});
            int index;
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            for (const auto& iter_node : ptr_tracking_info->map_tracking_node_) {
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
                        // which means that the space filling code of the node exceeds the maximum given by ull_max
                        LogManager::LogError("nodes is out of computational"
                         " domain in MpiManager::IniSendNReceiveTracking in "
                         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
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
                    vec_ptr_buffer.at(i_chunk) = IniSerializeTrackingNode(
                        vec_nodes_ranks.at(iter_rank).at(i_chunk), &buffer_size_send);
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
            TrackingNode tracking_node_instance = ptr_tracking_info->k0TrackNodeInstance_;
            MPI_Recv(&num_chunks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                int buffer_size_receive;
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
                MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                    i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                IniDeserializeTrackingNode(vec_ptr_buffer,
                    tracking_node_instance, &ptr_tracking_info->map_tracking_node_);
            }
        }
    }
    MPI_Type_free(&mpi_tracking_index_type);
}
/**
 * @brief function to serializes a set of tracking nodes into a buffer.
 * @param[in] interface_nodes space filling code for nodes on refinement interface to be serialized.
 * @param[out] buffer a pointer to a char array where the serialized data will be stored.
 * @return The size of the serialized data in bytes.
 */
int MpiManager::IniSerializeRefinementInterfaceNode(const DefMap<DefAmrUint>& interface_nodes,
        std::unique_ptr<char[]>& buffer) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes = 0;
    if  (sizeof(int) + interface_nodes.size() *key_size > 0x7FFFFFFF) {
        LogManager::LogError("size of the buffer is greater than"
         " the maximum of int in GridInfoInterface::SerializeTrackingNode in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        num_nodes = static_cast<int>(interface_nodes.size());
    }
    int buffer_size = sizeof(int) + key_size * num_nodes;
    // allocation buffer to store the serialized data
    buffer = std::make_unique<char[]>(buffer_size);
    char* ptr_buffer = buffer.get();
    int position = 0;
    std::memcpy(ptr_buffer + position, &num_nodes, sizeof(int));
    position += sizeof(int);
    // serialize data stored in nodes
    for (const auto& iter : interface_nodes) {
        // convert bitset into unsigned long long and save it to buffer
        std::memcpy(ptr_buffer + position, &iter.first, key_size);
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
void MpiManager::IniDeserializeRefinementInterfaceNode(const DefAmrUint flag0,
     const std::unique_ptr<char[]>& buffer, DefMap<DefAmrUint>* ptr_map_interface_layer) const {
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
        ptr_map_interface_layer->insert({key_net, flag0});
    }
}
/**
 * @brief  function to send and receive nodes in a given level of the interface layer.
 * @param[in] i_level  level of the interface layer to send/receive.
 * @param[in] flag0 flag for initialization.
 * @param[in] outmost_for_all_ranks nodes on the outmost fine to coarse interface for all ranks stored on rank 0. 
 * @param[out] ptr_map_interface_layer pointer to nodes on the refinement interface layer.
 */
void MpiManager::IniSendNReiveOneLayerRefinementInterface(
    const DefAmrUint flag0, const DefMap<std::set<int>>& outmost_for_all_ranks,
    DefMap<DefAmrUint>* const ptr_map_interface_layer) const {
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    std::vector<std::vector<DefMap<DefAmrUint>>> vec_nodes_ranks(num_ranks),
        vec_partition_nodes_ranks(num_ranks);
    std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0),
     i_chunk_partition_each_rank(num_ranks, -1), i_partition_counts(num_ranks, 0);
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefSFBitset) - 1;
    vec_nodes_ranks.at(0).push_back({});
    if (rank_id == 0) {
        std::vector<DefMap<DefAmrIndexUint>> partition_interface_current(num_ranks);
        std::vector<DefSFBitset> nodes_in_region;
        for (const auto& iter_node : *ptr_map_interface_layer) {
            if (outmost_for_all_ranks.find(iter_node.first) != outmost_for_all_ranks.end()) {
                for (const auto& iter_rank : outmost_for_all_ranks.at(iter_node.first)) {
                    if (iter_rank == 0) {  // data on rank 0 doesn't need partition
                        vec_nodes_ranks.at(0).at(0).insert({iter_node.first, flag0});
                    } else {
                        if (i_counts.at(iter_rank) == 0) {
                            vec_nodes_ranks.at(iter_rank).push_back({});
                            i_chunk_each_rank.at(iter_rank) +=1;
                            vec_nodes_ranks.at(iter_rank)
                                .at(i_chunk_each_rank.at(iter_rank)).insert({iter_node.first, flag0});
                            ++i_counts.at(iter_rank);
                        } else {
                            vec_nodes_ranks.at(iter_rank)
                                .at(i_chunk_each_rank.at(iter_rank)).insert({iter_node.first, flag0});
                            ++i_counts.at(iter_rank);
                            if (i_counts.at(iter_rank) == max_buffer) {
                                // check if size of send buffer exceeds limits of int
                                i_counts.at(iter_rank) = 0;
                            }
                        }
                    }
                }
            }
        }
        // keep nodes should be on rank 0 and delete others
        ptr_map_interface_layer->swap(vec_nodes_ranks.at(0).at(0));
        // send interface nodes to other ranks
        int buffer_size_send = 0;
        for (auto iter_rank = 1; iter_rank < num_ranks; ++iter_rank) {
            int num_chunks = static_cast<int>(vec_nodes_ranks.at(iter_rank).size());
            MPI_Send(&num_chunks, 1, MPI_INT, iter_rank, 0, MPI_COMM_WORLD);
            std::vector<std::unique_ptr<char[]>> vec_ptr_buffer(num_chunks);
            std::vector<MPI_Request> reqs_send(num_chunks);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                buffer_size_send = IniSerializeRefinementInterfaceNode(
                    vec_nodes_ranks.at(iter_rank).at(i_chunk), vec_ptr_buffer.at(i_chunk));
                MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                i_chunk, MPI_COMM_WORLD, &reqs_send[i_chunk]);
            }
            MPI_Waitall(num_chunks, reqs_send.data(), MPI_STATUS_IGNORE);
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
            IniDeserializeRefinementInterfaceNode(flag0, vec_ptr_buffer, ptr_map_interface_layer);
        }
    }
}
/**
 * @brief function to send/receive information of interface layers between grids at different refinement levels.
 * @param[in] dims dimension of the grid.
 * @param[in] i_level refinement level of the interface layer.
 * @param[in] num_of_layers_coarse2fine number of layers in the coarse to fine overlapping region.
 * @param[in] flag0 flag for initialization.
 * @param[in] outmost_for_all_ranks nodes on the outmost fine to coarse interface for all ranks stored on rank 0. 
 * @param[out] ptr_map_interface_info  pointer to refinement interface information.
 */
void MpiManager::IniSendNReceiveCoarse2Fine0Interface(const DefAmrIndexUint dims,
    const DefAmrIndexUint i_level, const DefAmrIndexUint num_of_layers_coarse2fine, const DefAmrUint flag0,
    const DefMap<std::set<int>>& outmost_for_all_ranks,
    std::map<std::pair<ECriterionType, DefAmrIndexUint>,
    std::shared_ptr<InterfaceLayerInfo>>* const ptr_map_interface_info) const {
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitCriterionIndexType(&mpi_interface_index_type);
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    DefAmrIndexUint num_interface = DefAmrIndexUint(ptr_map_interface_info->size());
    MPI_Bcast(&num_interface, 1, MPI_AMR_INDEX_UINT_TYPE, 0, MPI_COMM_WORLD);
    std::vector<std::pair<ECriterionType, DefAmrIndexUint>> vec_interface_indices;
    if (rank_id == 0) {
        // save interface indices stored in unordered_map to vector
        for (const auto& iter_interface : *ptr_map_interface_info) {
            vec_interface_indices.push_back(iter_interface.first);
        }
    }
    CriterionIndexForMpi interface_index_tmp;
    for (DefAmrIndexUint i_interface = 0; i_interface < num_interface; ++i_interface) {
        // broadcast the current index of the interfaces
        if (rank_id == 0) {
            interface_index_tmp.enum_criterion = vec_interface_indices[i_interface].first;
            interface_index_tmp.index_criterion =
                static_cast<uint64_t>(vec_interface_indices[i_interface].second);
            interface_index_tmp.index_creator = 0;
        }
        MPI_Bcast(&interface_index_tmp, 1, mpi_interface_index_type, 0, MPI_COMM_WORLD);
        std::pair<ECriterionType, DefAmrIndexUint> pair_interface = {
                interface_index_tmp.enum_criterion, interface_index_tmp.index_criterion};
        if (rank_id > 0) {  // create interface instance
            if (ptr_map_interface_info->find(pair_interface) == ptr_map_interface_info->end()) {
                ptr_map_interface_info->insert({pair_interface, std::make_shared<InterfaceLayerInfo>()});
            }
            ptr_map_interface_info->at(pair_interface)->k0ExtendInnerNeg_.resize(dims);
            ptr_map_interface_info->at(pair_interface)->k0ExtendInnerPos_.resize(dims);
            ptr_map_interface_info->at(pair_interface)->k0ExtendOuterNeg_.resize(dims);
            ptr_map_interface_info->at(pair_interface)->k0ExtendOuterPos_.resize(dims);
            ptr_map_interface_info->at(pair_interface)->vec_outer_coarse2fine_.resize(num_of_layers_coarse2fine);
            ptr_map_interface_info->at(pair_interface)->vec_inner_coarse2fine_.resize(num_of_layers_coarse2fine);
        } else {
#ifdef DEBUG_CHECK_GRID
            if (i_level > 0) {
                if (ptr_map_interface_info->at(pair_interface)->k0ExtendOuterNeg_.size() != dims) {
                    LogManager::LogError("size of k0ExtendOuterNeg_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendOuterNeg_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (ptr_map_interface_info->at(pair_interface)->k0ExtendOuterPos_.size() != dims) {
                    LogManager::LogError("size of k0ExtendOuterPos_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendOuterPos_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (ptr_map_interface_info->at(pair_interface)->k0ExtendInnerNeg_.size() != dims) {
                    LogManager::LogError("size of k0ExtendInnerNeg_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendInnerNeg_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (ptr_map_interface_info->at(pair_interface)->k0ExtendInnerPos_.size() != dims) {
                    LogManager::LogError("size of k0ExtendInnerPos_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendInnerPos_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (ptr_map_interface_info->at(pair_interface)->vec_outer_coarse2fine_.size()
                 != num_of_layers_coarse2fine) {
                    LogManager::LogError("size of vec_outer_coarse2fine_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->vec_outer_coarse2fine_.size())
                     + " is not equal to the number of layers " + std::to_string(num_of_layers_coarse2fine)
                     + " in IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (ptr_map_interface_info->at(pair_interface)->vec_inner_coarse2fine_.size()
                 != num_of_layers_coarse2fine) {
                    LogManager::LogError("size of vec_inner_coarse2fine_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->vec_inner_coarse2fine_.size())
                     + " is not equal to the number of layers " + std::to_string(num_of_layers_coarse2fine)
                     + " in IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
#endif  // DEBUG_CHECK_GRID
        }

        InterfaceLayerInfo* ptr_interface = ptr_map_interface_info->at(pair_interface).get();

        // broadcast number of extending layers
        if (i_level > 0) {
            MPI_Bcast(ptr_interface->k0ExtendInnerNeg_.data(),
             static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
            MPI_Bcast(ptr_interface->k0ExtendInnerPos_.data(),
             static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
            MPI_Bcast(ptr_interface->k0ExtendOuterNeg_.data(),
             static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
            MPI_Bcast(ptr_interface->k0ExtendOuterPos_.data(),
             static_cast<int>(dims), MPI_AMR_UINT_TYPE, 0, MPI_COMM_WORLD);
        }

        //  send and receive outer layers
        IniSendNReiveOneLayerRefinementInterface(flag0, outmost_for_all_ranks,
            &ptr_interface->vec_outer_coarse2fine_.at(num_of_layers_coarse2fine - 2));

        //  send and receive inner layers
        IniSendNReiveOneLayerRefinementInterface(flag0, outmost_for_all_ranks,
            &ptr_interface->vec_inner_coarse2fine_.at(num_of_layers_coarse2fine - 2));
    }

    MPI_Type_free(&mpi_interface_index_type);
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
