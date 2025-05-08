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
#include <numeric>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
#include "io/debug_write.h"
// struct to store the index of grid interfaces for mpi communication
struct InterfaceIndexForMpi {
    rootproject::DefInt i_level;
    rootproject::amrproject::ECriterionType enum_criterion;
    rootproject::DefInt criterion_count;
    int layer_indicator;
    int layer_count;

    bool operator==(const InterfaceIndexForMpi& other) const {
        return i_level == other.i_level &&
               enum_criterion == other.enum_criterion &&
               criterion_count == other.criterion_count &&
               layer_indicator == other.layer_indicator &&
               layer_count == other.layer_count;
    }
};
template <>
struct std::hash<InterfaceIndexForMpi> {
    std::size_t operator()(const InterfaceIndexForMpi& k) const {
        return ((std::hash<rootproject::DefInt>()(k.i_level) ^
                (std::hash<int>()(static_cast<int>(k.enum_criterion)) << 1)) >> 1) ^
               (std::hash<rootproject::DefInt>()(k.criterion_count) << 1) ^
               (std::hash<int>()(k.layer_indicator) << 2) ^
               (std::hash<int>()(k.layer_count) << 3);
    }
};
namespace rootproject {
namespace amrproject {
struct CriterionIndexForMpi {
    ECriterionType enum_criterion;
    DefInt index_criterion;
    DefInt index_creator;
};
/**
 * @brief function to create and commit MPI <ECriterionType, DefInt, DefInt> types.
 * @param[out] ptr_mpi_tracking_index_type pointer to the <ECriterionType, DefInt, DefInt> types.
 * @return if the MPI types were created and committed successfully.
 */
int MpiManager::CreateAndCommitCriterionIndexType(MPI_Datatype* ptr_mpi_tracking_index_type) const {
    MPI_Datatype mpi_tracking_enum_type;
    MPI_Type_contiguous(sizeof(ECriterionType), MPI_BYTE, &mpi_tracking_enum_type);
    MPI_Type_commit(&mpi_tracking_enum_type);
    // the last MPI_UNSIGNED store the index (k0IndexCreator) of the creator for tracking grid nodes
    MPI_Datatype mpi_create_type[3] = {mpi_tracking_enum_type, GetMpiIntType(), GetMpiIntType()};
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
    if (sizeof(int) + set_nodes.size() *key_size
        > static_cast<DefSizet>((std::numeric_limits<int>::max)())) {
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
    DefSFCodeToUint key_code;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        ptr_map_tracking->insert({static_cast<DefSFBitset>(key_code), tracking_node_instance});
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
void MpiManager::IniSendNReceiveTracking(const DefInt dims, const DefInt i_level,
    const std::vector<DefSFBitset>& bitset_max, const SFBitsetAuxInterface& bitset_aux,
    const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator,
    std::map<std::pair<ECriterionType, DefInt>, std::shared_ptr<TrackingGridInfoInterface>>*
    const ptr_map_tracking_info) const {
    MPI_Datatype mpi_tracking_index_type;
    CreateAndCommitCriterionIndexType(&mpi_tracking_index_type);
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (DefSizet i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_aux.SFBitsetToSFCode(bitset_max.at(i));
        }
    }

    int num_max = static_cast<int>(ull_max.size());
    DefInt num_tracking = DefInt(ptr_map_tracking_info->size());
    MPI_Bcast(&num_tracking, 1, GetMpiIntType(), 0, MPI_COMM_WORLD);
    CriterionIndexForMpi tracking_index_tmp;
    std::vector<std::pair<ECriterionType, DefInt>> vec_tracking_indices;
    if (rank_id == 0) {
        // save tracking indices stored in unordered_map to vector
        for (const auto& iter_tracking : *ptr_map_tracking_info) {
            vec_tracking_indices.emplace_back(iter_tracking.first);
        }
    }

    for (DefInt i_tracking = 0; i_tracking < num_tracking; ++i_tracking) {
        // broadcast the current index of the tracking grid
        if (rank_id == 0) {
            tracking_index_tmp.enum_criterion = vec_tracking_indices.at(i_tracking).first;
            tracking_index_tmp.index_criterion = vec_tracking_indices.at(i_tracking).second;
            tracking_index_tmp.index_creator =
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0IndexCreator;
#ifdef DEBUG_CHECK_GRID
            if (static_cast<DefInt>(ptr_map_tracking_info->at(
                vec_tracking_indices.at(i_tracking))->k0ExtendOuterNeg_.size()) != dims) {
                LogManager::LogError("size of k0ExtendOuterNeg_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterNeg_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            if (static_cast<DefInt>(ptr_map_tracking_info->at(
                vec_tracking_indices.at(i_tracking))->k0ExtendOuterPos_.size()) != dims) {
                LogManager::LogError("size of k0ExtendOuterPos_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendOuterPos_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            if (static_cast<DefInt>(ptr_map_tracking_info->at(
                vec_tracking_indices.at(i_tracking))->k0ExtendInnerNeg_.size()) != dims) {
                LogManager::LogError("size of k0ExtendInnerNeg_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerNeg_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            if (static_cast<DefInt>(ptr_map_tracking_info->at(
                vec_tracking_indices.at(i_tracking))->k0ExtendInnerPos_.size()) != dims) {
                LogManager::LogError("size of k0ExtendInnerPos_ " + std::to_string(
                ptr_map_tracking_info->at(vec_tracking_indices.at(i_tracking))->k0ExtendInnerPos_.size())
                + " is not equal to the dimension " + std::to_string(dims) + " in IniSendNReceiveTracking in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
#endif  // DEBUG_CHECK_GRID
        }
        MPI_Bcast(&tracking_index_tmp, 1, mpi_tracking_index_type, 0, MPI_COMM_WORLD);
        std::pair<ECriterionType, DefInt> pair_tracking = {
            tracking_index_tmp.enum_criterion, tracking_index_tmp.index_criterion};
        if (rank_id > 0) {
            DefInt index_creator = tracking_index_tmp.index_creator;
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
           static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendInnerPos_.data(),
           static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendOuterNeg_.data(),
           static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
        MPI_Bcast(ptr_tracking_info->k0ExtendOuterPos_.data(),
           static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);

        // broadcast type information
        DefInt uint_type[2] = {0, 0};
        if (rank_id == 0) {
            uint_type[0] = ptr_tracking_info->computational_cost_;
            uint_type[1] = static_cast<DefInt>(ptr_tracking_info->grid_extend_type_);
        }
        MPI_Bcast(uint_type, 2, GetMpiIntType(), 0, MPI_COMM_WORLD);
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
            std::vector<DefSFCodeToUint>::iterator iter_index;
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            for (const auto& iter_node : ptr_tracking_info->map_tracking_node_) {
                background_code = bitset_aux.SFBitsetToSFCode(
                    bitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first));
                iter_index = std::lower_bound(ull_max.begin(),
                    ull_max.end(), background_code);
                index = static_cast<int>(iter_index - ull_max.begin());
                key_ull = bitset_aux.SFBitsetToSFCode(iter_node.first);
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
 * @param[out] ptr_buffer_size pointer to size of the buffer in bytes.
 * @return unique pointer to a char array to store the serialized data.
 */
std::unique_ptr<char[]> MpiManager::IniSerializeRefinementInterfaceNode(
    const DefMap<DefInt>& interface_nodes, int* const ptr_buffer_size) const {
    int key_size = sizeof(DefSFBitset);
    int num_nodes = 0;
    if  (sizeof(int) + interface_nodes.size() *key_size
            > static_cast<DefSizet>((std::numeric_limits<int>::max)())) {
        LogManager::LogError("size of the buffer is greater than"
         " the maximum of int in GridInfoInterface::SerializeTrackingNode in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        num_nodes = static_cast<int>(interface_nodes.size());
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
    for (const auto& iter : interface_nodes) {
        // convert bitset into unsigned long long and save it to buffer
        std::memcpy(ptr_buffer + position, &iter.first, key_size);
        position += key_size;
    }
    return buffer;
}
/**
 * @brief function to deserialize tracking node data and insert it into the associated map.
 * @param[in] criterion_type type of criterion.
 * @param[in]  index_criterion index of the criterion.
 * @param[in]  buffer buffer containing the tracking node data.
 */
void MpiManager::IniDeserializeRefinementInterfaceNode(const DefInt flag0,
     const std::unique_ptr<char[]>& buffer, DefMap<DefInt>* ptr_map_interface_layer) const {
    char* ptr_buffer = buffer.get();
    int key_size = sizeof(DefSFBitset);
    // number of nodes
    int num_nodes;
    int position = 0;
    std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
    // deserialize data stored in buffer
    position += sizeof(int);
    DefSFCodeToUint key_code;
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
        std::memcpy(&key_code, ptr_buffer + position, key_size);
        position += key_size;
        ptr_map_interface_layer->insert({key_code, flag0});
    }
}
/**
 * @brief  function to send and receive nodes in a given level of the interface layer.
 * @param[in] i_level  level of the interface layer to send/receive.
 * @param[in] flag0 flag for initialization.
 * @param[in] outmost_for_all_ranks nodes on the outmost fine to coarse interface for all ranks stored on rank 0. 
 * @param[in, out] ptr_map_interface_layer pointer to nodes on the refinement interface layer.
 */
void MpiManager::IniSendNReiveOneLayerRefinementInterface(
    const DefInt flag0, const DefMap<std::set<int>>& outmost_for_all_ranks,
    DefMap<DefInt>* const ptr_map_interface_layer) const {
    int rank_id = rank_id_, num_ranks = num_of_ranks_;
    std::vector<std::vector<DefMap<DefInt>>> vec_nodes_ranks(num_ranks),
        vec_partition_nodes_ranks(num_ranks);
    std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0),
     i_chunk_partition_each_rank(num_ranks, -1), i_partition_counts(num_ranks, 0);
    int max_buffer = (std::numeric_limits<int>::max)() / sizeof(DefSFBitset) - 1;
    vec_nodes_ranks.at(0).push_back({});
    if (rank_id == 0) {
        std::vector<DefMap<DefInt>> partition_interface_current(num_ranks);
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
                vec_ptr_buffer.at(i_chunk) = IniSerializeRefinementInterfaceNode(
                    vec_nodes_ranks.at(iter_rank).at(i_chunk), &buffer_size_send);
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
void MpiManager::IniSendNReceiveCoarse2Fine0Interface(const DefInt dims,
    const DefInt i_level, const DefInt num_of_layers_coarse2fine, const DefInt flag0,
    const DefMap<std::set<int>>& outmost_for_all_ranks,
    std::map<std::pair<ECriterionType, DefInt>,
    std::shared_ptr<InterfaceLayerInfo>>* const ptr_map_interface_info) const {
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitCriterionIndexType(&mpi_interface_index_type);
    int rank_id = rank_id_;
    DefInt num_interface = DefInt(ptr_map_interface_info->size());
    MPI_Bcast(&num_interface, 1, GetMpiIntType(), 0, MPI_COMM_WORLD);
    std::vector<std::pair<ECriterionType, DefInt>> vec_interface_indices;
    if (rank_id == 0) {
        // save interface indices stored in unordered_map to vector
        for (const auto& iter_interface : *ptr_map_interface_info) {
            vec_interface_indices.emplace_back(iter_interface.first);
        }
    }
    CriterionIndexForMpi interface_index_tmp;
    for (DefInt i_interface = 0; i_interface < num_interface; ++i_interface) {
        // broadcast the current index of the interfaces
        if (rank_id == 0) {
            interface_index_tmp.enum_criterion = vec_interface_indices[i_interface].first;
            interface_index_tmp.index_criterion =
                static_cast<uint64_t>(vec_interface_indices[i_interface].second);
            interface_index_tmp.index_creator = 0;
        }
        MPI_Bcast(&interface_index_tmp, 1, mpi_interface_index_type, 0, MPI_COMM_WORLD);
        std::pair<ECriterionType, DefInt> pair_interface = {
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
                if (static_cast<DefInt>(ptr_map_interface_info->at(
                    pair_interface)->k0ExtendOuterNeg_.size()) != dims) {
                    LogManager::LogError("size of k0ExtendOuterNeg_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendOuterNeg_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in MpiManager::IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (static_cast<DefInt>(ptr_map_interface_info->at(
                    pair_interface)->k0ExtendOuterPos_.size()) != dims) {
                    LogManager::LogError("size of k0ExtendOuterPos_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendOuterPos_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in MpiManager::IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (static_cast<DefInt>(ptr_map_interface_info->at(
                    pair_interface)->k0ExtendInnerNeg_.size()) != dims) {
                    LogManager::LogError("size of k0ExtendInnerNeg_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendInnerNeg_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in MpiManager::IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (static_cast<DefInt>(ptr_map_interface_info->at(
                    pair_interface)->k0ExtendInnerPos_.size()) != dims) {
                    LogManager::LogError("size of k0ExtendInnerPos_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->k0ExtendInnerPos_.size())
                     + " is not equal to the dimension " + std::to_string(dims)
                     + " in MpiManager::IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (static_cast<DefInt>(ptr_map_interface_info->at(pair_interface)->vec_outer_coarse2fine_.size())
                    != num_of_layers_coarse2fine) {
                    LogManager::LogError("size of vec_outer_coarse2fine_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->vec_outer_coarse2fine_.size())
                     + " is not equal to the number of layers " + std::to_string(num_of_layers_coarse2fine)
                     + " in MpiManager::IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                if (static_cast<DefInt>(ptr_map_interface_info->at(pair_interface)->vec_inner_coarse2fine_.size())
                    != num_of_layers_coarse2fine) {
                    LogManager::LogError("size of vec_inner_coarse2fine_ " + std::to_string(
                     ptr_map_interface_info->at(pair_interface)->vec_inner_coarse2fine_.size())
                     + " is not equal to the number of layers " + std::to_string(num_of_layers_coarse2fine)
                     + " in MpiManager::IniSendNReceiveCoarse2Fine0Interface in "
                     + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
#endif  // DEBUG_CHECK_GRID
        }

        InterfaceLayerInfo* ptr_interface = ptr_map_interface_info->at(pair_interface).get();

        // broadcast number of extending layers
        if (i_level > 0) {
            MPI_Bcast(ptr_interface->k0ExtendInnerNeg_.data(),
             static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
            MPI_Bcast(ptr_interface->k0ExtendInnerPos_.data(),
             static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
            MPI_Bcast(ptr_interface->k0ExtendOuterNeg_.data(),
             static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
            MPI_Bcast(ptr_interface->k0ExtendOuterPos_.data(),
             static_cast<int>(dims), GetMpiIntType(), 0, MPI_COMM_WORLD);
        }

        //  send and receive outer layers
        IniSendNReiveOneLayerRefinementInterface(flag0, outmost_for_all_ranks,
            &ptr_interface->vec_outer_coarse2fine_.at(0));

        //  send and receive inner layers
        IniSendNReiveOneLayerRefinementInterface(flag0, outmost_for_all_ranks,
            &ptr_interface->vec_inner_coarse2fine_.at(0));
    }

    MPI_Type_free(&mpi_interface_index_type);
}
/**
 * @brief function to write grid nodes to a checkpoint file.
 * @param[in] file_name file name to write the file.
 * @param[in] num_of_nodes number of grid nodes at each level on current rank.
 * @param[in] map_node_level nodes at which levels on current rank.
 * @param[in] vec_grid_info grid information at all levels.
 */
void MpiManager::WriteCheckPointGridNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_nodes,
    const std::map<DefSFCodeToUint, BackgroundLoadData>& map_node_level,
    const std::vector<std::shared_ptr<GridInfoInterface>>& vec_grid_info) const {
    DefAmrLUint num_of_nodes_for_all_levels = std::accumulate(num_of_nodes.begin(), num_of_nodes.end(), 0);
    std::vector<DefAmrLUint> num_of_nodes_each_rank(num_of_ranks_, 0);
    MPI_Allgather(&num_of_nodes_for_all_levels, 1, GetMpiAmrLUintType(),
        num_of_nodes_each_rank.data(), 1, GetMpiAmrLUintType(), MPI_COMM_WORLD);

    MPI_File mpi_file;
    MPI_Status mpi_status;
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
    if (rank_id_ == 0) {
        DefAmrLUint total_nodes = std::accumulate(
            num_of_nodes_each_rank.begin(), num_of_nodes_each_rank.end(), 0);
        MPI_File_write_at_all(mpi_file, 0, &total_nodes, 1, GetMpiAmrLUintType(), &mpi_status);
    } else {
        MPI_File_write_at_all(mpi_file, 0, nullptr, 0, GetMpiAmrLUintType(), &mpi_status);
    }

    // compute offset for each rank
    MPI_Offset offset = sizeof(DefAmrLUint);
    int data_size = 0;
    for (const auto& iter_grid : vec_grid_info) {
        if (data_size < iter_grid->GetSizeOfGridNodeInfoForCheckPoint()) {
            data_size = iter_grid->GetSizeOfGridNodeInfoForCheckPoint();
        }
    }
    const int sfbitset_size = sizeof(DefSFBitset), int_size = sizeof(DefInt);
    const int node_size = sfbitset_size + int_size + data_size;
    for (int i_rank = 0; i_rank < rank_id_; ++i_rank) {
        if (offset > (std::numeric_limits<MPI_Offset>::max)()
            - node_size* static_cast<MPI_Offset>(num_of_nodes_each_rank.at(i_rank))) {
            LogManager::LogError("offset exceeds the maximum of MPI_Offset in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        offset += node_size * static_cast<MPI_Offset>(num_of_nodes_each_rank.at(i_rank));
    }

    DefInt max_level = vec_grid_info.back()->GetGridLevel();
    if (max_level != static_cast<DefInt>(num_of_nodes.size() - 1)) {
        LogManager::LogError("max_level is not equal to the side of input (num_of_nodes) which is"
            " number of nodes in each level in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    const DefSizet max_chunk_size = static_cast<DefSizet>((std::numeric_limits<int>::max)());
    std::vector<char> chunk_buffer;
    chunk_buffer.reserve(max_chunk_size);
    DefSizet chunk_index = 0;
    MPI_Offset current_offset = offset;
    std::vector<DefSFBitset> sfbitset_cell;
    DefInt node_load = 0;
    for (const auto& iter : map_node_level) {
        node_load = 0;
        for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
            if (iter.second.level_bitset_.Get(i_level)) {
                const auto& grid_info = *vec_grid_info.at(i_level).get();
                const SFBitsetAuxInterface& sfbitset_aux =  *grid_info.GetPtrSFBitsetAux();
                sfbitset_aux.SFBitsetHigherLevelInACell(i_level, iter.first, &sfbitset_cell);
                const DefInt node_cost = grid_info.GetComputationalCost();
                for (const auto& iter_cell : sfbitset_cell) {
                    if (grid_info.map_grid_node_.find(iter_cell) != grid_info.map_grid_node_.end()) {
                        node_load += node_cost;
                        // Check if there’s enough space left in the chunk buffer
                        if (chunk_index + node_size > max_chunk_size) {
                            MPI_File_write_at_all(mpi_file, current_offset, chunk_buffer.data(),
                                static_cast<int>(chunk_index), MPI_BYTE, &mpi_status);
                            current_offset += chunk_index;
                            chunk_index = 0;
                            chunk_buffer.clear();
                        }
                        // Append to buffer
                        chunk_buffer.resize(chunk_index + node_size);
                        std::memcpy(chunk_buffer.data() + chunk_index, &iter_cell, sfbitset_size);
                        chunk_index += sfbitset_size;
                        std::memcpy(chunk_buffer.data() + chunk_index, &i_level, int_size);
                        chunk_index += int_size;
                        grid_info.map_grid_node_.at(iter_cell)->CopyANodeToBufferForCheckpoint(
                            chunk_buffer.data() + chunk_index);
                        chunk_index += data_size;
                    }
                }
            }
        }
        if (node_load != iter.second.total_load_) {
            const SFBitsetAuxInterface& sfbitset_aux = *vec_grid_info.at(0)->GetPtrSFBitsetAux();
            DefInt dims = vec_grid_info.at(0)->GetPtrToParentGridManager()->k0GridDims_;
            std::vector<DefReal> coordinates(dims);
            sfbitset_aux.SFBitsetComputeCoordinateVir(
                iter.first, vec_grid_info.at(0)->GetGridSpace(), &coordinates);
            if (dims == 2) {
                LogManager::LogError("Computed load is not equal to the input load for node ("
                    + std::to_string(coordinates[kXIndex]) + ", "
                    + std::to_string(coordinates[kYIndex]) + ") in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            } else if (dims == 3) {
                LogManager::LogError("Computed load is not equal to the input load for node ("
                    + std::to_string(coordinates[kXIndex]) + ", "+ std::to_string(coordinates[kYIndex])
                    + ", "+ std::to_string(coordinates[kZIndex]) +  ") in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
        }
    }
    if (chunk_index > 0) {
        MPI_File_write_at_all(mpi_file, current_offset,
            chunk_buffer.data(), static_cast<int>(chunk_index), MPI_BYTE, &mpi_status);
    }
    MPI_File_close(&mpi_file);
}
/**
 * @brief function to read grid nodes from a checkpoint file.
 * @param[in] file_name file name to read the file.
 * @param[in] num_of_nodes number of nodes on each rank.
 * @param[in, out] ptr_vec_grid_info pointer to class storting grid information on current rank.
 */
void MpiManager::ReadCheckPointGridNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_nodes,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const {

    // compute offset for each rank
    MPI_Offset offset = sizeof(DefAmrLUint);
    int data_size = 0;
    for (const auto& iter_grid : *ptr_vec_grid_info) {
        if (data_size < iter_grid->GetSizeOfGridNodeInfoForCheckPoint()) {
            data_size = iter_grid->GetSizeOfGridNodeInfoForCheckPoint();
        }
    }
    const int sfbitset_size = sizeof(DefSFBitset), int_size = sizeof(DefInt);
    int node_size = sfbitset_size + int_size + data_size;
    for (int i_rank = 0; i_rank < rank_id_; ++i_rank) {
        if (offset > (std::numeric_limits<MPI_Offset>::max)()
            - node_size* static_cast<MPI_Offset>(num_of_nodes.at(i_rank))) {
            LogManager::LogError("offset exceeds the maximum of MPI_Offset in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        offset += node_size * static_cast<MPI_Offset>(num_of_nodes.at(i_rank));
    }
    const std::size_t max_chunk_size = static_cast<std::size_t>(INT_MAX);
    std::size_t total_size = num_of_nodes.at(rank_id_) * node_size;

    int max_chunks = 0;
    for (int i = 0; i < num_of_ranks_; ++i) {
        std::size_t size_i = static_cast<std::size_t>(num_of_nodes.at(i)) * node_size;
        int chunks_i = static_cast<int>((size_i + max_chunk_size - 1) / max_chunk_size);
        if (chunks_i > max_chunks) {
            max_chunks = chunks_i;
        }
    }

    MPI_Offset current_offset = offset;
    std::vector<char> chunk_buffer;
    chunk_buffer.resize((std::min)(max_chunk_size, total_size));
    std::size_t remaining = total_size;
    MPI_File mpi_file;
    MPI_Status mpi_status;
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
    for (int chunk_idx = 0; chunk_idx < max_chunks; ++chunk_idx) {
        std::size_t this_chunk_size = 0;
        void* buffer_ptr = nullptr;

        if (remaining > 0) {
            this_chunk_size = (std::min)(max_chunk_size, remaining);
            buffer_ptr = chunk_buffer.data();
        }

        int mpi_chunk_size = static_cast<int>(this_chunk_size);

        // All ranks call collectively, even if no data to read
        MPI_File_read_at_all(mpi_file, current_offset, buffer_ptr, mpi_chunk_size, MPI_BYTE, &mpi_status);

        if (mpi_chunk_size > 0) {
            std::size_t index = 0;
            while (index < this_chunk_size) {
                DefSFBitset sfbitset;
                std::memcpy(&sfbitset, chunk_buffer.data() + index, sfbitset_size);
                index += sfbitset_size;
                DefInt i_level_read;
                std::memcpy(&i_level_read, chunk_buffer.data() + index, int_size);
                index += int_size;
                if (static_cast<int>(ptr_vec_grid_info->size()) <= i_level_read) {
                    LogManager::LogError("Grid level exceeds available levels in " +
                        std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                } else {
                    auto& grid_info = ptr_vec_grid_info->at(i_level_read);
                    if (grid_info->GetPtrToParentGridManager()->InstantiateGridNode(sfbitset, grid_info.get())) {
                        auto& node = grid_info->map_grid_node_.at(sfbitset);
                        node->ReadANodeFromBufferForCheckpoint(chunk_buffer.data() + index);
                        index += data_size;
                    } else {
                        LogManager::LogError("Failed to instantiate grid node from checkpoint in " +
                            std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                    }
                }
            }

            current_offset += mpi_chunk_size;
            remaining -= mpi_chunk_size;
        }
    }

    MPI_File_close(&mpi_file);
}
/**
 * @brief function to create and commit MPI InterfaceIndexForMpi types.
 * @param[out] ptr_mpi_interface_index_type pointer to the InterfaceIndexForMpi types.
 * @return if the MPI types were created and committed successfully.
 */
int MpiManager::CreateAndCommitInterfaceIndexType(MPI_Datatype* ptr_mpi_interface_index_type) const {
    MPI_Datatype mpi_interface_enum_type;
    MPI_Type_contiguous(sizeof(ECriterionType), MPI_BYTE, &mpi_interface_enum_type);
    MPI_Type_commit(&mpi_interface_enum_type);
    // mpi type for: i_level(DefInt), criterion type(ECriterionType),
    // criterion count (DefInt), layer indicator(int), layer count (int)
    MPI_Datatype mpi_create_type[5] = {GetMpiIntType(), mpi_interface_enum_type, GetMpiIntType(), MPI_INT, MPI_INT};
    int mpi_block_length[5] = {1, 1, 1, 1, 1};
    MPI_Aint mpi_disp[5];
    mpi_disp[0] = offsetof(InterfaceIndexForMpi, i_level);
    mpi_disp[1] = offsetof(InterfaceIndexForMpi, enum_criterion);
    mpi_disp[2] = offsetof(InterfaceIndexForMpi, criterion_count);
    mpi_disp[3] = offsetof(InterfaceIndexForMpi, layer_indicator);
    mpi_disp[4] = offsetof(InterfaceIndexForMpi, layer_count);
    MPI_Type_create_struct(5, mpi_block_length, mpi_disp, mpi_create_type, ptr_mpi_interface_index_type);
    MPI_Type_commit(ptr_mpi_interface_index_type);
    MPI_Type_free(&mpi_interface_enum_type);
    return 0;
}
/**
 * @brief function to write interface nodes to a checkpoint file.
 * @param[in] file_name file name to write the file.
 * @param[in] num_of_nodes number of interface nodes at each level on current rank.
 * @param[in] map_node_level nodes at which levels on current rank.
 * @param[in] vec_grid_info grid information at all levels.
 */
void MpiManager::WriteCheckPointInterfaceNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_nodes,
    const std::map<DefSFCodeToUint, BackgroundLoadData>& map_node_level,
    const std::vector<std::shared_ptr<GridInfoInterface>>& vec_grid_info) const {
    DefAmrLUint num_of_nodes_for_all_levels = std::accumulate(num_of_nodes.begin(), num_of_nodes.end(), 0);
    std::vector<DefAmrLUint> num_of_nodes_each_rank(num_of_ranks_, 0);
    MPI_Allgather(&num_of_nodes_for_all_levels, 1, GetMpiAmrLUintType(),
        num_of_nodes_each_rank.data(), 1, GetMpiAmrLUintType(), MPI_COMM_WORLD);

    DefInt max_level = vec_grid_info.back()->GetGridLevel();
    if (max_level != static_cast<DefInt>(num_of_nodes.size() - 1)) {
        LogManager::LogError("max_level is not equal to the side of input (num_of_nodes) which is"
            " number of nodes in each level in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    MPI_File mpi_file;
    MPI_Status mpi_status;
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
    int num_indices = 0;
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitInterfaceIndexType(&mpi_interface_index_type);
    std::unordered_map<InterfaceIndexForMpi, int> index_map;
    int inner_fine2coarse = 1, outer_fine2coarse = 2, inner_coarse2fine = 3, outer_coarse2fine = 4;
    auto add_entries = [&num_indices, &index_map](const int layer_indicator,
        const DefInt i_level, const ECriterionType enum_criterion, const DefInt criterion_count,
        const std::vector<DefMap<DefInt>>& layer_vec,
        std::vector<InterfaceIndexForMpi>* const ptr_vec_interface_index) {
        const int layer_size = static_cast<int>(layer_vec.size());
        for (int i_layer = 0; i_layer < layer_size; ++i_layer) {
            InterfaceIndexForMpi entry;
            entry.i_level = i_level;
            entry.enum_criterion = enum_criterion;
            entry.criterion_count = criterion_count;
            entry.layer_indicator = layer_indicator;
            entry.layer_count = i_layer;
            if (ptr_vec_interface_index != nullptr) {
                ptr_vec_interface_index->push_back(entry);
            }
            index_map.insert({entry, num_indices});
            if (num_indices + 1 > (std::numeric_limits<int>::max)()) {
                LogManager::LogError("num_indices exceeds the maximum of int in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            ++num_indices;
        }
    };

    MPI_Offset offset = sizeof(int);
    if (rank_id_ == 0) {
        std::vector<InterfaceIndexForMpi> vec_interface_index;
        for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
            const auto& grid_info = *vec_grid_info.at(i_level).get();
            for (const auto& iter_interface : grid_info.map_ptr_interface_layer_info_) {
                const auto& info = *iter_interface.second;
                auto enum_criterion = iter_interface.first.first;
                auto criterion_count = iter_interface.first.second;

                add_entries(inner_fine2coarse, i_level, enum_criterion,
                    criterion_count, info.vec_inner_fine2coarse_, &vec_interface_index);
                add_entries(outer_fine2coarse, i_level, enum_criterion,
                    criterion_count, info.vec_outer_fine2coarse_, &vec_interface_index);
                add_entries(inner_coarse2fine, i_level, enum_criterion,
                    criterion_count, info.vec_inner_coarse2fine_, &vec_interface_index);
                add_entries(outer_coarse2fine, i_level, enum_criterion,
                    criterion_count, info.vec_outer_coarse2fine_, &vec_interface_index);
            }
        }

        MPI_File_write_at_all(mpi_file, 0, &num_indices, 1, MPI_INT, &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, vec_interface_index.data(),
            num_indices, mpi_interface_index_type, &mpi_status);
    } else {
        for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
            const auto& grid_info = *vec_grid_info.at(i_level).get();
            for (const auto& iter_interface : grid_info.map_ptr_interface_layer_info_) {
                const auto& info = *iter_interface.second;
                auto enum_criterion = iter_interface.first.first;
                auto criterion_count = iter_interface.first.second;

                add_entries(inner_fine2coarse, i_level, enum_criterion,
                    criterion_count, info.vec_inner_fine2coarse_, nullptr);
                add_entries(outer_fine2coarse, i_level, enum_criterion,
                    criterion_count, info.vec_outer_fine2coarse_, nullptr);
                add_entries(inner_coarse2fine, i_level, enum_criterion,
                    criterion_count, info.vec_inner_coarse2fine_, nullptr);
                add_entries(outer_coarse2fine, i_level, enum_criterion,
                    criterion_count, info.vec_outer_coarse2fine_, nullptr);
            }
        }
        MPI_File_write_at_all(mpi_file, 0, nullptr, 0, MPI_INT, &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, nullptr, 0, mpi_interface_index_type, &mpi_status);
    }

    // compute offset for each rank
    int mpi_type_size = 0;
    MPI_Type_size(mpi_interface_index_type, &mpi_type_size);
    offset += num_indices * mpi_type_size;
    MPI_Type_free(&mpi_interface_index_type);

    const int sfbitset_size = sizeof(DefSFBitset), int_size = sizeof(int);
    const int node_size = sfbitset_size + int_size;
    for (int i_rank = 0; i_rank < rank_id_; ++i_rank) {
        if (offset > (std::numeric_limits<MPI_Offset>::max)()
            - node_size* static_cast<MPI_Offset>(num_of_nodes_each_rank.at(i_rank))) {
            LogManager::LogError("offset exceeds the maximum of MPI_Offset in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        offset += node_size * static_cast<MPI_Offset>(num_of_nodes_each_rank.at(i_rank));
    }

    const DefSizet max_chunk_size = static_cast<DefSizet>((std::numeric_limits<int>::max)());
    std::vector<char> chunk_buffer;
    chunk_buffer.reserve(max_chunk_size);
    DefSizet chunk_index = 0;
    MPI_Offset current_offset = offset;
    std::vector<DefSFBitset> sfbitset_cell;
    auto add_node = [&mpi_file, &mpi_status, &current_offset,
        &index_map, &chunk_index](const int layer_indicator,
        const DefInt i_level, const ECriterionType enum_criterion, const DefInt criterion_count,
        const DefSFBitset& sfbitset_in, const std::vector<DefMap<DefInt>>& layer_vec,
        std::vector<char>* ptr_chunk_buffer) {
        const int layer_size = static_cast<int>(layer_vec.size());
        for (int i_layer = 0; i_layer < layer_size; ++i_layer) {
            if (layer_vec.at(i_layer).find(sfbitset_in) != layer_vec.at(i_layer).end()) {
                InterfaceIndexForMpi entry;
                entry.i_level = i_level;
                entry.enum_criterion = enum_criterion;
                entry.criterion_count = criterion_count;
                entry.layer_indicator = layer_indicator;
                entry.layer_count = i_layer;
                if (chunk_index + node_size > max_chunk_size) {
                    MPI_File_write_at_all(mpi_file, current_offset, ptr_chunk_buffer->data(),
                        static_cast<int>(chunk_index), MPI_BYTE, &mpi_status);
                    current_offset += chunk_index;
                    chunk_index = 0;
                    ptr_chunk_buffer->clear();
                }

                ptr_chunk_buffer->resize(chunk_index + node_size);
                std::memcpy(ptr_chunk_buffer->data() + chunk_index, &sfbitset_in, sfbitset_size);
                chunk_index += sfbitset_size;
                if (index_map.find(entry) == index_map.end()) {
                    LogManager::LogError("entry not found in interface indices in "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                } else {
                    int entry_index = index_map.at(entry);
                    std::memcpy(ptr_chunk_buffer->data() + chunk_index, &entry_index, int_size);
                }
                chunk_index += int_size;
                break;
            }
        }
    };
    for (const auto& iter : map_node_level) {
        for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
            if (iter.second.level_bitset_.Get(i_level)) {
                const auto& grid_info = *vec_grid_info.at(i_level).get();
                const SFBitsetAuxInterface& sfbitset_aux =  *grid_info.GetPtrSFBitsetAux();
                sfbitset_aux.SFBitsetHigherLevelInACell(i_level, iter.first, &sfbitset_cell);
                for (const auto& iter_cell : sfbitset_cell) {
                    const auto& grid_info = *vec_grid_info.at(i_level).get();
                    for (const auto& iter_interface : grid_info.map_ptr_interface_layer_info_) {
                        const auto& info = *iter_interface.second;
                        auto enum_criterion = iter_interface.first.first;
                        auto criterion_count = iter_interface.first.second;
                        add_node(inner_fine2coarse, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_inner_fine2coarse_, &chunk_buffer);
                        add_node(outer_fine2coarse, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_outer_fine2coarse_, &chunk_buffer);
                        add_node(inner_coarse2fine, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_inner_coarse2fine_, &chunk_buffer);
                        add_node(outer_coarse2fine, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_outer_coarse2fine_, &chunk_buffer);
                    }
                }
            }
        }
    }
    if (chunk_index > 0) {
        MPI_File_write_at_all(mpi_file, current_offset,
            chunk_buffer.data(), static_cast<int>(chunk_index), MPI_BYTE, &mpi_status);
    }
    MPI_File_close(&mpi_file);
}
/**
 * @brief function to read interface nodes from a checkpoint file.
 * @param[in] file_name file name to read the file.
 * @param[in] num_of_nodes number of interface nodes on each rank stored on rank 0.
 * @param[out] ptr_vec_grid_info pointer to class storting grid information on current rank.
 */
void MpiManager::ReadCheckPointInterfaceNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_nodes,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const {
    MPI_File mpi_file;
    MPI_Status mpi_status;
    int num_indices = 0;
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitInterfaceIndexType(&mpi_interface_index_type);

    int mpi_err = MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
    if (mpi_err != MPI_SUCCESS) {
        LogManager::LogError("Failed to open file for reading: " + file_name);
    }

    MPI_File_read_at_all(mpi_file, 0, &num_indices, 1, MPI_INT, &mpi_status);

    std::vector<InterfaceIndexForMpi> vec_interface_index(num_indices);
    MPI_Offset offset = sizeof(int);
    MPI_File_read_at_all(mpi_file, offset, vec_interface_index.data(),
        num_indices, mpi_interface_index_type, &mpi_status);

    int mpi_type_size = 0;
    MPI_Type_size(mpi_interface_index_type, &mpi_type_size);
    MPI_Type_free(&mpi_interface_index_type);
    offset += num_indices * mpi_type_size;
    const int sfbitset_size = sizeof(DefSFBitset), int_size = sizeof(DefInt);
    int node_size = sfbitset_size + int_size;
    for (int i_rank = 0; i_rank < rank_id_; ++i_rank) {
        if (offset > (std::numeric_limits<MPI_Offset>::max)()
            - node_size* static_cast<MPI_Offset>(num_of_nodes.at(i_rank))) {
            LogManager::LogError("offset exceeds the maximum of MPI_Offset in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        offset += node_size * static_cast<MPI_Offset>(num_of_nodes.at(i_rank));
    }
    const std::size_t max_chunk_size = static_cast<std::size_t>(INT_MAX);
    std::size_t total_size = num_of_nodes.at(rank_id_) * node_size;

    int max_chunks = 0;
    for (int i = 0; i < num_of_ranks_; ++i) {
        std::size_t size_i = static_cast<std::size_t>(num_of_nodes.at(i)) * node_size;
        int chunks_i = static_cast<int>((size_i + max_chunk_size - 1) / max_chunk_size);
        if (chunks_i > max_chunks) {
            max_chunks = chunks_i;
        }
    }
    MPI_Offset current_offset = offset;
    std::vector<char> chunk_buffer;
    chunk_buffer.resize((std::min)(max_chunk_size, total_size));
    std::size_t remaining = total_size;
    int inner_fine2coarse = 1, outer_fine2coarse = 2, inner_coarse2fine = 3, outer_coarse2fine = 4;

    for (int chunk_idx = 0; chunk_idx < max_chunks; ++chunk_idx) {
        std::size_t this_chunk_size = 0;
        void* buffer_ptr = nullptr;
        if (remaining > 0) {
            this_chunk_size = (std::min)(max_chunk_size, remaining);
            buffer_ptr = chunk_buffer.data();
        }

        int mpi_chunk_size = static_cast<int>(this_chunk_size);

        // All ranks call collectively, even if no data to read
        MPI_File_read_at_all(mpi_file, current_offset, buffer_ptr, mpi_chunk_size, MPI_BYTE, &mpi_status);

        if (mpi_chunk_size > 0) {
            std::size_t index = 0;
            while (index < this_chunk_size) {
                DefSFBitset sfbitset;
                std::memcpy(&sfbitset, chunk_buffer.data() + index, sfbitset_size);
                index += sfbitset_size;
                DefInt index_read;
                std::memcpy(&index_read, chunk_buffer.data() + index, int_size);
                index += int_size;
                const DefInt i_level = vec_interface_index.at(index_read).i_level;
                const ECriterionType enum_criterion = vec_interface_index.at(index_read).enum_criterion;
                const DefInt criterion_count = vec_interface_index.at(index_read).criterion_count;
                const int layer_indicator = vec_interface_index.at(index_read).layer_indicator;
                const int layer_count = vec_interface_index.at(index_read).layer_count;
                if (static_cast<DefInt>(ptr_vec_grid_info->size()) <= i_level) {
                    LogManager::LogError("Grid level exceeds available levels in " +
                        std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                } else {
                    auto& grid_info = *ptr_vec_grid_info->at(i_level).get();
                    const std::pair<ECriterionType, DefInt> pair_interface(enum_criterion, criterion_count);
                    if (grid_info.map_ptr_interface_layer_info_.find(pair_interface)
                        == grid_info.map_ptr_interface_layer_info_.end()) {
                        grid_info.map_ptr_interface_layer_info_.insert(
                            std::make_pair(pair_interface, std::make_shared<InterfaceLayerInfo>()));
                    }
                    auto& ptr_interface = grid_info.map_ptr_interface_layer_info_.at(pair_interface);
                    if (layer_indicator == inner_fine2coarse) {
                        if (static_cast<int>(ptr_interface->vec_inner_fine2coarse_.size()) <= layer_count) {
                            ptr_interface->vec_inner_fine2coarse_.resize(layer_count + 1);
                        }
                        ptr_interface->vec_inner_fine2coarse_.at(layer_count).insert({sfbitset, 0});
                    }
                    if (layer_indicator == outer_fine2coarse) {
                        if (static_cast<int>(ptr_interface->vec_outer_fine2coarse_.size()) <= layer_count) {
                            ptr_interface->vec_outer_fine2coarse_.resize(layer_count + 1);
                        }
                        ptr_interface->vec_outer_fine2coarse_.at(layer_count).insert({sfbitset, 0});
                    }
                    if (layer_indicator == inner_coarse2fine) {
                        if (static_cast<int>(ptr_interface->vec_inner_coarse2fine_.size()) <= layer_count) {
                            ptr_interface->vec_inner_coarse2fine_.resize(layer_count + 1);
                        }
                        ptr_interface->vec_inner_coarse2fine_.at(layer_count).insert({sfbitset, 0});
                    }
                    if (layer_indicator == outer_coarse2fine) {
                        if (static_cast<int>(ptr_interface->vec_outer_coarse2fine_.size()) <= layer_count) {
                            ptr_interface->vec_outer_coarse2fine_.resize(layer_count + 1);
                        }
                        ptr_interface->vec_outer_coarse2fine_.at(layer_count).insert({sfbitset, 0});
                    }
                }
            }
            current_offset += mpi_chunk_size;
            remaining -= mpi_chunk_size;
        }
    }
    MPI_File_close(&mpi_file);
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
