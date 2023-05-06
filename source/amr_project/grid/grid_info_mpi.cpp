//  Copyright (c) 2022, Zhengliang Liu
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
struct TrackingIndexForMpi {
    ECriterionType enum_criterion;
    uint64_t index_criterion;
    unsigned int index_creator;
};
/**
 * @brief function to create and commit MPI std::tuple<ECriterionType, unit64_t, unsigned int> types.
 * @param[out] ptr_mpi_tracking_index_type pointer to the std::tuple<ECriterionType, unit64_t, unsigned int> types.
 * @return if the MPI types were created and committed successfully.
 */
int GridInfoInterface::CreateAndCommitTrackingIndexType(MPI_Datatype* ptr_mpi_tracking_index_type) {
    MPI_Datatype mpi_tracking_enum_type;
    MPI_Type_contiguous(sizeof(ECriterionType), MPI_BYTE, &mpi_tracking_enum_type);
    MPI_Type_commit(&mpi_tracking_enum_type);
    int ss; MPI_Type_size(mpi_tracking_enum_type, &ss);
    std::tuple<ECriterionType, uint64_t, unsigned int> tuple_tmp;
    // the last MPI_UNSIGNED store the index (k0IndexCreator) of the creator for tracking grid nodes
    MPI_Datatype mpi_create_type[3] = {mpi_tracking_enum_type, MPI_UINT64_T, MPI_UNSIGNED};
    int mpi_block_length[3] = {1, 1, 1};
    MPI_Aint mpi_disp[3];
    mpi_disp[0] = offsetof(TrackingIndexForMpi, enum_criterion);
    mpi_disp[1] = offsetof(TrackingIndexForMpi, index_criterion);
    mpi_disp[2] = offsetof(TrackingIndexForMpi, index_creator);
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
int GridInfoInterface::SerializeTrackingNode(const std::set<DefSFCodeToUint>& set_nodes,
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

void GridInfoInterface::DeserializeTrackingNode(const ECriterionType& criterion_type,
    DefSizet index_criterion, const std::unique_ptr<char[]>& buffer) const {
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
        map_ptr_tracking_grid_info_.at({criterion_type, index_criterion})->map_tracking_node_.insert(
            {static_cast<DefSFBitset>(key_net),
            map_ptr_tracking_grid_info_.at({criterion_type, index_criterion})->k0TrackNodeInstance_});
    }
}
void GridInfoInterface::IniSendNReceiveTracking(const std::vector<DefSFBitset>& bitset_max,
     const MpiManager& mpi_manager, const SFBitsetAuxInterface& bistset_aux,
    const std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>>& vec_tracking_info_creator) {
    MPI_Datatype mpi_tracking_index_type;
    CreateAndCommitTrackingIndexType(&mpi_tracking_index_type);
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    std::vector<DefSFCodeToUint> ull_max(bitset_max.size());
    if (rank_id == 0) {
        for (auto i = 0; i < bitset_max.size(); ++i) {
            ull_max.at(i) = bitset_max.at(i).to_ullong();
        }
    }

    DefSizet num_max = ull_max.size();
    uint64_t num_tracking = static_cast<uint64_t>(map_ptr_tracking_grid_info_.size());
    MPI_Bcast(&num_tracking, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    TrackingIndexForMpi tracking_index_tmp;
    std::vector<std::pair<ECriterionType, DefSizet>> vec_tracking_indices;
    if (rank_id == 0) {
        // save tracking indices stored in unordered_map to vector
        for (const auto& iter_tracking : map_ptr_tracking_grid_info_) {
            vec_tracking_indices.push_back(iter_tracking.first);
        }
    }
    for (uint64_t i = 0; i < num_tracking; ++i) {
        // broad cast the current index of the tracking grid
        if (rank_id == 0) {
            tracking_index_tmp.enum_criterion = vec_tracking_indices.at(i).first;
            tracking_index_tmp.index_criterion = static_cast<uint64_t>(vec_tracking_indices.at(i).second);
            tracking_index_tmp.index_creator =
                map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i))->k0IndexCreator;
        }
        MPI_Bcast(&tracking_index_tmp, 1, mpi_tracking_index_type, 0, MPI_COMM_WORLD);

        if (rank_id == 0 && num_ranks > 1) {
            //   send tracking grid information
            DefSFCodeToUint background_code, key_ull;
            TrackingGridInfoInterface* ptr_tracking_info =
                map_ptr_tracking_grid_info_.at(vec_tracking_indices.at(i)).get();
            int max_buffer = (std::numeric_limits<int>::max)() / ptr_tracking_info->SizeOfEachNodeForMpi() - 1;
            std::vector<std::vector<std::set<DefSFCodeToUint>>> vec_nodes_ranks(num_ranks);
            vec_nodes_ranks.at(0).push_back({});
            int index;
            std::vector<int> i_chunk_each_rank(num_ranks, -1), i_counts(num_ranks, 0);
            for (const auto& iter_node : ptr_tracking_info->map_tracking_node_) {
                background_code = (bistset_aux.SFBitsetToNLowerLevelVir(i_level_, iter_node.first)).to_ullong();
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
                std::vector<std::unique_ptr<MPI_Request>> reqs_send(num_chunks);
                std::vector<std::unique_ptr<MPI_Status>> stats_send(num_chunks);
                for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                    buffer_size_send = SerializeTrackingNode(vec_nodes_ranks.at(iter_rank).at(i_chunk),
                    vec_ptr_buffer.at(i_chunk));
                    reqs_send[i_chunk].reset(new MPI_Request);
                    stats_send[i_chunk].reset(new MPI_Status);
                    MPI_Send(&buffer_size_send, 1, MPI_INT, iter_rank, i_chunk, MPI_COMM_WORLD);
                    MPI_Isend(vec_ptr_buffer.at(i_chunk).get(), buffer_size_send, MPI_BYTE, iter_rank,
                    i_chunk, MPI_COMM_WORLD, reqs_send[i_chunk].get());
                }
                MPI_Waitall(num_chunks, reinterpret_cast<MPI_Request*>(reqs_send.data()),
                reinterpret_cast<MPI_Status*>(stats_send.data()));
            }
            for (const auto& iter_key : vec_nodes_ranks.at(0).at(0)) {
                ptr_tracking_info->map_tracking_node_.erase(static_cast<DefSFBitset>(iter_key));
            }
        } else if (rank_id > 0) {
            std::pair<ECriterionType, DefSizet> pair_tracking = {
                tracking_index_tmp.enum_criterion, tracking_index_tmp.index_criterion};
            unsigned int index_creator = tracking_index_tmp.index_creator;
            if (map_ptr_tracking_grid_info_.find(pair_tracking) == map_ptr_tracking_grid_info_.end()) {
                map_ptr_tracking_grid_info_.insert({pair_tracking,
                    vec_tracking_info_creator.at(index_creator).get()->CreateTrackingGridInfo()});
                map_ptr_tracking_grid_info_.at(pair_tracking)->k0IndexCreator = index_creator;
            }
            // receive tracking grid information
            int num_chunks = 0;
            MPI_Recv(&num_chunks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk) {
                int buffer_size_receive;
                MPI_Recv(&buffer_size_receive, sizeof(int), MPI_INT, 0, i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::unique_ptr<char[]> vec_ptr_buffer = std::make_unique<char[]>(buffer_size_receive);
                MPI_Recv(vec_ptr_buffer.get(), buffer_size_receive, MPI_BYTE, 0,
                    i_chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                DeserializeTrackingNode(pair_tracking.first, pair_tracking.second, vec_ptr_buffer);
            }
        }
    }

    MPI_Type_free(&mpi_tracking_index_type);
}
void GridInfoInterface::IniSendNReceiveInterface(const std::vector<DefSFBitset>& ptr_bitset_max,
    const MpiManager& mpi_manager) {
    int rank_id = mpi_manager.rank_id_, num_ranks = mpi_manager.num_of_ranks_;
    if (rank_id == 0) {
        uint64_t uint64_num;
        for (auto iter_tracking : map_ptr_tracking_grid_info_) {
            // define custom datatype for enum class
            MPI_Datatype e_criterion_type;
            MPI_Type_contiguous(sizeof(ECriterionType), MPI_BYTE, &e_criterion_type);
            MPI_Type_commit(&e_criterion_type);

            uint64_num = static_cast<uint64_t>(iter_tracking.first.second);
            MPI_Datatype pair_type;    // define datatype for pair<ECriterionType, uint64_t>
            MPI_Datatype mpi_create_type[2] = {e_criterion_type, MPI_UINT64_T};
            int mpi_block_length[2] = {1, 1};
            MPI_Aint mpi_disp[2] = {0, sizeof(ECriterionType)};
            MPI_Type_create_struct(2, mpi_block_length, mpi_disp, mpi_create_type, &pair_type);
            MPI_Type_commit(&pair_type);
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
