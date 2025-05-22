//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_checkpoint.cpp
* @author Zhengliang Liu
* @brief functions used to write and write checkpoint.
*/
#include <stddef.h>
#include <utility>
#include <numeric>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
#include "io/debug_write.h"
// struct to store the index of grid interfaces for mpi communication
template <>
struct std::hash<rootproject::amrproject::InterfaceIndexForMpi> {
    std::size_t operator()(const rootproject::amrproject::InterfaceIndexForMpi& k) const {
        return ((std::hash<rootproject::DefInt>()(k.i_level) ^
                (std::hash<int>()(static_cast<int>(k.enum_criterion)) << 1)) >> 1) ^
               (std::hash<rootproject::DefInt>()(k.criterion_count) << 1) ^
               (std::hash<int>()(k.layer_indicator) << 2) ^
               (std::hash<int>()(k.layer_count) << 3);
    }
};
namespace rootproject {
namespace amrproject {
/** 
 * @brief function to compute computational load on each rank.
 * @param[in] max_level maximum level of all grids.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[in] vec_grid_info grid information at all levels.
 * @param[out] ptr_node_level pointer to nodes at which levels on current rank.
 * @param[out] ptr_num_of_nodes pointer to number of nodes at each level on current rank.
 * @param[out] ptr_num_of_interface_nodes pointer to number of interface nodes at each level on current rank.
 * @return total computational load on current rank.
 */
DefAmrLUint MpiManager::ComputeComputationalLoadOnEachRank(
    const DefInt max_level, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<std::shared_ptr<GridInfoInterface>>& vec_grid_info,
    std::map<DefSFCodeToUint, BackgroundLoadData>* const ptr_node_level,
    std::vector<DefAmrLUint>* const ptr_num_of_nodes,
    std::vector<DefAmrLUint>* const ptr_num_of_interface_nodes) const {
    DefAmrLUint load_sum = 0;
    ptr_node_level->clear();
    ptr_num_of_nodes->resize(max_level + 1, 0);
    ptr_num_of_interface_nodes->resize(max_level + 1, 0);
    if (max_level >= vec_grid_info.size()) {
        LogManager::LogError("max_level is smaller than the number of grids in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefSFCodeToUint code_background;
    BackgroundLoadData background_load_data(max_level + 1);

    for (const auto iter_grid : vec_grid_info) {
        DefInt node_cost = iter_grid->GetComputationalCost();
        DefInt i_level = iter_grid->GetGridLevel();
        for (const auto& iter_node : iter_grid->map_grid_node_) {
            if (!(iter_node.second->flag_status_ & NodeBitStatus::kNodeStatusMpiPartitionOuter_)) {
                if (load_sum + node_cost > (std::numeric_limits<DefAmrLUint>::max)()) {
                    LogManager::LogError("computational load exceeds the maximum of DefAmrLUint in "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
                load_sum += node_cost;
                ptr_num_of_nodes->at(i_level) += 1;
                code_background = sfbitset_aux.SFBitsetToSFCode(
                    sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first));
                if (ptr_node_level->find(code_background) == ptr_node_level->end()) {
                    ptr_node_level->insert({code_background, background_load_data});
                }
                ptr_node_level->at(code_background).level_bitset_.Set(i_level, true);
                ptr_node_level->at(code_background).num_of_grid_nodes_ += 1;
                ptr_node_level->at(code_background).total_load_ += node_cost;
            }
        }
        auto func_add_node = [i_level, &sfbitset_aux,
            &ptr_node_level, &iter_grid, &ptr_num_of_interface_nodes](
            const std::vector<rootproject::DefMap<rootproject::DefInt>>& vec_interface) {
            DefSFCodeToUint code_background;
            for (const auto& iter_layer : vec_interface) {
                for (const auto& iter_node : iter_layer) {
                    if (iter_grid->map_grid_node_.find(iter_node.first)
                        != iter_grid->map_grid_node_.end()) {
                        if (!(iter_grid->map_grid_node_.at(iter_node.first)->flag_status_
                            & NodeBitStatus::kNodeStatusMpiPartitionOuter_)) {
                            code_background = sfbitset_aux.SFBitsetToSFCode(
                                sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first));
                            ptr_node_level->at(code_background).num_of_interface_nodes_ += 1;
                            ptr_num_of_interface_nodes->at(i_level) += 1;
                        }
                    } else {
                        LogManager::LogError("grid node not exist but interface node exist in "
                            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                    }
                }
            }
        };
        for (const auto& iter_interface : iter_grid->map_ptr_interface_layer_info_) {
            func_add_node(iter_interface.second->vec_inner_coarse2fine_);
            func_add_node(iter_interface.second->vec_inner_fine2coarse_);
            func_add_node(iter_interface.second->vec_outer_coarse2fine_);
            func_add_node(iter_interface.second->vec_outer_fine2coarse_);
        }
    }
    return load_sum;
}
/**
 * @brief function to write computational load to file.
 * @param[in] max_level maximum level of all grids.
 * @param[in] load_sum total computational load on current rank.
 * @param[in] file_name file name to store the file.
 * @param[in] map_node_level nodes at which levels on current rank.
 */
void MpiManager::WriteCheckPointNodesAtWhichLevels(const DefInt max_level,
    const DefAmrLUint load_sum, const std::string& file_name,
    const std::map<DefSFCodeToUint, BackgroundLoadData>& map_node_level) const {
    DefAmrLUint load_all;
    std::vector<DefAmrLUint> vec_load_size(num_of_ranks_, 0),
    vec_load_size_current(num_of_ranks_, static_cast<DefAmrLUint>(map_node_level.size()));
    MPI_Reduce(&load_sum, &load_all, 1, GetMpiAmrLUintType(), MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Alltoall(vec_load_size_current.data(), 1, GetMpiAmrLUintType(), vec_load_size.data(), 1,
        GetMpiAmrLUintType(), MPI_COMM_WORLD);

    MPI_File mpi_file;
    MPI_Status mpi_status;

    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
    DefAmrLUint total_nodes = std::accumulate(vec_load_size.begin(), vec_load_size.end(), 0);
    MPI_Offset offset = sizeof(DefAmrLUint);
    if (rank_id_ == 0) {
        MPI_File_write_at_all(mpi_file, 0, &load_all, 1, GetMpiAmrLUintType(), &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, &total_nodes, 1, GetMpiAmrLUintType(), &mpi_status);
        offset += sizeof(DefAmrLUint);
        MPI_File_write_at_all(mpi_file, offset, &max_level, 1, GetMpiIntType(), &mpi_status);
    } else {
        MPI_File_write_at_all(mpi_file, 0, nullptr, 0, GetMpiAmrLUintType(), &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, nullptr, 0, GetMpiAmrLUintType(), &mpi_status);
        offset += sizeof(DefAmrLUint);
        MPI_File_write_at_all(mpi_file, offset, nullptr, 0, GetMpiIntType(), &mpi_status);
    }
    offset += sizeof(DefInt);
    const int num_size = sizeof(int32_t);
    const int load_size = sizeof(DefInt);
    int level_size = 0, data_size = 0;
    if (!map_node_level.empty()) {
        level_size = static_cast<int>(map_node_level.begin()->second.level_bitset_.GetNumBytes());
        data_size += level_size + 2 *num_size + load_size;
    }
    const int node_size = sizeof(DefSFCodeToUint) + data_size;

    for (int i_rank = 0; i_rank < rank_id_; ++i_rank) {
        if (offset > (std::numeric_limits<MPI_Offset>::max)()
            - node_size* static_cast<MPI_Offset>(vec_load_size.at(i_rank))) {
            LogManager::LogError("offset exceeds the maximum of MPI_Offset in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        offset += node_size * vec_load_size.at(i_rank);
    }

    if (map_node_level.size() * node_size > (std::numeric_limits<int>::max)()) {
        LogManager::LogError("buffer size exceeds the maximum of int in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    std::vector<char> buffer(map_node_level.size() * node_size);
    size_t index = 0;
    const int key_size = sizeof(DefSFCodeToUint);
    for (const auto& iter : map_node_level) {
        std::memcpy(buffer.data() + index, &iter.first, key_size);
        index += key_size;
        std::memcpy(buffer.data() + index, iter.second.level_bitset_.Data(), level_size);
        index += level_size;
        std::memcpy(buffer.data() + index, &iter.second.num_of_grid_nodes_, num_size);
        index += num_size;
        std::memcpy(buffer.data() + index, &iter.second.num_of_interface_nodes_, num_size);
        index += num_size;
        std::memcpy(buffer.data() + index, &iter.second.total_load_, load_size);
        index += load_size;
    }

    MPI_File_write_at_all(mpi_file, offset, buffer.data(),
        static_cast<int>(map_node_level.size()) * node_size, MPI_BYTE, &mpi_status);
    MPI_File_close(&mpi_file);
}
/**
 * @brief function to read computational load from file.
 * @param[in] file_name file name to store the file.
 * @param[in] ptr_map_node_level pointer to nodes at which levels on current rank.
 * @return total computational load for all ranks.
 */
DefAmrLUint MpiManager::ReadCheckPointNodesAtWhichLevels(
    const std::string& file_name, std::map<DefSFCodeToUint, BackgroundLoadData>* const ptr_map_node_level) const {
    ptr_map_node_level->clear();
    MPI_File mpi_file;
    MPI_Status mpi_status;
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);

    MPI_Offset offset = sizeof(DefAmrLUint);
    DefAmrLUint total_nodes = 0, load_read = 0;
    DefInt max_level = 0;
    MPI_File_read_at_all(mpi_file, 0, &load_read, 1, GetMpiAmrLUintType(), &mpi_status);
    MPI_File_read_at_all(mpi_file, offset, &total_nodes, 1, GetMpiAmrLUintType(), &mpi_status);
    offset += sizeof(DefAmrLUint);
    MPI_File_read_at_all(mpi_file, offset, &max_level, 1, GetMpiIntType(), &mpi_status);
    offset += sizeof(DefInt);

    const int num_size = sizeof(int32_t);
    const int load_size = sizeof(DefInt);
    const int code_size = sizeof(DefSFCodeToUint);
    DefSFCodeToUint code_read;
    BackgroundLoadData value(max_level + 1);
    int level_size = 0, data_size = 0;
    level_size = static_cast<int>(value.level_bitset_.GetNumBytes());
    data_size += level_size + 2 *num_size + load_size;

    const std::size_t max_chunk_size = static_cast<std::size_t>((std::numeric_limits<int>::max)());
    std::size_t total_size = total_nodes * (code_size + data_size);

    MPI_Offset current_offset = offset;
    std::vector<char> chunk_buffer;
    chunk_buffer.resize((std::min)(max_chunk_size, total_size));
    std::size_t remaining = total_size;
    std::size_t size_i = static_cast<std::size_t>(total_nodes) * (code_size + data_size);
    int max_chunks = static_cast<int>((size_i + max_chunk_size - 1) / max_chunk_size);

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
                std::memcpy(&code_read, chunk_buffer.data() + index, code_size);
                index += code_size;
                std::memcpy(value.level_bitset_.Data(), chunk_buffer.data() + index, level_size);
                index += level_size;
                std::memcpy(&value.num_of_grid_nodes_, chunk_buffer.data() + index, num_size);
                index += num_size;
                std::memcpy(&value.num_of_interface_nodes_, chunk_buffer.data() + index, num_size);
                index += num_size;
                std::memcpy(&value.total_load_, chunk_buffer.data() + index, load_size);
                index += load_size;

                ptr_map_node_level->insert({code_read, value});
            }

            current_offset += mpi_chunk_size;
            remaining -= mpi_chunk_size;
        }
    }

    MPI_File_close(&mpi_file);

    return load_read;
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

    int32_t data_size = 0;
    for (const auto& iter_grid : vec_grid_info) {
        if (data_size < iter_grid->GetSizeOfGridNodeInfoForCheckPoint()) {
            data_size = iter_grid->GetSizeOfGridNodeInfoForCheckPoint();
        }
    }
    if (rank_id_ == 0) {
        DefAmrLUint total_nodes = std::accumulate(
            num_of_nodes_each_rank.begin(), num_of_nodes_each_rank.end(), 0);
        MPI_File_write_at_all(mpi_file, 0, &total_nodes, 1, GetMpiAmrLUintType(), &mpi_status);
        MPI_File_write_at_all(mpi_file, sizeof(DefAmrLUint), &data_size, 1, MPI_INT32_T, &mpi_status);
    } else {
        MPI_File_write_at_all(mpi_file, 0, nullptr, 0, GetMpiAmrLUintType(), &mpi_status);
        MPI_File_write_at_all(mpi_file, sizeof(DefAmrLUint), nullptr, 0, MPI_INT, &mpi_status);
    }

    MPI_Offset offset = sizeof(DefAmrLUint) + sizeof(int32_t);
    const int32_t sfbitset_size = sizeof(DefSFBitset), int_size = sizeof(DefInt);
    const int32_t node_size = sfbitset_size + int_size + data_size;
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
 * @param[in] mpi_interface_background nodes on mpi interface at background level.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[out] ptr_mpi_interface_each_level pointer to nodes on mpi interface at each level.
 * @param[in, out] ptr_grid_manager pointer to class manage grids.
 */
void MpiManager::ReadCheckPointGridNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_nodes,
    const DefMap<DefInt>& mpi_interface_background,
    const SFBitsetAuxInterface& sfbitset_aux,
    std::vector<DefMap<DefInt>>* const ptr_mpi_interface_each_level,
    GridManagerInterface* const ptr_grid_manager) const {

    DefInt dims = ptr_grid_manager->k0GridDims_;
    std::vector<bool> periodic_min(dims, false), periodic_max(dims, false);
    ptr_grid_manager->vec_ptr_grid_info_.at(0)->CheckIfPeriodicDomainRequired(dims, &periodic_min, &periodic_max);
    std::vector<DefAmrLUint> indices_min = ptr_grid_manager->GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = ptr_grid_manager->GetMaxIndexOfBackgroundNodeArrAsVec();
    DefInt max_level = ptr_grid_manager->GetMaxLevel();
    ptr_mpi_interface_each_level->resize(max_level + 1);
    std::vector<std::vector<DefSFBitset>> domain_min_level(max_level + 1), domain_max_level(max_level + 1);
    std::vector<DefSFBitset> sfbitset_neighbors;
    DefSFCodeToUint code_tmp, code_min = sfbitset_aux.SFBitsetToSFCode(sfbitset_min_current_rank_),
        code_max = sfbitset_aux.SFBitsetToSFCode(sfbitset_max_current_rank_);
    DefInt remove_mpi_flag =  ~(NodeBitStatus::kNodeStatusMpiPartitionOuter_|
        NodeBitStatus::kNodeStatusMpiPartitionInner_|NodeBitStatus::kNodeStatusMpiInterpInner_);

    MPI_File mpi_file;
    MPI_Status mpi_status;
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
    DefAmrLUint total_node = 0;
    MPI_File_read_at_all(mpi_file, 0, &total_node, 1, GetMpiAmrLUintType(), &mpi_status);
    int32_t data_size = 0;
    MPI_File_read_at_all(mpi_file, sizeof(DefAmrLUint), &data_size, 1, MPI_INT32_T, &mpi_status);

    MPI_Offset offset = sizeof(DefAmrLUint) + sizeof(int32_t);
    for (const auto& iter_grid : ptr_grid_manager->vec_ptr_grid_info_) {
        DefInt i_level = iter_grid->GetGridLevel();
        domain_min_level.at(i_level).resize(dims);
        domain_max_level.at(i_level).resize(dims);
        sfbitset_aux.GetMinAtGivenLevel(i_level, indices_min, &domain_min_level.at(i_level));
        sfbitset_aux.GetMaxAtGivenLevel(i_level, indices_max, &domain_max_level.at(i_level));
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
    const std::size_t max_chunk_size = static_cast<std::size_t>((std::numeric_limits<int>::max)());
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
    DefAmrLUint node_count = 0;
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
                ++node_count;
                DefSFBitset sfbitset;
                std::memcpy(&sfbitset, chunk_buffer.data() + index, sfbitset_size);
                index += sfbitset_size;
                DefInt i_level_read;
                std::memcpy(&i_level_read, chunk_buffer.data() + index, int_size);
                index += int_size;
                if (ptr_grid_manager->vec_ptr_grid_info_.size() <= i_level_read) {
                    LogManager::LogError("Grid level exceeds available levels in " +
                        std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                } else {
                    auto& grid_info = ptr_grid_manager->vec_ptr_grid_info_.at(i_level_read);
                    DefSFBitset sfbitset_background = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level_read, sfbitset);
                    if (mpi_interface_background.find(sfbitset_background) != mpi_interface_background.end()) {
                        sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(sfbitset, periodic_min, periodic_max,
                            domain_min_level.at(i_level_read), domain_max_level.at(i_level_read), &sfbitset_neighbors);
                        for (const auto& iter : sfbitset_neighbors) {
                            code_tmp = sfbitset_aux.SFBitsetToSFCode(
                                sfbitset_aux.SFBitsetToNLowerLevelVir(i_level_read, iter));
                            if (code_tmp < code_min || code_tmp > code_max) {
                                ptr_mpi_interface_each_level->at(i_level_read).insert({sfbitset, 0});
                                break;
                            }
                        }
                    }

                    if (grid_info->GetPtrToParentGridManager()->InstantiateGridNode(sfbitset, grid_info.get())) {
                        auto& node = grid_info->map_grid_node_.at(sfbitset);
                        node->ReadANodeFromBufferForCheckpoint(chunk_buffer.data() + index);
                        node->flag_status_ &= remove_mpi_flag;
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

    if (num_of_nodes.at(rank_id_) != node_count) {
        LogManager::LogError("Failed to read all nodes from checkpoint in " +
            std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    MPI_File_close(&mpi_file);
}
/**
 * @brief function to compute minimum and maximum space filling codes for each rank.
 * @param[in] load_all total computational load for all ranks.
 * @param[in] vec_cost computational cost at each level.
 * @param[in] map_node_level nodes at which levels on current rank.
 * @param[out] ptr_num_of_grid_nodes pointer to number of grid nodes at each level on current rank.
 * @param[out] ptr_num_of_interface_nodes pointer to number of grid interface nodes at each level on current rank.
 */
void MpiManager::ComputeMinNMaxSFbitsetForEachRank(const DefAmrLUint load_all, const std::vector<DefInt>& vec_cost,
    const std::map<DefSFCodeToUint, BackgroundLoadData>& map_node_level,
    std::vector<DefAmrLUint>* const ptr_num_of_grid_nodes,
    std::vector<DefAmrLUint>* const ptr_num_of_interface_nodes) {
    const int num_ranks = num_of_ranks_;
    DefAmrLUint ave_load = static_cast<DefAmrLUint>(load_all / num_ranks) + 1;
    DefAmrLUint load_rank0 = load_all - (num_ranks - 1) * ave_load;
    std::vector<DefAmrLUint> rank_load(num_ranks, ave_load);
    rank_load.at(0) = load_rank0;   // load at rank 0, assuming lower than other ranks
    // traverse background nodes
    vec_sfcode_min_all_ranks_.resize(num_ranks);
    vec_sfcode_max_all_ranks_.resize(num_ranks);
    DefAmrLUint load_count = 0, load_count_all = 0;
    int i_rank = 0;
    vec_sfcode_min_all_ranks_.at(i_rank) = 0;
    vec_sfcode_max_all_ranks_.back() = std::prev(map_node_level.end())->first;
    DefInt node_cost = 0;
    DefInt max_level = static_cast<DefInt>(vec_cost.size() - 1);
    ptr_num_of_grid_nodes->clear();
    ptr_num_of_grid_nodes->resize(num_of_ranks_, 0);
    ptr_num_of_interface_nodes->clear();
    ptr_num_of_interface_nodes->resize(num_of_ranks_, 0);
    for (const auto& iter : map_node_level) {
        if (load_count >= rank_load.at(i_rank)) {
            load_count_all += load_count;
            load_count = 0;
            ++i_rank;
            vec_sfcode_min_all_ranks_.at(i_rank) = iter.first;
        }
        node_cost = iter.second.total_load_;
        ptr_num_of_grid_nodes->at(i_rank) += iter.second.num_of_grid_nodes_;
        ptr_num_of_interface_nodes->at(i_rank) += iter.second.num_of_interface_nodes_;
        if (load_count + node_cost >= rank_load.at(i_rank)) {
            vec_sfcode_max_all_ranks_.at(i_rank) = iter.first;
        }
        load_count += node_cost;
    }
    sfbitset_min_current_rank_ = vec_sfcode_min_all_ranks_.at(rank_id_);
    sfbitset_max_current_rank_ = vec_sfcode_max_all_ranks_.at(rank_id_);
    load_count_all += load_count;
    if (load_count_all != load_all) {
        LogManager::LogError("computed load (" + std::to_string(load_count_all)
            + ") is not equal to the input load(" + std::to_string(load_all) + ") in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
 * @brief function to search for outer MPI nodes based on the interface.
 * @param[in] mpi_interface_nodes_each_level nodes on mpi interface at each level.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[in] grid_manager class manage grids.
 * @param[out] ptr_outer_mpi_nodes pointer to nodes on mpi outer layers for each rank.
 */
void MpiManager::SearchForMpiOuterLayerBasedOnInterface(
    const std::vector<DefMap<DefInt>>& mpi_interface_nodes_each_level,
    const SFBitsetAuxInterface& sfbitset_aux, const GridManagerInterface& grid_manager,
    std::vector<std::map<int, DefMap<DefInt>>>* const ptr_outer_mpi_nodes) const {
    DefInt dims = grid_manager.k0GridDims_;
    std::vector<DefAmrLUint> indices_min = grid_manager.GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = grid_manager.GetMaxIndexOfBackgroundNodeArrAsVec();
    std::vector<DefSFBitset> domain_min_n_level(dims), domain_max_n_level(dims);
    std::vector<DefInt> search_length(dims, k0NumPartitionOuterLayers_);
    std::vector<DefSFBitset> sfbitset_region;
    DefSFCodeToUint code_min_background_level = sfbitset_aux.SFBitsetToSFCode(sfbitset_min_current_rank_),
        code_max_background_level = sfbitset_aux.SFBitsetToSFCode(sfbitset_max_current_rank_);

    int current_rank = rank_id_;
    const std::vector<DefSFCodeToUint>& ull_min = vec_sfcode_min_all_ranks_,
        ull_max = vec_sfcode_max_all_ranks_;
    DefInt max_level = grid_manager.GetMaxLevel();
    if (static_cast<DefInt>(mpi_interface_nodes_each_level.size()) != max_level + 1) {
         LogManager::LogError("numer of interface layers is not equal to refinement levels in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    ptr_outer_mpi_nodes->resize(max_level + 1);
    for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
        std::vector<bool> periodic_min(dims, false), periodic_max(dims, false);
        grid_manager.vec_ptr_grid_info_.at(i_level)->CheckIfPeriodicDomainRequired(dims, &periodic_min, &periodic_max);
        sfbitset_aux.GetMinAtGivenLevel(i_level, indices_min, &domain_min_n_level);
        sfbitset_aux.GetMaxAtGivenLevel(i_level, indices_max, &domain_max_n_level);
        auto func_set_inner = [&sfbitset_aux, &ull_min, &ull_max,
            &mpi_interface_nodes_each_level, &periodic_min, &periodic_max,
            &domain_min_n_level, &domain_max_n_level, &search_length, &ptr_outer_mpi_nodes,
            current_rank, i_level, code_min_background_level, code_max_background_level](
            const DefSFBitset& sfbitset_in) {
            std::vector<DefSFBitset> vec_in_region;
            std::vector<std::pair<DefAmrLUint, DefSFBitset>> vec_overlap;
            DefSFCodeToUint code_tmp;
            sfbitset_aux.FindNodesInPeriodicRegionCenterOverlap(sfbitset_in, search_length, search_length,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region, &vec_overlap);
            if (!vec_overlap.empty()) {
                for (const auto& iter_overlap : vec_overlap) {
                    vec_in_region.push_back(iter_overlap.second);
                }
            }
            std::vector<DefSFCodeToUint>::const_iterator iter_index;
            int node_rank;
            for (const auto& iter_node : vec_in_region) {
                if (iter_node ==  SFBitsetAuxInterface::kInvalidSFbitset) {
                    continue;
                }
                code_tmp = sfbitset_aux.SFBitsetToSFCode(
                    sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node));
                if ((code_tmp < code_min_background_level || code_tmp > code_max_background_level)) {
                    iter_index = std::lower_bound(ull_max.cbegin(), ull_max.cend(), code_tmp);
                    node_rank = static_cast<int>(iter_index - ull_max.cbegin());
                    if (node_rank != current_rank) {
                        if (ptr_outer_mpi_nodes->at(i_level).find(node_rank)
                            == ptr_outer_mpi_nodes->at(i_level).end()) {
                            ptr_outer_mpi_nodes->at(i_level).insert({node_rank, DefMap<DefInt>()});
                        }
                        ptr_outer_mpi_nodes->at(i_level).at(node_rank).insert({iter_node, i_level});
                    }
                }
            }
        };

        for (const auto& iter_node : mpi_interface_nodes_each_level.at(i_level)) {
            func_set_inner(iter_node.first);
        }
    }
}
/**
 * @brief function to send and receive nodes on outer mpi layers for checkpointing.
 * @param[in] outer_mpi_nodes nodes on outer mpi layers for each rank at all levels.
 * @param[out] ptr_grid_manager pointer to class manage grids.
 * */
void MpiManager::SendNReceiveCheckPointMpiLayers(const std::vector<std::map<int, DefMap<DefInt>>>& outer_mpi_nodes,
    GridManagerInterface* const ptr_grid_manager) {
    int node_size = sizeof(DefInt);
    DefInt max_level = ptr_grid_manager->GetMaxLevel();
    if (max_level + 1 < static_cast<DefInt>(outer_mpi_nodes.size())) {
        LogManager::LogError("number of input outer mpi layers is greater than number of refinement levels in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
        std::vector<BufferSizeInfo> send_buffer_info, receive_buffer_info;
        SendNReceiveGridNodeBufferSize(node_size, 0,
            outer_mpi_nodes.at(i_level), &send_buffer_info, &receive_buffer_info);

        std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
        std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send, vec_ptr_buffer_receive(num_of_ranks_);
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
        for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
            if (send_buffer_info.at(i_rank).bool_exist_) {
                vec_vec_reqs_send.push_back({});
                if (outer_mpi_nodes.at(i_level).find(i_rank) != outer_mpi_nodes.at(i_level).end()) {
                    const int& num_chunks = send_buffer_info.at(i_rank).num_chunks_;
                    vec_vec_reqs_send.back().resize(num_chunks);
                    int size_tmp;
                    vec_ptr_buffer_send.emplace_back(SerializeNodeSFBitset(
                        outer_mpi_nodes.at(i_level).at(i_rank), &size_tmp));
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
        DefInt num_levels = static_cast<DefInt>(ptr_grid_manager->vec_ptr_grid_info_.size());
        mpi_communication_inner_layers_.resize(num_levels);
        for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
            if (receive_buffer_info.at(i_rank).bool_exist_) {
                MPI_Waitall(static_cast<int>(vec_vec_reqs_receive.at(i_rev).size()),
                    vec_vec_reqs_receive.at(i_rev).data(), MPI_STATUSES_IGNORE);
                char* ptr_buffer = vec_ptr_buffer_receive.at(i_rank).get();
                int key_size = sizeof(DefSFBitset);
                int num_nodes;
                int position = 0;
                std::memcpy(&num_nodes, ptr_buffer, sizeof(int));
                // deserialize data stored in buffer
                position += sizeof(int);
                DefSFBitset key_host;
                for (int i_node = 0; i_node < num_nodes; ++i_node) {
                    std::memcpy(&key_host, ptr_buffer + position, key_size);
                    position += key_size;
                    auto& grid_info = ptr_grid_manager->vec_ptr_grid_info_.at(i_level);
                    if (grid_info->map_grid_node_.find(key_host) != grid_info->map_grid_node_.end()) {
                        grid_info->map_grid_node_.at(key_host)->flag_status_ |=
                            NodeBitStatus::kNodeStatusMpiPartitionInner_;
                        if (mpi_communication_inner_layers_.at(i_level).find(i_rank)
                            == mpi_communication_inner_layers_.at(i_level).end()) {
                            mpi_communication_inner_layers_.at(i_level).insert({i_rank, {}});
                        }
                        mpi_communication_inner_layers_.at(i_level).at(i_rank).insert({key_host, 0});
                    }
                    // else { node requested by other ranks does not exist}
                }
                ++i_rev;
            }
        }
        int i_send = 0;
        for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
            if (send_buffer_info.at(i_rank).bool_exist_) {
                MPI_Waitall(static_cast<int>(vec_vec_reqs_send.at(i_send).size()),
                    vec_vec_reqs_send.at(i_send).data(), MPI_STATUSES_IGNORE);
                ++i_send;
            }
        }
    }
}
/**
 * @brief function to wait for mpi communication and instantiate grid node information from received buffer.
 * @param[in] send_buffer_info  buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] vec_ptr_buffer_receive  buffer storing received grid node information.
 * @param[in] func_read_a_node_from_buffer function to read grid node infomation from buffer.
 * @param[out] ptr_vec_vec_reqs_send  pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive  pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 */
void MpiManager::ReadNInstantiateCheckPointGridNodesInMpiLayers(
    const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
    const std::function<void(const char*,  GridNode* const)>& func_read_a_node_from_buffer,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    GridInfoInterface* const ptr_grid_info) {
    int i_rev = 0;
    DefInt i_level = ptr_grid_info->GetGridLevel();
    GridManagerInterface* ptr_grid_manager = ptr_grid_info->GetPtrToParentGridManager();
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

            char* ptr_buffer = vec_ptr_buffer_receive.at(i_rev).get();
            int key_size = sizeof(DefSFBitset);
            int node_info_size = ptr_grid_info->GetSizeOfGridNodeInfoForMpiCommunication();
            DefSizet num_nodes = buffer_size/(key_size + node_info_size);
            DefSizet position = 0;
            DefSFBitset key_code;
            for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
                std::memcpy(&key_code, ptr_buffer + position, key_size);
                position += key_size;
                if (ptr_grid_info->map_grid_node_.find(key_code) == ptr_grid_info->map_grid_node_.end()) {
                    ptr_grid_manager->InstantiateGridNode(key_code, ptr_grid_info);
                }
                GridNode* ptr_node = ptr_grid_info->map_grid_node_.at(key_code).get();
                func_read_a_node_from_buffer(ptr_buffer + position, ptr_node);
                ptr_node->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                mpi_communication_outer_layers_.at(i_level).insert({key_code, 0});
                position += node_info_size;
                if (position > buffer_size) {
                    LogManager::LogError("Buffer to store node information overflows, please check"
                        " size of node info for mpi communication in "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
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
 * @brief function to send and receive grid nodes on mpi layers.
 * @param[in] send_buffer_info  buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] vec_vec_reqs_send  vector storing mpi requests information of sending.
 * @param[in] vec_vec_reqs_receive  vector storing mpi requests information of receiving.
 * @param[out] vec_ptr_buffer_send  pointer to buffer storing sending grid node information.
 * @param[out] vec_ptr_buffer_receive  pointer to buffer storing receiving grid node information.
 * @param[out] ptr_grid_info pointer to class storting grid information.
 */
void MpiManager::CommunicateCheckPointMpiLayers(GridManagerInterface* const ptr_grid_manager) {
    DefInt max_level = ptr_grid_manager->GetMaxLevel();
    std::function<void(const char*,  amrproject::GridNode* const)> func_read_a_node_from_buffer =
        [](const char* const ptr_node_buffer, amrproject::GridNode* const ptr_node) {
        ptr_node->ReadANodeFromBufferForMpi(ptr_node_buffer);
    };
    mpi_communication_outer_layers_.resize(max_level + 1);
    // send and receive grid nodes on mpi layers
    for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
        auto& grid_info = ptr_grid_manager->vec_ptr_grid_info_.at(i_level);
        std::vector<amrproject::MpiManager::BufferSizeInfo> send_buffer_info, receive_buffer_info;
        std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
        std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send, vec_ptr_buffer_receive;
        SendAndReceiveGridNodesOnMpiLayers(&send_buffer_info, &receive_buffer_info,
            &vec_vec_reqs_send, &vec_vec_reqs_receive,
            &vec_ptr_buffer_send, &vec_ptr_buffer_receive, grid_info.get());
        ReadNInstantiateCheckPointGridNodesInMpiLayers(send_buffer_info,
            receive_buffer_info, vec_ptr_buffer_receive, func_read_a_node_from_buffer,
            &vec_vec_reqs_send, &vec_vec_reqs_receive, grid_info.get());
    }
}
/**
 * @brief function to send, receive and instantiate nodes on both grid interface and mpi layers.
 * @param[in] vec_interface_index interface indices.
 * @param[in] interface_to_send interface nodes need to be sent.
 * @param[out] ptr_grid_manager pointer to class manage grids.
 */
void MpiManager::CommunicateCheckPointInterfaceLayers(const std::vector<InterfaceIndexForMpi>& vec_interface_index,
    const std::map<int, std::pair<int, DefMap<std::vector<InterfaceIndexForMpi>>>>& interface_to_send,
    GridManagerInterface* const ptr_grid_manager) {
    std::vector<BufferSizeInfo> send_buffer_info(num_of_ranks_), receive_buffer_info(num_of_ranks_);
    // send buffer information of interface nodes
    const int key_size = sizeof(DefSFBitset);
    const int node_buffer_size = key_size + sizeof(DefInt);
    std::vector<MPI_Request> reqs_receive(num_of_ranks_);
    for (int i = 1; i < num_of_ranks_; ++i) {
        // send and receive number of chunks if needed
        int i_rank_send = (rank_id_ + i) % num_of_ranks_;
        int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
        bool bool_send_to_i_rank = (interface_to_send.find(i_rank_send)
            != interface_to_send.end());
        int num_node_all = 0, last_num_nodes = 0, other_num_nodes = 0;
        if (bool_send_to_i_rank) {
            num_node_all = interface_to_send.at(i_rank_send).first;
            int num_chunks = static_cast<int>(node_buffer_size * num_node_all
                /((std::numeric_limits<int>::max)() - sizeof(int)))  + 1;
            send_buffer_info.at(i_rank_send).bool_exist_ = true;
            send_buffer_info.at(i_rank_send).num_chunks_ = num_chunks;
            if (num_chunks < 2) {  // buffer is enough to store information of all nodes
                other_num_nodes = 0;
                last_num_nodes = num_node_all;
            } else {
                last_num_nodes = num_node_all%(num_chunks - 1);
                other_num_nodes = (num_node_all - last_num_nodes)/(num_chunks - 1);
            }
            send_buffer_info.at(i_rank_send).num_chunks_ = num_chunks;
            MPI_Send(&num_chunks, 1, MPI_INT, i_rank_send, 0, MPI_COMM_WORLD);
        } else {
            send_buffer_info.at(i_rank_send).bool_exist_ = false;
            send_buffer_info.at(i_rank_send).num_chunks_ = 0;
            MPI_Send(&send_buffer_info.at(i_rank_send).num_chunks_,
                1, MPI_INT, i_rank_send, 0, MPI_COMM_WORLD);
        }

        MPI_Recv(&receive_buffer_info.at(i_rank_receive).num_chunks_,
            1, MPI_INT, i_rank_receive, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receive_buffer_info.at(i_rank_receive).num_chunks_ == 0) {
            receive_buffer_info.at(i_rank_receive).bool_exist_ = false;
        } else {
            receive_buffer_info.at(i_rank_receive).bool_exist_ = true;
        }

        // send and receive buffer size of each chunk if needed
        if (bool_send_to_i_rank) {
            DefSizet buffer_size_send = node_buffer_size * other_num_nodes;
            if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                send_buffer_info.at(i_rank_send).array_buffer_size_.at(0) = static_cast<int>(buffer_size_send);
            } else {
                LogManager::LogError("size of the buffer is greater than the maximum value of int");
            }
            buffer_size_send = node_buffer_size * last_num_nodes;
            if (CheckBufferSizeNotExceedMax(buffer_size_send)) {
                send_buffer_info.at(i_rank_send).array_buffer_size_.at(1) = static_cast<int>(buffer_size_send);
            } else {
                LogManager::LogError("size of the buffer is greater than the maximum value of int");
            }
            MPI_Send(send_buffer_info.at(i_rank_send).array_buffer_size_.data(),
                2, MPI_INT, i_rank_send, 0, MPI_COMM_WORLD);
        }
        if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
            MPI_Recv(receive_buffer_info.at(i_rank_receive).array_buffer_size_.data(),
                2, MPI_INT, i_rank_receive, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_send, vec_ptr_buffer_receive(num_of_ranks_);
    SendNReceiveMpiInterfaceIndices(vec_interface_index,
        send_buffer_info, receive_buffer_info, interface_to_send,
        &vec_vec_reqs_send, &vec_vec_reqs_receive, &vec_ptr_buffer_send,
        &vec_ptr_buffer_receive);

    WaitAndReadMpiInterfaceIndicesFromBuffer(vec_interface_index,
        send_buffer_info, receive_buffer_info, vec_ptr_buffer_receive,
        &vec_vec_reqs_send, &vec_vec_reqs_receive, &ptr_grid_manager->vec_ptr_grid_info_);
}
/**
 * @brief function to send and receive interface indices.
 * @param[in] vec_interface_index interface indices.
 * @param[in] send_buffer_info buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] interface_to_send interface nodes need to be sent.
 * @param[out] ptr_vec_vec_reqs_send pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_vec_ptr_buffer_send pointer to buffer storing sending grid node information.
 * @param[out] ptr_vec_ptr_buffer_receive pointer to buffer storing receiving grid node information.
 */
void MpiManager::SendNReceiveMpiInterfaceIndices(
    const std::vector<InterfaceIndexForMpi>& vec_interface_index,
    const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::map<int, std::pair<int, DefMap<std::vector<InterfaceIndexForMpi>>>>& interface_to_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_send,
    std::vector<std::unique_ptr<char[]>>* const ptr_vec_ptr_buffer_receive) const {
    ptr_vec_vec_reqs_send->clear();
    ptr_vec_ptr_buffer_send->clear();
    ptr_vec_vec_reqs_receive->clear();
    ptr_vec_ptr_buffer_receive->clear();

    const int key_size = sizeof(DefSFBitset);
    const int node_buffer_size = key_size + sizeof(DefInt);

    for (int i = 1; i < num_of_ranks_; ++i) {
        int i_rank_receive = (rank_id_ - i + num_of_ranks_)% num_of_ranks_;
        if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
            ptr_vec_vec_reqs_receive->push_back({});
            ptr_vec_ptr_buffer_receive->emplace_back(ReceiveCharBufferFromOtherRanks(i_rank_receive,
                receive_buffer_info.at(i_rank_receive), &ptr_vec_vec_reqs_receive->back()));
        }
    }
    std::unordered_map<InterfaceIndexForMpi, int> index_map;
    DefInt num_indices = 0;
    for (const auto& iter : vec_interface_index) {
        index_map.insert({iter, num_indices});
        ++num_indices;
    }

    for (int i = 1; i < num_of_ranks_; ++i) {
        int i_rank_send = (rank_id_ + i) % num_of_ranks_;
        if (send_buffer_info.at(i_rank_send).bool_exist_) {
            ptr_vec_vec_reqs_send->push_back({});
            ptr_vec_vec_reqs_send->push_back({});
            ptr_vec_ptr_buffer_send->emplace_back(std::make_unique<char[]>(
                node_buffer_size*interface_to_send.at(i_rank_send).first));

            const int& num_chunks = send_buffer_info.at(i_rank_send).num_chunks_;
            ptr_vec_vec_reqs_send->back().resize(num_chunks);
            // copy interface information to buffer
            int position = 0;
            char* ptr_buffer_send = ptr_vec_ptr_buffer_send->back().get();
            for (const auto& iter : interface_to_send.at(i_rank_send).second) {
                for (const auto& iter_index : iter.second) {
                    if (index_map.find(iter_index) == index_map.end()) {
                        LogManager::LogError("interface index not found in the map in "
                            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                    } else {
                        std::memcpy(ptr_buffer_send + position, &(iter.first), key_size);
                        position += key_size;
                        std::memcpy(ptr_buffer_send + position, &(index_map.at(iter_index)), sizeof(DefInt));
                        position += sizeof(DefInt);
                    }
                }
            }

            // send grid node information chunk by chunk via non-blocking communication
            const int& buffer_size_send = send_buffer_info.at(i_rank_send).array_buffer_size_.at(0);
            for (int i_chunk = 0; i_chunk < num_chunks - 1; ++i_chunk) {
                MPI_Isend(ptr_buffer_send+i_chunk*buffer_size_send,
                    send_buffer_info.at(i_rank_send).array_buffer_size_.at(0),
                    MPI_BYTE, i_rank_send, i_chunk, MPI_COMM_WORLD, &ptr_vec_vec_reqs_send->back().at(i_chunk));
            }
            int i_chunk_last = num_chunks - 1;
            MPI_Isend(ptr_buffer_send+i_chunk_last*buffer_size_send,
                send_buffer_info.at(i_rank_send).array_buffer_size_.at(1), MPI_BYTE, i_rank_send, i_chunk_last,
                MPI_COMM_WORLD, &ptr_vec_vec_reqs_send->back().at(i_chunk_last));
        }
    }
}
/**
 * @brief function to instantiate grid interface based on index.
 * @param[in] index index of the interface.
 * @param[in] sfbitset_in input space filling code.
 * @param[in] vec_interface_index vector of interface indices.
 * @param[out] ptr_vec_grid_info pointer to vector storing grid information.
 */
void MpiManager::InstantiateGridInterfaceBaseOnIndex(const DefInt index,
    const DefSFBitset sfbitset_in, const std::vector<InterfaceIndexForMpi>& vec_interface_index,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) const {
    const DefInt& i_level = vec_interface_index.at(index).i_level;
    const ECriterionType& enum_criterion = vec_interface_index.at(index).enum_criterion;
    const DefInt& criterion_count = vec_interface_index.at(index).criterion_count;
    const int& layer_indicator = vec_interface_index.at(index).layer_indicator;
    const int& layer_count = vec_interface_index.at(index).layer_count;
    if (ptr_vec_grid_info->size() <= i_level) {
        LogManager::LogError("Grid level exceeds available levels in " +
            std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        auto& grid_info = *ptr_vec_grid_info->at(i_level).get();
        const std::pair<ECriterionType, DefInt> pair_interface(enum_criterion, criterion_count);
        if (grid_info.map_ptr_interface_layer_info_.find(pair_interface)
            == grid_info.map_ptr_interface_layer_info_.end()) {
            LogManager::LogError("Grid interface not found in " + std::string(__FILE__) + " at line ");
        }
        auto& ptr_interface = grid_info.map_ptr_interface_layer_info_.at(pair_interface);
        if (layer_indicator == kIndictorInnerFine2Coarse_) {
            if (ptr_interface->vec_inner_fine2coarse_.size() <= layer_count) {
                ptr_interface->vec_inner_fine2coarse_.resize(layer_count + 1);
            }
            ptr_interface->vec_inner_fine2coarse_.at(layer_count).insert({sfbitset_in, 0});
        }
        if (layer_indicator == kIndictorOuterFine2Coarse_) {
            if (ptr_interface->vec_outer_fine2coarse_.size() <= layer_count) {
                ptr_interface->vec_outer_fine2coarse_.resize(layer_count + 1);
            }
            ptr_interface->vec_outer_fine2coarse_.at(layer_count).insert({sfbitset_in, 0});
        }
        if (layer_indicator == kIndictorInnerCoarse2Fine_) {
            if (ptr_interface->vec_inner_coarse2fine_.size() <= layer_count) {
                ptr_interface->vec_inner_coarse2fine_.resize(layer_count + 1);
            }
            ptr_interface->vec_inner_coarse2fine_.at(layer_count).insert({sfbitset_in, 0});
        }
        if (layer_indicator == kIndictorOuterCoarse2Fine_) {
            if (ptr_interface->vec_outer_coarse2fine_.size() <= layer_count) {
                ptr_interface->vec_outer_coarse2fine_.resize(layer_count + 1);
            }
            ptr_interface->vec_outer_coarse2fine_.at(layer_count).insert({sfbitset_in, 0});
        }
    }
}
/**
 * @brief function to wait for mpi communication and instantiate grid interface based on index.
 * @param[in] vec_interface_index vector of interface indices.
 * @param[in] send_buffer_info buffer size information of sending.
 * @param[in] receive_buffer_info buffer size information of receiving.
 * @param[in] vec_ptr_buffer_receive buffer storing received grid node information.
 * @param[out] ptr_vec_vec_reqs_send pointer to vector storing mpi requests information of sending.
 * @param[out] ptr_vec_vec_reqs_receive pointer to vector storing mpi requests information of receiving.
 * @param[out] ptr_vec_grid_info pointer to vector storing grid information.
 */
void MpiManager::WaitAndReadMpiInterfaceIndicesFromBuffer(
    const std::vector<InterfaceIndexForMpi>& vec_interface_index,
    const std::vector<BufferSizeInfo>& send_buffer_info,
    const std::vector<BufferSizeInfo>& receive_buffer_info,
    const std::vector<std::unique_ptr<char[]>>& vec_ptr_buffer_receive,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_send,
    std::vector<std::vector<MPI_Request>>* const ptr_vec_vec_reqs_receive,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info) {
    const int key_size = sizeof(DefSFBitset);
    const int node_info_size = sizeof(DefInt);
    const int node_buffer_size = key_size + node_info_size;
    int i_rev = 0;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        int i_rank_receive = (rank_id_ - i_rank + num_of_ranks_)% num_of_ranks_;
        if (receive_buffer_info.at(i_rank_receive).bool_exist_) {
            MPI_Waitall(static_cast<int>(ptr_vec_vec_reqs_receive->at(i_rev).size()),
                ptr_vec_vec_reqs_receive->at(i_rev).data(), MPI_STATUSES_IGNORE);
            const DefSizet buffer_size = receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(1)
                + (receive_buffer_info.at(i_rank_receive).num_chunks_ - 1)
                *receive_buffer_info.at(i_rank_receive).array_buffer_size_.at(0);
            const DefSizet num_nodes = buffer_size/node_buffer_size;
            DefSizet position = 0;
            DefSFBitset key_code;
            const char* ptr_buffer = vec_ptr_buffer_receive.at(i_rev).get();
            for (DefSizet i_node = 0; i_node < num_nodes; ++i_node) {
                std::memcpy(&key_code, ptr_buffer + position, key_size);
                position += key_size;
                DefInt interface_index;
                std::memcpy(&interface_index, ptr_buffer + position, node_info_size);
                position += node_info_size;
                const DefInt i_level = vec_interface_index.at(interface_index).i_level;
                if (ptr_vec_grid_info->size() <= i_level) {
                    LogManager::LogError("Grid level exceeds available levels in " +
                        std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                } else {
                    auto& grid_info = *ptr_vec_grid_info->at(i_level).get();
                    InstantiateGridInterfaceBaseOnIndex(interface_index, key_code,
                        vec_interface_index, ptr_vec_grid_info);
                    const DefInt& num_f2c_ghost_layer = grid_info.GetNumFine2CoarseGhostLayer();
                    const DefInt& num_c2f_ghost_layer = grid_info.GetNumCoarse2FineGhostLayer();
                    const DefInt& num_f2c_max_layer = grid_info.GetNumFine2CoarseLayer();
                    const DefInt& num_c2f_max_layer = grid_info.GetNumCoarse2FineLayer();
                    if (grid_info.map_grid_node_.find(key_code) == grid_info.map_grid_node_.end()) {
                        LogManager::LogError("Grid node does not exist in " +
                            std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                    }
                    if (mpi_communication_outer_layers_.at(i_level).find(key_code)
                        == mpi_communication_outer_layers_.at(i_level).end()) {
                        mpi_communication_outer_layers_.at(i_level).insert({key_code, 0});
                    }
                    auto& ptr_node = grid_info.map_grid_node_.at(key_code);
                    ptr_node->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                    const int& indicator = vec_interface_index.at(interface_index).layer_indicator;
                    const int& ilayer = vec_interface_index.at(interface_index).layer_count;
                    if (indicator == kIndictorOuterFine2Coarse_ || indicator == kIndictorInnerFine2Coarse_) {
                        if (ilayer ==  num_f2c_max_layer - 1) {
                            ptr_node->flag_status_ |= NodeBitStatus::kNodeStatusFine2Coarse0_;
                        }
                        if (ilayer >= num_f2c_max_layer - num_f2c_ghost_layer) {
                            ptr_node->flag_status_ |=  NodeBitStatus::kNodeStatusFine2CoarseGhost_;
                        }
                    }
                    if (indicator == kIndictorOuterCoarse2Fine_ || indicator == kIndictorInnerCoarse2Fine_) {
                        if (ilayer == num_c2f_max_layer - 1) {
                            ptr_node->flag_status_ |= NodeBitStatus::kNodeStatusCoarse2Fine0_;
                        }
                        if (ilayer >= num_c2f_max_layer - num_c2f_ghost_layer) {
                            ptr_node->flag_status_ |= NodeBitStatus::kNodeStatusCoarse2FineGhost_;
                        }
                    }
                }

                if (position > buffer_size) {
                    LogManager::LogError("Buffer to store node information overflows, please check"
                        " size of node info for mpi communication in " + std::string(__FILE__) + " at line "
                        + std::to_string(__LINE__));
                }
            }
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
 * @brief function to write grid interface nodes to a checkpoint file.
 * @param[in] file_name file name to write the file.
 * @param[in] num_of_interface_nodes number of interface nodes at each level on current rank.
 * @param[in] map_node_level nodes at which levels on current rank.
 * @param[in] vec_grid_info grid information at all levels.
 */
void MpiManager::WriteCheckPointInterfaceNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_interface_nodes,
    const std::map<DefSFCodeToUint, BackgroundLoadData>& map_node_level,
    const std::vector<std::shared_ptr<GridInfoInterface>>& vec_grid_info) const {
    DefAmrLUint num_of_nodes_for_all_levels =
        std::accumulate(num_of_interface_nodes.begin(), num_of_interface_nodes.end(), 0);
    std::vector<DefAmrLUint> num_of_nodes_each_rank(num_of_ranks_, 0);
    MPI_Allgather(&num_of_nodes_for_all_levels, 1, GetMpiAmrLUintType(),
        num_of_nodes_each_rank.data(), 1, GetMpiAmrLUintType(), MPI_COMM_WORLD);

    DefInt max_level = vec_grid_info.back()->GetGridLevel();
    if (max_level != static_cast<DefInt>(num_of_interface_nodes.size() - 1)) {
        LogManager::LogError("max_level is not equal to the side of input (num_of_nodes) which is"
            " number of nodes in each level in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    MPI_File mpi_file;
    MPI_Status mpi_status;
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
    int32_t num_indices = 0;
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitInterfaceIndexType(&mpi_interface_index_type);
    std::unordered_map<InterfaceIndexForMpi, int32_t> index_map;
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

    MPI_Offset offset = sizeof(int32_t);
    if (rank_id_ == 0) {
        std::vector<InterfaceIndexForMpi> vec_interface_index;
        for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
            const auto& grid_info = *vec_grid_info.at(i_level).get();
            for (const auto& iter_interface : grid_info.map_ptr_interface_layer_info_) {
                const auto& info = *iter_interface.second;
                auto enum_criterion = iter_interface.first.first;
                auto criterion_count = iter_interface.first.second;

                add_entries(kIndictorInnerFine2Coarse_, i_level, enum_criterion,
                    criterion_count, info.vec_inner_fine2coarse_, &vec_interface_index);
                add_entries(kIndictorOuterFine2Coarse_, i_level, enum_criterion,
                    criterion_count, info.vec_outer_fine2coarse_, &vec_interface_index);
                add_entries(kIndictorInnerCoarse2Fine_ , i_level, enum_criterion,
                    criterion_count, info.vec_inner_coarse2fine_, &vec_interface_index);
                add_entries(kIndictorOuterCoarse2Fine_, i_level, enum_criterion,
                    criterion_count, info.vec_outer_coarse2fine_, &vec_interface_index);
            }
        }

        MPI_File_write_at_all(mpi_file, 0, &num_indices, 1, MPI_INT32_T, &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, vec_interface_index.data(),
            num_indices, mpi_interface_index_type, &mpi_status);
    } else {
        for (DefInt i_level = 0; i_level <= max_level; ++i_level) {
            const auto& grid_info = *vec_grid_info.at(i_level).get();
            for (const auto& iter_interface : grid_info.map_ptr_interface_layer_info_) {
                const auto& info = *iter_interface.second;
                auto enum_criterion = iter_interface.first.first;
                auto criterion_count = iter_interface.first.second;

                add_entries(kIndictorInnerFine2Coarse_, i_level, enum_criterion,
                    criterion_count, info.vec_inner_fine2coarse_, nullptr);
                add_entries(kIndictorOuterFine2Coarse_, i_level, enum_criterion,
                    criterion_count, info.vec_outer_fine2coarse_, nullptr);
                add_entries(kIndictorInnerCoarse2Fine_, i_level, enum_criterion,
                    criterion_count, info.vec_inner_coarse2fine_, nullptr);
                add_entries(kIndictorOuterCoarse2Fine_, i_level, enum_criterion,
                    criterion_count, info.vec_outer_coarse2fine_, nullptr);
            }
        }
        MPI_File_write_at_all(mpi_file, 0, nullptr, 0, MPI_INT32_T, &mpi_status);
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
    int max_chunks = 0;
    for (int i_rank = 0; i_rank < num_of_ranks_; ++i_rank) {
        std::size_t size_i = static_cast<std::size_t>(num_of_nodes_each_rank.at(i_rank)) * node_size;
        int chunks_i = static_cast<int>((size_i + max_chunk_size - 1) / max_chunk_size);
        if (chunks_i > max_chunks) {
            max_chunks = chunks_i;
        }
    }
    DefSizet chunk_index = 0;
    MPI_Offset current_offset = offset;
    std::vector<DefSFBitset> sfbitset_cell;
    int local_chunk_count = 0;
    auto add_node = [sfbitset_size, int_size, &mpi_file, &mpi_status, &current_offset,
        &index_map, &chunk_index, &local_chunk_count](const int layer_indicator,
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
                    ++local_chunk_count;
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
                        add_node(kIndictorInnerFine2Coarse_, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_inner_fine2coarse_, &chunk_buffer);
                        add_node(kIndictorOuterFine2Coarse_, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_outer_fine2coarse_, &chunk_buffer);
                        add_node(kIndictorInnerCoarse2Fine_, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_inner_coarse2fine_, &chunk_buffer);
                        add_node(kIndictorOuterCoarse2Fine_, i_level, enum_criterion,
                            criterion_count, iter_cell, info.vec_outer_coarse2fine_, &chunk_buffer);
                    }
                }
            }
        }
    }
    if (chunk_index > 0) {
        MPI_File_write_at_all(mpi_file, current_offset,
            chunk_buffer.data(), static_cast<int>(chunk_index), MPI_BYTE, &mpi_status);
        ++local_chunk_count;
    }
    for (; local_chunk_count < max_chunks; ++local_chunk_count) {
        MPI_File_write_at_all(mpi_file, current_offset, nullptr, 0, MPI_BYTE, &mpi_status);
    }
    MPI_File_close(&mpi_file);
}
/**
 * @brief function to read interface nodes from a checkpoint file.
 * @param[in] file_name file name to read the file.
 * @param[in] num_of_nodes number of interface nodes on each rank.
 * @param[out] ptr_vec_interface_index pointer to vector storing interface indices.
 * @param[out] ptr_vec_grid_info pointer to class storting grid information on current rank.
 * @param[out] ptr_interface_to_send pointer to interface nodes need to send to other ranks.
 */
void MpiManager::ReadCheckPointInterfaceNodes(const std::string& file_name,
    const std::vector<DefAmrLUint>& num_of_nodes,
    std::vector<InterfaceIndexForMpi>* const ptr_vec_interface_index,
    std::vector<std::shared_ptr<GridInfoInterface>>* const ptr_vec_grid_info,
    std::map<int, std::pair<int, DefMap<std::vector<InterfaceIndexForMpi>>>>* const ptr_interface_to_send) const {
    // musted by called after all nodes excluding those on mpi outer layers are instantiated
    MPI_File mpi_file;
    MPI_Status mpi_status;
    int32_t num_indices = 0;
    MPI_Datatype mpi_interface_index_type;
    CreateAndCommitInterfaceIndexType(&mpi_interface_index_type);

    int mpi_err = MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
    if (mpi_err != MPI_SUCCESS) {
        LogManager::LogError("Failed to open file for reading: " + file_name);
    }

    MPI_File_read_at_all(mpi_file, 0, &num_indices, 1, MPI_INT32_T, &mpi_status);

    ptr_vec_interface_index->resize(num_indices);
    MPI_Offset offset = sizeof(MPI_INT32_T);
    MPI_File_read_at_all(mpi_file, offset, ptr_vec_interface_index->data(),
        num_indices, mpi_interface_index_type, &mpi_status);
    for (const auto& iter : *ptr_vec_interface_index) {
        if (ptr_vec_grid_info->size() <=  iter.i_level) {
            LogManager::LogError("Grid level exceeds available levels in " +
                std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        } else {
            auto& grid_info = *ptr_vec_grid_info->at(iter.i_level).get();
            const std::pair<ECriterionType, DefInt> pair_interface(iter.enum_criterion,  iter.criterion_count);
            if (grid_info.map_ptr_interface_layer_info_.find(pair_interface)
                == grid_info.map_ptr_interface_layer_info_.end()) {
                grid_info.map_ptr_interface_layer_info_.insert(
                    std::make_pair(pair_interface, std::make_shared<InterfaceLayerInfo>()));
            }
            grid_info.map_ptr_interface_layer_info_.at(pair_interface)
                ->vec_inner_coarse2fine_.resize(grid_info.GetNumCoarse2FineLayer());
            grid_info.map_ptr_interface_layer_info_.at(pair_interface)
                ->vec_inner_fine2coarse_.resize(grid_info.GetNumFine2CoarseLayer());
            grid_info.map_ptr_interface_layer_info_.at(pair_interface)
                ->vec_outer_coarse2fine_.resize(grid_info.GetNumCoarse2FineLayer());
            grid_info.map_ptr_interface_layer_info_.at(pair_interface)
                ->vec_outer_fine2coarse_.resize(grid_info.GetNumFine2CoarseLayer());
        }
    }

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
    const std::size_t max_chunk_size = static_cast<std::size_t>((std::numeric_limits<int>::max)());
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
                const DefInt i_level = ptr_vec_interface_index->at(index_read).i_level;
                if (ptr_vec_grid_info->size() <= i_level) {
                    LogManager::LogError("Grid level exceeds available levels in " +
                        std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                } else {
                    auto& grid_info = *ptr_vec_grid_info->at(i_level).get();
                    if (grid_info.map_grid_node_.find(sfbitset) == grid_info.map_grid_node_.end()) {
                        LogManager::LogError("Grid node does not exist in " +
                            std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                    }
                    if (grid_info.map_grid_node_.at(sfbitset)->flag_status_&
                        NodeBitStatus::kNodeStatusMpiPartitionInner_) {
                        for (auto& iter_inner_layer : mpi_communication_inner_layers_.at(i_level)) {
                            if (iter_inner_layer.second.find(sfbitset) != iter_inner_layer.second.end()) {
                                if (ptr_interface_to_send->find(iter_inner_layer.first)
                                    == ptr_interface_to_send->end()) {
                                    ptr_interface_to_send->insert({iter_inner_layer.first, {0, {}}});
                                }
                                ptr_interface_to_send->at(iter_inner_layer.first).first++;
                                if (ptr_interface_to_send->at(iter_inner_layer.first).second.find(sfbitset)
                                    == ptr_interface_to_send->at(iter_inner_layer.first).second.end()) {
                                    ptr_interface_to_send->at(iter_inner_layer.first).second.insert(
                                        {sfbitset, {}});
                                }
                                ptr_interface_to_send->at(iter_inner_layer.first).second.at(sfbitset)
                                    .push_back(ptr_vec_interface_index->at(index_read));
                            }
                        }
                    }
                    InstantiateGridInterfaceBaseOnIndex(index_read, sfbitset,
                        *ptr_vec_interface_index, ptr_vec_grid_info);
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
