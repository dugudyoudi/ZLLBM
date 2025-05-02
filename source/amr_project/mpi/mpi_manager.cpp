//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage mpi processes.
* @date  2023-4-16
* @note .
*/
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <limits>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
static MPI_Datatype MPI_REAL_DATA_TYPE;
static MPI_Datatype MPI_INT_DATA_TYPE;
static MPI_Datatype MPI_UINT_DATA_TYPE;
static MPI_Datatype MPI_AMR_LUINT_TYPE;
static MPI_Datatype MPI_CODE_UINT_TYPE;
static std::once_flag mpi_init_flag;
void SetupMpiTypes() {
    MPI_Type_match_size(MPI_TYPECLASS_REAL, sizeof(DefReal), &MPI_REAL_DATA_TYPE);
    MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(DefInt), &MPI_INT_DATA_TYPE);
    MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(DefUint), &MPI_UINT_DATA_TYPE);
    MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(DefAmrLUint), &MPI_AMR_LUINT_TYPE);
    MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(DefSFCodeToUint), &MPI_CODE_UINT_TYPE);
}
MPI_Datatype GetMpiRealType() {
    std::call_once(mpi_init_flag, SetupMpiTypes);
    return MPI_REAL_DATA_TYPE;
}
MPI_Datatype GetMpiIntType() {
    std::call_once(mpi_init_flag, SetupMpiTypes);
    return MPI_INT_DATA_TYPE;
}
MPI_Datatype GetMpiUintType() {
    std::call_once(mpi_init_flag, SetupMpiTypes);
    return MPI_UINT_DATA_TYPE;
}
MPI_Datatype GetMpiAmrLUintType() {
    std::call_once(mpi_init_flag, SetupMpiTypes);
    return MPI_AMR_LUINT_TYPE;
}
MPI_Datatype GetMpiCodeUintType() {
    std::call_once(mpi_init_flag, SetupMpiTypes);
    return MPI_CODE_UINT_TYPE;
}
void MpiManager::SetUpMpi() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_ranks_);

    SetupMpiTypes();

    int data_size;
    MPI_Type_size(GetMpiRealType(), &data_size);
    if (sizeof(DefReal) != data_size) {
        LogManager::LogError("size of DefReal is not equal to size of MPI_REAL_DATA_TYPE) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    MPI_Type_size(GetMpiIntType(), &data_size);
    if (sizeof(DefInt) != data_size) {
        LogManager::LogError("size of DefInt is not equal to size of MPI_INT_DATA_TYPE) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    MPI_Type_size(GetMpiUintType(), &data_size);
    if (sizeof(DefInt) != data_size) {
        LogManager::LogError("size of DefInt is not equal to size of MPI_UINT_DATA_TYPE) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    MPI_Type_size(GetMpiAmrLUintType(), &data_size);
    if (sizeof(DefAmrLUint) != data_size) {
        LogManager::LogError("size of DefAmrLUint is not equal to size of MPI_AMR_LUINT_TYPE) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    MPI_Type_size(GetMpiCodeUintType(), &data_size);
    if (sizeof(DefSFCodeToUint) != data_size) {
        LogManager::LogError("size of DefSFCodeToUint is not equal to size of MPI_CODE_UINT_TYPE) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
void MpiManager::FinalizeMpi() {
    MPI_Finalize();
}
/**
 * @brief function to broadcast grid bounds of all ranks on rank 0 to other ranks.
 * @param[in] ptr_bitset_bounds pointer to bounds.
 */
void MpiManager::IniBroadcastSFBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds) {
    if (rank_id_ == 0) {  // bitset_bounds on rank has be calculated
        int bit_size = static_cast<int>(ptr_bitset_bounds->size()) * sizeof(DefSFBitset);
        MPI_Bcast(ptr_bitset_bounds->data(), bit_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
}
/**
 * @brief function to send and receive nodes for interpolation.
 * @param[in] sfbitset_aux   class manage space filling curves.
 * @param[in] grid_info_lower class storting grid node information at one lower level on current rank.
 * @param[in, out] ptr_grid_info pointer to class storting grid node information on current rank
 */
void MpiManager::MpiCommunicationForInterpolation(
    const SFBitsetAuxInterface& sfbitset_aux, const GridInfoInterface& grid_info_lower,
    GridInfoInterface* const ptr_grid_info) const {
    std::vector<BufferSizeInfo> send_buffer_info, receive_buffer_info;
    std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
    std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_receive, vec_ptr_buffer_send;
    SendGhostNodeForInterpolation(sfbitset_aux, ptr_grid_info->vec_num_interp_nodes_receive_,
        grid_info_lower, ptr_grid_info, &send_buffer_info, &receive_buffer_info, &vec_vec_reqs_send,
        &vec_vec_reqs_receive, &vec_ptr_buffer_send, &vec_ptr_buffer_receive);

    WaitAndReadGhostNodeForInterpolation(send_buffer_info, receive_buffer_info,
        vec_ptr_buffer_receive, &vec_vec_reqs_send, &vec_vec_reqs_receive, ptr_grid_info);
}
/** 
 * @brief function to compute computational load on each rank.
 * @param[in] max_level maximum level of all grids.
 * @param[in] sfbitset_aux class manage space filling curves.
 * @param[in] vec_grid_info grid information at all levels.
 * @param[out] ptr_node_level pointer to nodes at which levels on current rank.
 * @param[out] ptr_num_of_nodes pointer to number of nodes at each level on current rank.
 * @return total computational load on current rank.
 */
DefAmrLUint MpiManager::ComputeComputationalLoadOnEachRank(
    const DefInt max_level, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<std::shared_ptr<GridInfoInterface>>& vec_grid_info,
    std::map<DefSFCodeToUint, BackgroundLoadData>* const ptr_node_level,
    std::vector<DefAmrLUint>* const ptr_num_of_nodes) const {
    DefAmrLUint load_sum = 0;
    ptr_node_level->clear();
    ptr_num_of_nodes->resize(num_of_ranks_, 0);
    if (ptr_num_of_nodes->size() < vec_grid_info.size()) {
        LogManager::LogError("max_level is smaller than the number of grids in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefSFCodeToUint code_background;
    BackgroundLoadData background_load_data(max_level + 1);
    auto func_set_inner = [this, &sfbitset_aux](const DefSFBitset& sfbitset_in) {
    };

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
                ptr_num_of_nodes->at(node_cost) += 1;
                code_background = sfbitset_aux.SFBitsetoSFCode(
                    sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first));
                if (ptr_node_level->find(code_background) == ptr_node_level->end()) {
                    ptr_node_level->insert({code_background, background_load_data});
                }
                ptr_node_level->at(code_background).level_bitset_.Set(i_level, true);
                ptr_node_level->at(code_background).num_of_grid_nodes_ += 1;
                ptr_node_level->at(code_background).total_load_ += node_cost;
            }
        }
        auto func_add_node = [i_level, &sfbitset_aux, &ptr_node_level, &iter_grid](
            const std::vector<rootproject::DefMap<rootproject::DefInt>>& vec_interface) {
            DefSFCodeToUint code_background;
            for (const auto& iter_layer : vec_interface) {
                for (const auto& iter_node : iter_layer) {
                    if (iter_grid->map_grid_node_.find(iter_node.first)
                        != iter_grid->map_grid_node_.end()) {
                        if (!(iter_grid->map_grid_node_.at(iter_node.first)->flag_status_
                            & NodeBitStatus::kNodeStatusMpiPartitionOuter_)) {
                            code_background = sfbitset_aux.SFBitsetoSFCode(
                                sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first));
                            ptr_node_level->at(code_background).num_of_interface_nodes_ += 1;
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
            func_add_node(iter_interface.second->vec_outer_fine2coarse_);
            func_add_node(iter_interface.second->vec_inner_coarse2fine_);
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
    MPI_Offset offset = sizeof(DefAmrLUint);
    if (rank_id_ == 0) {
        MPI_File_write_at_all(mpi_file, 0, &load_all, 1, GetMpiAmrLUintType(), &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, &max_level, 1, GetMpiIntType(), &mpi_status);
    } else {
        MPI_File_write_at_all(mpi_file, 0, nullptr, 0, GetMpiAmrLUintType(), &mpi_status);
        MPI_File_write_at_all(mpi_file, offset, nullptr, 0, GetMpiIntType(), &mpi_status);
    }
    offset += sizeof(DefInt);
    const int num_size = sizeof(int);
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

    if (map_node_level.size() * node_size > (std::numeric_limits<int>::max)()) {
        LogManager::LogError("buffer size exceeds the maximum of int in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
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
    std::ifstream file_data(file_name, std::ios::binary);
    if (!file_data.is_open()) {
        LogManager::LogError("Failed to open file: " + file_name + " in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return 0;
    }
    DefAmrLUint load_read = 0;
    file_data.read(reinterpret_cast<char*>(&load_read), sizeof(DefAmrLUint));
    DefInt max_level = 0;
    file_data.read(reinterpret_cast<char*>(&max_level), sizeof(DefInt));
    ptr_map_node_level->clear();
    DefSFCodeToUint key;
    BackgroundLoadData value(max_level + 1);
    const int level_size = static_cast<int>(value.level_bitset_.GetNumBytes());
    const int num_size = sizeof(int), sfbitset_size = sizeof(DefSFBitset);
    const int load_size = sizeof(int);

    while (file_data) {
        file_data.read(reinterpret_cast<char*>(&key), sfbitset_size);
        if (!file_data) break;

        file_data.read(reinterpret_cast<char*>(value.level_bitset_.Data()), level_size);
        if (!file_data) break;

        file_data.read(reinterpret_cast<char*>(&value.num_of_grid_nodes_), num_size);
        if (!file_data) break;

        file_data.read(reinterpret_cast<char*>(&value.num_of_interface_nodes_), num_size);
        if (!file_data) break;

        file_data.read(reinterpret_cast<char*>(&value.total_load_), load_size);
        if (!file_data) break;

        ptr_map_node_level->insert({key, value});
    }
    file_data.close();

    return load_read;
}
/**
 * @brief function to compute minimum and maximum space filling codes for each rank.
 * @param[in] load_all total computational load for all ranks.
 * @param[in] vec_cost computational cost at each level.
 * @param[in] map_node_level nodes at which levels on current rank.
 * @param[out] ptr_bitset_min pointer to minimum space filling codes for each rank.
 * @param[out] ptr_bitset_max pointer to maximum space filling codes for each rank.
 * @param[out] ptr_num_of_grid_nodes pointer to number of grid nodes at each level on current rank.
 * @param[out] ptr_num_of_interface_nodes pointer to number of grid nodes at each level on current rank.
 */
void MpiManager::ComputeMinNMaxSFbitsetForEachRank(const DefAmrLUint load_all, const std::vector<DefInt>& vec_cost,
    const std::map<DefSFCodeToUint, BackgroundLoadData>& map_node_level,
    std::vector<DefSFBitset>* const ptr_bitset_min, std::vector<DefSFBitset>* const ptr_bitset_max,
    std::vector<DefAmrLUint>* const ptr_num_of_grid_nodes, std::vector<DefAmrLUint>* const ptr_interface_of_nodes) {
    const int num_ranks = num_of_ranks_;
    DefAmrLUint ave_load = static_cast<DefAmrLUint>(load_all / num_ranks) + 1;
    DefAmrLUint load_rank0 = load_all - (num_ranks - 1) * ave_load;
    std::vector<DefAmrLUint> rank_load(num_ranks, ave_load);
    rank_load.at(0) = load_rank0;   // load at rank 0, assuming lower than other ranks
    // traverse background nodes
    ptr_bitset_min->resize(num_ranks);
    ptr_bitset_max->resize(num_ranks);
    DefAmrLUint load_count = 0;
    int i_rank = 0;
    ptr_bitset_min->at(i_rank) = static_cast<DefSFBitset>(std::prev(map_node_level.begin())->first);
    ptr_bitset_max->back() = static_cast<DefSFBitset>(std::prev(map_node_level.end())->first);
    DefInt node_cost = 0;
    DefInt max_level = static_cast<DefInt>(vec_cost.size() - 1);
    ptr_num_of_grid_nodes->clear();
    ptr_num_of_grid_nodes->resize(num_of_ranks_, 0);
    ptr_interface_of_nodes->clear();
    ptr_interface_of_nodes->resize(num_of_ranks_, 0);
    for (const auto& iter : map_node_level) {
        if (load_count >= rank_load.at(i_rank)) {
            load_count = 0;
            ++i_rank;
            ptr_bitset_min->at(i_rank) = static_cast<DefSFBitset>(iter.first);
        }
        node_cost = iter.second.total_load_;
        ptr_num_of_grid_nodes->at(i_rank) += iter.second.num_of_grid_nodes_;
        ptr_interface_of_nodes->at(i_rank) += iter.second.num_of_interface_nodes_;
        if (load_count + node_cost >= rank_load.at(i_rank)) {
            ptr_bitset_max->at(i_rank) = static_cast<DefSFBitset>(iter.first);
        }
        load_count += node_cost;
    }
    if (load_count != load_all) {
        LogManager::LogError("computed load is not equal to the input load in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
