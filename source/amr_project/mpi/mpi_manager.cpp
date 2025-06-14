//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
