//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2023-4-16
* @note .
*/
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "grid/grid_manager.h"
#include "io/log_write.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
void MpiManager::StartupMpi(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_ranks_);

    if (sizeof(DefReal) != sizeof(MPI_REAL_DATA_TYPE)) {
        LogError("size of DefReal is not equal to size of MPI_REAL_DATA_TYPE)");
    } else if (sizeof(DefInt) != sizeof(MPI_INT_DATA_TYPE)) {
        LogError("size of DefInt is not equal to size of MPI_INT_DATA_TYPE)");
    } else if (sizeof(DefLInt) != sizeof(MPI_LINT_DATA_TYPE)) {
        LogError("size of DefLInt is not equal to size of MPI_LINT_DATA_TYPE)");
    } else if (sizeof(DefUint) != sizeof(MPI_UINT_DATA_TYPE)) {
        LogError("size of DefUint is not equal to size of MPI_UINT_DATA_TYPE)");
    } else if (sizeof(DefLUint) != sizeof(MPI_LUINT_DATA_TYPE)) {
        LogError("size of DefLUint is not equal to size of MPI_LUINT_DATA_TYPE)");
    } else if (sizeof(DefTypeUint) != sizeof(MPI_TYPEUINT_DATA_TYPE)) {
        LogError("size of DefTypeUint is not equal to size of MPI_TYPEUINT_DATA_TYPE)");
    } else if (sizeof(DefSFCodeToUint) != sizeof(MPI_CODEUINT_DATA_TYPE)) {
        LogError("size of DefSFCodeToUint is not equal to size of MPI_CODEUINT_DATA_TYPE)");
    } else if (sizeof(DefSizet) != sizeof(MPI_SIZET_DATA_TYPE)) {
        LogError("size of DefSizet is not equal to size of MPI_SIZET_DATA_TYPE)");
    }
}
void MpiManager::FinalizeMpi() {
    MPI_Finalize();
}
void MpiManager::SetMpiParameters() {
}
/**
 * @brief function to broadcast grid bounds of all ranks on rank 0 to other ranks.
 * @param[in] ptr_bitset_bounds pointer to bounds.
 */
void MpiManager::IniBroadcastBitsetBounds(std::vector<DefSFBitset>* const ptr_bitset_bounds) {
    if (rank_id_ == 0) {  // bitset_bounds on rank has be calculated
        int bit_size = static_cast<int>(ptr_bitset_bounds->size()) * sizeof(DefSFBitset);
        MPI_Bcast(ptr_bitset_bounds->data(), bit_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ENABLE_MPI
