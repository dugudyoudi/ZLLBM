//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include "mpi/mpi_manager.h"
#ifdef ENABLE_MPI
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
void MpiManager::StartupMpi(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &numb_of_ranks_);
    if (rank_id_ == 0) {
        LogInfo(
            "number of processes is: " + std::to_string(numb_of_ranks_));
}
void MpiManager::PartiteBackground() {
    num_of_nodes_per_mpiblock_ =
        static_cast<DefSizet>(1) << level_for_mpiblock_;

}
void MpiManager::SetMpiParameters() {
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif ENABLE_MPI
