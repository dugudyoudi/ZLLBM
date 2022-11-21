//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#ifdef ENABLE_MPI
#include "mpi/mpi_manager.h"
#include "grid/grid_manager.h"
#include "io/io_manager.h"
namespace rootproject {
namespace amrproject {
void MpiManager::StartupMpi(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id_);
    MPI_Comm_size(MPI_COMM_WORLD, &numb_of_ranks_);
}
void MpiManager::PartiteBackground() {
    num_of_nodes_per_mpiblock_ =
        static_cast<DefSizet>(1) << level_for_mpiblock_;

}
void MpiManager::SetMpiParameters() {
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif
