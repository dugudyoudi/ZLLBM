//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file mpi_manager_debug.cpp
* @author Zhengliang Liu
* @brief functions used to debug mpi communication.
*/
#include "./amr_manager.h"
#ifdef DEBUG_CHECK_GRID
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
void AmrManager::CheckMeshAfterInitialization() const {
#ifdef ENABLE_MPI
    ptr_mpi_manager_->DebugMpiForAllGrids(*ptr_grid_manager_.get());
#endif  // ENABLE_MPI
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_CHECK_GRID
