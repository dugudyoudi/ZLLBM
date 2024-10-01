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
    for (const auto& iter : ptr_grid_manager_->vec_ptr_grid_info_) {
        if (iter->i_level_ == 1) {
            ptr_mpi_manager_->CheckMpiNodesCorrespondence(*iter.get());
            ptr_mpi_manager_->CheckMpiPeriodicCorrespondence(*iter.get());
        }
    }
#endif  // ENABLE_MPI
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_CHECK_GRID
