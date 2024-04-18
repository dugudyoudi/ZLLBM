//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_info_mpi.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid informaiton for mpi communication.
* @date  2022-6-7
*/
#include <string>
#include "grid/grid_info_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to set mpi communication region for periodic boundary conditions.
 * @param[in] dims  dimensions of the grid.
 * @param[in] sfbitset_aux class to manage space filling coe.
 * @param[in] bool_periodic_min booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] bool_periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[out] ptr_partition_interface pointer to container storing nodes on partition information.
 */
void GridInfoInterface::SetPeriodicBoundaryAsPartitionInterface(DefAmrIndexUint dims,
    const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
    DefMap<DefAmrIndexUint>* const ptr_partition_interface) {
    DefAmrUint flag0 = ~0;
    DefAmrIndexUint i_level = i_level_;
    if (bool_periodic_min.size() != dims) {
        LogManager::LogError("size of bool_periodic_min " + std::to_string(bool_periodic_min.size())
            + " is not equal to the grid dimension " + std::to_string(dims) + "in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (bool_periodic_max.size() != dims) {
        LogManager::LogError("size of bool_periodic_max " + std::to_string(bool_periodic_min.size())
            + " is not equal to the grid dimension " + std::to_string(dims) + "in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    for (DefAmrIndexUint i_dims = 0; i_dims < dims; ++i_dims) {
        if (bool_periodic_min.at(i_dims)) {
            ptr_partition_interface->insert(
                domain_boundary_min_.at(i_dims).begin(), domain_boundary_min_.at(i_dims).end());
        }
        if (bool_periodic_max.at(i_dims)) {
            ptr_partition_interface->insert(
                domain_boundary_max_.at(i_dims).begin(), domain_boundary_max_.at(i_dims).end());
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
