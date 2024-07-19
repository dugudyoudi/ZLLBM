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
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to set mpi communication region for periodic boundary conditions.
 * @param[in] dims  dimensions of the grid.
 * @param[in] sfbitset_aux class to manage space filling codes.
 * @param[in] bool_periodic_min booleans indicating if the boundary is periodic at minimum domain boundaries.
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
/**
 * @brief function to set mpi communication region for periodic boundary conditions.
 * @param[in] bool_periodic_min booleans indicating if the boundary is periodic at minimum domain boundaries.
 * @param[in] bool_periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] sfbitset_aux class to manage space filling codes.
 * @param[in] refinement_interface nodes on the refinement interface.
 * @param[in] map_nodes_lower grid nodes at one lower refinement level.
 */
int GridInfoInterface::AddGhostNodesForInterpolation(const std::vector<bool>& bool_periodic_min,
    const std::vector<bool>& bool_periodic_max, const SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<DefAmrUint>& refinement_interface, const DefMap<std::unique_ptr<GridNode>>& map_nodes_lower) {
    DefSFBitset sfbitset_lower, sfbitset_current;
    std::vector<DefSFBitset> vec_in_region;
    size_t dims = k0VecBitsetDomainMin_.size();
    std::vector<DefSFBitset> domain_min(dims), domain_max(dims);
    for (DefAmrIndexUint i_dims = 0; i_dims < dims; ++i_dims) {
        domain_min.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(1, k0VecBitsetDomainMin_.at(i_dims));
        domain_max.at(i_dims) = sfbitset_aux.SFBitsetToNLowerLevelVir(1, k0VecBitsetDomainMax_.at(i_dims));
    }
    for (auto& iter_node : refinement_interface) {
        sfbitset_lower = sfbitset_aux.SFBitsetToNLowerLevelVir(1, iter_node.first);
        sfbitset_aux.FindNodesInPeriodicReginOfGivenLength(sfbitset_lower,
            max_interp_length_, bool_periodic_min, bool_periodic_max,
            domain_min, domain_max, &vec_in_region);
        for (const auto& iter_node_region : vec_in_region) {
            if (map_nodes_lower.find(iter_node_region) ==  map_nodes_lower.end()) {
                sfbitset_current = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_node_region);
                if (map_grid_node_.find(iter_node_region) == map_grid_node_.end()) {
                    interp_nodes_outer_layer_.insert({iter_node_region, GridNodeCreator()});
                }
            }
        }
    }
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject
