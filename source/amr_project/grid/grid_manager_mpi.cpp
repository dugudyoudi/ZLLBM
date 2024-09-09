//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager_mpi.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is used.
* @date  2022-6-7
*/
#include <mpi.h>
#include <string>
#include <set>
#include <functional>
#include "./auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to instantiate grid nodes on refinement interface and mpi communication layer.
 * @param[in] i_level refinement level.
 * @param[in] num_partition_outer_layer number of outer layers for mpi communication.
 * @param[in] flag_outmost_refinement_n_outer_mpi  flag indicates the node is on both
 *                  outmost layer of coarse2fine or fine2coarse and mpi outer layer.
 * @param[in] code_min  minimum space filling code of current rank at background level.
 * @param[in] code_max  maximum space filling code of current rank at background level.
 * @param[in] sfbitset_aux  class manage space filling curves.
 * @param[in] sfbitset_one_lower_level space filling codes at one lower refinement level.
 * @param[in] sfbitset_partition_interface_background nodes on the background partitioned interface of background level on current rank.
 * @param[out] ptr_outer_layer_current_level pointer to nodes on the outer layer at current refinement level.
 * @param[out] ptr_outer_layer_lower_level pointer to nodes on the outer layer at one lower refinement level.
 */
void GridManagerInterface::InstantiateOverlapLayerOfRefinementInterfaceMpi(
    const DefAmrIndexUint i_level, const DefAmrIndexUint num_partition_outer_layer,
    const DefAmrUint flag_outmost_refinement_n_outer_mpi,
    const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefAmrIndexUint>& map_sfbitset_one_lower_level,
    const DefMap<DefAmrIndexUint>& sfbitset_partition_interface_background,
    DefMap<DefAmrIndexUint>* const ptr_outer_layer_current_level,
    DefMap<DefAmrIndexUint>* const ptr_outer_layer_lower_level) {
    InterfaceLayerInfo* ptr_interface_info = nullptr;
    InterfaceLayerInfo* ptr_interface_info_lower = nullptr;
    DefAmrUint flag_tmp, flag_0 = NodeBitStatus::kNodeStatus0_;
    int maxlayer;
    if (i_level == 0) {
        LogManager::LogError("input refinement level should be greater than 0 in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefAmrIndexUint i_lower_level = i_level - 1;
    GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
    GridInfoInterface& grid_info_lower = *(vec_ptr_grid_info_.at(i_lower_level));
    DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
    DefMap<std::unique_ptr<GridNode>>& map_grid_lower = grid_info_lower.map_grid_node_;
    DefSFCodeToUint code_background;
#ifdef DEBUG_CHECK_GRID
    if (grid_info_lower.k0NumCoarse2FineLayer_ < 2) {
        LogManager::LogError("number of coarse to fine layers at level "
            + std::to_string(i_level - 1) + " is at least 1"
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (grid_info.k0NumFine2CoarseLayer_ < 3) {
        LogManager::LogError("number of fine to coarse layers at level "
            + std::to_string(i_level) + " is at least 3"
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID
    // find interface node at current level
    DefAmrIndexUint layer_coarse_outer = 1,
        layer_coarse_innermost = 0,
        layer0 = grid_info.k0NumFine2CoarseLayer_ - 1,
        layer_m1 = layer0 - 1,
        layer_m2 = layer_m1 - 1;
    std::vector<DefSFBitset> sfbitset_neighbors;
    for (auto& iter_interface : grid_info_lower.map_ptr_interface_layer_info_) {
        ptr_interface_info_lower = iter_interface.second.get();
        if (grid_info.map_ptr_interface_layer_info_.find(iter_interface.first)
            == grid_info.map_ptr_interface_layer_info_.end()) {
            grid_info.map_ptr_interface_layer_info_.insert(
                { iter_interface.first, std::make_shared<InterfaceLayerInfo>() });
        }
        ptr_interface_info = grid_info.map_ptr_interface_layer_info_
            .at(iter_interface.first).get();
        // interface inside the geometry
        if (ptr_interface_info_lower->vec_inner_coarse2fine_.size() > 0) {
            ptr_interface_info->vec_inner_fine2coarse_.resize(grid_info.k0NumFine2CoarseLayer_);
            FindOverlappingLayersBasedOnOutermostCoarse(
                ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_innermost),
                map_sfbitset_one_lower_level,
                &ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_outer),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
            // insert node instance
            maxlayer = static_cast<int>(ptr_interface_info->vec_inner_fine2coarse_.size());
            for (int ilayer = 0; ilayer < maxlayer; ++ilayer) {
                flag_tmp = flag_0;
                if (ilayer == maxlayer - 1) {
                    flag_tmp |= NodeBitStatus::kNodeStatusFine2Coarse0_;
                }
                if (ilayer >= maxlayer - grid_info.k0NumFine2CoarseGhostLayer_) {
                    flag_tmp |=  NodeBitStatus::kNodeStatusFine2CoarseGhost_;
                }
                for (const auto& iter_layer_node : ptr_interface_info
                    ->vec_inner_fine2coarse_.at(ilayer)) {
                    if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                            map_grid.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                        }
                    } else {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                    }
                    // check if node is on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                        if (ptr_outer_layer_current_level->find(iter_layer_node.first)
                            == ptr_outer_layer_current_level->end()) {
                            ptr_outer_layer_current_level->insert(
                                {iter_layer_node.first, flag_outmost_refinement_n_outer_mpi});
                        } else {
                            ptr_outer_layer_current_level->at(iter_layer_node.first) =
                                flag_outmost_refinement_n_outer_mpi;
                        }
                    }
                }
            }
            maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                flag_tmp = flag_0;
                if (ilayer == maxlayer - 1) {
                    flag_tmp |= NodeBitStatus::kNodeStatusCoarse2Fine0_;
                }
                if (ilayer >= maxlayer - grid_info_lower.k0NumCoarse2FineGhostLayer_) {
                    flag_tmp |= NodeBitStatus::kNodeStatusCoarse2FineGhost_;
                }
                for (const auto& iter_layer_node : ptr_interface_info_lower
                    ->vec_inner_coarse2fine_.at(ilayer)) {
                    if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                        }
                    } else {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_lower_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                        if (ptr_outer_layer_lower_level->find(iter_layer_node.first)
                            == ptr_outer_layer_lower_level->end()) {
                            ptr_outer_layer_lower_level->insert({iter_layer_node.first,
                                flag_outmost_refinement_n_outer_mpi});
                        } else {
                            ptr_outer_layer_lower_level->at(iter_layer_node.first) =
                                flag_outmost_refinement_n_outer_mpi;
                        }
                    }
                }
            }
        }

        // interface outside the geometry
        if (ptr_interface_info_lower->vec_outer_coarse2fine_.size() > 0) {
            ptr_interface_info->vec_outer_fine2coarse_.resize(
                grid_info.k0NumFine2CoarseLayer_);
            FindOverlappingLayersBasedOnOutermostCoarse(
                ptr_interface_info_lower->vec_outer_coarse2fine_.at(layer_coarse_innermost),
                map_sfbitset_one_lower_level,
                &ptr_interface_info_lower->vec_outer_coarse2fine_.at(layer_coarse_outer),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer0),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m1),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m2));

            // insert node instance
            DefAmrIndexUint maxlayer = DefAmrIndexUint(ptr_interface_info->vec_outer_fine2coarse_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                flag_tmp = flag_0;
                if (ilayer == maxlayer - 1) {
                    flag_tmp |= NodeBitStatus::kNodeStatusFine2Coarse0_;
                }
                if (ilayer >= maxlayer - grid_info.k0NumFine2CoarseGhostLayer_) {
                    flag_tmp |=  NodeBitStatus::kNodeStatusFine2CoarseGhost_;
                }
                for (const auto& iter_layer_node : ptr_interface_info->vec_outer_fine2coarse_.at(ilayer)) {
                    if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                            map_grid.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                        }
                    } else {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                    }

                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                        if (ptr_outer_layer_current_level->find(iter_layer_node.first)
                            == ptr_outer_layer_current_level->end()) {
                            ptr_outer_layer_current_level->insert(
                                {iter_layer_node.first, flag_outmost_refinement_n_outer_mpi});
                        } else {
                            ptr_outer_layer_current_level->at(iter_layer_node.first) =
                                flag_outmost_refinement_n_outer_mpi;
                        }
                    }
                }
            }
            maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_outer_coarse2fine_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                flag_tmp = flag_0;
                if (ilayer == maxlayer - 1) {
                    flag_tmp |= NodeBitStatus::kNodeStatusCoarse2Fine0_;
                }
                if (ilayer >= maxlayer - grid_info_lower.k0NumCoarse2FineGhostLayer_) {
                    flag_tmp |= NodeBitStatus::kNodeStatusCoarse2FineGhost_;
                }
                for (const auto& iter_layer_node : ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(ilayer)) {
                    if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                        }
                    } else {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_lower_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                        if (ptr_outer_layer_lower_level->find(iter_layer_node.first)
                            == ptr_outer_layer_lower_level->end()) {
                            ptr_outer_layer_lower_level->insert(
                                {iter_layer_node.first, flag_outmost_refinement_n_outer_mpi});
                        } else {
                            ptr_outer_layer_lower_level->at(iter_layer_node.first) =
                                flag_outmost_refinement_n_outer_mpi;
                        }
                    }
                }
            }
        }
    }
}
/**
 * @brief function to instantiate grid noes for all refinement levels and find layers for mpi communication.
 * @param[in] i_rank current rank id
 * @param[in] num_partition_inner_layer number of inner layers for mpi communication
 * @param[in] num_partition_outer_layer number of outer layers for mpi communication
 * @param[in] vec_sfcode_min  minimum space filling code of all ranks.
 * @param[in] vec_sfcode_max  maximum space filling code of all ranks.
 * @param[in] sfbitset_aux  class manage space filling curves.
 * @param[in] sfbitset_one_lower_level space filling codes at one lower refinement level.
 * @param[in] sfbitset_ghost_one_lower_level space filling codes of mpi communication node near coarse to fine refinement interface.
 * @param[in] sfbitset_partition_interface_0  nodes on the background partitioned interface of background level on current rank.
 * @param[out] ptr_mpi_inner_layer pointer to nodes on the inner layer for mpi communication (sending).
 * @param[out] ptr_mpi_outer_layer pointer to nodes on the outer layer for mpi communication (receiving). 
 */
void GridManagerInterface::InstantiateGridNodeAllLevelMpi(const int i_rank,
    const DefAmrIndexUint num_partition_inner_layer, const DefAmrIndexUint num_partition_outer_layer,
    const std::vector<DefSFCodeToUint>& vec_sfcode_min, const std::vector<DefSFCodeToUint>& vec_sfcode_max,
    const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level,
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_ghost_one_lower_level,
    const DefMap<DefAmrIndexUint>& sfbitset_partition_interface_0,
    std::vector<std::map<int, DefMap<DefAmrIndexUint>>>* const ptr_mpi_inner_layer,
    std::vector<DefMap<DefAmrIndexUint>>* const ptr_mpi_outer_layer) {
    DefAmrIndexUint flag_normal_outer_node = ~0, flag_on_coarse2fine = ~0 - 1,
        flag_outmost_refinement_n_outer_mpi = flag_on_coarse2fine, flag_on_periodic = flag_on_coarse2fine;
    DefMap<DefAmrIndexUint> background_occupied;  // background nodes coincide with those at higher levels
    ptr_mpi_inner_layer->resize(k0MaxLevel_ + 1);
    ptr_mpi_outer_layer->resize(k0MaxLevel_ + 1);
    DefSFCodeToUint code_min_background_level = vec_sfcode_min.at(i_rank),
        code_max_background_level = vec_sfcode_max.at(i_rank);
    std::vector<DefSFBitset> vec_in_region;
    std::vector<DefSFCodeToUint>::iterator iter_index;

    // initialize grid information, will be used for instantiate grid nodes
    for (DefAmrIndexUint i_level = 0; i_level <= k0MaxLevel_; ++i_level) {
        vec_ptr_grid_info_.at(i_level)->InitialGridInfo();
    }
    std::vector<DefSFCodeToUint> ull_max(vec_sfcode_max);
    std::vector<DefAmrIndexLUint> indices_min = GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = GetMaxIndexOfBackgroundNodeArrAsVec();

    std::vector<DefMap<DefAmrIndexUint>> sfbitset_extra_lower_level(k0MaxLevel_ + 1);
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
        DefAmrIndexUint i_level_lower = i_level - 1;
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        std::vector<bool> periodic_min(k0GridDims_, false), periodic_max(k0GridDims_, false);
        grid_info.CheckIfPeriodicDomainRequired(k0GridDims_, &periodic_min, &periodic_max);
        // initialize grid information
        DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
        std::vector<DefSFBitset> bitset_cell_lower, bitset_all, bitset_neighbors;
        DefSFBitset bitset_background;
        std::vector<DefSFBitset> domain_min_n_level(k0GridDims_), domain_max_n_level(k0GridDims_);
        sfbitset_aux.GetMinAtGivenLevel(i_level, indices_min, &domain_min_n_level);
        sfbitset_aux.GetMaxAtGivenLevel(i_level, indices_max, &domain_max_n_level);

        if (grid_info.ptr_solver_ == nullptr) {
            LogManager::LogError("solver has not been assigned to grid info"
            " of refinement level " + std::to_string(i_level)
            + " in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }

        // instantiate nodes on outer mpi communication layer and refinement interface
        // The spatial step of node on the coarse to fine interface at i_level
        // is half of those in sfbitset_one_lower_level used for instantiation.
        // Thus nodes on coarse to fine interfaces are sent from rank0 individually
        for (const auto& iter_ghost : sfbitset_ghost_one_lower_level.at(i_level)) {
            if (InstantiateGridNode(iter_ghost.first, vec_ptr_grid_info_.at(i_level_lower).get())) {
                vec_ptr_grid_info_.at(i_level_lower)->map_grid_node_.at(iter_ghost.first)->flag_status_
                    |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                ptr_mpi_outer_layer->at(i_level_lower).insert({iter_ghost.first, flag_on_coarse2fine});
            }
        }

        // instantiate refinement interfaces
        InstantiateOverlapLayerOfRefinementInterfaceMpi(i_level, num_partition_outer_layer,
            flag_outmost_refinement_n_outer_mpi, code_min_background_level,
            code_max_background_level, sfbitset_aux, sfbitset_one_lower_level.at(i_level),
            sfbitset_partition_interface_0, &ptr_mpi_outer_layer->at(i_level), &ptr_mpi_outer_layer->at(i_level_lower));

        // instantiate nodes stored in sfbitset_one_lower_level
        DefMap<DefAmrIndexUint> partition_interface_level;   // nodes on partition interface at current level
        DefMap<DefAmrIndexUint> outer_layer_tmp;
        int flag_node;
        DefSFCodeToUint code_tmp, code_neighbor_tmp;
        DefSFBitset sfbitset_current_level;
        bool bool_partition_outside, bool_partition_inside;

        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            if (CheckCoincideBackground(i_level_lower, iter_low.first, &bitset_background)) {
                if (background_occupied.find(bitset_background) == background_occupied.end()) {
                    background_occupied.insert({ bitset_background, kFlagSize0_ });
                }
            }
            sfbitset_current_level = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_low.first);
            flag_node = grid_info.CheckIfNodeOutsideCubicDomain(k0GridDims_, sfbitset_current_level, sfbitset_aux);
            (this->*ptr_func_insert_domain_boundary_)(flag_node, sfbitset_current_level, &grid_info);

            if (NodesBelongToOneCell(iter_low.first,
                sfbitset_one_lower_level.at(i_level), &bitset_cell_lower)) {
                bool_partition_outside = false;
                bool_partition_inside = false;
                FindAllNodesInACellAtOneLevelLower(bitset_cell_lower, &bitset_all);
                for (const auto& iter_node : bitset_all) {
                    code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                    if (code_tmp < code_min_background_level || code_tmp > code_max_background_level) {
                        bool_partition_outside = true;
                        if (outer_layer_tmp.find(iter_node) == outer_layer_tmp.end()) {
                            outer_layer_tmp.insert({iter_node, flag_normal_outer_node});
                        }
                    } else {
                        bool_partition_inside = true;
                        if (map_grid.find(iter_node) == map_grid.end()) {
                            InstantiateGridNode(iter_node, &grid_info);
                        }
                    }
                }
                // one or more nodes in the cell should be on the partition interface.
                // check if one node is on the partition interface by finding if it is inside
                // and at least one of its neighbors is outside the partition block
                if (bool_partition_outside && bool_partition_inside) {
                    for (const auto& iter_node : bitset_all) {
                        code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                        if (code_tmp >= code_min_background_level && code_tmp <= code_max_background_level
                            && partition_interface_level.find(iter_node) == partition_interface_level.end()) {
                            sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_node,
                                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                            for (const auto& iter_neighbor : bitset_neighbors) {
                                code_neighbor_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(
                                    i_level, iter_neighbor).to_ullong();
                                if (code_neighbor_tmp < code_min_background_level
                                    || code_neighbor_tmp > code_max_background_level) {
                                    partition_interface_level.insert({iter_node, 0});
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        DefMap<DefAmrIndexUint> inner_layer_tmp;
        InstantiateDomainBoundaryForMpi(num_partition_outer_layer - 1, code_min_background_level,
            code_max_background_level, periodic_min, periodic_max, sfbitset_one_lower_level.at(i_level),
            sfbitset_aux, vec_ptr_grid_info_.at(i_level).get(), &inner_layer_tmp, &ptr_mpi_outer_layer->at(i_level));

        // Extend mpi interface nodes for a given number of inner communication layers.
        // Noting that nodes on the partitioned interface is considered in the inner layer.
        // Use temporal container (inner_layer_tmp) at this stage since some of them inserted
        // here are on the periodic boundaries but not on the inner layer. In addition, some
        // inner nodes should have duplicates since they will be send to different ranks (overlap region).
        for (const auto& iter_interface : partition_interface_level) {
            sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_interface.first, num_partition_inner_layer-1,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                if ((code_tmp >= code_min_background_level && code_tmp <= code_max_background_level)
                    && inner_layer_tmp.find(iter_node) == inner_layer_tmp.end()
                    && map_grid.find(iter_node) != map_grid.end()) {
                    inner_layer_tmp.insert({iter_node, kFlag0_});
                }
            }
            sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_interface.first, num_partition_outer_layer,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                if ((code_tmp < code_min_background_level || code_tmp > code_max_background_level)
                    &&outer_layer_tmp.find(iter_node) != outer_layer_tmp.end()) {
                    // instantiate node in mpi outer layer only
                    ptr_mpi_outer_layer->at(i_level).insert({iter_node, flag_normal_outer_node});
                    InstantiateGridNode(iter_node, &grid_info);
                    map_grid.at(iter_node)->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                }
            }
        }

        // extend one mpi layer when necessary since mesh generation is based on one lower level
        DefMap<DefAmrIndexUint> map_extra_expand_outer;
        bool bool_interface_upper_extra = ((num_partition_outer_layer%2) == 0);
        DefSFBitset sfbitset_tmp;
        DefSFCodeToUint code_interface;
        std::vector<DefSFBitset> vec_periodic;
        bool bool_has_periodic_boundary =
            grid_info.CheckIfPeriodicDomainRequired(k0GridDims_, &periodic_min, &periodic_max);
        for (const auto& iter_interface : ptr_mpi_outer_layer->at(i_level)) {
            // set ranks that should be sent by mpi inner layer nodes
            sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_interface.first, num_partition_inner_layer,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            code_interface = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_interface.first).to_ullong();
            SetUpRanksSentForMpiInnerLayer(i_rank, code_interface, vec_in_region,
                ull_max, inner_layer_tmp, &ptr_mpi_inner_layer->at(i_level));
            // extend the mpi outer layers when necessary
            if ((bool_interface_upper_extra && code_interface > code_max_background_level)
                && iter_interface.second == flag_normal_outer_node) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                    periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                for (const auto& iter_neighbor : bitset_neighbors) {
                    if (map_extra_expand_outer.find(iter_neighbor) == map_extra_expand_outer.end()) {
                        if (map_grid.find(iter_neighbor) == map_grid.end()) {
                            InstantiateGridNode(iter_neighbor, &grid_info);
                            map_grid.at(iter_neighbor)->flag_status_
                                |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                            map_extra_expand_outer.insert({iter_neighbor, 0});
                        } else {
                            code_neighbor_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(
                                i_level, iter_neighbor).to_ullong();
                            if (code_neighbor_tmp > code_max_background_level) {
                                map_grid.at(iter_neighbor)->flag_status_
                                    |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                                map_extra_expand_outer.insert({iter_neighbor, 0});
                            }
                        }
                    }
                }
            } else if ((!bool_interface_upper_extra && code_interface < code_min_background_level)
                && iter_interface.second == flag_normal_outer_node) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                    periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                for (const auto& iter_neighbor : bitset_neighbors) {
                    if (map_extra_expand_outer.find(iter_neighbor) == map_extra_expand_outer.end()) {
                        if (map_grid.find(iter_neighbor) == map_grid.end()) {
                            InstantiateGridNode(iter_neighbor, &grid_info);
                            map_grid.at(iter_neighbor)->flag_status_
                                |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                            map_extra_expand_outer.insert({iter_neighbor, 0});
                        } else {
                            code_neighbor_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(
                                i_level, iter_neighbor).to_ullong();
                            if (code_neighbor_tmp < code_min_background_level) {
                                map_grid.at(iter_neighbor)->flag_status_
                                    |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                                map_extra_expand_outer.insert({iter_neighbor, 0});
                            }
                        }
                    }
                }
            }
        }

        // add the outmost layer to mpi outer communication layer
        std::vector<DefSFBitset> vec_coarse;
        for (const auto& iter_node : map_extra_expand_outer) {
            ptr_mpi_outer_layer->at(i_level).insert({iter_node.first, flag_normal_outer_node});
            sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_node.first, num_partition_inner_layer,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            code_interface = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first).to_ullong();
            SetUpRanksSentForMpiInnerLayer(i_rank, code_interface, vec_in_region,
                ull_max, inner_layer_tmp, &ptr_mpi_inner_layer->at(i_level));
            sfbitset_aux.FindNeighboringCoarseFromFine(iter_node.first, &vec_coarse);
            for (const auto& iter_coarse : vec_coarse) {
                if (sfbitset_extra_lower_level.at(i_level).find(iter_coarse)
                    == sfbitset_extra_lower_level.at(i_level).end()) {
                    sfbitset_extra_lower_level.at(i_level).insert({iter_coarse, 0});
                }
            }
        }

        // setup node flag for mpi inner layers
        for (const auto& iter_layer : ptr_mpi_inner_layer->at(i_level)) {
            for (const auto& iter_node : iter_layer.second) {
                map_grid.at(iter_node.first)->flag_status_
                    |= amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_;
            }
        }
        if (i_level < k0MaxLevel_) {
            SetUpMpiInnerLayerForHigherLevel(i_rank, i_level, num_partition_outer_layer,
                code_min_background_level, code_max_background_level,
                periodic_min, periodic_max, domain_min_n_level,
                domain_max_n_level, ull_max, partition_interface_level,
                ptr_mpi_outer_layer->at(i_level), vec_ptr_grid_info_.at(i_level + 1)->map_grid_node_, sfbitset_aux,
                &ptr_mpi_inner_layer->at(i_level), &vec_ptr_grid_info_.at(i_level)->map_grid_node_);
        }
    }

    // find overlapping node for refinement levels of 0 and 1
    for (const auto& iter_interfaces : vec_ptr_grid_info_.at(0)->map_ptr_interface_layer_info_) {
        for (const auto& iter_coarse2fine : iter_interfaces.second->vec_inner_coarse2fine_) {
            for (const auto& iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
        }
        for (const auto& iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
            for (const auto& iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
        }
    }
    if (k0MaxLevel_ > 0) {
        for (const auto& iter_ghost : sfbitset_ghost_one_lower_level.at(1)) {
            background_occupied.erase(iter_ghost.first);
        }
    }

    InstantiateBackgroundGrid(code_min_background_level, code_max_background_level, background_occupied);

    // find background nodes on mpi communication interface
    DefMap<std::unique_ptr<GridNode>>& map_grid_level_0 = vec_ptr_grid_info_.at(0)->map_grid_node_;
    DefMap<DefAmrIndexUint> map_outmost_current, map_outmost_pre;
    std::vector<DefSFBitset> vec_neighbors;
    DefMap<DefAmrIndexUint> partition_interface_level(sfbitset_partition_interface_0);
    std::vector<bool> periodic_min(k0GridDims_, false), periodic_max(k0GridDims_, false);;
    vec_ptr_grid_info_.at(0)->CheckIfPeriodicDomainRequired(k0GridDims_, &periodic_min, &periodic_max);
    for (const auto& iter_interface : partition_interface_level) {
        if (background_occupied.find(iter_interface.first) == background_occupied.end()) {
            map_outmost_current.insert({iter_interface.first, 0});
        }
    }

    DefMap<DefAmrIndexUint> inner_layer_tmp;
    std::vector<DefSFBitset> domain_min_n_level(k0GridDims_), domain_max_n_level(k0GridDims_);
    sfbitset_aux.GetMinAtGivenLevel(0, indices_min, &domain_min_n_level);
    sfbitset_aux.GetMaxAtGivenLevel(0, indices_max, &domain_max_n_level);

    vec_ptr_grid_info_.at(0)->SetPeriodicBoundaryAsPartitionInterface(
        k0GridDims_, sfbitset_aux, periodic_min, periodic_max, &map_outmost_current);
    // find nodes on inner communication layers at level 0
    DefSFCodeToUint code_tmp;
    for (const auto& iter_interface : map_outmost_current) {
        if (map_grid_level_0.find(iter_interface.first) != map_grid_level_0.end()) {
            sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_interface.first, num_partition_inner_layer-1,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                if (iter_node != SFBitsetAuxInterface::kInvalidSFbitset) {
                    code_tmp = iter_node.to_ullong();
                    if ((code_tmp >= code_min_background_level && code_tmp <= code_max_background_level)
                        && inner_layer_tmp.find(iter_node) == inner_layer_tmp.end()
                        && background_occupied.find(iter_node) == background_occupied.end()) {
                        inner_layer_tmp.insert({iter_node, 0});
                    }
                }
            }
        }
    }

    // add ghost layers to outer layer
    if (k0MaxLevel_ > 0) {
        for (const auto& iter_node : sfbitset_ghost_one_lower_level.at(1)) {
            if (map_grid_level_0.find(iter_node.first) != map_grid_level_0.end()) {
                code_tmp = iter_node.first.to_ullong();
                if (code_tmp < code_min_background_level || code_tmp > code_max_background_level) {
                    ptr_mpi_outer_layer->at(0).insert({iter_node.first, flag_normal_outer_node});
                    map_grid_level_0.at(iter_node.first)->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                }
            }
        }
    }

    // find and instantiate nodes on outer communication layers layer by layer
    for (auto i_layer = 0; i_layer < num_partition_outer_layer; ++i_layer) {
        for (const auto& iter_interface : map_outmost_current) {
            if (map_grid_level_0.find(iter_interface.first) != map_grid_level_0.end()) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                    periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_neighbors);
                for (const auto& iter_node : vec_neighbors) {
                    code_tmp = iter_node.to_ullong();
                    if ((code_tmp < code_min_background_level || code_tmp > code_max_background_level)
                        && background_occupied.find(iter_node) == background_occupied.end()) {
                        if (ptr_mpi_outer_layer->at(0).find(iter_node) == ptr_mpi_outer_layer->at(0).end()) {
                            ptr_mpi_outer_layer->at(0).insert({iter_node, flag_normal_outer_node});
                            if (map_grid_level_0.find(iter_node) == map_grid_level_0.end()) {
                                if (InstantiateGridNode(iter_node, vec_ptr_grid_info_.at(0).get())) {
                                    map_grid_level_0.at(iter_node)->flag_status_
                                        |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
                                }
                            }
                            map_outmost_pre.insert({iter_node, 0});
                        } else if (ptr_mpi_outer_layer->at(0).at(iter_node) != flag_outmost_refinement_n_outer_mpi) {
                            ptr_mpi_outer_layer->at(0).at(iter_node) = flag_normal_outer_node;
                            map_outmost_pre.insert({iter_node, 0});
                        }
                    }
                }
            }
        }
        map_outmost_current = map_outmost_pre;
        map_outmost_pre.clear();
    }

    int node_rank;
    int flag_node;
    // set ranks that should be sent by mpi inner layer nodes (level 0)
    for (const auto& iter_interface : ptr_mpi_outer_layer->at(0)) {
        sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_interface.first, num_partition_inner_layer,
            periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
        for (const auto& iter_node : vec_in_region) {
            code_tmp = iter_node.to_ullong();
            if ((code_tmp >= code_min_background_level && code_tmp <= code_max_background_level)
                && inner_layer_tmp.find(iter_node) != inner_layer_tmp.end()) {
                iter_index = std::lower_bound(ull_max.begin(), ull_max.end(), iter_interface.first.to_ullong());
                node_rank = static_cast<int>(iter_index - ull_max.begin());
                if (node_rank != i_rank) {
                    if (map_grid_level_0.find(iter_node) != map_grid_level_0.end()) {
                        if (ptr_mpi_inner_layer->at(0).find(node_rank)
                            == ptr_mpi_inner_layer->at(0).end()) {
                            ptr_mpi_inner_layer->at(0).insert({node_rank, {}});
                        }
                        if ((ptr_mpi_inner_layer->at(0).at(node_rank).find(iter_node)
                            == ptr_mpi_inner_layer->at(0).at(node_rank).end())) {
                            ptr_mpi_inner_layer->at(0).at(node_rank).insert({iter_node, 0});
                            map_grid_level_0.at(iter_node)->flag_status_
                                |= amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_;
                        }
                    }
                }
            }
        }
        flag_node = vec_ptr_grid_info_.at(0)->CheckIfNodeOutsideCubicDomain(
            k0GridDims_, iter_interface.first, sfbitset_aux);
        (this->*ptr_func_insert_domain_boundary_)(
            flag_node, iter_interface.first, vec_ptr_grid_info_.at(0).get());
    }

    if (k0MaxLevel_ > 0) {
        SetUpMpiInnerLayerForHigherLevel(i_rank, 0, num_partition_outer_layer, code_min_background_level,
            code_max_background_level, periodic_min, periodic_max,
            domain_min_n_level, domain_max_n_level, ull_max, partition_interface_level,
            ptr_mpi_outer_layer->at(0), vec_ptr_grid_info_.at(1)->map_grid_node_, sfbitset_aux,
            &ptr_mpi_inner_layer->at(0), &vec_ptr_grid_info_.at(0)->map_grid_node_);
    }

    for (DefAmrIndexUint i_level = 1; i_level <= k0MaxLevel_; ++i_level) {
        // setup flags for refinement interface
        MarkRefinementInterface(i_level, sfbitset_one_lower_level.at(i_level), sfbitset_extra_lower_level.at(i_level));
        // add nodes for interpolation
        DefAmrIndexUint maxlayer = vec_ptr_grid_info_.at(i_level)->k0NumFine2CoarseLayer_;
        for (auto& iter_interface : vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_) {
            for (DefAmrIndexUint i_layer = maxlayer - vec_ptr_grid_info_.at(i_level)->k0NumFine2CoarseGhostLayer_;
                i_layer < maxlayer; ++i_layer) {
                vec_ptr_grid_info_.at(i_level)->AddGhostNodesForInterpolation(periodic_min, periodic_max,
                    sfbitset_aux, iter_interface.second->vec_inner_fine2coarse_.at(i_layer),
                    vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_);
                vec_ptr_grid_info_.at(i_level)->AddGhostNodesForInterpolation(periodic_min, periodic_max,
                    sfbitset_aux, iter_interface.second->vec_outer_fine2coarse_.at(i_layer),
                    vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_);
            }
        }
    }
}
/**
 * @brief function to setup flags for refinement interface.
 * @param[in] i_level refinement level.
 * @param[in] sfbitset_one_lower_level_current space filling codes at one lower refinement level for current refinement level
 * @param[in] sfbitset_extra space filling codes in extended mpi layer.
 */
int GridManagerInterface::MarkRefinementInterface(const DefAmrIndexUint i_level,
    const DefMap<DefAmrIndexUint>& sfbitset_one_lower_level_current,
    const DefMap<DefAmrIndexUint>& sfbitset_extra) {
    for (auto& iter_interface : vec_ptr_grid_info_.at(i_level - 1)->map_ptr_interface_layer_info_) {
        DefAmrIndexUint layer_fine = vec_ptr_grid_info_.at(i_level)->k0NumFine2CoarseLayer_ - 1;
        std::vector<DefSFBitset> corner_bitsets;
        // find overlap nodes based on the outmost fine to coarse refinement interface
        int layer0 = 0;
        layer_fine = vec_ptr_grid_info_.at(i_level)->k0NumFine2CoarseLayer_ - 1;
        DefMap<DefAmrIndexUint> sfbitset_one_lower_level_tmp(sfbitset_one_lower_level_current);
        sfbitset_one_lower_level_tmp.insert(sfbitset_extra.begin(), sfbitset_extra.end());

        for (const auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(layer0)) {
            sfbitset_one_lower_level_tmp.insert({iter_node.first, kFlag0_});
        }
        for (const auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(layer0)) {
            FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
            for (const auto& iter_conner : corner_bitsets) {
                IdentifyInnermostInterfaceForACell(iter_conner,
                    iter_interface.second->vec_outer_coarse2fine_.at(layer0),
                    vec_ptr_grid_info_.at(i_level)->map_grid_node_, sfbitset_one_lower_level_tmp,
                    &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .at(iter_interface.first)->vec_outer_fine2coarse_.at(layer_fine - 2),
                    &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .at(iter_interface.first)->vec_outer_fine2coarse_.at(layer_fine - 1),
                    &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .at(iter_interface.first)->vec_outer_fine2coarse_.at(layer_fine),
                    &iter_interface.second->vec_outer_coarse2fine_.at(layer0 + 1));
            }
        }
        for (const auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(layer0)) {
            sfbitset_one_lower_level_tmp.insert({iter_node.first, kFlag0_});
        }
        for (const auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(layer0)) {
            FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
            for (auto& iter_conner : corner_bitsets) {
                IdentifyInnermostInterfaceForACell(iter_conner,
                    iter_interface.second->vec_inner_coarse2fine_.at(layer0),
                    vec_ptr_grid_info_.at(i_level)->map_grid_node_, sfbitset_one_lower_level_tmp,
                    &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .at(iter_interface.first)->vec_inner_fine2coarse_.at(layer_fine - 2),
                    &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .at(iter_interface.first)->vec_inner_fine2coarse_.at(layer_fine - 1),
                    &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .at(iter_interface.first)->vec_inner_fine2coarse_.at(layer_fine),
                    &iter_interface.second->vec_inner_coarse2fine_.at(layer0 + 1));
            }
        }

        // find overlap nodes based on layers other than the outmost fine to coarse refinement interface
        DefAmrUint flag_tmp = NodeBitStatus::kNodeStatus0_;
        DefAmrIndexUint maxlayer = vec_ptr_grid_info_.at(i_level - 1)->k0NumCoarse2FineLayer_ - 1;
        for (DefAmrIndexUint i_layer = 1;
            i_layer < vec_ptr_grid_info_.at(i_level - 1)->k0NumCoarse2FineLayer_ - 1; ++i_layer) {
            layer_fine = (vec_ptr_grid_info_.at(i_level - 1)->k0NumCoarse2FineLayer_ - i_layer - 1)*2;
            flag_tmp = NodeBitStatus::kNodeStatus0_;
            if (i_layer >= maxlayer - vec_ptr_grid_info_.at(i_level - 1)->k0NumCoarse2FineGhostLayer_) {
                flag_tmp |= NodeBitStatus::kNodeStatusCoarse2FineGhost_;
            }
            for (auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(i_layer)) {
                if (vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.find(iter_node.first)
                    != vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.end()) {
                    vec_ptr_grid_info_.at(i_level-1)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
                } else {
                    LogManager::LogError("node on coarse to fine interface does not exist in map_grid_node_");
                }
                FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
                for (auto& iter_conner : corner_bitsets) {
                    IdentifyInterfaceForACell(iter_conner,
                        iter_interface.second->vec_outer_coarse2fine_.at(i_layer - 1),
                        iter_interface.second->vec_outer_coarse2fine_.at(i_layer),
                        vec_ptr_grid_info_.at(i_level)->map_grid_node_,
                        vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_,
                        &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                        .at(iter_interface.first)->vec_outer_fine2coarse_.at(layer_fine),
                        &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                        .at(iter_interface.first)->vec_outer_fine2coarse_.at(layer_fine + 1),
                        &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                        .at(iter_interface.first)->vec_outer_fine2coarse_.at(layer_fine + 2),
                        &iter_interface.second->vec_outer_coarse2fine_.at(i_layer + 1));
                }
            }
            for (auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(i_layer)) {
                if (vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.find(iter_node.first)
                    != vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.end()) {
                    vec_ptr_grid_info_.at(i_level-1)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
                } else {
                    LogManager::LogError("node on coarse to fine interface does not exist in map_grid_node_");
                }
                FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
                for (auto& iter_conner : corner_bitsets) {
                    IdentifyInterfaceForACell(iter_conner,
                        iter_interface.second->vec_outer_coarse2fine_.at(i_layer - 1),
                        iter_interface.second->vec_outer_coarse2fine_.at(i_layer),
                        vec_ptr_grid_info_.at(i_level)->map_grid_node_,
                        vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_,
                        &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                        .at(iter_interface.first)->vec_inner_fine2coarse_.at(layer_fine),
                        &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                        .at(iter_interface.first)->vec_inner_fine2coarse_.at(layer_fine + 1),
                        &vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                        .at(iter_interface.first)->vec_inner_fine2coarse_.at(layer_fine + 2),
                        &iter_interface.second->vec_inner_coarse2fine_.at(i_layer + 1));
                }
            }
        }
        // set flags for refinement interface
        flag_tmp = NodeBitStatus::kNodeStatusCoarse2Fine0_|NodeBitStatus::kNodeStatusCoarse2FineGhost_;
        for (auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(maxlayer)) {
            if (vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.find(iter_node.first)
                != vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.end()) {
                vec_ptr_grid_info_.at(i_level-1)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
            } else {
                LogManager::LogError("node on coarse to fine interface does not exist in map_grid_node_");
            }
        }
        for (auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(maxlayer)) {
            if (vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.find(iter_node.first)
                != vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_.end()) {
                vec_ptr_grid_info_.at(i_level-1)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
            } else {
                LogManager::LogError("node on coarse to fine interface does not exist in map_grid_node_");
            }
        }
        maxlayer = vec_ptr_grid_info_.at(i_level)->k0NumFine2CoarseLayer_;
        flag_tmp = NodeBitStatus::kNodeStatusFine2CoarseGhost_;
        for (DefAmrIndexUint i_layer = maxlayer - vec_ptr_grid_info_.at(i_level)->k0NumFine2CoarseGhostLayer_;
            i_layer < maxlayer - 1; ++i_layer) {
            for (auto& iter_node : vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                .at(iter_interface.first)->vec_outer_fine2coarse_.at(i_layer)) {
                if (vec_ptr_grid_info_.at(i_level)->map_grid_node_.find(iter_node.first)
                    != vec_ptr_grid_info_.at(i_level)->map_grid_node_.end()) {
                    vec_ptr_grid_info_.at(i_level)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
                } else {
                    LogManager::LogError("node on fine to coarse interface does not exist in map_grid_node_");
                }
            }
            for (auto& iter_node : vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                .at(iter_interface.first)->vec_inner_fine2coarse_.at(i_layer)) {
                if (vec_ptr_grid_info_.at(i_level)->map_grid_node_.find(iter_node.first)
                    != vec_ptr_grid_info_.at(i_level)->map_grid_node_.end()) {
                    vec_ptr_grid_info_.at(i_level)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
                } else {
                    LogManager::LogError("node on fine to coarse interface does not exist in map_grid_node_");
                }
            }
        }
        flag_tmp = NodeBitStatus::kNodeStatusFine2Coarse0_|NodeBitStatus::kNodeStatusFine2CoarseGhost_;
        for (auto& iter_node : vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
            .at(iter_interface.first)->vec_outer_fine2coarse_.at(maxlayer - 1)) {
            if (vec_ptr_grid_info_.at(i_level)->map_grid_node_.find(iter_node.first)
                != vec_ptr_grid_info_.at(i_level)->map_grid_node_.end()) {
                vec_ptr_grid_info_.at(i_level)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
            } else {
                LogManager::LogError("node on fine to coarse interface does not exist in map_grid_node_");
            }
        }
        for (auto& iter_node : vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
            .at(iter_interface.first)->vec_inner_fine2coarse_.at(maxlayer - 1)) {
            if (vec_ptr_grid_info_.at(i_level)->map_grid_node_.find(iter_node.first)
                != vec_ptr_grid_info_.at(i_level)->map_grid_node_.end()) {
                vec_ptr_grid_info_.at(i_level)->map_grid_node_.at(iter_node.first)->flag_status_ |= flag_tmp;
            } else {
                LogManager::LogError("node on fine to coarse interface does not exist in map_grid_node_");
            }
        }
    }
    return 0;
}
/**
 * @brief function to instantiate nodes on periodic boundary.
 * @param[in] bool_min if true, it is minimum domain boundary.
 * @param[in] i_dim corresponding dimension of the boundary.
 * @param[in] i_level refinement level.
 * @param[in] boundary_min the minium value of the indicator for current boundary.
 * @param[in] sfbitset_in space filling codes at current refinement level.
 * @param[in] code_min_background_level the minimum space filling codes at background level.
 * @param[in] sfbitset_aux  class manage space filling curves.
 * @param[out] ptr_grid_info pointer to class store grid information.
 * @param[out] ptr_inner_layer pointer to nodes on the inner layer for mpi communication.
 * @param[out] ptr_outer_layer pointer to nodes on the outer layer for mpi communication.
 * @param[out] ptr_boundary_corner_y_min pointer to nodes at the conner of y min boundary.
 * @param[out] ptr_boundary_corner_y_max pointer to nodes at the conner of y max boundary.
 * @param[out] ptr_boundary_corner_z_min pointer to nodes at the conner of z min boundary.
 * @param[out] ptr_boundary_corner_z_max pointer to nodes at the conner of z max boundary. 
 */
void GridManagerInterface::InstantiatePeriodicNodes(const bool bool_min,
    const DefAmrIndexUint i_dim, const DefAmrIndexUint i_level,
    const DefAmrIndexUint num_ghost_layer, const int boundary_min, const DefSFBitset sfbitset_in,
    const DefSFCodeToUint code_min_background_level, const DefSFCodeToUint code_max_background_level,
    const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info,
    DefMap<DefAmrIndexUint>* const ptr_inner_layer, DefMap<DefAmrIndexUint>* const ptr_outer_layer,
    DefMap<DefAmrIndexUint>* const ptr_boundary_corner_y_min,
    DefMap<DefAmrIndexUint>* const ptr_boundary_corner_y_max,
    DefMap<DefAmrIndexUint>* const ptr_boundary_corner_z_min,
    DefMap<DefAmrIndexUint>* const ptr_boundary_corner_z_max) {
    DefSFBitset sfbitset_tmp = sfbitset_in;
    DefMap<std::unique_ptr<GridNode>>& map_grid = ptr_grid_info->map_grid_node_;
    int boundary_flag;
    DefSFCodeToUint code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, sfbitset_tmp).to_ullong();
    if (code_tmp >= code_min_background_level && code_tmp <= code_max_background_level) {
        ptr_inner_layer->insert({sfbitset_tmp, kFlag0_});
    }
    DefSFBitset sfbitset_tmp_pre = sfbitset_tmp;
    for (DefAmrIndexUint i = 0; i < num_ghost_layer; ++i) {
        if (bool_min) {
            sfbitset_tmp = sfbitset_aux.FindNodeInPos(i_dim, sfbitset_tmp);
        } else {
            sfbitset_tmp = sfbitset_aux.FindNodeInNeg(i_dim, sfbitset_tmp);
        }
        if (map_grid.find(sfbitset_tmp) == map_grid.end()) {
            InstantiateGridNode(sfbitset_tmp, ptr_grid_info);
        } else if ((map_grid.at(sfbitset_tmp)->flag_status_&NodeBitStatus::kNodeStatusFine2Coarse0_) &&
            ((map_grid.at(sfbitset_tmp_pre)->flag_status_&NodeBitStatus::kNodeStatusFine2Coarse0_) == 0)) {
            // break when the node status changes
            break;
        }

        code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, sfbitset_tmp).to_ullong();
        if (code_tmp < code_min_background_level || code_tmp > code_max_background_level) {
            ptr_outer_layer->insert({sfbitset_tmp, kFlag0_});
            map_grid.at(sfbitset_tmp)->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
            boundary_flag = CheckNodeOnDomainBoundary(i_level, sfbitset_tmp);
            if (boundary_flag > boundary_min) {
                // find edge or corner nodes which should on periodic boundary
                // to avoid missing nodes which are not on this mpi rank
                if (boundary_flag&4) {
                    ptr_boundary_corner_y_min->insert({sfbitset_tmp, kFlag0_});
                }
                if (boundary_flag&8) {
                    ptr_boundary_corner_y_max->insert({sfbitset_tmp, kFlag0_});
                }
                if (boundary_flag&16) {
                    ptr_boundary_corner_z_min->insert({sfbitset_tmp, kFlag0_});
                }
                if (boundary_flag&32) {
                    ptr_boundary_corner_z_max->insert({sfbitset_tmp, kFlag0_});
                }
            }
        } else {
            ptr_inner_layer->insert({sfbitset_tmp, kFlag0_});
        }
        sfbitset_tmp_pre = sfbitset_tmp;
    }
}
/**
 * @brief function to instantiate nodes on periodic boundary and outer mpi layers.
 * @param[in] num_ghost_layer number of mpi communication layers extended interior from the domain boundary.
 * @param[in] code_min_background_level minimum space filling code at background level for current rank.
 * @param[in] code_max_background_level maximum space filling code at background level for current rank.
 * @param[in] periodic_min booleans indicating if the boundary is periodic at minimum domain boundaries.
 * @param[in] periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] sfbitset_one_lower_level  space filling codes at one lower refinement level.
 * @param[in] sfbitset_aux  class manage space filling curves.
 * @param[out] ptr_grid_info pointer to class store grid information.
 * @param[out] ptr_inner_layer pointer to nodes on the inner layer for mpi communication.
 * @param[out] ptr_outer_layer pointer to nodes on the outer layer for mpi communication.
 * @note num_ghost_layer only valid for boundaries where periodic_min or periodic_max is true
 */
void  GridManagerInterface::InstantiateDomainBoundaryForMpi(const DefAmrIndexLUint num_ghost_layer,
    const DefSFCodeToUint code_min_background_level, const DefSFCodeToUint code_max_background_level,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const DefMap<DefAmrIndexUint>& sfbitset_one_lower_level,
    const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info,
    DefMap<DefAmrIndexUint>* const ptr_inner_layer, DefMap<DefAmrIndexUint>* const ptr_outer_layer) {
    DefAmrIndexUint i_level = ptr_grid_info->i_level_;
    DefSFBitset sfbitset_tmp, sfbitset_tmp_pre;
    DefMap<std::unique_ptr<GridNode>>& map_grid = ptr_grid_info->map_grid_node_;
    std::vector<DefSFBitset> nodes_on_surface;
    DefSFCodeToUint code_tmp;
    DefMap<DefAmrIndexUint> boundary_corner_y_min, boundary_corner_y_max,
        boundary_corner_z_min, boundary_corner_z_max;
    for (DefAmrIndexUint i_dim = 0; i_dim < k0GridDims_; ++i_dim) {
        DefMap<DefAmrIndexUint> boundary_min_tmp;
        int boundary_min = (1 << (i_dim + 2)) - 1;
        for (const auto& iter_node : ptr_grid_info->domain_boundary_min_.at(i_dim)) {
            NodesBelongToOneSurfAtHigherLevel(sfbitset_aux.SFBitsetToNLowerLevelVir(
                1, iter_node.first), i_dim, sfbitset_one_lower_level, &nodes_on_surface);
            for (const auto& iter_surface : nodes_on_surface) {
                boundary_min_tmp.insert({iter_surface, kFlag0_});
            }
        }
        for (const auto& iter_node : boundary_min_tmp) {
            ptr_grid_info->domain_boundary_min_.at(i_dim).insert({iter_node.first, kFlag0_});
            if (map_grid.find(iter_node.first) == map_grid.end()) {
                InstantiateGridNode(iter_node.first, ptr_grid_info);
            }
            code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first).to_ullong();
            if (code_tmp < code_min_background_level || code_tmp > code_max_background_level) {
                ptr_outer_layer->insert({iter_node.first, kFlag0_});
                map_grid.at(iter_node.first)->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
            }

            if (periodic_min.at(i_dim)) {
                InstantiatePeriodicNodes(true, i_dim, i_level, num_ghost_layer, boundary_min, iter_node.first,
                    code_min_background_level, code_max_background_level, sfbitset_aux, ptr_grid_info,
                    ptr_inner_layer, ptr_outer_layer, &boundary_corner_y_min, &boundary_corner_y_max,
                    &boundary_corner_z_min, &boundary_corner_z_max);
            }
        }

        DefMap<DefAmrIndexUint> boundary_max_tmp;
        for (const auto& iter_node : ptr_grid_info->domain_boundary_max_.at(i_dim)) {
            NodesBelongToOneSurfAtHigherLevel(sfbitset_aux.SFBitsetToNLowerLevelVir(
                1, iter_node.first), i_dim, sfbitset_one_lower_level, &nodes_on_surface);
            for (const auto& iter_surface : nodes_on_surface) {
                boundary_max_tmp.insert({iter_surface, kFlag0_});
            }
        }
        for (const auto& iter_node : boundary_max_tmp) {
            ptr_grid_info->domain_boundary_max_.at(i_dim).insert({iter_node.first, kFlag0_});
            if (map_grid.find(iter_node.first) == map_grid.end()) {
                InstantiateGridNode(iter_node.first, ptr_grid_info);
            }
            code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node.first).to_ullong();
            if (code_tmp < code_min_background_level || code_tmp > code_max_background_level) {
                ptr_outer_layer->insert({iter_node.first, kFlag0_});
                map_grid.at(iter_node.first)->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOuter_;
            }
            if (periodic_max.at(i_dim)) {
                InstantiatePeriodicNodes(false, i_dim, i_level, num_ghost_layer, boundary_min, iter_node.first,
                    code_min_background_level, code_max_background_level, sfbitset_aux, ptr_grid_info,
                    ptr_inner_layer, ptr_outer_layer, &boundary_corner_y_min, &boundary_corner_y_max,
                    &boundary_corner_z_min, &boundary_corner_z_max);
            }
        }
    }

    if (periodic_min.at(kYIndex)) {
        for (const auto& iter_node : boundary_corner_y_min) {
            InstantiatePeriodicNodes(true, kYIndex, i_level, num_ghost_layer, 15, iter_node.first,
                code_min_background_level, code_max_background_level, sfbitset_aux, ptr_grid_info,
                ptr_inner_layer, ptr_outer_layer, &boundary_corner_y_max, &boundary_corner_y_max,
                &boundary_corner_z_min, &boundary_corner_z_max);
        }
    }
    if (periodic_max.at(kYIndex)) {
        for (const auto& iter_node : boundary_corner_y_max) {
            InstantiatePeriodicNodes(false, kYIndex, i_level, num_ghost_layer, 15, iter_node.first,
                code_min_background_level, code_max_background_level, sfbitset_aux, ptr_grid_info,
                ptr_inner_layer, ptr_outer_layer, &boundary_corner_y_min, &boundary_corner_y_min,
                &boundary_corner_z_min, &boundary_corner_z_max);
        }
    }
    if (k0GridDims_ == 3) {
        if (periodic_min.at(kZIndex)) {
            for (const auto& iter_node : boundary_corner_z_min) {
                InstantiatePeriodicNodes(true, kZIndex, i_level, num_ghost_layer, 63, iter_node.first,
                    code_min_background_level, code_max_background_level, sfbitset_aux, ptr_grid_info,
                    ptr_inner_layer, ptr_outer_layer, &boundary_corner_y_min, &boundary_corner_y_max,
                    &boundary_corner_z_max, &boundary_corner_z_max);
            }
        }
        if (periodic_max.at(kZIndex)) {
            for (const auto& iter_node : boundary_corner_z_max) {
                InstantiatePeriodicNodes(false, kZIndex, i_level, num_ghost_layer, 63, iter_node.first,
                    code_min_background_level, code_max_background_level, sfbitset_aux, ptr_grid_info,
                    ptr_inner_layer, ptr_outer_layer, &boundary_corner_y_min, &boundary_corner_y_max,
                    &boundary_corner_z_min, &boundary_corner_z_min);
            }
        }
    }
}
/**
 * @brief function to set indices for nodes on inner communication layers.
 * @param[in] i_rank current rank id
 * @param[in] interface_code_background_level space filling code of node in outer mpi layer at background level.
 * @param[in] vec_nodes_in_region nodes in a region of given length.
 * @param[in] ull_max  the maximum space filling code for each partition rank.
 * @param[in] mpi_inner_layer_tmp container temporally storing nodes are possible in inner layer.
 * @param[in] ptr_mpi_inner_layer pointer to nodes on the inner layer for mpi communication (sending).
 */
void GridManagerInterface::SetUpRanksSentForMpiInnerLayer(int i_rank,
    const DefSFCodeToUint interface_code_background_level, const std::vector<DefSFBitset>& vec_nodes_in_region,
    const std::vector<DefSFCodeToUint>& ull_max, const DefMap<DefAmrIndexUint>& mpi_inner_layer_tmp,
    std::map<int, DefMap<DefAmrIndexUint>>* const ptr_mpi_inner_layer) {
    int node_rank;
    std::vector<DefSFCodeToUint>::const_iterator iter_index;
    for (const auto& iter_node : vec_nodes_in_region) {
        if ((iter_node != SFBitsetAuxInterface::kInvalidSFbitset)
            &&(mpi_inner_layer_tmp.find(iter_node) != mpi_inner_layer_tmp.end())) {
            iter_index = std::lower_bound(ull_max.cbegin(), ull_max.cend(), interface_code_background_level);
            node_rank = static_cast<int>(iter_index - ull_max.cbegin());
            if (node_rank != i_rank) {
                if (ptr_mpi_inner_layer->find(node_rank)
                    == ptr_mpi_inner_layer->end()) {
                    ptr_mpi_inner_layer->insert({node_rank, {}});
                }
                if (ptr_mpi_inner_layer->at(node_rank).find(iter_node)
                    == ptr_mpi_inner_layer->at(node_rank).end()) {
                    ptr_mpi_inner_layer->at(node_rank).insert({iter_node, 0});
                }
            }
        }
    }
}
/**
 * @brief function to find nodes in inner mpi communication layers where their counterparts
 *    in outer mpi communication layers are on grid of one higher refinement level .
 * @param[in] i_rank rank of the current process.
 * @param[in] i_level refinement level of target grid.
 * @param[in] num_partition_outer_layer number of outer mpi communication layers.
 * @param[in] code_min_background_level minimum space filling code at background level for current rank.
 * @param[in] code_max_background_level maximum space filling code at background level for current rank.
 * @param[in] periodic_min booleans indicating if the boundary is periodic at minimum domain boundaries.
 * @param[in] periodic_max booleans indicating if the boundary is periodic at maximum domain boundaries.
 * @param[in] domain_min_n_level space filling codes representing the minimum domain at each level.
 * @param[in] domain_max_n_level space filling codes representing the maximum domain at each level.
 * @param[in] ull_max the maximum space filling code for each partition rank.
 * @param[in] partition_interface_level nodes on partition interface at current refinement level.
 * @param[in] mpi_outer_layer  nodes on the outer layer for mpi communication (receiving).
 * @param[in] map_grid_node_higher_level grid nodes at one higher refinement level .
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[out] ptr_mpi_inner_layer pointer to nodes on the inner layer for mpi communication (sending).
 * @param[out] ptr_map_grid_node the pointer to target grid nodes at given refinement level.
 */
void GridManagerInterface::SetUpMpiInnerLayerForHigherLevel(int i_rank,
    const DefAmrIndexUint i_level, const DefAmrIndexLUint num_partition_outer_layer,
    const DefSFCodeToUint code_min_background_level, const DefSFCodeToUint code_max_background_level,
    const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
    const std::vector<DefSFBitset>& domain_min_n_level, const std::vector<DefSFBitset>& domain_max_n_level,
    const std::vector<DefSFCodeToUint>& ull_max, const DefMap<DefAmrIndexUint>& partition_interface_level,
    const DefMap<DefAmrIndexUint>& mpi_outer_layer, const DefMap<std::unique_ptr<GridNode>>& map_grid_node_higher_level,
    const SFBitsetAuxInterface& sfbitset_aux, std::map<int, DefMap<DefAmrIndexUint>>* const ptr_mpi_inner_layer,
    DefMap<std::unique_ptr<GridNode>> *const ptr_map_grid_node) {
    // find nodes on outer communication layer but not in ptr_mpi_outer_layer,
    // such as those in layers across grids of different refinement levels
    DefSFCodeToUint code_background;
    std::vector<DefSFBitset> vec_outer, vec_inner;
    DefSFBitset sfbitset_higher_level;
    std::vector<DefSFCodeToUint>::const_iterator iter_index;
    int node_rank;
    for (const auto& iter_interface : partition_interface_level) {
        sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_interface.first, num_partition_outer_layer,
            periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_outer);
        for (const auto& iter_outer : vec_outer) {
            if (iter_outer != SFBitsetAuxInterface::kInvalidSFbitset) {
                code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_outer).to_ullong();
                if (code_background < code_min_background_level || code_background > code_max_background_level) {
                    if (mpi_outer_layer.find(iter_outer) == mpi_outer_layer.end()) {
                        sfbitset_higher_level = sfbitset_aux.SFBitsetToNHigherLevelVir(1, iter_outer);
                        // nodes in mpi outer layers at higher level exist indicating that
                        // informaiton of inner nodes need to be sent to other ranks even if
                        // there is no corresponding outer nodes at the current refinement level
                        if (map_grid_node_higher_level.find(sfbitset_higher_level)
                            != map_grid_node_higher_level.end()) {
                            iter_index = std::lower_bound(ull_max.cbegin(), ull_max.cend(), code_background);
                            node_rank = static_cast<int>(iter_index - ull_max.begin());
                            sfbitset_aux.FindNodesInPeriodicRegionCenter(iter_outer, num_partition_outer_layer,
                                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_inner);
                            for (const auto& iter_inner : vec_inner) {
                                code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                                    i_level, iter_inner).to_ullong();
                                if (code_background >= code_min_background_level
                                    && code_background <= code_max_background_level) {
                                    if (ptr_map_grid_node->find(iter_inner)
                                        != ptr_map_grid_node->end()) {
                                        ptr_map_grid_node->at(iter_inner)->flag_status_ |=
                                            NodeBitStatus::kNodeStatusMpiPartitionInner_;
                                        if (ptr_mpi_inner_layer->find(node_rank)
                                            == ptr_mpi_inner_layer->end()) {
                                            ptr_mpi_inner_layer->insert({node_rank, {}});
                                        }
                                        if (ptr_mpi_inner_layer->at(node_rank).find(iter_inner)
                                            == ptr_mpi_inner_layer->at(node_rank).end()) {
                                            ptr_mpi_inner_layer->at(node_rank).insert({iter_inner, 0});
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
