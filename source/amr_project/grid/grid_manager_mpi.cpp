//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager_mpi.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes when mpi is used.
* @date  2022-6-7
*/
#include <string>
#include <set>
#include <mpi.h>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to instantiate grid noes on refinement interface and mpi communication layer.
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
    DefAmrUint flag_temp, flag_0 = 0;
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
    DefAmrIndexUint layer_coarse0 = grid_info_lower.k0NumCoarse2FineLayer_ - 1,
        layer_coarse_m1 = layer_coarse0 - 1,
        layer0 = grid_info.k0NumFine2CoarseLayer_ - 1,
        layer_m1 = layer0 - 1,
        layer_m2 = layer_m1 - 1;
    std::vector<DefSFBitset> sfbitset_neighbors;
    for (auto& iter_interface : grid_info_lower.map_ptr_interface_layer_info_) {
        ptr_interface_info_lower = iter_interface.second.get();
        if (grid_info.map_ptr_interface_layer_info_
            .find(iter_interface.first)
            == grid_info.map_ptr_interface_layer_info_.end()) {
            grid_info.map_ptr_interface_layer_info_.insert(
                { iter_interface.first,
                std::make_shared<InterfaceLayerInfo>() });
        }
        ptr_interface_info = grid_info.map_ptr_interface_layer_info_
            .at(iter_interface.first).get();
        // interface inside the geometry
        if (ptr_interface_info_lower->vec_inner_coarse2fine_.size() > 0) {
            ptr_interface_info->vec_inner_fine2coarse_.resize(
                grid_info.k0NumFine2CoarseLayer_);
            FindOverlappingLayersBasedOnOutermostCoarse(
                ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_m1),
                map_sfbitset_one_lower_level,
                &ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse0),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
            // insert node instance
            maxlayer = static_cast<int>(ptr_interface_info->vec_inner_fine2coarse_.size());
            for (int ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusFine2Coarse0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusFine2CoarseM1_;
                } else {
                    flag_temp = flag_0;
                }
                for (const auto& iter_layer_node : ptr_interface_info
                    ->vec_inner_fine2coarse_.at(ilayer)) {
                    if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                            map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        }
                    } else {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node is on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
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
            maxlayer = static_cast<int>(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusCoarse2Fine0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusCoarse2FineM1_;
                } else {
                    flag_temp = flag_0;
                }
                for (const auto& iter_layer_node : ptr_interface_info_lower
                    ->vec_inner_coarse2fine_.at(ilayer)) {
                    if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        }
                    } else {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_lower_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
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
                ptr_interface_info_lower->vec_outer_coarse2fine_.at(layer_coarse_m1),
                map_sfbitset_one_lower_level,
                &ptr_interface_info_lower->vec_outer_coarse2fine_.at(layer_coarse0),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer0),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m1),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m2));
            // insert node instance
            DefAmrIndexUint maxlayer = DefAmrIndexUint(ptr_interface_info->vec_outer_fine2coarse_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusFine2Coarse0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusFine2CoarseM1_;
                } else {
                    flag_temp = flag_0;
                }
                for (const auto& iter_layer_node : ptr_interface_info->vec_outer_fine2coarse_.at(ilayer)) {
                    if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                            map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        }
                    } else {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }

                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
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
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusCoarse2Fine0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_0 | NodeBitStatus::kNodeStatusCoarse2FineM1_;
                } else {
                    flag_temp = flag_0;
                }
                for (const auto& iter_layer_node : ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(ilayer)) {
                    if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                        if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        }
                    } else {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_lower_level, iter_layer_node.first).to_ullong();
                    if ((ilayer == maxlayer - 1) && (code_background < code_min || code_background > code_max)) {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
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
 * @param[out] ptr_mpi_outer_layer pointer to nodes on the outer layer for mpi communication (sending). 
 */
void GridManagerInterface::InstantiateGridNodeAllLevelMpi(const int i_rank,
    const DefAmrIndexUint num_partition_inner_layer, const DefAmrIndexUint num_partition_outer_layer,
    const std::vector<DefSFCodeToUint>& vec_sfcode_min, const std::vector<DefSFCodeToUint>& vec_sfcode_max,
    const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level,
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_ghost_one_lower_level,
    const DefMap<DefAmrIndexUint>& sfbitset_partition_interface_0,
    std::vector<std::map<int, DefMap<DefAmrIndexUint>>>* const ptr_mpi_inner_layer,
    std::vector<DefMap<DefAmrIndexUint>>* const ptr_mpi_outer_layer) {
    DefAmrIndexUint flag_normal_outer_node = ~0, flag_on_coarse2fine = ~0 - 1,
        flag_outmost_refinement_n_outer_mpi = flag_on_coarse2fine;
    DefMap<DefAmrIndexUint> background_occupied;  // background nodes coincide with those at higher levels
    ptr_mpi_inner_layer->resize(k0MaxLevel_ + 1);
    ptr_mpi_outer_layer->resize(k0MaxLevel_ + 1);
    DefSFCodeToUint code_min_background_level = vec_sfcode_min.at(i_rank),
        code_max_background_level = vec_sfcode_max.at(i_rank);
    DefSFCodeToUint code_tmp, code_neighbor_tmp;
    bool bool_partition_outside, bool_partition_inside;
    std::vector<DefSFBitset> vec_in_region;
    std::vector<DefSFCodeToUint>::iterator iter_index;
    int node_rank;
    // initialize grid information, will be used for instantiate grid nodes
    for (DefAmrIndexUint i_level = 0; i_level <= k0MaxLevel_; ++i_level) {
        vec_ptr_grid_info_.at(i_level)->InitialGridInfo();
    }
    std::vector<DefSFCodeToUint> ull_max(vec_sfcode_max);
    std::vector<DefAmrIndexLUint> indices_min = GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = GetMaxIndexOfBackgroundNodeArrAsVec();
    int flag_node;
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
                    |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
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
        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            if (CheckCoincideBackground(i_level_lower, iter_low.first, &bitset_background)) {
                if (background_occupied.find(bitset_background)
                 == background_occupied.end()) {
                    background_occupied.insert({ bitset_background, kFlagSize0_ });
                }
            }
            if (NodesBelongToOneCell(iter_low.first,
                sfbitset_one_lower_level.at(i_level), &bitset_cell_lower)) {
                bool_partition_outside = false;
                bool_partition_inside = false;
                FindAllNodesInACellAtOneLevelLower(bitset_cell_lower, &bitset_all);
                for (const auto& iter_node : bitset_all) {
                    if (map_grid.find(iter_node) == map_grid.end()) {
                        InstantiateGridNode(iter_node, &grid_info);
                    }
                    code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                    if (code_tmp < code_min_background_level || code_tmp > code_max_background_level) {
                        map_grid.at(iter_node)->flag_status_ |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
                        if (ptr_mpi_outer_layer->at(i_level).find(iter_node)
                            == ptr_mpi_outer_layer->at(i_level).end()) {
                            ptr_mpi_outer_layer->at(i_level).insert({iter_node, flag_normal_outer_node});
                        }
                        bool_partition_outside = true;
                    } else {
                        bool_partition_inside = true;
                        flag_node = grid_info.CheckIfNodeOutsideCubicDomain(k0GridDims_, iter_node, sfbitset_aux);
                        (this->*ptr_func_insert_domain_boundary_)(flag_node, iter_node, &grid_info);
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
        vec_ptr_grid_info_.at(i_level)->SetPeriodicBoundaryAsPartitionInterface(
            k0GridDims_, sfbitset_aux, periodic_min, periodic_max, &partition_interface_level);
        DefMap<DefAmrIndexUint> inner_layer_tmp;
        // Extend mpi interface nodes for a given number of inner communication layers.
        // Noting that nodes ont the partitioned interface is considered in the inner layer.
        // Use temporal container (inner_layer_tmp) at this stager since some of them inserted
        // here are on the periodic boundaries but not on the inner layer. In addition, some
        // inner nodes should have duplicate since they will be send to different ranks (overlap region).
        for (const auto& iter_interface : partition_interface_level) {
            sfbitset_aux.FindNodesInPeriodicReginNearby(iter_interface.first, num_partition_inner_layer-1,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                if ((code_tmp >= code_min_background_level && code_tmp <= code_max_background_level)
                    && inner_layer_tmp.find(iter_node) == inner_layer_tmp.end()
                    && map_grid.find(iter_node) != map_grid.end()) {
                    inner_layer_tmp.insert({iter_node, 0});
                }
            }
        }

        DefMap<DefAmrIndexUint> map_extra_expand_outer;
        bool bool_interface_upper_extra = ((num_partition_outer_layer%2) == 0);
        // extend mpi communication layers for periodic boundaries if necessary
        DefSFBitset sfbitset_tmp;
        for (DefAmrIndexUint i_dim = 0; i_dim < k0GridDims_; ++i_dim) {
            if (periodic_min.at(i_dim) && (bool_interface_upper_extra)) {
                for (const auto iter_interface : grid_info.domain_boundary_min_.at(i_dim)) {
                    sfbitset_aux.FindNodesInPeriodicReginOfGivenLength(iter_interface.first,
                        num_partition_outer_layer, periodic_min, periodic_max,
                        domain_min_n_level, domain_max_n_level, &vec_in_region);
                    for (const auto& iter_node : vec_in_region) {
                        code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                        if (code_tmp > code_max_background_level
                            &&ptr_mpi_outer_layer->at(i_level).find(iter_node)
                            != ptr_mpi_outer_layer->at(i_level).end()
                            &&ptr_mpi_outer_layer->at(i_level).at(iter_node) == flag_normal_outer_node) {
                            sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_node,
                                periodic_min, periodic_max, domain_min_n_level,
                                domain_max_n_level, &bitset_neighbors);
                            for (const auto& iter_neighbor : bitset_neighbors) {
                                if (map_extra_expand_outer.find(iter_neighbor) == map_extra_expand_outer.end()
                                    && map_grid.find(iter_neighbor) == map_grid.end()) {
                                    if (InstantiateGridNode(iter_neighbor, &grid_info)) {
                                        map_grid.at(iter_neighbor)->flag_status_
                                            |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
                                        map_extra_expand_outer.insert({iter_neighbor, 0});
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (periodic_max.at(i_dim) && (bool_interface_upper_extra)) {
                for (const auto iter_interface : grid_info.domain_boundary_max_.at(i_dim)) {
                    sfbitset_aux.FindNodesInPeriodicReginOfGivenLength(iter_interface.first,
                        num_partition_outer_layer, periodic_min, periodic_max,
                        domain_min_n_level, domain_max_n_level, &vec_in_region);
                    for (const auto& iter_node : vec_in_region) {
                        code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                        if (code_tmp > code_max_background_level
                            &&ptr_mpi_outer_layer->at(i_level).find(iter_node)
                            != ptr_mpi_outer_layer->at(i_level).end()
                            &&ptr_mpi_outer_layer->at(i_level).at(iter_node) == flag_normal_outer_node) {
                            sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_node,
                                periodic_min, periodic_max, domain_min_n_level,
                                    domain_max_n_level, &bitset_neighbors);
                            for (const auto& iter_neighbor : bitset_neighbors) {
                                if (map_extra_expand_outer.find(iter_neighbor) == map_extra_expand_outer.end()
                                    && map_grid.find(iter_neighbor) == map_grid.end()) {
                                    if (InstantiateGridNode(iter_neighbor, &grid_info)) {
                                        map_grid.at(iter_neighbor)->flag_status_ |=
                                            NodeBitStatus::kNodeStatusMpiPartitionOutside_;
                                        map_extra_expand_outer.insert({iter_neighbor, 0});
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        DefSFCodeToUint code_in_region, code_interface;
        std::vector<DefSFBitset> vec_periodic;
        bool bool_has_periodic_boundary =
            grid_info.CheckIfPeriodicDomainRequired(k0GridDims_, &periodic_min, &periodic_max);
        for (const auto& iter_interface : ptr_mpi_outer_layer->at(i_level)) {
            // set ranks that should be sent by mpi inner layer nodes
            if (bool_has_periodic_boundary) {
                // counterparts of nodes on periodic boundary wil be used here as what have been done
                // in grid generation since only nodes can construct a cell at one lower level will be instantiated
                grid_info.CheckNodesOnCubicPeriodicBoundary(k0GridDims_, iter_interface.first,
                    periodic_min, periodic_max, sfbitset_aux, &vec_periodic);
                for (const auto& iter_periodic : vec_periodic) {
                    sfbitset_aux.FindNodesInPeriodicReginNearby(iter_periodic, num_partition_inner_layer,
                        periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
                    for (const auto& iter_node : vec_in_region) {
                        code_in_region = sfbitset_aux.SFBitsetToNLowerLevelVir(
                                i_level, iter_node).to_ullong();
                        if (code_in_region >= code_min_background_level
                            &&code_in_region <= code_max_background_level
                            &&map_grid.find(iter_node) != map_grid.end()) {
                            iter_index = std::lower_bound(
                                ull_max.begin(), ull_max.end(), code_interface);
                            node_rank = static_cast<int>(iter_index - ull_max.begin());
                            if (node_rank != i_rank) {
                                if (ptr_mpi_inner_layer->at(i_level).find(node_rank)
                                    == ptr_mpi_inner_layer->at(i_level).end()) {
                                    ptr_mpi_inner_layer->at(i_level).insert({node_rank, {}});
                                }
                                if (ptr_mpi_inner_layer->at(i_level).at(node_rank).find(iter_node)
                                    == ptr_mpi_inner_layer->at(i_level).at(node_rank).end()) {
                                    ptr_mpi_inner_layer->at(i_level).at(node_rank).insert({iter_node, 0});
                                }
                            }
                        }
                    }
                }
            }
            sfbitset_aux.FindNodesInPeriodicReginNearby(iter_interface.first, num_partition_inner_layer,
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
                    if (map_extra_expand_outer.find(iter_neighbor) == map_extra_expand_outer.end()
                        && map_grid.find(iter_neighbor) == map_grid.end()) {
                        if (InstantiateGridNode(iter_neighbor, &grid_info)) {
                            map_grid.at(iter_neighbor)->flag_status_
                                |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
                            map_extra_expand_outer.insert({iter_neighbor, 0});
                        }

                    }
                }
            } else if ((!bool_interface_upper_extra && code_interface < code_min_background_level)
                && iter_interface.second == flag_normal_outer_node) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                    periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                for (const auto& iter_neighbor : bitset_neighbors) {
                    if (map_extra_expand_outer.find(iter_neighbor) == map_extra_expand_outer.end()
                        && map_grid.find(iter_neighbor) == map_grid.end()) {
                        if (InstantiateGridNode(iter_neighbor, &grid_info)) {
                            map_grid.at(iter_neighbor)->flag_status_
                                |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
                            map_extra_expand_outer.insert({iter_neighbor, 0});
                        }
                    }
                }
            }
        }
        // add the outmost layer to mpi outer communication layer
        for (const auto& iter_interface : map_extra_expand_outer) {
            ptr_mpi_outer_layer->at(i_level).insert({iter_interface.first, flag_normal_outer_node});
            sfbitset_aux.FindNodesInPeriodicReginNearby(iter_interface.first, num_partition_inner_layer,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            code_interface = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_interface.first).to_ullong();
            SetUpRanksSentForMpiInnerLayer(i_rank, code_interface, vec_in_region,
                ull_max, inner_layer_tmp, &ptr_mpi_inner_layer->at(i_level));
        }
    }
    // find overlapping node for refinement levels of 0 and 1
    for (const auto& iter_interfaces : vec_ptr_grid_info_.at(0)->map_ptr_interface_layer_info_) {
        for (const auto& iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
            for (const auto& iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
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
    vec_ptr_grid_info_.at(0)->SetPeriodicBoundaryAsPartitionInterface(
        k0GridDims_, sfbitset_aux, periodic_min, periodic_max, &partition_interface_level);
    for (const auto& iter_interface : partition_interface_level) {
        if (background_occupied.find(iter_interface.first) == background_occupied.end()) {
            map_outmost_current.insert({iter_interface.first, 0});
        }
    }
    DefMap<DefAmrIndexUint> inner_layer_tmp;
    std::vector<DefSFBitset> domain_min_n_level(k0GridDims_), domain_max_n_level(k0GridDims_);
    sfbitset_aux.GetMinAtGivenLevel(0, indices_min, &domain_min_n_level);
    sfbitset_aux.GetMaxAtGivenLevel(0, indices_max, &domain_max_n_level);
    // find nodes on inner communication layers at level 0
    for (const auto& iter_interface : map_outmost_current) {
        if (map_grid_level_0.find(iter_interface.first) != map_grid_level_0.end()) {
            sfbitset_aux.FindNodesInPeriodicReginNearby(iter_interface.first, num_partition_inner_layer-1,
                periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = iter_node.to_ullong();
                if ((code_tmp >= code_min_background_level && code_tmp <= code_max_background_level)
                    && inner_layer_tmp.find(iter_node) == inner_layer_tmp.end()
                    && background_occupied.find(iter_node) == background_occupied.end()) {
                    inner_layer_tmp.insert({iter_node, 0});
                }
            }
        }
    }
    // find nodes on outer communication layers layer by layer
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
                                        |= NodeBitStatus::kNodeStatusMpiPartitionOutside_;
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

    // set ranks that should be sent by mpi inner layer nodes (level 0)
    for (const auto& iter_interface : ptr_mpi_outer_layer->at(0)) {
        sfbitset_aux.FindNodesInPeriodicReginNearby(iter_interface.first, num_partition_inner_layer,
            periodic_min, periodic_max, domain_min_n_level, domain_max_n_level, &vec_in_region);
        for (const auto& iter_node : vec_in_region) {
            code_tmp = iter_node.to_ullong();
            if ((code_tmp >= code_min_background_level && code_tmp <= code_max_background_level)
                && inner_layer_tmp.find(iter_node) != inner_layer_tmp.end()) {
                iter_index = std::lower_bound(ull_max.begin(), ull_max.end(), iter_interface.first.to_ullong());
                node_rank = static_cast<int>(iter_index - ull_max.begin());
                if (node_rank != i_rank) {
                    if (ptr_mpi_inner_layer->at(0).find(node_rank)
                        == ptr_mpi_inner_layer->at(0).end()) {
                        ptr_mpi_inner_layer->at(0).insert({node_rank, {}});
                    }
                    if (ptr_mpi_inner_layer->at(0).at(node_rank).find(iter_node)
                        == ptr_mpi_inner_layer->at(0).at(node_rank).end()) {
                        ptr_mpi_inner_layer->at(0).at(node_rank).insert({iter_node, 0});
                    }
                }
            }
        }
    }
}
/**
 * @brief function to set indices for nodes on inner communication layers.
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
 * @param[out] ptr_mpi_outer_layer pointer to nodes on the outer layer for mpi communication (sending). 
 */
void GridManagerInterface::SetUpRanksSentForMpiInnerLayer(int i_rank,
    const DefSFCodeToUint interface_code_background_level, const std::vector<DefSFBitset>& vec_nodes_in_region,
    const std::vector<DefSFCodeToUint>& ull_max, const DefMap<DefAmrIndexUint>& mpi_inner_layer_tmp,
    std::map<int, DefMap<DefAmrIndexUint>>* const ptr_mpi_inner_layer) {
    int node_rank;
    std::vector<DefSFCodeToUint>::const_iterator iter_index;
    for (const auto& iter_node : vec_nodes_in_region) {
        if (mpi_inner_layer_tmp.find(iter_node) != mpi_inner_layer_tmp.end()) {
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
}  // end namespace amrproject
}  // end namespace rootproject
