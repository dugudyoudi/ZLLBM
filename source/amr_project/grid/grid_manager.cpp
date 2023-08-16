//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
*/
#include <string>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to create instance based on the type of solver.
*/
void GridManagerInterface::CreateSameGridInstanceForAllLevel(
    GridInfoCreatorInterface* ptr_grid_creator) {
    for (DefAmrIndexUint i_level = 0; i_level < k0MaxLevel_ + 1; ++i_level) {
        vec_ptr_grid_info_.emplace_back(ptr_grid_creator->CreateGridInfo());
        GridInfoInterface& grid_ref = *(vec_ptr_grid_info_).back();
        grid_ref.i_level_ = i_level;
        grid_ref.grid_space_ = std::vector<DefReal>(k0GridDims_, 0.);
        // set computational cost for each node 2^i_level
        grid_ref.computational_cost_ = static_cast<DefAmrUint>(TwoPowerN(i_level));
        std::vector<DefReal> domain_dx_ = GetDomainDxArrAsVec();
        for (DefAmrIndexUint idim = 0; idim < k0GridDims_; ++idim) {
            grid_ref.grid_space_.at(idim) =
                domain_dx_.at(idim) /
                static_cast<DefReal>(TwoPowerN(i_level));
        }
        grid_ref.k0GridNodeInstance_.flag_status_ = kNodeStatusExist_;
        grid_ref.SetNumberOfVecElements();
    }
}
/**
* @brief function to setup default grid related parameters.
*/
void GridManagerInterface::DefaultInitialization(const DefAmrIndexUint max_level) {
    k0MaxLevel_ = max_level;
}
/**
* @brief   function to check if the geometry is out of domain.
* @param[in]  coordinate_min  minimum coordinates.
* @param[in]  coordinate_min  maximum coordinates.
* @param[in]  domain_offset  offset of the background grid.
* @param[in]  domain_size  size of the computational domain.
* @return  0: coordinates does not exceed the computational domain; 1, 4:
*          coordinates are less than the lower bound; 2, 8: coordinates
*          greater than the upper bound;
*/
int GridManagerInterface::CheckIfPointOutsideDomain(
    const std::array<DefReal, 2>& coordinate_min,
    const std::array<DefReal, 2>& coordinate_max,
    const std::array<DefReal, 2>& domain_offset,
    const std::array<DefReal, 2>& domain_size) const {
    int status = 0;
    if (coordinate_min[kXIndex] < domain_offset[kXIndex]) {
        status |= 1;
    }
    if (coordinate_max[kXIndex] > domain_size[kXIndex]
        + domain_offset[kXIndex]) {
        status |= 2;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        status |= 4;
    }
    if (coordinate_max[kYIndex] > domain_size[kYIndex]
        + domain_offset[kYIndex]) {
        status |= 8;
    }
    return status;
}
/**
* @brief   function to check if the geometry is out of domain.
* @param[in]  coordinate_min  minimum coordinates.
* @param[in]  coordinate_min  maximum coordinates.
* @param[in]  domain_offset  offset of the background grid.
* @param[in]  domain_size  size of the computational domain.
* @return  0: coordinates does not exceed the computational domain; 1, 4, 16:
*          coordinates are less than the lower bound; 2, 8, 32: coordinates
*          greater than the upper bound;
*/
int GridManagerInterface::CheckIfPointOutsideDomain(
    const std::array<DefReal, 3>& coordinate_min,
    const std::array<DefReal, 3>& coordinate_max,
    const std::array<DefReal, 3>& domain_offset,
    const std::array<DefReal, 3>& domain_size) const {
    int status = 0;
    if (coordinate_min[kXIndex] < domain_offset[kXIndex]) {
        status |= 1;
    }
    if (coordinate_max[kXIndex] > domain_size[kXIndex]
        + domain_offset[kXIndex]) {
        status |= 2;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        status |= 4;
    }
    if (coordinate_max[kYIndex] > domain_size[kYIndex]
        + domain_offset[kYIndex]) {
        status |= 8;
    }
    if (coordinate_min[kZIndex] < domain_offset[kZIndex]) {
        status |= 16;
    }
    if (coordinate_max[kZIndex] > domain_size[kZIndex]
        + domain_offset[kZIndex]) {
        status |= 32;
    }
    return status;
}
/**
* @brief   function to set extened layer of grid based on geometry information.
* @param[in]  i_level  refinement level.
* @param[in]  geo_info  geometry information.
* @param[out]  ptr_inner_layer_neg  number of extended layer inside geometry
*               in negative directions.
* @param[out]  ptr_inner_layer_pos  number of extended layer inside geometry
*               in positive directions.
* @param[out]  ptr_outer_layer_neg  number of extended layer outside geometry
*               in negative directions.
* @param[out]  ptr_outer_layer_pos  number of extended layer outside geometry
*               in positive directions.
*/
void GridManagerInterface::SetNumberOfExtendLayerForGrid(const DefAmrIndexUint i_level,
    const GeometryInfoInterface& geo_info,
    std::vector<DefAmrIndexLUint>* const ptr_inner_layer_neg,
    std::vector<DefAmrIndexLUint>* const ptr_inner_layer_pos,
    std::vector<DefAmrIndexLUint>* const ptr_outer_layer_neg,
    std::vector<DefAmrIndexLUint>* const ptr_outer_layer_pos) {
    std::vector<DefAmrIndexLUint> layer_min(k0GridDims_, k0IntExtendMin_);
    ptr_inner_layer_neg->assign(layer_min.begin(), layer_min.end());
    ptr_inner_layer_pos->assign(layer_min.begin(), layer_min.end());
    ptr_outer_layer_neg->assign(layer_min.begin(), layer_min.end());
    ptr_outer_layer_pos->assign(layer_min.begin(), layer_min.end());
    if ((geo_info.k0IntInnerExtend_.size() > i_level)
        && (geo_info.k0IntInnerExtend_.at(i_level) > k0IntExtendMin_)) {
        ptr_inner_layer_neg->at(kXIndex) =
            geo_info.k0IntInnerExtend_.at(i_level);
        ptr_inner_layer_neg->at(kYIndex) =
            geo_info.k0IntInnerExtend_.at(i_level);
        ptr_inner_layer_pos->at(kXIndex) =
            geo_info.k0IntInnerExtend_.at(i_level);
        ptr_inner_layer_pos->at(kYIndex) =
            geo_info.k0IntInnerExtend_.at(i_level);
        if (k0GridDims_ == 3) {
            ptr_inner_layer_neg->at(kZIndex) =
                geo_info.k0IntInnerExtend_.at(i_level);
            ptr_inner_layer_pos->at(kZIndex) =
                geo_info.k0IntInnerExtend_.at(i_level);
        }
    } else {
        ptr_inner_layer_neg->at(kXIndex) = k0IntExtendMin_;
        ptr_inner_layer_neg->at(kYIndex) = k0IntExtendMin_;
        ptr_inner_layer_pos->at(kXIndex) = k0IntExtendMin_;
        ptr_inner_layer_pos->at(kYIndex) = k0IntExtendMin_;
        if (k0GridDims_ == 3) {
            ptr_inner_layer_neg->at(kZIndex) = k0IntExtendMin_;
            ptr_inner_layer_pos->at(kZIndex) = k0IntExtendMin_;
        }
    }
    if ((geo_info.k0XIntExtendNegative_.size() > i_level)
        && (geo_info.k0XIntExtendNegative_.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_neg->at(kXIndex) =
            geo_info.k0XIntExtendNegative_.at(i_level);
    } else {
        ptr_outer_layer_neg->at(kXIndex) = k0IntExtendMin_;
    }
    if ((geo_info.k0XIntExtendPositive_.size() > i_level)
        && (geo_info.k0XIntExtendPositive_.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_pos->at(kXIndex) =
            geo_info.k0XIntExtendPositive_.at(i_level);
    } else {
        ptr_outer_layer_pos->at(kXIndex) = k0IntExtendMin_;
    }
    if ((geo_info.k0YIntExtendNegative_.size() > i_level)
        && (geo_info.k0YIntExtendNegative_.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_neg->at(kYIndex) =
            geo_info.k0YIntExtendNegative_.at(i_level);
    } else {
        ptr_outer_layer_neg->at(kYIndex) = k0IntExtendMin_;
    }
    if ((geo_info.k0YIntExtendPositive_.size() > i_level)
        && (geo_info.k0YIntExtendPositive_.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_pos->at(kYIndex) =
            geo_info.k0YIntExtendPositive_.at(i_level);
    } else {
        ptr_outer_layer_pos->at(kYIndex) = k0IntExtendMin_;
    }
    if (k0GridDims_ == 3) {
        if ((geo_info.k0ZIntExtendNegative_.size() > i_level)
            && (geo_info.k0ZIntExtendNegative_.at(i_level) >
                k0IntExtendMin_)) {
            ptr_outer_layer_neg->at(kZIndex) =
                geo_info.k0ZIntExtendNegative_.at(i_level);
        } else {
            ptr_outer_layer_neg->at(kZIndex) = k0IntExtendMin_;
        }
        if ((geo_info.k0ZIntExtendPositive_.size() > i_level)
            && (geo_info.k0ZIntExtendPositive_.at(i_level) >
                k0IntExtendMin_)) {
            ptr_outer_layer_pos->at(kZIndex) =
                geo_info.k0ZIntExtendPositive_.at(i_level);
        } else {
            ptr_outer_layer_pos->at(kZIndex) = k0IntExtendMin_;
        }
    }
}
/**
 * @brief function to find overlapping layers between grid of adjacent refinement levels based on the outermost coarse layer.
 * @param[in]  layer_coarse_0  nodes on the outermost coarse layer.
 * @param[in]  sfbitset_exist   existent nodes.
 * @param[out] ptr_layer_coarse_m1 pointer to nodes in the second outermost coarse layer.
 * @param[out] ptr_layer_fine_0 pointer to nodes in the outermost fine layer.
 * @param[out] ptr_layer_fine_m1 pointer to nodes in the second outermost fine layer.
 * @param[out] ptr_layer_fine_m2 pointer to nodes in the third outermost fine layer.
 */
void GridManagerInterface::FindOverlappingLayersBasedOnOutermostCoarse(
    const DefMap<DefAmrUint>& layer_coarse_0, const DefMap<DefAmrIndexUint>& sfbitset_exist,
    DefMap<DefAmrUint>* const ptr_layer_coarse_m1, DefMap<DefAmrUint>* const ptr_layer_fine_0,
    DefMap<DefAmrUint>* const ptr_layer_fine_m1, DefMap<DefAmrUint>* const ptr_layer_fine_m2) {
#ifdef DEBUG_CHECK_GRID
    if (&layer_coarse_0 == ptr_layer_coarse_m1) {
        LogManager::LogError("input (layer_coarse_0)"
         " should not be the same as output (ptr_layer_coarse_m1)");
    }
#endif  // DEBUG_CHECK_GRID

    std::vector<DefSFBitset> corner_bitsets;
    for (const auto& iter_node : layer_coarse_0) {
        FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
        for (auto& iter_conner : corner_bitsets) {
            IdentifyInterfaceForACell(iter_conner, layer_coarse_0, sfbitset_exist,
             ptr_layer_fine_m2, ptr_layer_fine_m1, ptr_layer_fine_0);
        }
    }
    // coarse node on the layer_coarse_m1 layer, which overlap with
    // layer0 layer at one level higher
#ifdef DEBUG_CHECK_GRID
    if (ptr_layer_fine_0 == ptr_layer_coarse_m1) {
        LogManager::LogError("input (ptr_layer_fine_0)"
         " should not be the same as output (ptr_layer_coarse_m1)");
    }
#endif  // DEBUG_CHECK_GRID
    OverlapLayerFromHighToLow(*ptr_layer_fine_m2, ptr_layer_coarse_m1);
}
/**
 * @brief function to find and instantiate overlapping layers of refinement interfaces.
 * @param[in] sfbitset_one_lower_level  space filling codes at one lower refinement level.
 */
void GridManagerInterface::InstantiateOverlapLayerOfRefinementInterface(
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level) {
    InterfaceLayerInfo* ptr_interface_info = nullptr;
    InterfaceLayerInfo* ptr_interface_info_lower = nullptr;
    DefAmrIndexUint layer_coarse0, layer_coarse_m1, layer0, layer_m1, layer_m2;
    DefAmrUint flag_temp, flag_refined = kNodeStatusExist_;
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        GridInfoInterface& grid_info_lower =
            *(vec_ptr_grid_info_.at(i_level - 1));
        DefMap<GridNode>& map_grid = grid_info.map_grid_node_;
        DefMap<GridNode>& map_grid_lower = grid_info_lower.map_grid_node_;
        GridNode node_instance(grid_info.k0GridNodeInstance_);
#ifdef DEBUG_CHECK_GRID
        if (grid_info_lower.k0NumCoarse2FineLayer_ < 2) {
            LogManager::LogError("number of coarse to fine layers at level "
                + std::to_string(i_level - 1) + " is at least 1");
        }
        if (grid_info.k0NumFine2CoarseLayer_ < 3) {
            LogManager::LogError("number of fine to coarse layers at level "
                + std::to_string(i_level) + " is at least 3");
        }
#endif  // DEBUG_CHECK_GRID
        // find interface node at current level
        layer_coarse0 = grid_info_lower.k0NumCoarse2FineLayer_ - 1;
        layer_coarse_m1 = layer_coarse0 - 1;
        layer0 = grid_info.k0NumFine2CoarseLayer_ - 1;
        layer_m1 = layer0 - 1;
        layer_m2 = layer_m1 - 1;
        for (auto& iter_interface :
            grid_info_lower.map_ptr_interface_layer_info_) {
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
                    ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse0),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_m1),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
                // insert node instance
                DefAmrIndexUint maxlayer = DefAmrIndexUint(ptr_interface_info->vec_inner_fine2coarse_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_inner_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first)
                            == map_grid.end()) {
                            map_grid.insert({ iter_layer_node.first, node_instance });
                            map_grid.at(iter_layer_node.first).flag_status_ =
                                flag_temp;
                        } else {
                            map_grid.at(iter_layer_node.first).flag_status_ |=
                                flag_temp;
                        }
                    }
                }
                maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_inner_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first)
                            == map_grid_lower.end()) {
                            map_grid_lower.insert({ iter_layer_node.first, node_instance});
                            map_grid_lower.at(iter_layer_node.first).flag_status_ =
                                flag_temp;
                        } else {
                            map_grid_lower.at(iter_layer_node.first).flag_status_ |=
                                flag_temp;
                        }
                    }
                }
            }
            // interface outside the geometry
            if (ptr_interface_info_lower->vec_outer_coarse2fine_.size() > 0) {
                ptr_interface_info->vec_outer_fine2coarse_.resize(
                    grid_info.k0NumFine2CoarseLayer_);
                FindOverlappingLayersBasedOnOutermostCoarse(
                    ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(layer_coarse0),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(layer_coarse_m1),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m2));
                // insert node instance
                DefAmrIndexUint maxlayer = DefAmrIndexUint(ptr_interface_info->vec_outer_fine2coarse_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_outer_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first)
                            == map_grid.end()) {
                            map_grid.insert({ iter_layer_node.first, node_instance});
                            map_grid.at(iter_layer_node.first).flag_status_ =
                                flag_temp;
                        } else {
                            map_grid.at(iter_layer_node.first).flag_status_ |=
                                flag_temp;
                        }
                    }
                }
                maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_outer_coarse2fine_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_outer_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first)
                            == map_grid_lower.end()) {
                            map_grid_lower.insert({ iter_layer_node.first, node_instance});
                            map_grid_lower.at(iter_layer_node.first).flag_status_ =
                                flag_temp;
                        } else {
                            map_grid_lower.at(iter_layer_node.first).flag_status_ |=
                                flag_temp;
                        }
                    }
                }
            }
        }
    }
}
/**
 * @brief function to instantiate grid noes for all refinement levels.
 * @param[in] sfbitset_min  minimum space filling code of the current rank.
 * @param[in] sfbitset_max  maximum space filling code of the current rank.
 * @param[in] partitioned_interface_background interface of partitioned block of current rank.
 * @param[in] space filling codes at one lower refinement level.
 */
void GridManagerInterface::InstantiateGridNodeAllLevel(const DefSFBitset sfbitset_min,
    const DefSFBitset sfbitset_max, const DefMap<DefAmrIndexUint>& partitioned_interface_background,
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level) {
    DefMap<DefAmrIndexUint> background_occupied;
    InstantiateOverlapLayerOfRefinementInterface(sfbitset_one_lower_level);
    DefSFCodeToUint code_min = sfbitset_min.to_ullong(), code_max = sfbitset_max.to_ullong();
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
#ifdef ENABLE_MPI
        std::vector<DefSFBitset> corresponding_ones(k0GridDims_),
            domain_min_m1_n_level(k0GridDims_), domain_max_p1_n_level(k0GridDims_);
        std::vector<DefSFBitset> vec_ghost_for_mpi_communication;
        GetNLevelCorrespondingOnes(i_level, &corresponding_ones);
        GetMinM1AtGivenLevel(i_level, &domain_min_m1_n_level);
        GetMaxP1AtGivenLevel(i_level, &domain_max_p1_n_level);
#endif  // ENABLE_MPI
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        DefMap<GridNode>& map_grid = grid_info.map_grid_node_;
        GridNode node_instance(grid_info.k0GridNodeInstance_);
        node_instance.flag_status_ = kNodeStatusExist_;
        std::vector<DefSFBitset> bitset_cell_lower, bitset_all;
        DefSFBitset bitset_background;
        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            if (CheckCoincideBackground(i_level - 1,
                iter_low.first, &bitset_background)) {
                if (background_occupied.find(bitset_background)
                 == background_occupied.end()) {
                    background_occupied.insert({ bitset_background, kFlagSize0_ });
                }
            }
            if (NodesBelongToOneCell(iter_low.first,
                sfbitset_one_lower_level.at(i_level), &bitset_cell_lower)) {
                FindAllNodesInACellAtLowerLevel(bitset_cell_lower, &bitset_all);
                for (const auto& iter_node : bitset_all) {
                    if (map_grid.find(iter_node) == map_grid.end()) {
                        map_grid.insert({ iter_node, node_instance });
                    }
#ifdef ENABLE_MPI
                    // node is on the partitioned interface
                    if (CheckNodeOnOuterBoundaryOfBackgroundCell(i_level, code_min, code_max,
                     iter_node, domain_min_m1_n_level, domain_max_p1_n_level,
                     corresponding_ones, partitioned_interface_background)) {
                        SearchForGhostLayerForMinNMax(iter_node, grid_info.num_of_ghost_layer_,
                         code_min, code_max, domain_min_m1_n_level, domain_max_p1_n_level,
                         &vec_ghost_for_mpi_communication);
                        for (const auto& iter_ghost : vec_ghost_for_mpi_communication) {
                            map_grid.insert({ iter_ghost, node_instance });
                            map_grid.at(iter_ghost).flag_status_ |= kNodeStatusGhostCommunication_;
                        }
                    }
#endif  // ENABLE_MPI
                }
            }
        }
    }
    InstantiateBackgroundGrid(code_min, code_max, background_occupied);
}
}  // end namespace amrproject
}  // end namespace rootproject
