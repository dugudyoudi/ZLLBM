//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved
/**
* @file grid_generation_serial.cpp
* @author Zhengliang Liu
* @brief functions used to generate grid serially.
* @date  2022-11-24
* @note  functions from geometry_manager will be called.
*/
#include "./auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "io/log_write.h"
#ifdef ENABLE_MPI
#include "mpi/mpi_manager.h"
#endif  // ENABLE_MPI
namespace rootproject {
namespace amrproject {
void GridManagerInterface::GenerateInitialMeshBasedOnGeoSerial(
    const std::vector<std::shared_ptr<GeometryInfoInterface>> &vec_geo_info) {
    std::vector<DefMap<DefUint>> sfbitset_one_lower_level;
    GenerateGridFromHighToLowLevelSerial(
        vec_geo_info, &sfbitset_one_lower_level);
    InstantiateGridNodeAllLevelSerial(sfbitset_one_lower_level);
}
/**
* @brief   function to generate grid for all levels of refinement.
* @param[in]  vec_geo_info vector containing vertices
*             to instance storing geometry information.
* @param[out]  ptr_sfbitset_one_lower_level   maps to store 
*            nodes at each refinement level
*/
void GridManagerInterface::GenerateGridFromHighToLowLevelSerial(
    const std::vector<std::shared_ptr<GeometryInfoInterface>>&vec_geo_info,
    std::vector<DefMap<DefUint>>* const ptr_sfbitset_one_lower_level) {
    ptr_sfbitset_one_lower_level->resize(k0MaxLevel_ + 1);
    SFBitsetAuxInterface* sfbitset_aux_ptr = this->GetSFBitsetAuxPtr();

    if (vec_ptr_grid_info_.size() != k0MaxLevel_ + 1) {
        LogError("Number of grid refinement level " + std::to_string(
            vec_ptr_grid_info_.size() - 1) + " is different from k0MaxLevel_ "
            + std::to_string(k0MaxLevel_) + ", need to create grid instance "
            + "before calling the function: GenerateGridFromHighToLowLevelSerial.");
    }

    // generate tracking and ghost nodes based on geometries
    DefSizet i_geo = 0;
    for (auto& iter : vec_geo_info) {
        iter->i_geo_ = i_geo;
        iter->FindTrackingNodeBasedOnGeo(sfbitset_aux_ptr,
            vec_ptr_grid_info_.at(iter->i_level_).get());
        ++i_geo;
    }

    std::vector<DefReal> flood_fill_origin(k0GridDims_);
    std::map<std::pair<ECriterionType, DefSizet>, DefMap<DefUint>>
        innermost_layer, outermost_layer;
    GeometryInfoInterface* ptr_geo_info = nullptr;
    struct NumExtendLayer {
        std::vector<DefLUint> neg;
        std::vector<DefLUint> pos;
    } num_extend_layer;
    num_extend_layer.neg = std::vector<DefLUint>(k0GridDims_, 0);
    num_extend_layer.pos = std::vector<DefLUint>(k0GridDims_, 0);
    /** @todo extended outer layers based on each outmost layer could be 
    stored in individual map, which will avoid interfering between different
    layers but is more time consuming since additional nodes for every 
    outmost layers are need to identify inside and outside region. Here
    nodes are inserted in one map (ptr_sfbitset_one_lower_level) rather
    than several ones, thus processes of extending inner and outer
    are the same*/
    for (DefSizet i_level = k0MaxLevel_; i_level > 0; --i_level) {
        std::map<std::pair<ECriterionType, DefSizet>, NumExtendLayer>
            map_num_extend_inner_layer, map_num_extend_outer_layer;
        std::map<std::pair<ECriterionType, DefSizet>, DefMap<DefUint>>
            innermost_layer_current, outermost_layer_current;
        // add nodes based on tracking nodes
        for (auto& iter_tracking_grid_info : vec_ptr_grid_info_.at(i_level)
            ->map_ptr_tracking_grid_info_) {
            DefMap<DefUint> node_near_tracking;
            // set number of extended layers
            SetNumberOfExtendLayerForGrid(i_level,
                *(vec_geo_info.at(iter_tracking_grid_info.first.second)),
                &iter_tracking_grid_info.second->k0ExtendInnerNeg,
                &iter_tracking_grid_info.second->k0ExtendInnerPos,
                &iter_tracking_grid_info.second->k0ExtendOuterNeg,
                &iter_tracking_grid_info.second->k0ExtendOuterPos);

            for (DefUint idims = 0; idims < k0GridDims_; ++idims) {
                num_extend_layer.neg[idims] = iter_tracking_grid_info.second->
                    k0ExtendOuterNeg[idims] - 1;
                num_extend_layer.pos[idims] = iter_tracking_grid_info.second->
                    k0ExtendOuterPos[idims] - 1;
            }
            map_num_extend_outer_layer.insert(
                { iter_tracking_grid_info.first, num_extend_layer });
            // tracking nodes need to identify inside and outside
            if (iter_tracking_grid_info.second->grid_extend_type_
                == EGridExtendType::kInAndOut) {
                GenerateGridNodeNearTrackingNode(
                    i_level, iter_tracking_grid_info.first, &node_near_tracking);
                for (DefUint idims = 0; idims < k0GridDims_; ++idims) {
                    num_extend_layer.neg[idims] = iter_tracking_grid_info.second->
                        k0ExtendInnerNeg[idims] - 1;
                    num_extend_layer.pos[idims] = iter_tracking_grid_info.second->
                        k0ExtendInnerPos[idims] - 1;
                }
                map_num_extend_inner_layer.insert(
                    { iter_tracking_grid_info.first, num_extend_layer });
                DefMap<DefUint> node_inside, node_outside;
                IdentifyTypeOfLayerByFloodFill(i_level - 1,
                    iter_tracking_grid_info.first.second,
                    vec_geo_info[iter_tracking_grid_info.first.second]
                    ->GetFloodFillOriginArrAsVec(),
                    node_near_tracking, &node_outside, &node_inside);
                innermost_layer.insert(
                    { iter_tracking_grid_info.first, node_inside });
                outermost_layer.insert(
                    { iter_tracking_grid_info.first, node_outside });
            } else {
                outermost_layer.insert(
                    { iter_tracking_grid_info.first, node_near_tracking });
            }
            ptr_sfbitset_one_lower_level->at(i_level).insert(
                node_near_tracking.begin(), node_near_tracking.end());
        }
        // extend the inner layer (number of extended layer
        // can be identified for each tracking grid)
        InterfaceLayerInfo* ptr_interface_info = nullptr;
        for (auto& iter_inner : innermost_layer) {
            if (map_num_extend_inner_layer.find(iter_inner.first)
                == map_num_extend_inner_layer.end()) {
                vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .insert({ iter_inner.first,
                        std::make_shared<InterfaceLayerInfo>() });
                ptr_interface_info = vec_ptr_grid_info_.at(i_level)
                    ->map_ptr_interface_layer_info_.at(iter_inner.first).get();
                SetNumberOfExtendLayerForGrid(i_level,
                    *(vec_geo_info.at(iter_inner.first.second)),
                    &ptr_interface_info->k0ExtendInnerNeg,
                    &ptr_interface_info->k0ExtendInnerPos,
                    &ptr_interface_info->k0ExtendOuterNeg,
                    &ptr_interface_info->k0ExtendOuterPos);
                num_extend_layer.neg =
                    ptr_interface_info->k0ExtendInnerNeg;
                num_extend_layer.pos =
                    ptr_interface_info->k0ExtendInnerPos;
                map_num_extend_inner_layer.insert(
                    { iter_inner.first, num_extend_layer });
                num_extend_layer.neg =
                    ptr_interface_info->k0ExtendOuterNeg;
                num_extend_layer.pos =
                    ptr_interface_info->k0ExtendOuterPos;
                map_num_extend_outer_layer.insert(
                    { iter_inner.first, num_extend_layer });
            }
            if (innermost_layer_current.find(iter_inner.first)
                == innermost_layer_current.end()) {
                innermost_layer_current.insert({ iter_inner.first, {} });
                ExtendGivenNumbOfLayer(i_level,
                map_num_extend_inner_layer.at(iter_inner.first).neg,
                map_num_extend_inner_layer.at(iter_inner.first).pos,
                iter_inner.second,
                &ptr_sfbitset_one_lower_level->at(i_level),
                &innermost_layer_current.at(iter_inner.first));
            } else {
                innermost_layer_current.clear();
            }
        }

        // extend the outer layer (number of extended layer
        // can be identified for each tracking grid)
        for (auto& iter_outer : outermost_layer) {
            if (map_num_extend_outer_layer.find(iter_outer.first)
                == map_num_extend_outer_layer.end()) {
                vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .insert({ iter_outer.first,
                        std::make_shared<InterfaceLayerInfo>() });
                ptr_interface_info = vec_ptr_grid_info_.at(i_level)
                    ->map_ptr_interface_layer_info_.at(iter_outer.first).get();
                SetNumberOfExtendLayerForGrid(i_level,
                    *(vec_geo_info.at(iter_outer.first.second)),
                    &ptr_interface_info->k0ExtendInnerNeg,
                    &ptr_interface_info->k0ExtendInnerPos,
                    &ptr_interface_info->k0ExtendOuterNeg,
                    &ptr_interface_info->k0ExtendOuterPos);
                num_extend_layer.neg =
                    ptr_interface_info->k0ExtendInnerNeg;
                num_extend_layer.pos =
                    ptr_interface_info->k0ExtendInnerPos;
                map_num_extend_inner_layer.insert(
                    { iter_outer.first, num_extend_layer });
                num_extend_layer.neg =
                    ptr_interface_info->k0ExtendOuterNeg;
                num_extend_layer.pos =
                    ptr_interface_info->k0ExtendOuterPos;
                map_num_extend_outer_layer.insert(
                    { iter_outer.first, num_extend_layer });
            }
            if (outermost_layer_current.find(iter_outer.first)
                == outermost_layer_current.end()) {
                outermost_layer_current.insert({ iter_outer.first, {} });
                ExtendGivenNumbOfLayer(i_level,
                map_num_extend_outer_layer.at(iter_outer.first).neg,
                map_num_extend_outer_layer.at(iter_outer.first).pos,
                iter_outer.second,
                &ptr_sfbitset_one_lower_level->at(i_level),
                &outermost_layer_current.at(iter_outer.first));
            } else {
                outermost_layer_current.clear();
            }
        }
        // find interface between different grid
        DefSizet level_low = i_level - 1;
        for (auto& iter_inner : innermost_layer_current) {
            innermost_layer.at(iter_inner.first).clear();
            if (vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                .find(iter_inner.first) == vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.end()) {
                vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                    .insert({ iter_inner.first,
                        std::make_shared<InterfaceLayerInfo>() });
            }
            ptr_interface_info = vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.at(iter_inner.first).get();
            ptr_interface_info->vec_inner_coarse2fine_.resize(
                vec_ptr_grid_info_.at(level_low)->k0NumCoarse2FineLayer_);
            FindInterfaceBetweenGrid(i_level, iter_inner.second,
                &ptr_sfbitset_one_lower_level->at(i_level),
                &ptr_interface_info->vec_inner_coarse2fine_.back(),
                &ptr_sfbitset_one_lower_level->at(level_low),
                &innermost_layer.at(iter_inner.first));
        }
        for (auto& iter_outer : outermost_layer_current) {
            outermost_layer.at(iter_outer.first).clear();
            if (vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                .find(iter_outer.first) == vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.end()) {
                vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                    .insert({ iter_outer.first,
                        std::make_shared<InterfaceLayerInfo>() });
            }
            ptr_interface_info = vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.at(iter_outer.first).get();
            ptr_interface_info->vec_outer_coarse2fine_.resize(
                vec_ptr_grid_info_.at(level_low)->k0NumCoarse2FineLayer_);
            FindInterfaceBetweenGrid(i_level, iter_outer.second,
                &ptr_sfbitset_one_lower_level->at(i_level),
                &ptr_interface_info->vec_outer_coarse2fine_.back(),
                &ptr_sfbitset_one_lower_level->at(level_low),
                &outermost_layer.at(iter_outer.first));
        }
    }
}
void GridManagerInterface::InstantiateGridNodeAllLevelSerial(
    const std::vector<DefMap<DefUint>>& sfbitset_one_lower_level) {
    InterfaceLayerInfo* ptr_interface_info = nullptr;
    InterfaceLayerInfo* ptr_interface_info_lower = nullptr;
    DefSizet layer_coarse0, layer_coarse_m1, layer0, layer_m1, layer_m2;

#ifdef ENABLE_MPI

#endif  // ENABLE_MPI

    DefMap<DefUint> background_occupied;
    DefUint flag_temp, flag_refined = kFlagExist_;
    for (DefSizet i_level = k0MaxLevel_; i_level > 0; --i_level) {
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        GridInfoInterface& grid_info_lower =
            *(vec_ptr_grid_info_.at(i_level - 1));
        DefMap<GridNode>& map_grid = grid_info.map_grid_node_;
        DefMap<GridNode>& map_grid_lower = grid_info_lower.map_grid_node_;
#ifdef DEBUG_CHECK_GRID
        if (grid_info_lower.k0NumCoarse2FineLayer_ < 2) {
            LogError("number of coarse to fine layers at level "
                + std::to_string(i_level - 1) + " is at least 1");
        }
        if (grid_info.k0NumFine2CoarseLayer_ < 3) {
            LogError("number of fine to coarse layers at level "
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
                    ptr_interface_info_lower
                    ->vec_inner_coarse2fine_.at(layer_coarse0),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower
                    ->vec_inner_coarse2fine_.at(layer_coarse_m1),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
                // insert node instance
                DefSizet maxlayer = ptr_interface_info
                    ->vec_inner_fine2coarse_.size();
                for (DefSizet ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kFlagFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kFlagFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_inner_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first)
                            == map_grid.end()) {
                            map_grid.insert({ iter_layer_node.first,
                                grid_info.k0GridNodeInstance_ });
                            map_grid.at(iter_layer_node.first).flag_status_ =
                                flag_temp;
                        } else {
                            map_grid.at(iter_layer_node.first).flag_status_ |=
                                flag_temp;
                        }
                    }
                }
                maxlayer = ptr_interface_info_lower
                    ->vec_inner_coarse2fine_.size();
                for (DefSizet ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kFlagCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kFlagCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_inner_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first)
                            == map_grid_lower.end()) {
                            map_grid_lower.insert({ iter_layer_node.first,
                                grid_info.k0GridNodeInstance_ });
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
                DefSizet maxlayer = ptr_interface_info
                    ->vec_outer_fine2coarse_.size();
                for (DefSizet ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kFlagFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kFlagFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_outer_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first)
                            == map_grid.end()) {
                            map_grid.insert({ iter_layer_node.first,
                                grid_info.k0GridNodeInstance_ });
                            map_grid.at(iter_layer_node.first).flag_status_ =
                                flag_temp;
                        } else {
                            map_grid.at(iter_layer_node.first).flag_status_ |=
                                flag_temp;
                        }
                    }
                }
                maxlayer = ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.size();
                for (DefSizet ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kFlagCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kFlagCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_outer_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first)
                            == map_grid_lower.end()) {
                            map_grid_lower.insert({ iter_layer_node.first,
                                grid_info.k0GridNodeInstance_ });
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
        std::vector<DefSFBitset> bitset_cell_lower, bitset_all;

        DefSFBitset bitset_background;
        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            if (CheckCoincideBackground(i_level - 1,
                iter_low.first, &bitset_background)) {
                if (background_occupied.find(bitset_background) ==
                    background_occupied.end()) {
                    background_occupied.insert({ bitset_background, kFlag0_ });
                }
            }
            if (NodesBelongToOneCell(iter_low.first,
                sfbitset_one_lower_level.at(i_level), &bitset_cell_lower)) {
                FindAllNodesInACellAtLowerLevel(
                    bitset_cell_lower, &bitset_all);
                for (const auto& iter_node : bitset_all) {
                    if (map_grid.find(iter_node) == map_grid.end()) {
                        map_grid.insert({ iter_node,
                                grid_info.k0GridNodeInstance_ });
                    }
                }
            }
        }
    }

    // instantiate background nodes
    GenerateBackgroundGrid(background_occupied);
}
void GridManagerInterface::FindOverlappingLayersBasedOnOutermostCoarse(
    const DefMap<DefUint>& layer_coarse_0,
    const DefMap<DefUint>& sfbitset_exist,
    DefMap<DefUint>* const ptr_layer_coarse_m1,
    DefMap<DefUint>* const ptr_layer_fine_0,
    DefMap<DefUint>* const ptr_layer_fine_m1,
    DefMap<DefUint>* const ptr_layer_fine_m2) {
#ifdef DEBUG_CHECK_GRID
    if (&layer_coarse_0 == ptr_layer_coarse_m1) {
        LogError("input (layer_coarse_0)"
            " should not be the same as output (ptr_layer_coarse_m1)");
    }
#endif  // DEBUG_CHECK_GRID

    std::vector<DefSFBitset> corner_bitsets;
    for (const auto& iter_node : layer_coarse_0) {
        FindCornersForNeighbourCells(
            iter_node.first, &corner_bitsets);
        for (auto& iter_conner : corner_bitsets) {
            IdentifyInterfaceForACell(kFlagGridInterfaceOutermost_,
                iter_conner, sfbitset_exist,
                ptr_layer_fine_m2, ptr_layer_fine_m1, ptr_layer_fine_0);
        }
    }
    // coarse node on the layer_coarse_m1 layer, which overlap with
    // layer0 layer at one level higher
#ifdef DEBUG_CHECK_GRID
    if (ptr_layer_fine_0 == ptr_layer_coarse_m1) {
        LogError("input (ptr_layer_fine_0)"
            " should not be the same as output (ptr_layer_coarse_m1)");
    }
#endif  // DEBUG_CHECK_GRID
    OverlapLayerFromHighToLow(*ptr_layer_fine_m2, ptr_layer_coarse_m1);
}
/**
* @brief   function to add a given number of layers.
* @param[in]  i_level   level of refinement.
* @param[in]  flag_extend type of extending grid which determines
*                  the number of layers need to be extended in all directions.
* @param[in]  the based layer for extension.
* @param[out]  ptr_map_exist existing nodes.
* @param[out]  ptr_map_outmost nodes on the outmost layer.
* @note  add the number of layer around the geometry in all the directions.
*        Used to set number of layers inside the geometry different from that
*        outsite the geometry.
*/
void GridManagerInterface::ExtendGivenNumbOfLayer(
    const DefSizet i_level, const std::vector<DefLUint> num_extend_neg,
    const std::vector<DefLUint> num_extend_pos,
    const DefMap<DefUint>& map_start_layer,
    DefMap<DefUint>* const ptr_map_exist,
    DefMap<DefUint>* const ptr_map_outmost) const {
    std::vector<DefSFBitset> vec_bitset_min(k0GridDims_, 0),
        vec_bitset_max(k0GridDims_, 0);
    // space filling code is at one level lower (ilevel -1)
    ComputeSFBitsetOnBoundaryAtGivenLevel(
        i_level - 1, &vec_bitset_min, &vec_bitset_max);
    // calculate layers need to be extended in each direction
    if (num_extend_neg.size() != k0GridDims_) {
        LogError("Dimension of num_extend_neg at" + std::to_string(i_level)
            + " should be " + std::to_string(k0GridDims_)
            + " rather than " + std::to_string(num_extend_neg.size())
            + "in GridManagerInterface::ExtendGivenNumbOfLayer.");
    }
    if (num_extend_pos.size() != k0GridDims_) {
        LogError("Dimension of num_extend_pos at" + std::to_string(i_level)
            + " should be " + std::to_string(k0GridDims_)
            + " rather than " + std::to_string(num_extend_pos.size())
            + "in GridManagerInterface::ExtendGivenNumbOfLayer.");
    }
    // find the maximum layers need to be extended
    DefLUint extend_layer_max = 0;
    for (DefUint i = 0; i < k0GridDims_; ++i) {
        if (num_extend_neg[i] > extend_layer_max) {
            extend_layer_max = num_extend_neg[i];
        }
        if (num_extend_pos[i] > extend_layer_max) {
            extend_layer_max = num_extend_pos[i];
        }
    }
    std::vector<bool> bool_extend_neg(k0GridDims_, true),
        bool_extend_pos(k0GridDims_, true);
    DefMap<DefUint> map_input_layer(map_start_layer);
    std::vector<DefMap<DefUint>> vec_boundary_min, vec_boundary_max;
    std::vector<DefLUint> num_extend_neg_temp(num_extend_neg),
        num_extend_pos_temp(num_extend_pos);
    // extend grid layer by layer
    for (DefLUint i_layer = 0; i_layer < extend_layer_max; ++i_layer) {
        DefMap<DefUint> map_output_layer;
        ExtendOneLayerGrid(map_input_layer, vec_bitset_min, vec_bitset_max,
            bool_extend_neg, bool_extend_pos, &map_output_layer, ptr_map_exist,
            ptr_map_outmost, &vec_boundary_min, &vec_boundary_max);
        map_input_layer = std::move(map_output_layer);
        for (DefUint i = 0; i < k0GridDims_; ++i) {
            if (bool_extend_neg[i]) {
                num_extend_neg_temp[i] -= 1;
                if (num_extend_neg_temp[i] == 0) {
                    bool_extend_neg[i] = false;
                }
            }
            if (bool_extend_pos[i]) {
                num_extend_pos_temp[i] -= 1;
                if (num_extend_pos_temp[i] == 0) {
                    bool_extend_pos[i] = false;
                }
            }
        }
    }
    // added the last extended layer to the outmost one
    for (const auto& iter : map_input_layer) {
        ptr_map_outmost->insert(iter);
    }
}
/**
* @brief   function to do flood fill on existing nodes around geometry points
*          (only two layers)
* @param[in]  sfbitset_inside  space filling code of a node inside geometry
* @param[in]  map_nodes_exist nodes exist for flood fill.
* @param[out]  ptr_map_nodes_outside
*                  nodes outside the geometry.
* @param[out]  ptr_map_nodes_inside
*                  nodes inside the geometry (colored by flood fill method).
*/
int GridManagerInterface::FloodFillForInAndOut(
    const DefSFBitset& sfbitset_inside,
    const DefMap<DefUint>&  map_nodes_exist,
    DefMap<DefUint>* const ptr_map_nodes_outside,
    DefMap<DefUint>* const ptr_map_nodes_inside) const {

    std::vector<DefSFBitset> bitset_neigh;
    if (map_nodes_exist.find(sfbitset_inside) != map_nodes_exist.end()) {
        // node inside the geometry is in map_nodes_exist, then
        // the colored nodes are indentified as the outside nodes
        std::vector<DefSFBitset> bitset_neigh;
        for (auto iter = map_nodes_exist.begin();
            iter != map_nodes_exist.end(); ++iter) {
            FindAllNeighboursSFBitset(iter->first, &bitset_neigh);
            for (const auto& iter_neigh : bitset_neigh) {
                if (map_nodes_exist.find(iter_neigh)
                    == map_nodes_exist.end()) {
                    ptr_map_nodes_outside->insert({ iter->first, kFlag0_ });
                    break;
                }
            }
        }
        return 1;
    } else {
        // create a copy of map_nodes_exist and add one layer for flood fill
        DefUint flag_bit_colored = 1 << (sizeof(DefUint) * 8 - 1);
        DefUint flag_bit_exist = flag_bit_colored | (sizeof(DefUint) * 8 - 2);
        DefMap<DefUint> map_nodes_temp(map_nodes_exist);
        // generate one more layer around the input grid
        for (auto iter = map_nodes_exist.begin();
            iter != map_nodes_exist.end(); ++iter) {
            map_nodes_temp.at(iter->first) = flag_bit_exist;
            FindAllNeighboursSFBitset(iter->first, &bitset_neigh);
            for (const auto& iter_neigbour : bitset_neigh) {
                if (map_nodes_temp.find(iter_neigbour)
                    == map_nodes_temp.end()) {
                    map_nodes_temp.insert({ iter_neigbour, kFlag0_ });
                }
            }
        }
        // flood fill
        std::vector<DefSFBitset> vec_sfbitset_stk, vec_bitset(k0NumNeighbours_);
        DefUint i = 0;
        DefSFBitset sfbitset_seed;
        vec_sfbitset_stk.push_back(sfbitset_inside);
        while (!vec_sfbitset_stk.empty() && i < imax_flood_fill) {
            sfbitset_seed = vec_sfbitset_stk.back();
            vec_sfbitset_stk.pop_back();
            ++i;
            if (map_nodes_temp.find(sfbitset_seed) == map_nodes_temp.end()) {
            } else if ((map_nodes_temp.at(sfbitset_seed) & flag_bit_exist)
                == flag_bit_exist) {
                ptr_map_nodes_inside->insert({ sfbitset_seed, kFlag0_ });
            } else if ((map_nodes_temp.at(sfbitset_seed) & flag_bit_colored)
                != flag_bit_colored) {
                // color the node
                map_nodes_temp.at(sfbitset_seed) |= flag_bit_colored;
                // add neighbouring nodes to seed
                PushBackSFBitsetInFloodFill(sfbitset_seed, &vec_sfbitset_stk);
            }
        }
        if (i == imax_flood_fill) {
            return 2;
        }
        std::vector<DefSFBitset> bitset_neigh;
        for (auto iter = map_nodes_exist.begin();
            iter != map_nodes_exist.end(); ++iter) {
            if (ptr_map_nodes_inside->find(iter->first)
                == ptr_map_nodes_inside->end()) {
                FindAllNeighboursSFBitset(iter->first, &bitset_neigh);
                for (const auto& iter_neigh : bitset_neigh) {
                    if (map_nodes_exist.find(iter_neigh)
                        == map_nodes_exist.end()) {
                        ptr_map_nodes_outside->insert(
                            { iter->first, kFlag0_ });
                        break;
                    }
                }
            }
        }
        return 0;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
