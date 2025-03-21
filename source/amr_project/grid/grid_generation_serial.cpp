//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
namespace rootproject {
namespace amrproject {
/**
* @brief   function to generate grid for all levels of refinement.
* @param[in]  vec_geo_info vector containing vertices to instance storing geometry information.
* @param[out]  ptr_sfbitset_one_lower_level   maps to store nodes at each refinement level
*/
void GridManagerInterface::GenerateGridFromHighToLowLevelSerial(
    const std::vector<std::shared_ptr<GeometryInfoInterface>>&vec_geo_info,
    std::vector<DefMap<DefInt>>* const ptr_sfbitset_one_lower_level) {
    SFBitsetAuxInterface* ptr_sfbitset_aux = this->GetPtrToSFBitsetAux();
    if (static_cast<DefInt>(vec_ptr_grid_info_.size()) != k0MaxLevel_ + 1) {
        LogManager::LogError("Number of grid refinement level " + std::to_string(
            vec_ptr_grid_info_.size() - 1) + " is different from k0MaxLevel_ "
            + std::to_string(k0MaxLevel_) + ", need to create grid instance "
            + "before calling the function: GenerateGridFromHighToLowLevelSerial"
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    // generate tracking and ghost nodes based on geometries
    DefInt i_geo = 0;
    for (auto& iter : vec_geo_info) {
        iter->FindTrackingNodeBasedOnGeo(*ptr_sfbitset_aux, vec_ptr_grid_info_.at(iter->GetLevel()).get());
        ++i_geo;
    }

    std::vector<DefReal> flood_fill_origin(k0GridDims_);
    std::map<std::pair<ECriterionType, DefInt>, DefMap<DefInt>> innermost_layer, outermost_layer;
    struct NumExtendLayer {
        std::vector<DefAmrLUint> neg;
        std::vector<DefAmrLUint> pos;
    } num_extend_layer;
    num_extend_layer.neg = std::vector<DefAmrLUint>(k0GridDims_, 0);
    num_extend_layer.pos = std::vector<DefAmrLUint>(k0GridDims_, 0);
    /** @todo extended outer layers based on each outmost layer could be 
    stored in individual map, which will avoid interfering between different
    layers but is more time consuming since additional nodes for every 
    outmost layers are needed to identify inside and outside region. Here
    nodes are inserted in one map (ptr_sfbitset_one_lower_level) rather
    than several ones, thus processes of extending inner and outer
    are the same*/
    for (DefInt i_level = k0MaxLevel_; i_level > 0; --i_level) {
        std::map<std::pair<ECriterionType, DefInt>, NumExtendLayer>
            map_num_extend_inner_layer, map_num_extend_outer_layer;
        std::map<std::pair<ECriterionType, DefInt>, DefMap<DefInt>> innermost_layer_current, outermost_layer_current;
        // add nodes based on tracking nodes
        for (auto& iter_tracking_grid_info : vec_ptr_grid_info_.at(i_level)
            ->map_ptr_tracking_grid_info_) {
            DefMap<DefInt> node_near_tracking;
            // set number of extended layers
            SetNumberOfExtendLayerForGrid(i_level,
                *(vec_geo_info.at(iter_tracking_grid_info.first.second)),
                &iter_tracking_grid_info.second->k0ExtendInnerNeg_,
                &iter_tracking_grid_info.second->k0ExtendInnerPos_,
                &iter_tracking_grid_info.second->k0ExtendOuterNeg_,
                &iter_tracking_grid_info.second->k0ExtendOuterPos_);

            for (DefInt idims = 0; idims < k0GridDims_; ++idims) {
                num_extend_layer.neg[idims] = iter_tracking_grid_info.second->
                    k0ExtendOuterNeg_[idims] - 1;
                num_extend_layer.pos[idims] = iter_tracking_grid_info.second->
                    k0ExtendOuterPos_[idims] - 1;
            }
            map_num_extend_outer_layer.insert(
                { iter_tracking_grid_info.first, num_extend_layer });
            // tracking nodes need to identify inside and outside

            GenerateGridNodeNearTrackingNode(
                i_level, iter_tracking_grid_info.first, &node_near_tracking);
            if (iter_tracking_grid_info.second->grid_extend_type_
                == EGridExtendType::kInAndOut) {
                for (DefInt idims = 0; idims < k0GridDims_; ++idims) {
                    num_extend_layer.neg[idims] = iter_tracking_grid_info.second->
                        k0ExtendInnerNeg_[idims] - 1;
                    num_extend_layer.pos[idims] = iter_tracking_grid_info.second->
                        k0ExtendInnerPos_[idims] - 1;
                }
                map_num_extend_inner_layer.insert(
                    { iter_tracking_grid_info.first, num_extend_layer });
                DefMap<DefInt> node_inside, node_outside;
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
            for (const auto& iter_tracking : node_near_tracking) {
                ptr_sfbitset_one_lower_level->at(i_level).insert({iter_tracking.first, kFlagSize0_});
            }
        }
        // extend the inner layer (number of extended layer
        // can be identified for each tracking grid)
        InterfaceLayerInfo* ptr_interface_info = nullptr;
        for (auto& iter_inner : innermost_layer) {
            if (map_num_extend_inner_layer.find(iter_inner.first)
                == map_num_extend_inner_layer.end()) {
                vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_
                    .insert({ iter_inner.first, std::make_shared<InterfaceLayerInfo>() });
                ptr_interface_info = vec_ptr_grid_info_.at(i_level)
                    ->map_ptr_interface_layer_info_.at(iter_inner.first).get();
                SetNumberOfExtendLayerForGrid(i_level,
                    *(vec_geo_info.at(iter_inner.first.second)),
                    &ptr_interface_info->k0ExtendInnerNeg_,
                    &ptr_interface_info->k0ExtendInnerPos_,
                    &ptr_interface_info->k0ExtendOuterNeg_,
                    &ptr_interface_info->k0ExtendOuterPos_);
                num_extend_layer.neg =
                    ptr_interface_info->k0ExtendInnerNeg_;
                num_extend_layer.pos =
                    ptr_interface_info->k0ExtendInnerPos_;
                map_num_extend_inner_layer.insert(
                    { iter_inner.first, num_extend_layer });
                num_extend_layer.neg =
                    ptr_interface_info->k0ExtendOuterNeg_;
                num_extend_layer.pos =
                    ptr_interface_info->k0ExtendOuterPos_;
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
                    &ptr_interface_info->k0ExtendInnerNeg_,
                    &ptr_interface_info->k0ExtendInnerPos_,
                    &ptr_interface_info->k0ExtendOuterNeg_,
                    &ptr_interface_info->k0ExtendOuterPos_);
                num_extend_layer.neg = ptr_interface_info->k0ExtendInnerNeg_;
                num_extend_layer.pos = ptr_interface_info->k0ExtendInnerPos_;
                map_num_extend_inner_layer.insert({ iter_outer.first, num_extend_layer });
                num_extend_layer.neg = ptr_interface_info->k0ExtendOuterNeg_;
                num_extend_layer.pos = ptr_interface_info->k0ExtendOuterPos_;
                map_num_extend_outer_layer.insert({ iter_outer.first, num_extend_layer });
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

        // find interface between different grids
        DefInt level_low = i_level - 1;
        const DefInt num_coarse2fine_layer = vec_ptr_grid_info_.at(level_low)->GetNumCoarse2FineLayer();
        if (num_coarse2fine_layer < 2) {
        LogManager::LogError("Number of coarse to fine layers should be greater than 2, which is"
            + std::to_string(num_coarse2fine_layer) + "at refinement level " + std::to_string(level_low));
        }
        for (auto& iter_inner : innermost_layer_current) {
            innermost_layer.at(iter_inner.first).clear();
            if (vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                .find(iter_inner.first) == vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.end()) {
                vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                    .insert({ iter_inner.first, std::make_shared<InterfaceLayerInfo>() });
            }
            ptr_interface_info = vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.at(iter_inner.first).get();
            ptr_interface_info->vec_inner_coarse2fine_.resize(num_coarse2fine_layer);
            ptr_interface_info->vec_outer_coarse2fine_.resize(num_coarse2fine_layer);
            FindOutmostLayerForFineGrid(i_level, iter_inner.second,
                &ptr_sfbitset_one_lower_level->at(i_level),
                &ptr_interface_info->vec_inner_coarse2fine_.at(0),
                &ptr_sfbitset_one_lower_level->at(level_low),
                &innermost_layer.at(iter_inner.first));
        }
        for (auto& iter_outer : outermost_layer_current) {
            outermost_layer.at(iter_outer.first).clear();
            if (vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                .find(iter_outer.first) == vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.end()) {
                vec_ptr_grid_info_.at(level_low)->map_ptr_interface_layer_info_
                    .insert({ iter_outer.first, std::make_shared<InterfaceLayerInfo>() });
            }
            ptr_interface_info = vec_ptr_grid_info_.at(level_low)
                ->map_ptr_interface_layer_info_.at(iter_outer.first).get();
            ptr_interface_info->vec_inner_coarse2fine_.resize(num_coarse2fine_layer);
            ptr_interface_info->vec_outer_coarse2fine_.resize(num_coarse2fine_layer);
            FindOutmostLayerForFineGrid(i_level, iter_outer.second,
                &ptr_sfbitset_one_lower_level->at(i_level),
                &ptr_interface_info->vec_outer_coarse2fine_.at(0),
                &ptr_sfbitset_one_lower_level->at(level_low),
                &outermost_layer.at(iter_outer.first));
        }
    }
}

/**
* @brief   function to add a given number of layers.
* @param[in]  i_level   level of refinement.
* @param[in]  num_extend_neg number of extended layers in negative direction.
* @param[in]  num_extend_pos number of extended layers in negative direction.
* @param[in]  map_start_layer the layer to be extended.
* @param[out]  ptr_map_exist pointer to all the nodes.
* @param[out]  ptr_map_outmost pointer to nodes on the outmost layer.
*/
void GridManagerInterface::ExtendGivenNumbOfLayer(
    const DefInt i_level, const std::vector<DefAmrLUint> num_extend_neg,
    const std::vector<DefAmrLUint> num_extend_pos,
    const DefMap<DefInt>& map_start_layer,
    DefMap<DefInt>* const ptr_map_exist,
    DefMap<DefInt>* const ptr_map_outmost) const {
    std::vector<DefSFBitset> vec_bitset_min(k0GridDims_, 0),
        vec_bitset_max(k0GridDims_, 0);
    // space filling code is at one level lower (ilevel -1)
    ComputeSFBitsetOnBoundaryAtGivenLevel(
        i_level - 1, &vec_bitset_min, &vec_bitset_max);
    // calculate layers need to be extended in each direction
    if (static_cast<DefInt>(num_extend_neg.size()) != k0GridDims_) {
        LogManager::LogError("Dimension of num_extend_neg at" + std::to_string(i_level)
            + " should be " + std::to_string(k0GridDims_)
            + " rather than " + std::to_string(num_extend_neg.size())
            + " in GridManagerInterface::ExtendGivenNumbOfLayer in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (static_cast<DefInt>(num_extend_pos.size()) != k0GridDims_) {
        LogManager::LogError("Dimension of num_extend_pos at" + std::to_string(i_level)
            + " should be " + std::to_string(k0GridDims_)
            + " rather than " + std::to_string(num_extend_pos.size())
            + " in GridManagerInterface::ExtendGivenNumbOfLayer in t"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    // find the maximum layers need to be extended
    DefAmrLUint extend_layer_max = 0;
    for (DefInt i = 0; i < k0GridDims_; ++i) {
        if (num_extend_neg[i] > extend_layer_max) {
            extend_layer_max = num_extend_neg[i];
        }
        if (num_extend_pos[i] > extend_layer_max) {
            extend_layer_max = num_extend_pos[i];
        }
    }
    std::vector<bool> bool_extend_neg(k0GridDims_, true),
        bool_extend_pos(k0GridDims_, true);
    DefMap<DefInt> map_input_layer(map_start_layer);
    std::vector<DefMap<DefInt>> vec_boundary_min, vec_boundary_max;
    std::vector<DefAmrLUint> num_extend_neg_tmp(num_extend_neg),
        num_extend_pos_tmp(num_extend_pos);
    // extend grid layer by layer
    for (DefAmrLUint i_layer = 0; i_layer < extend_layer_max; ++i_layer) {
        DefMap<DefInt> map_output_layer;
        ExtendOneLayerGrid(map_input_layer, vec_bitset_min, vec_bitset_max,
            bool_extend_neg, bool_extend_pos, &map_output_layer, ptr_map_exist,
            ptr_map_outmost, &vec_boundary_min, &vec_boundary_max);
        map_input_layer = std::move(map_output_layer);
        for (DefInt i = 0; i < k0GridDims_; ++i) {
            if (bool_extend_neg[i]) {
                num_extend_neg_tmp[i] -= 1;
                if (num_extend_neg_tmp[i] == 0) {
                    bool_extend_neg[i] = false;
                }
            }
            if (bool_extend_pos[i]) {
                num_extend_pos_tmp[i] -= 1;
                if (num_extend_pos_tmp[i] == 0) {
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
* @brief   function to do flood fill on existing nodes around geometry points (only two layers)
* @param[in]  sfbitset_inside  space filling code of a node inside geometry
* @param[in]  map_nodes_exist nodes exist for flood fill.
* @param[out]  ptr_map_nodes_outside  nodes outside the geometry.
* @param[out]  ptr_map_nodes_inside nodes inside the geometry (colored by flood fill method).
*/
int GridManagerInterface::FloodFillForInAndOut(
    const DefSFBitset& sfbitset_inside,
    const DefMap<DefInt>&  map_nodes_exist,
    DefMap<DefInt>* const ptr_map_nodes_outside,
    DefMap<DefInt>* const ptr_map_nodes_inside) const {

    std::vector<DefSFBitset> bitset_neigh;
    if (map_nodes_exist.find(sfbitset_inside) != map_nodes_exist.end()) {
        // node inside the geometry is in map_nodes_exist, then
        // nodes which is lack one or more neighbors are identified as outside nodes
        std::vector<DefSFBitset> bitset_neigh;
        for (auto iter = map_nodes_exist.begin();
            iter != map_nodes_exist.end(); ++iter) {
            GridFindAllNeighborsVir(iter->first, &bitset_neigh);
            for (const auto& iter_neigh : bitset_neigh) {
                if (map_nodes_exist.find(iter_neigh)
                    == map_nodes_exist.end()) {
                    ptr_map_nodes_outside->insert({ iter->first, kFlag0_ });
                    break;
                }
            }
        }
        return -1;
    } else {
        // create a copy of map_nodes_exist and add one layer for flood fill
        DefInt flag_bit_colored = DefInt(1 << (sizeof(DefInt) * 8 - 1));
        DefInt flag_bit_exist = flag_bit_colored | (sizeof(DefInt) * 8 - 2);
        DefMap<DefInt> map_nodes_tmp(map_nodes_exist);
        // generate one more layer around the input grid
        for (auto iter = map_nodes_exist.begin();
            iter != map_nodes_exist.end(); ++iter) {
            map_nodes_tmp.at(iter->first) = flag_bit_exist;
            GridFindAllNeighborsVir(iter->first, &bitset_neigh);
            for (const auto& iter_neighbour : bitset_neigh) {
                if (map_nodes_tmp.find(iter_neighbour)
                    == map_nodes_tmp.end()) {
                    map_nodes_tmp.insert({ iter_neighbour, kFlagSize0_ });
                }
            }
        }
        // flood fill
        std::vector<DefSFBitset> vec_sfbitset_stk, vec_bitset(k0NumNeighbors_);
        DefInt i = 0;
        DefSFBitset sfbitset_seed;
        vec_sfbitset_stk.emplace_back(sfbitset_inside);
        while (!vec_sfbitset_stk.empty() && i < K0IMaxFloodFill_) {
            sfbitset_seed = vec_sfbitset_stk.back();
            vec_sfbitset_stk.pop_back();
            ++i;
            if (map_nodes_tmp.find(sfbitset_seed) == map_nodes_tmp.end()) {
            } else if ((map_nodes_tmp.at(sfbitset_seed) & flag_bit_exist)
                == flag_bit_exist) {
                ptr_map_nodes_inside->insert({ sfbitset_seed, kFlagSize0_ });
            } else if ((map_nodes_tmp.at(sfbitset_seed) & flag_bit_colored)
                != flag_bit_colored) {
                // color the node
                map_nodes_tmp.at(sfbitset_seed) |= flag_bit_colored;
                // add neighboring nodes to seed
                PushBackSFBitsetInFloodFill(sfbitset_seed, &vec_sfbitset_stk);
            }
        }
        if (i == K0IMaxFloodFill_) {
            return 2;
        }
        std::vector<DefSFBitset> bitset_neigh;
        for (auto iter = map_nodes_exist.begin();
            iter != map_nodes_exist.end(); ++iter) {
            if (ptr_map_nodes_inside->find(iter->first)
                == ptr_map_nodes_inside->end()) {
                GridFindAllNeighborsVir(iter->first, &bitset_neigh);
                for (const auto& iter_neigh : bitset_neigh) {
                    if (map_nodes_exist.find(iter_neigh)
                        == map_nodes_exist.end()) {
                        ptr_map_nodes_outside->insert({ iter->first, kFlagSize0_ });
                        break;
                    }
                }
            }
        }
        return 0;
    }
}
/**
* @brief   function to reset background grid as mpi outer layer near coarse to fine interface.
* @param[in]  num_mpi_layers  number of mpi outer layers.
* @param[in]  background_grid_info class of grid information for background grid.
* @param[in]  level1_grid  nodes in level 1 grid.
* @param[out]  ptr_background_grid pointer to nodes in level 0 grid (background).
* @note nodes in grid of level 1 (on fine to coarse interface) is not included in background grid.
*/
void GridManagerInterface::ResetBackgroundGridAsMpiLayer(const DefInt num_mpi_layers,
    const GridInfoInterface& background_grid_info, const DefMap<DefInt>& level1_grid,
    DefMap<DefInt>* const ptr_background_grid) {
    SFBitsetAuxInterface* ptr_sfbitset_aux = GetPtrToSFBitsetAux();
    if (k0MaxLevel_ > 0) {
        ptr_background_grid->clear();
        std::vector<DefSFBitset> domain_min_m1_n_level(k0GridDims_), domain_max_p1_n_level(k0GridDims_), vec_region;
        ptr_sfbitset_aux->GetMinM1AtGivenLevel(
            0, ptr_sfbitset_aux->GetMinBackgroundIndices(), &domain_min_m1_n_level);
        ptr_sfbitset_aux->GetMaxP1AtGivenLevel(
            0, ptr_sfbitset_aux->GetMaxBackgroundIndices(), &domain_max_p1_n_level);
        for (const auto& iter_interface : background_grid_info.map_ptr_interface_layer_info_) {
            for (const auto& iter_node : iter_interface.second->vec_inner_coarse2fine_.at(0)) {
                ptr_sfbitset_aux->FindNodesInReginOfGivenLength(iter_node.first, num_mpi_layers,
                    domain_min_m1_n_level, domain_max_p1_n_level, &vec_region);
                for (const auto& iter_region : vec_region) {
                    if (level1_grid.find(iter_region) == level1_grid.end()) {
                        ptr_background_grid->insert({iter_region, kFlag0_});
                    }
                }
            }
            for (const auto& iter_node : iter_interface.second->vec_outer_coarse2fine_.at(0)) {
                ptr_sfbitset_aux->FindNodesInReginOfGivenLength(iter_node.first, num_mpi_layers,
                    domain_min_m1_n_level, domain_max_p1_n_level, &vec_region);
                for (const auto& iter_region : vec_region) {
                    if (level1_grid.find(iter_region) == level1_grid.end()) {
                        ptr_background_grid->insert({iter_region, kFlag0_});
                    }
                }
            }
        }
    }

}
}  // end namespace amrproject
}  // end namespace rootproject
