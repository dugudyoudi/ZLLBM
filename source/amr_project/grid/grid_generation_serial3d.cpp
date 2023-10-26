//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_generation_serial3d.cpp
* @author Zhengliang Liu
* @brief functions used to generate grid serially.
* @date  2022-11-28
* @note  functions from geometry_manager will be called.
*/
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to add one layer near the tracking grid.
* @param[in] i_level refinement level.
* @param[in] tracking_grid_key key of the tracking grid used for
*                  generating one layer of grid.
* @param[out] ptr_map_node_temp  nodes near the tracking grid.
*/
void GridManager3D::GenerateGridNodeNearTrackingNode(const DefAmrIndexUint i_level,
    const std::pair<ECriterionType, DefAmrIndexUint>& tracking_grid_key,
    DefMap<DefAmrUint>* const ptr_map_node_temp) const {
    DefSFBitset bitset_lower_level;
    std::array<DefSFBitset, 27> bitset_of_a_cell;
    for (const auto& iter : this->vec_ptr_grid_info_.at(i_level)
        ->map_ptr_tracking_grid_info_.at(tracking_grid_key)
        ->map_tracking_node_) {
        bitset_lower_level = SFBitsetToOneLowerLevel(iter.first);
        SFBitsetFindAllNeighbors(bitset_lower_level, &bitset_of_a_cell);
        for (const auto& iter_bitset : bitset_of_a_cell) {
            if (ptr_map_node_temp->find(iter_bitset)
                == ptr_map_node_temp->end()) {
                ptr_map_node_temp->insert({ iter_bitset, kFlag0_ });
            }
        }
    }
}
/**
* @brief   function to add layers near tracking grids.
* @param[in] i_level refinement level.
* @param[in]  i_geo index of the geometry (only for write log).
* @param[in]  flood_fill_start_point starting point to do flood fill.
* @param[in]  map_nodes_exist   existing nodes.
* @param[out] ptr_map_nodes_outside   nodes haven't been colored.
* @param[out] ptr_map_nodes_inside   nodes have been colored.
*/
void GridManager3D::IdentifyTypeOfLayerByFloodFill(
    const DefAmrIndexUint i_level, const DefAmrIndexUint i_geo,
    const std::vector<DefReal> flood_fill_start_point,
    const DefMap<DefAmrUint>& map_nodes_exist,
    DefMap<DefAmrUint>* const ptr_map_nodes_outside,
    DefMap<DefAmrUint>* const ptr_map_nodes_inside) const {
    // step 1: find start point for flood fill
    if (flood_fill_start_point.size() != k0GridDims_) {
        LogManager::LogWarning("Dimension of flood_fill_start_point is different from k0GridDims_.");
    }
    std::array<DefReal, 3> flood_fill_origin =
    { flood_fill_start_point[kXIndex], flood_fill_start_point[kYIndex],
        flood_fill_start_point[kZIndex] };

    // calculate bounds of searching step based on domain boundary
    DefAmrIndexLUint scale_i_level = static_cast<DefAmrIndexLUint>(TwoPowerN(i_level));
    DefAmrIndexLUint x_index = static_cast<DefAmrIndexLUint>(flood_fill_origin.at(kXIndex)
        / (k0DomainDx_[kXIndex] / scale_i_level) + kEps);
    DefAmrIndexLUint y_index = static_cast<DefAmrIndexLUint>(flood_fill_origin.at(kYIndex)
        / (k0DomainDx_[kYIndex] / scale_i_level) + kEps);
    DefAmrIndexLUint z_index = static_cast<DefAmrIndexLUint>(flood_fill_origin.at(kZIndex)
        / (k0DomainDx_[kZIndex] / scale_i_level) + kEps);
    DefSizet x_index_max =
        k0MaxIndexOfBackgroundNode_[kXIndex] * scale_i_level - x_index;
    DefSizet y_index_max =
        k0MaxIndexOfBackgroundNode_[kYIndex] * scale_i_level - y_index;
    DefSizet z_index_max =
        k0MaxIndexOfBackgroundNode_[kZIndex] * scale_i_level - z_index;

    bool bool_find_node_for_flood_fill = false;
    DefSizet i_count = 0, count_sum = 0;;
    DefSFBitset sfbitset_origin_vertex = SFBitsetEncoding(
        std::array<DefAmrIndexLUint, 3>({ x_index , y_index, z_index }));
    // search in -x direction from the vec_origin until meet the first vertex
    // in map_nodes_exist
    DefSFBitset sfbitset_temp = sfbitset_origin_vertex, sfbitset_start_vertex;
    while (i_count < x_index) {
        sfbitset_temp = FindXNeg(sfbitset_temp);
        if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
            sfbitset_start_vertex = FindXPos(sfbitset_temp);
            bool_find_node_for_flood_fill = true;
            break;
        }
        ++i_count;
    }
    count_sum += i_count;
    // search in -y direction from the vec_origin until meet the first vertex
    // in map_nodes_exist
    if (!bool_find_node_for_flood_fill) {
        sfbitset_temp = sfbitset_origin_vertex;
        i_count = 0;
        while (i_count < y_index) {
            sfbitset_temp = FindYNeg(sfbitset_temp);
            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
                sfbitset_start_vertex = FindYPos(sfbitset_temp);
                bool_find_node_for_flood_fill = true;
                break;
            }
            ++i_count;
        }
    }
    count_sum += i_count;
    // search in -z direction from the vec_origin until meet the first vertex
    // in map_nodes_exist
    if (!bool_find_node_for_flood_fill) {
        sfbitset_temp = sfbitset_origin_vertex;
        i_count = 0;
        while (i_count < z_index) {
            sfbitset_temp = FindZNeg(sfbitset_temp);
            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
                sfbitset_start_vertex = FindZPos(sfbitset_temp);
                bool_find_node_for_flood_fill = true;
                break;
            }
            ++i_count;
        }
    }
    count_sum += i_count;
    // search in +x direction from the vec_origin until meet the first vertex
    // in map_nodes_exist
    if (!bool_find_node_for_flood_fill) {
        sfbitset_temp = sfbitset_origin_vertex;
        i_count = 0;
        while (i_count < x_index_max) {
            sfbitset_temp = FindXPos(sfbitset_temp);
            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
                sfbitset_start_vertex = FindXNeg(sfbitset_temp);
                bool_find_node_for_flood_fill = true;
                break;
            }
            ++i_count;
        }
    }
    count_sum += i_count;
    // search in +y direction from the vec_origin until meet the first vertex
    // in map_nodes_exist
    if (!bool_find_node_for_flood_fill) {
        sfbitset_temp = sfbitset_origin_vertex;
        i_count = 0;
        while (i_count < y_index_max) {
            sfbitset_temp = FindYPos(sfbitset_temp);
            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
                sfbitset_start_vertex = FindYNeg(sfbitset_temp);
                bool_find_node_for_flood_fill = true;
                break;
            }
            ++i_count;
        }
    }
    count_sum += i_count;
    // search in +z direction from the vec_origin until meet the first vertex
    // in map_nodes_exist
    if (!bool_find_node_for_flood_fill) {
        sfbitset_temp = sfbitset_origin_vertex;
        i_count = 0;
        while (i_count < z_index_max) {
            sfbitset_temp = FindZPos(sfbitset_temp);
            if (map_nodes_exist.find(sfbitset_temp) != map_nodes_exist.end()) {
                sfbitset_start_vertex = FindZNeg(sfbitset_temp);
                bool_find_node_for_flood_fill = true;
                break;
            }
            ++i_count;
        }
    }
    count_sum += i_count;
    if (map_nodes_exist.find(sfbitset_origin_vertex) != map_nodes_exist.end()) {
        sfbitset_start_vertex = sfbitset_origin_vertex;
        LogManager::LogInfo("input node for flood fill coincideS with existing nodes near geometry("
         + std::to_string(i_geo) + "), may not distinguish inside and outside");
    }

    if (bool_find_node_for_flood_fill) {
        int flag_floodfill;
        flag_floodfill = FloodFillForInAndOut(sfbitset_start_vertex,
            map_nodes_exist, ptr_map_nodes_outside, ptr_map_nodes_inside);
       if (flag_floodfill == 2) {
            LogManager::LogWarning("Iteration of flood fill exceed preset limits for."
            " geometry (" + std::to_string(i_geo) + ").");
        }
    } else {
        LogManager::LogError("Can't find starting node for food fill in geometry: "
            + std::to_string(i_geo) + " in IdentifyTypeOfLayerByFloodFill(2D)"
            + " after " + std::to_string(count_sum) + " iterations"
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief   function to add layers near tracking grids.
* @param[in] sfbitset_in space filling code of at the center.
* @param[out] ptr_vec_stk adjecent nodes in positive and negative directions.
*/
void GridManager3D::PushBackSFBitsetInFloodFill(const DefSFBitset& sfbitset_in,
    std::vector<DefSFBitset>* const ptr_vec_stk) const {
    std::array<DefSFBitset, 27> array_neighbors;
    SFBitsetFindAllNeighbors(sfbitset_in, &array_neighbors);
    std::vector<DefSFBitset> vec_temp(27);
    memcpy(vec_temp.data(), array_neighbors.data(),
        27 * sizeof(DefSFBitset));
    ptr_vec_stk->insert(ptr_vec_stk->end(), vec_temp.begin() + 1,
        vec_temp.end());
}
/**
* @brief   function to extend the grid by one layer.
* @param[in] map_start_layer the based layer for extension.
* @param[in] vec_bitset_min space filling code of the minimum boundary.
* @param[in] vec_bitset_max space filling code of the maximum boundary.
* @param[in] bool_extend_neg indicators of extending in negative directions.
* @param[in] bool_extend_pos indicators of extending in positive directions.
* @param[out] ptr_map_output_layer the extended layer.
* @param[out] ptr_map_exist existing nodes.
* @param[out] ptr_outmost_layer the outmost layer where bool_extend_xxx is false.
* @param[out] ptr_vector_boundary_neg nodes on the minimum boundary.
* @param[out] ptr_vector_boundary_pos nodes on the maximum boundary.
*/
void GridManager3D::ExtendOneLayerGrid(
    const DefMap<DefAmrUint>& map_start_layer,
    const std::vector<DefSFBitset>& vec_bitset_min,
    const std::vector<DefSFBitset>& vec_bitset_max,
    const std::vector<bool>& bool_extend_neg,
    const std::vector<bool>& bool_extend_pos,
    DefMap<DefAmrUint>* const ptr_map_output_layer,
    DefMap<DefAmrIndexUint>* const ptr_map_exist,
    DefMap<DefAmrUint>* const ptr_outmost_layer,
    std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_min,
    std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_max) const {
    std::array<bool, 3> bool_neg, bool_pos,
        bool_neg_boundary, bool_pos_boundary;
    const DefMap<DefAmrIndexUint> map_exist_temp(*ptr_map_exist);
    ptr_vector_boundary_min->resize(k0GridDims_);
    ptr_vector_boundary_max->resize(k0GridDims_);
    DefAmrUint flag_node_boundary;
    std::vector<DefSFBitset> vec_neighbors;
    DefMap<DefAmrUint> outmost_temp;
    for (const auto& iter : map_start_layer) {
        bool_neg_boundary[kXIndex]
            = (iter.first & k0SFBitsetTakeXRef_[kRefCurrent_])
            == vec_bitset_min[kXIndex];
        bool_pos_boundary[kXIndex]
            = (iter.first & k0SFBitsetTakeXRef_[kRefCurrent_])
            == vec_bitset_max[kXIndex];
        bool_neg_boundary[kYIndex]
            = (iter.first & k0SFBitsetTakeYRef_[kRefCurrent_])
            == vec_bitset_min[kYIndex];
        bool_pos_boundary[kYIndex]
            = (iter.first & k0SFBitsetTakeYRef_[kRefCurrent_])
            == vec_bitset_max[kYIndex];
        bool_neg_boundary[kZIndex]
            = (iter.first & k0SFBitsetTakeZRef_[kRefCurrent_])
            == vec_bitset_min[kZIndex];
        bool_pos_boundary[kZIndex]
            = (iter.first & k0SFBitsetTakeZRef_[kRefCurrent_])
            == vec_bitset_max[kZIndex];
        bool_neg[kXIndex] = (!bool_neg_boundary[kXIndex])
            && bool_extend_neg[kXIndex];
        bool_pos[kXIndex] = (!bool_pos_boundary[kXIndex])
            && bool_extend_pos[kXIndex];
        bool_neg[kYIndex] = (!bool_neg_boundary[kYIndex])
            && bool_extend_neg[kYIndex];
        bool_pos[kYIndex] = (!bool_pos_boundary[kYIndex])
            && bool_extend_pos[kYIndex];
        bool_neg[kZIndex] = (!bool_neg_boundary[kZIndex])
            && bool_extend_neg[kZIndex];
        bool_pos[kZIndex] = (!bool_pos_boundary[kZIndex])
            && bool_extend_pos[kZIndex];
        flag_node_boundary = FindAllNeighborsWithSpecifiedDirection(
            iter.first, bool_neg, bool_pos, &vec_neighbors);
        for (const auto& iter_neighbour : vec_neighbors) {
            if (ptr_map_exist->find(iter_neighbour) == ptr_map_exist->end()) {
                ptr_map_output_layer->insert({ iter_neighbour, kFlag0_ });
                ptr_map_exist->insert({ iter_neighbour, kFlagSize0_ });
            }
        }
        if ((flag_node_boundary & kFlagCurrentNodeXNeg_)
            == kFlagCurrentNodeXNeg_) {
            if (bool_neg_boundary[kXIndex]) {  // on the minimum x boundary
                ptr_vector_boundary_min->at(kXIndex)
                    .insert({ iter.first, kFlag0_ });
            } else {
                outmost_temp.insert({ iter.first, flag_node_boundary });
            }
        }
        if ((flag_node_boundary & kFlagCurrentNodeXPos_)
            == kFlagCurrentNodeXPos_) {
            if (bool_pos_boundary[kXIndex]) {  // on the maximum x boundary
                ptr_vector_boundary_max->at(kXIndex)
                    .insert({ iter.first, kFlag0_ });
            } else {
                outmost_temp.insert({ iter.first, flag_node_boundary });
            }
        }
        if ((flag_node_boundary & kFlagCurrentNodeYNeg_)
            == kFlagCurrentNodeYNeg_) {
            if (bool_neg_boundary[kYIndex]) {  // on the minimum y boundary
                ptr_vector_boundary_min->at(kYIndex)
                    .insert({ iter.first, kFlag0_ });
            } else {
                outmost_temp.insert({ iter.first, flag_node_boundary });
            }
        }
        if ((flag_node_boundary & kFlagCurrentNodeYPos_)
            == kFlagCurrentNodeYPos_) {
            if (bool_pos_boundary[kYIndex]) {  // on the maximum y boundary
                ptr_vector_boundary_max->at(kYIndex)
                    .insert({ iter.first, kFlag0_ });
            } else {
                outmost_temp.insert({ iter.first, flag_node_boundary });
            }
        }
        if ((flag_node_boundary & kFlagCurrentNodeZNeg_)
            == kFlagCurrentNodeZNeg_) {
            if (bool_neg_boundary[kZIndex]) {  // on the maximum z boundary
                ptr_vector_boundary_min->at(kZIndex)
                    .insert({ iter.first, kFlag0_ });
            } else {
                outmost_temp.insert({ iter.first, flag_node_boundary });
            }
        }
        if ((flag_node_boundary & kFlagCurrentNodeZPos_)
            == kFlagCurrentNodeZPos_) {
            if (bool_pos_boundary[kZIndex]) {  // on the maximum z boundary
                ptr_vector_boundary_max->at(kZIndex)
                    .insert({ iter.first, kFlag0_ });
            } else {
                outmost_temp.insert({ iter.first, flag_node_boundary });
            }
        }
    }
    // check if nodes on the outmost layer
    for (const auto iter : outmost_temp) {
        bool_neg[kXIndex]
            = !((iter.first & k0SFBitsetTakeXRef_[kRefCurrent_])
                == vec_bitset_min[kXIndex]);
        bool_pos[kXIndex]
            = !((iter.first & k0SFBitsetTakeXRef_[kRefCurrent_])
                == vec_bitset_max[kXIndex]);
        bool_neg[kYIndex]
            = !((iter.first & k0SFBitsetTakeYRef_[kRefCurrent_])
                == vec_bitset_min[kYIndex]);
        bool_pos[kYIndex]
            = !((iter.first & k0SFBitsetTakeYRef_[kRefCurrent_])
                == vec_bitset_max[kYIndex]);
        bool_neg[kZIndex]
            = !((iter.first & k0SFBitsetTakeZRef_[kRefCurrent_])
                == vec_bitset_min[kZIndex]);
        bool_pos[kZIndex]
            = !((iter.first & k0SFBitsetTakeZRef_[kRefCurrent_])
                == vec_bitset_max[kZIndex]);
        flag_node_boundary = FindAllNeighborsWithSpecifiedDirection(
            iter.first, bool_neg, bool_pos, &vec_neighbors);
        for (const auto& iter_neighbour : vec_neighbors) {
            if (ptr_map_exist->find(iter_neighbour) == ptr_map_exist->end()) {
                ptr_outmost_layer->insert({ iter.first, kFlag0_ });
                break;
            }
        }
    }
}
/**
* @brief   function to compute space filling code of domain boundary at given level from the background level.
* @param[in] i_level refinement level.
* @param[out] ptr_vec_bitset_min pointer to spacing filling code of the minimum domain boundary.
* @param[out] ptr_vec_bitset_max pointer to spacing filling code of the minimum maximum boundary.
*/
void GridManager3D::ComputeSFBitsetOnBoundaryAtGivenLevel(
    const DefAmrIndexUint i_level,
    std::vector<DefSFBitset>* const ptr_vec_bitset_min,
    std::vector<DefSFBitset>* const ptr_vec_bitset_max)  const {
    ptr_vec_bitset_min->at(kXIndex) = SFBitsetToNHigherLevel(
        i_level, SFBitsetMin_[kXIndex]);
    ptr_vec_bitset_min->at(kYIndex) = SFBitsetToNHigherLevel(
        i_level, SFBitsetMin_[kYIndex]);
    ptr_vec_bitset_min->at(kZIndex) = SFBitsetToNHigherLevel(
        i_level, SFBitsetMin_[kZIndex]);
    ptr_vec_bitset_max->at(kXIndex) = SFBitsetToNHigherLevel(
        i_level, SFBitsetMax_[kXIndex]);
    ptr_vec_bitset_max->at(kYIndex) = SFBitsetToNHigherLevel(
        i_level, SFBitsetMax_[kYIndex]);
    ptr_vec_bitset_max->at(kZIndex) = SFBitsetToNHigherLevel(
        i_level, SFBitsetMax_[kZIndex]);
}
/**
* @brief   function to compute space filling code of neighbors in given directions.
* @param[in] bitset_in space filling code of the center node.
* @param[in] bool_neg indicators of extending in negative directions.
* @param[in] bool_pos indicators of extending in positive directions.
* @param[out] ptr_vec_neigbours space filling codes of neighbors.
* @return indicators of current node on boundaries.
*/
DefAmrUint GridManager3D::FindAllNeighborsWithSpecifiedDirection(
    const DefSFBitset bitset_in,
    const std::array<bool, 3>& bool_neg, const std::array<bool, 3>& bool_pos,
    std::vector <DefSFBitset>* const ptr_vec_neigbours) const {
    ptr_vec_neigbours->clear();
    DefAmrUint flag_current_node = 0;
    DefSFBitset bitset_temp, bitset_temp1, bitset_temp2;
    if (bool_neg[kXIndex]) {  // (-x, 0, 0)
        bitset_temp = FindXNeg(bitset_in);
        ptr_vec_neigbours->push_back(bitset_temp);
        if (bool_neg[kYIndex]) {  // (-x, -y, 0)
            bitset_temp1 = FindYNeg(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
            if (bool_neg[kZIndex]) {  // (-x, -y, -z)
                bitset_temp2 = FindZNeg(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
            if (bool_pos[kZIndex]) {  // (-x, -y, +z)
                bitset_temp2 = FindZPos(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
        }
        if (bool_pos[kYIndex]) {  // (-x, +y, 0)
            bitset_temp1 = FindYPos(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
            if (bool_neg[kZIndex]) {  // (-x, +y, -z)
                bitset_temp2 = FindZNeg(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
            if (bool_pos[kZIndex]) {  // (-x, +y, +z)
                bitset_temp2 = FindZPos(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
        }
        if (bool_neg[kZIndex]) {  // (-x, 0, -z)
            bitset_temp1 = FindZNeg(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
        if (bool_pos[kZIndex]) {  // (-x, 0, +z)
            bitset_temp1 = FindZPos(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
    } else {
        flag_current_node |= kFlagCurrentNodeXNeg_;
    }
    if (bool_pos[kXIndex]) {  // (+x, 0, 0)
        bitset_temp = FindXPos(bitset_in);
        ptr_vec_neigbours->push_back(bitset_temp);
        if (bool_neg[kYIndex]) {  // (+x, -y, 0)
            bitset_temp1 = FindYNeg(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
            if (bool_neg[kZIndex]) {  // (+x, -y, -z)
                bitset_temp2 = FindZNeg(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
            if (bool_pos[kZIndex]) {  // (+x, -y, +z)
                bitset_temp2 = FindZPos(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
        }
        if (bool_pos[kYIndex]) {  // (+x, +y, 0)
            bitset_temp1 = FindYPos(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
            if (bool_neg[kZIndex]) {  // (+x, +y, -z)
                bitset_temp2 = FindZNeg(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
            if (bool_pos[kZIndex]) {  // (+x, +y, +z)
                bitset_temp2 = FindZPos(bitset_temp1);
                ptr_vec_neigbours->push_back(bitset_temp2);
            }
        }
        if (bool_neg[kZIndex]) {  // (+x, 0, -z)
            bitset_temp1 = FindZNeg(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
        if (bool_pos[kZIndex]) {  // (+x, 0, +z)
            bitset_temp1 = FindZPos(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
    } else {
        flag_current_node |= kFlagCurrentNodeXPos_;
    }
    if (bool_neg[kYIndex]) {  // (0, -y, 0)
        bitset_temp = FindYNeg(bitset_in);
        ptr_vec_neigbours->push_back(bitset_temp);
        if (bool_neg[kZIndex]) {  // (0, -y, -z)
            bitset_temp1 = FindZNeg(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
        if (bool_pos[kZIndex]) {  // (0, -y, +z)
            bitset_temp1 = FindZPos(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
    } else {
        flag_current_node |= kFlagCurrentNodeYNeg_;
    }
    if (bool_pos[kYIndex]) {  // (0, +y, 0)
        bitset_temp = FindYPos(bitset_in);
        ptr_vec_neigbours->push_back(bitset_temp);
        if (bool_neg[kZIndex]) {  // (0, +y, -z)
            bitset_temp1 = FindZNeg(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
        if (bool_pos[kZIndex]) {  // (0, +y, +z)
            bitset_temp1 = FindZPos(bitset_temp);
            ptr_vec_neigbours->push_back(bitset_temp1);
        }
    } else {
        flag_current_node |= kFlagCurrentNodeYPos_;
    }
    if (bool_neg[kZIndex]) {  // (0, 0, -z)
        bitset_temp = FindZNeg(bitset_in);
        ptr_vec_neigbours->push_back(bitset_temp);
    } else {
        flag_current_node |= kFlagCurrentNodeZNeg_;
    }
    if (bool_pos[kZIndex]) {  // (0, 0, +z)
        bitset_temp = FindZPos(bitset_in);
        ptr_vec_neigbours->push_back(bitset_temp);
    } else {
        flag_current_node |= kFlagCurrentNodeZPos_;
    }
    return flag_current_node;
}
/**
* @brief   function to find interface between grids of different
*          refinement levels.
* @param[in] i_level higher refinement level.
* @param[in] map_outmost_layer outmost layer of the grid at the higher
*            refinement level.
* @param[out] ptr_map_exist existing nodes at the higher refinement level.
* @param[out] ptr_interface_outmost nodes on the outmost layer.
* @param[out] ptr_layer_lower_level nodes at the lower refinement level.
* @param[out] ptr_layer_lower_level_outer outer nodes at the lower
*             refinement level.
*/
void  GridManager3D::FindInterfaceBetweenGrid(
    const DefAmrIndexUint i_level, const DefMap<DefAmrUint>& map_outmost_layer,
    DefMap<DefAmrIndexUint>* const map_exist,
    DefMap<DefAmrUint>* const ptr_interface_outmost,
    DefMap<DefAmrIndexUint>* const ptr_layer_lower_level,
    DefMap<DefAmrUint>* const ptr_layer_lower_level_outer) {
#ifdef DEBUG_CHECK_GRID
    if (&map_outmost_layer == ptr_interface_outmost) {
        LogManager::LogError("input (map_outmost_layer)"
            " should not be the same as output (ptr_interface_outmost) in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (&map_outmost_layer == ptr_layer_lower_level_outer) {
        LogManager::LogError("input (map_outmost_layer)"
            " should not be the same as output (ptr_layer_lower_level_outer) in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID
    std::vector<DefSFBitset> vec_bitset_min(k0GridDims_, 0),
        vec_bitset_max(k0GridDims_, 0);
    ComputeSFBitsetOnBoundaryAtGivenLevel(
        i_level - 1, &vec_bitset_min, &vec_bitset_max);
    std::vector<DefSFBitset> vec_bitset_min_lower(k0GridDims_, 0),
        vec_bitset_max_lower(k0GridDims_, 0);
    for (DefAmrIndexUint id = 0; id < k0GridDims_; ++id) {
        vec_bitset_min_lower[id] = SFBitsetToOneLowerLevel(vec_bitset_min[id]);
        vec_bitset_max_lower[id] = SFBitsetToOneLowerLevel(vec_bitset_max[id]);
    }
    std::array<bool, 3> bool_neg_not_boundary, bool_pos_not_boundary;
    DefAmrUint flag_node_boundary;
    std::vector<DefSFBitset> vec_neighbors, vec_lower_neighbors;
    DefSFBitset bitset_lower_level, bitset_temp, bitset_neighbor;
    for (const auto& iter : map_outmost_layer) {
        bool_neg_not_boundary[kXIndex]
            = !((iter.first & k0SFBitsetTakeXRef_[kRefCurrent_])
                == vec_bitset_min[kXIndex]);
        bool_pos_not_boundary[kXIndex]
            = !((iter.first & k0SFBitsetTakeXRef_[kRefCurrent_])
                == vec_bitset_max[kXIndex]);
        bool_neg_not_boundary[kYIndex]
            = !((iter.first & k0SFBitsetTakeYRef_[kRefCurrent_])
                == vec_bitset_min[kYIndex]);
        bool_pos_not_boundary[kYIndex]
            = !((iter.first & k0SFBitsetTakeYRef_[kRefCurrent_])
                == vec_bitset_max[kYIndex]);
        bool_neg_not_boundary[kZIndex]
            = !((iter.first & k0SFBitsetTakeZRef_[kRefCurrent_])
                == vec_bitset_min[kZIndex]);
        bool_pos_not_boundary[kZIndex]
            = !((iter.first & k0SFBitsetTakeZRef_[kRefCurrent_])
                == vec_bitset_max[kZIndex]);
        flag_node_boundary = FindAllNeighborsWithSpecifiedDirection(
            iter.first, bool_neg_not_boundary,
            bool_pos_not_boundary, &vec_neighbors);
        for (const auto& iter_neighour : vec_neighbors) {
            if (map_exist->find(iter_neighour) == map_exist->end()) {
                ptr_interface_outmost->insert({ iter.first, kFlag0_ });
                map_exist->at(iter.first) |= kNodeStatusCoarse2Fine0_;
                // find interface at lower level
                bitset_lower_level = SFBitsetToOneLowerLevel(iter.first);
                if (ptr_layer_lower_level->find(bitset_lower_level)
                    == ptr_layer_lower_level->end()) {
                    ptr_layer_lower_level->insert({ bitset_lower_level, kFlagSize0_ });
                }
                bool_neg_not_boundary[kXIndex]
                    = !((bitset_lower_level
                        & k0SFBitsetTakeXRef_[kRefCurrent_])
                        == vec_bitset_min_lower[kXIndex]);
                bool_pos_not_boundary[kXIndex]
                    = !((bitset_lower_level
                        & k0SFBitsetTakeXRef_[kRefCurrent_])
                        == vec_bitset_max_lower[kXIndex]);
                bool_neg_not_boundary[kYIndex]
                    = !((bitset_lower_level
                        & k0SFBitsetTakeYRef_[kRefCurrent_])
                        == vec_bitset_min_lower[kYIndex]);
                bool_pos_not_boundary[kYIndex]
                    = !((bitset_lower_level
                        & k0SFBitsetTakeYRef_[kRefCurrent_])
                        == vec_bitset_max_lower[kYIndex]);
                bool_neg_not_boundary[kZIndex]
                    = !((bitset_lower_level
                        & k0SFBitsetTakeZRef_[kRefCurrent_])
                        == vec_bitset_min_lower[kZIndex]);
                bool_pos_not_boundary[kZIndex]
                    = !((bitset_lower_level
                        & k0SFBitsetTakeZRef_[kRefCurrent_])
                        == vec_bitset_max_lower[kZIndex]);
                FindAllNeighborsWithSpecifiedDirection(
                    bitset_lower_level, bool_neg_not_boundary,
                    bool_pos_not_boundary, &vec_neighbors);
                bitset_temp = SFBitsetToOneHigherLevel(bitset_lower_level);
                if (map_exist->find(bitset_temp) == map_exist->end()) {
                    ptr_layer_lower_level_outer->insert({
                        bitset_lower_level, kFlag0_ });
                    for (const auto& iter_lower : vec_neighbors) {
                        bitset_neighbor = SFBitsetToOneHigherLevel(iter_lower);
                        if (map_exist->find(bitset_neighbor) != map_exist->end()) {
                            ptr_layer_lower_level->insert({ iter_lower, kFlagSize0_ });
                        }
                    }
                } else {
                    for (const auto& iter_lower : vec_neighbors) {
                        bitset_neighbor = SFBitsetToOneHigherLevel(iter_lower);
                        if (map_exist->find(bitset_neighbor) == map_exist->end()) {
                            ptr_layer_lower_level->insert({ iter_lower, kFlagSize0_ });
                            ptr_layer_lower_level_outer->insert({ iter_lower, kFlag0_ });
                        }
                    }
                }
                break;
            }
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
