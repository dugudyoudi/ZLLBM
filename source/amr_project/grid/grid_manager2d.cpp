//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
* @note  functions from other managers will be called.
*/
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
#include <string>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to setup and check grid related parameters.
*/
void GridManager2D::SetGridParameters() {
    SFBitsetMinAndMaxGlobal(
        k0IntOffset_, k0MaxIndexOfBackgroundNode_);
    SFBitsetMinAndMaxCoordinates(k0MaxLevel_,
        k0IntOffset_, k0MaxIndexOfBackgroundNode_);

    // check if length of computational domain is given
    if (k0DomainSize_.at(kXIndex) < kEps) {
        LogManager::LogError("Domain length in x direction (k0DomainSize_[0])"
            " should be a positive value");
    } else if (k0DomainSize_.at(kYIndex) < kEps) {
        LogManager::LogError("Domain length in x direction (k0DomainSize_[1])"
            " should be a positive value");
    }
    // check if grid space is given
    if (k0DomainDx_.at(kXIndex) < kEps
        && k0DomainDx_.at(kYIndex) < kEps) {
        LogManager::LogError("Grid space of x or y(k0DomainDx_)"
            " should be positive values");
    }

    // set grid space if not all grid spaces are given
    if (k0DomainDx_.at(kXIndex) < kEps) {
        k0DomainDx_.at(kXIndex) = k0DomainDx_.at(kYIndex);
    }
    if (k0DomainDx_.at(kYIndex) < kEps) {
        k0DomainDx_.at(kYIndex) = k0DomainDx_.at(kXIndex);
    }
    k0SpaceBackground_ = { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex] };

    // calculate number of background nodes in each direction
    DefAmrIndexLUint x_max = static_cast<DefAmrIndexLUint>(k0DomainSize_[kXIndex]
        / k0DomainDx_[kXIndex] + kEps);

    k0MaxIndexOfBackgroundNode_ = {
        static_cast<DefAmrIndexLUint>(k0DomainSize_[kXIndex]
        / k0DomainDx_[kXIndex] + kEps) + k0IntOffset_[kXIndex],
        static_cast<DefAmrIndexLUint>(k0DomainSize_[kYIndex]
        / k0DomainDx_[kYIndex] + kEps) + k0IntOffset_[kYIndex] };
    // set offsets
    k0RealOffset_[kXIndex] = k0IntOffset_[kXIndex] * k0DomainDx_[kXIndex];
    k0RealOffset_[kYIndex] = k0IntOffset_[kYIndex] * k0DomainDx_[kYIndex];

    k0SFBitsetDomainMin_ = SFBitsetEncoding({k0IntOffset_[kXIndex], k0IntOffset_[kYIndex]});
    k0SFBitsetDomainMax_ = SFBitsetEncoding({k0MaxIndexOfBackgroundNode_[kXIndex],
         k0MaxIndexOfBackgroundNode_[kYIndex]});

    // check if domain size may exceed range of morton code
    /* the criterion is the maximum index for
    *  (domain size + offset distance)/(minimum grid spacing).
    *  Actually, this is a strict requirement considering high resolution
    *  might be needed near the computational domain. If high resolution
    *  region is far from the boundary, maximum index for
    *  (domain size + offset distance)/(maximum grid spacing)
    *  could be used for larger domain
    */
    // number of bits available for background mesh in one dimension
    DefAmrIndexUint bit_max = kSFBitsetBit / k0GridDims_ - k0MaxLevel_;
    DefSizet index_max = TwoPowerN(bit_max);
    DefAmrIndexLUint scale_i_level = static_cast<DefAmrIndexLUint>(TwoPowerN(k0MaxLevel_));

    if (k0MaxIndexOfBackgroundNode_.at(kXIndex) > index_max) {
        LogManager::LogError("Domain size exceeds the limits of space filling code in"
            " x direction, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h");
    }
    if (k0MaxIndexOfBackgroundNode_.at(kYIndex) > index_max) {
        LogManager::LogError("Domain size exceeds the limits of space filling code in"
            " y direction, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h");
    }
}
/**
* @brief   function to write grid information in log file.
*/
void GridManager2D::PrintGridInfo(void) const {
    // print information of grid parameters
    LogManager::LogInfo("Dimension is: " + std::to_string(k0GridDims_));
    LogManager::LogInfo("Maximum refinement level is: "
        + std::to_string(k0MaxLevel_));
    LogManager::LogInfo("Domain size is: "
        + std::to_string(k0DomainSize_.at(kXIndex)) + " X "
        + std::to_string(k0DomainSize_.at(kYIndex)));
    LogManager::LogInfo("Grid space dx is: "
        + std::to_string(k0DomainDx_.at(kXIndex)) + ", and dy is: "
        + std::to_string(k0DomainDx_.at(kYIndex))); 
}
/**
* @brief   function to reset number of extended layers.
* @param[in]  bitset_in space filling code of the center node.
* @param[out] ptr_vec_neighbors space filling codes of neighbors
*            of the center node.
*/
void GridManager2D::FindAllNeighborsSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>*const ptr_vec_neighbors) const {
    std::array<DefSFBitset, 9> array_neighbors;
    SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
    ptr_vec_neighbors->resize(9);
    memcpy(ptr_vec_neighbors->data(), array_neighbors.data(),
        9 * sizeof(DefSFBitset));
}
/**
* @brief   function to reset number of extended layers.
* @param[in]  i_level refinement level.
* @param[in]  sfbitset_in space filling code of a given node.
* @param[out] ptr_vec_extend_neg number of extended layers
*            in negative directions.
* @param[out] ptr_vec_extend_neg number of extended layers
*            in positive directions
*/
void GridManager2D::ResetExtendLayerBasedOnDomainSize(
    const DefAmrIndexUint i_level, const DefSFBitset& sfbitset_in,
    std::vector<DefAmrIndexLUint>* const ptr_vec_extend_neg,
    std::vector<DefAmrIndexLUint>* const ptr_vec_extend_pos) const {

    // extended layer at the background refinement level
    DefAmrIndexLUint two_power_i_level = static_cast<DefAmrIndexLUint>(TwoPowerN(i_level));
    DefAmrIndexLUint index_xmin = k0IntOffset_[kXIndex] * two_power_i_level;
    DefAmrIndexLUint index_xmax = k0MaxIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrIndexLUint index_ymin = k0IntOffset_[kYIndex] * two_power_i_level;
    DefAmrIndexLUint index_ymax = k0MaxIndexOfBackgroundNode_[kYIndex] * two_power_i_level;

    std::array<DefAmrIndexLUint, 2> indices;
    SFBitsetComputeIndices(sfbitset_in, &indices);

    // reset extended layer
    if ((indices[kXIndex] - index_xmin) < ptr_vec_extend_neg->at(kXIndex)) {
        ptr_vec_extend_neg->at(kXIndex) = (indices[kXIndex] - index_xmin);
    }
    if ((index_xmax - indices[kXIndex])
        < ptr_vec_extend_pos->at(kXIndex)) {
        ptr_vec_extend_pos->at(kXIndex) = (index_xmax - indices[kXIndex]);
    }
    if ((indices[kYIndex] - index_ymin) < ptr_vec_extend_neg->at(kYIndex)) {
        ptr_vec_extend_neg->at(kYIndex) = (indices[kYIndex] - index_ymin);
    }
    if ((index_ymax - indices[kYIndex])
        < ptr_vec_extend_pos->at(kYIndex)) {
        ptr_vec_extend_pos->at(kYIndex) = (index_ymax - indices[kYIndex]);
    }
}
/**
* @brief   function to find space filling code of nodes on a cell (2D)
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  node_exist grid containing nodes at the same level
* @param[out] ptr_bitsets nodes of the cell
*/
bool GridManager2D::NodesBelongToOneCell(const DefSFBitset bitset_in,
    const DefMap<DefAmrIndexUint>& node_exist,
    std::vector<DefSFBitset>* const ptr_bitsets) {
    bool bool_cell;
    ptr_bitsets->clear();
    std::array<DefSFBitset, 4> bitset_cell;
    bool_cell = SFBitsetBelongToOneCell(bitset_in, node_exist, &bitset_cell);
    ptr_bitsets->assign(bitset_cell.begin(), bitset_cell.end());
    return bool_cell;
}
/**
* @brief   function to find conners of cells including the given node
* @param[in]  sfbitset_in   space filling code of the given node
* @param[out] ptr_bitsets conners of neighboring cells
*/
void GridManager2D::FindCornersForNeighbourCells(const DefSFBitset bitset_in,
    std::vector<DefSFBitset>* const ptr_bitsets)  const {
    ptr_bitsets->resize(4);
    ptr_bitsets->at(0) = bitset_in;
    ptr_bitsets->at(1) = FindXNeg(bitset_in);
    ptr_bitsets->at(2) = FindYNeg(bitset_in);
    ptr_bitsets->at(3) = FindXNeg(ptr_bitsets->at(2));
}
/**
* @brief   function to Identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node (lower level)
* @param[in] node_coarse_outmost nodes on the outmost layer of coarser grid
* @param[in] node_exist   existing nodes
* @param[out]  ptr_inner_layer map store nodes on the inner layer
* @param[out]  ptr_mid_layer map store nodes on the middle layer
* @param[out]  ptr_outer_layer map store nodes on the outer layer
* @note
*  o     o     o     o coarse grid \n
*
*  o  x  o  x  o  x  o outer layer \n
*  x  x  x  x  x  x  x mid layer \n
*  o  x  o  x  o  x  o inner layer (outmost coarse layer) \n
*  x  x  x  x  x  x  x fine grid \n
*  x is node at higher level and o is node at lower level
*/
void GridManager2D::IdentifyInterfaceForACell(
    const DefSFBitset bitset_in,
    const DefMap<DefAmrUint>& node_coarse_outmost,
    const DefMap<DefAmrIndexUint>& node_exist,
    DefMap<DefAmrUint>* const ptr_inner_layer,
    DefMap<DefAmrUint>* const ptr_mid_layer,
    DefMap<DefAmrUint>* const ptr_outer_layer) {
    // bitset_neighbors[0]:(0, 0); bitset_neighbors[1]:(+x, 0);
    // bitset_neighbors[2]:(0, +y); bitset_neighbors[3]:(+x, +y);
    std::array<DefSFBitset, 4> bitset_neighbors;
    DefSFBitset bitset_mid_higher;
    std::array<DefMap<DefAmrUint>* const, 3> arr_ptr_layer = {
    ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    if (SFBitsetBelongToOneCell<DefAmrIndexUint>(
        bitset_in, node_exist, &bitset_neighbors)) {
        // mid (0 + dx/2, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[1] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid (0, 0 + dy/2)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[2] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid (0 + dx/2, y)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[2], bitset_neighbors[3] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid (x, 0 + dy/2)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[1], bitset_neighbors[3] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (0 + x/2, 0 + dy/2)
        ptr_mid_layer->insert({ FindXNeg(bitset_mid_higher), kFlag0_ });
    }
}
/**
* @brief   function to identify types of interface nodes
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_outmost   existing nodes on the outermost layer of coarse grid
* @param[out]  arr_ptr_layer pointer to map store interface layers
*/
void GridManager2D::IdentifyInterfaceNodeOnEdge(
    const std::array<DefSFBitset, 2>& arr_bitset_lower,
    const DefSFBitset bitset_mid_higher,
    const DefMap<DefAmrUint>& node_outmost,
    const std::array<DefMap<DefAmrUint>* const, 3>& arr_ptr_layer) {
    bool node0_flag = node_outmost.find(arr_bitset_lower[0]) != node_outmost.end(),
        node1_flag = node_outmost.find(arr_bitset_lower[1]) != node_outmost.end();
    if (node0_flag == node1_flag) {
        if (node0_flag) {
            arr_ptr_layer[2]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[2]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[2]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[1]), kFlag0_ });
        } else {
            arr_ptr_layer[0]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[0]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[0]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[1]), kFlag0_ });
        }
    } else {
        if (node0_flag) {
            arr_ptr_layer[2]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[1]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[0]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[1]), kFlag0_ });
        } else {
            arr_ptr_layer[0]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[1]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[2]->insert({
                SFBitsetToOneHigherLevel(arr_bitset_lower[1]), kFlag0_ });
        }
    }
}
/**
* @brief   function to check if the node at a given level coincides
*          with the background node
* @param[in]  i_level   refinement level
* @param[in] bitset_higher   space filling code at the given refinement level
* @param[out] ptr_bitset   space filling code at the background level (level 0)
*/
bool GridManager2D::CheckCoincideBackground(const DefAmrIndexUint i_level,
    const DefSFBitset& bitset_higher, DefSFBitset* const ptr_bitset) const {
    DefSFBitset bitset_refine = SFBitsetBitsForRefinement(i_level);
    if ((bitset_higher & bitset_refine) == 0) {
        // a node at the current level whose coordinates
        // are the same as those of a background node
        *ptr_bitset = SFBitsetToNLowerLevel(i_level, bitset_higher);
        return true;
    } else {
        return false;
    }
}
/**
* @brief   function to find all nodes in a cell at one lower refinement level
* @param[in]  bitset_cell   spacing fill codes of nodes belong to a cell
* @param[out] ptr_bitset_all   space filling codes of nodes at higher
*                              refinement level in the given cell at 
*                              one lower level
*/
void GridManager2D::FindAllNodesInACellAtLowerLevel(
    const std::vector<DefSFBitset> bitset_cell,
    std::vector<DefSFBitset>* const ptr_bitset_all) const {
    ptr_bitset_all->resize(9);
    // bitset_cell[0]:(0, 0); bitset_cell[1]:(+x, 0);
    // bitset_cell[2]:(0, +y); bitset_cell[3]:(+x, +y);
    DefSFBitset bitset_temp;
    // bottom
    bitset_temp = SFBitsetToOneHigherLevel(bitset_cell.at(0));
    ptr_bitset_all->at(0) = bitset_temp;
    ptr_bitset_all->at(1) = FindXPos(bitset_temp);
    ptr_bitset_all->at(2) = SFBitsetToOneHigherLevel(bitset_cell.at(1));
    // middle
    ptr_bitset_all->at(3) = FindYPos(ptr_bitset_all->at(0));
    ptr_bitset_all->at(4) = FindYPos(ptr_bitset_all->at(1));
    ptr_bitset_all->at(5) = FindYPos(ptr_bitset_all->at(2));
    // top
    bitset_temp = SFBitsetToOneHigherLevel(bitset_cell.at(2));
    ptr_bitset_all->at(6) = bitset_temp;
    ptr_bitset_all->at(7) = FindXPos(bitset_temp);
    ptr_bitset_all->at(8) = SFBitsetToOneHigherLevel(bitset_cell.at(3));
}
/**
* @brief function to space filling code to n level lower
* @param[in]  n_level  number of levels need to be shrink.
* @param[in]  bitset_in   input space filling code.
* @return   space filling code at lower levels.
*/
DefSFBitset GridManager2D::NodeAtNLowerLevel(
    const DefAmrIndexUint n_level, const DefSFBitset& bitset_in) const {
    return SFBitsetToNLowerLevel(n_level, bitset_in);
}
/**
* @brief   function to find layer in the overlapping region based on the layer
*          at one level higher
* @param[in]  layer_high_level   layer in the overlapping region at high level
* @param[out] ptr_layer_low_level layer in the overlapping region at low level
*/
void GridManager2D::OverlapLayerFromHighToLow(
    const DefMap<DefAmrUint>& layer_high_level,
    DefMap<DefAmrUint>* const ptr_layer_low_level) {
    for (const auto& iter : layer_high_level) {
        if ((iter.first & k0SfBitsetCurrentLevelBits_) == 0) {
            ptr_layer_low_level->insert({
                SFBitsetToOneLowerLevel(iter.first), kFlag0_ });
        }
    }
}
/**
* @brief   function to find background node in the offset region.
* @param[in] bitset_in space filling code of the node
* @return if node is in the offset region
*/
bool GridManager2D::CheckBackgroundOffset(const DefSFBitset& bitset_in) const {
    std::array<DefAmrIndexLUint, 2> indices;
    SFBitsetComputeIndices(bitset_in, &indices);
    if ((indices[kXIndex] < k0IntOffset_[kXIndex])
        || (indices[kYIndex] < k0IntOffset_[kYIndex])) {
        return true;
    } else {
        return false;
    }
}
/**
* @brief   function to instantiate background node.
* @param[in] map_occupied node exist in grid at high refinement level.
*/
void GridManager2D::InstantiateBackgroundGrid(const DefSFBitset bitset_min,
    const DefSFBitset bitset_max, const DefMap<DefAmrIndexUint>& map_occupied) {
    DefSFBitset bitset_temp;
    GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(0));

    GridNode node_instance(grid_info.k0GridNodeInstance_);
    DefSFCodeToUint i_code = bitset_min.to_ullong(), code_max = bitset_max.to_ullong();
    while (i_code <= code_max) {
        ResetIndicesExceedingDomain(k0IntOffset_, k0MaxIndexOfBackgroundNode_, &i_code, &bitset_temp);
        if (map_occupied.find(bitset_temp) == map_occupied.end()) {
                grid_info.map_grid_node_.insert({ bitset_temp, node_instance});
        }
        ++i_code;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
