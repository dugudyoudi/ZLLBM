//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
* @note  functions from other managers will be called.
*/
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
#include <string>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to setup and check grid related parameters.
*/
void GridManager3D::SetGridParameters() {
    // check if length of computational domain is given
    if (k0DomainSize_.at(kXIndex) < kEps) {
        LogManager::LogError("Domain length in x direction (k0DomainSize_[0])"
            " should be a positive value in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else if (k0DomainSize_.at(kYIndex) < kEps) {
        LogManager::LogError("Domain length in y direction (k0DomainSize_[1])"
            " should be a positive value in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else if (k0DomainSize_.at(kZIndex) < kEps) {
        LogManager::LogError("Domain length in z direction (k0DomainSize_[2])"
            " should be a positive value in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    // check if grid space is given
    if (k0DomainDx_.at(kXIndex) < kEps
        && k0DomainDx_.at(kYIndex) < kEps
        && k0DomainDx_.at(kZIndex < kEps)) {
        LogManager::LogError("Grid space of x, y, or z (k0DomainDx_)"
            " should be positive values in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    // set grid space if not all grid spaces are given
    if (k0DomainDx_.at(kXIndex) < kEps) {
        if (k0DomainDx_.at(kYIndex) > kEps) {
            k0DomainDx_.at(kXIndex) = k0DomainDx_.at(kYIndex);
        } else if (k0DomainDx_.at(kZIndex) > kEps) {
            k0DomainDx_.at(kXIndex) = k0DomainDx_.at(kZIndex);
        }
    }
    if (k0DomainDx_.at(kYIndex) < kEps) {
        if (k0DomainDx_.at(kXIndex) > kEps) {
            k0DomainDx_.at(kYIndex) = k0DomainDx_.at(kXIndex);
        } else if (k0DomainDx_.at(kZIndex) > kEps) {
            k0DomainDx_.at(kYIndex) = k0DomainDx_.at(kZIndex);
        }
    }
    if (k0DomainDx_.at(kZIndex) < kEps) {
        if (k0DomainDx_.at(kXIndex) > kEps) {
            k0DomainDx_.at(kZIndex) = k0DomainDx_.at(kXIndex);
        } else if (k0DomainDx_.at(kYIndex) > kEps) {
            k0DomainDx_.at(kZIndex) = k0DomainDx_.at(kYIndex);
        }
    }
    k0SpaceBackground_ = { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex],
    k0DomainDx_[kZIndex] };

    // calculate number of background nodes in each direction
    k0MaxIndexOfBackgroundNode_ = {
        static_cast<DefAmrIndexLUint>(k0DomainSize_[kXIndex]
        / k0DomainDx_[kXIndex] + kEps) + k0MinIndexOfBackgroundNode_[kXIndex],
        static_cast<DefAmrIndexLUint>(k0DomainSize_[kYIndex]
        / k0DomainDx_[kYIndex] + kEps) + k0MinIndexOfBackgroundNode_[kYIndex],
        static_cast<DefAmrIndexLUint>(k0DomainSize_[kZIndex]
        / k0DomainDx_[kZIndex] + kEps) + k0MinIndexOfBackgroundNode_[kZIndex]};

    SFBitsetMinAndMaxGlobal(
        k0MinIndexOfBackgroundNode_, k0MaxIndexOfBackgroundNode_);
    SFBitsetMinAndMaxCoordinates(k0MaxLevel_,
        k0MinIndexOfBackgroundNode_, k0MaxIndexOfBackgroundNode_);

    // set offsets
    k0RealMin_[kXIndex] = k0MinIndexOfBackgroundNode_[kXIndex] * k0DomainDx_[kXIndex];
    k0RealMin_[kYIndex] = k0MinIndexOfBackgroundNode_[kYIndex] * k0DomainDx_[kYIndex];
    k0RealMin_[kZIndex] = k0MinIndexOfBackgroundNode_[kZIndex] * k0DomainDx_[kZIndex];

    k0SFBitsetDomainMin_ = SFBitsetEncoding({k0MinIndexOfBackgroundNode_[kXIndex], k0MinIndexOfBackgroundNode_[kYIndex], k0MinIndexOfBackgroundNode_[kZIndex]});
    k0SFBitsetDomainMax_ = SFBitsetEncoding({k0MaxIndexOfBackgroundNode_[kXIndex],
        k0MaxIndexOfBackgroundNode_[kYIndex], k0MaxIndexOfBackgroundNode_[kZIndex]});

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
            " storing space filling code (kSFBitsetBit) in defs_libs.h. Error in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (k0MaxIndexOfBackgroundNode_.at(kYIndex) > index_max) {
        LogManager::LogError("Domain size exceeds the limits of space filling code  in"
            " y direction, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h. Error in"
            +  std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (k0MaxIndexOfBackgroundNode_.at(kZIndex) > index_max) {
        LogManager::LogError("Domain size exceeds the limits of space filling code  in"
            " z direction, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h. Error in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief   function to write grid information in log file.
*/
void GridManager3D::PrintGridInfo(void) const {
    LogManager::LogInfo("Dimension is: " + std::to_string(k0GridDims_));
    LogManager::LogInfo("Maximum refinement level is: "
        + std::to_string(k0MaxLevel_));
    if (k0GridDims_ == 2) {
        LogManager::LogInfo("Domain size is: "
            + std::to_string(k0DomainSize_.at(kXIndex)) + " X "
            + std::to_string(k0DomainSize_.at(kYIndex)));
        LogManager::LogInfo("Grid space dx is: "
            + std::to_string(k0DomainDx_.at(kXIndex)) + ", and dy is: "
            + std::to_string(k0DomainDx_.at(kYIndex)));
    } else if (k0GridDims_ == 3) {
        LogManager::LogInfo("Domain size is: "
            + std::to_string(k0DomainSize_.at(kXIndex)) + " X "
            + std::to_string(k0DomainSize_.at(kYIndex)) + " X "
            + std::to_string(k0DomainSize_.at(kZIndex)));
        LogManager::LogInfo("Grid space dx is: "
            + std::to_string(k0DomainDx_.at(kXIndex)) + " , dy is: "
            + std::to_string(k0DomainDx_.at(kYIndex)) + " , and dz is: "
            + std::to_string(k0DomainDx_.at(kZIndex)));
    }
}
/**
* @brief   function to reset number of extended layers.
* @param[in]  bitset_in space filling code of the center node.
* @param[out] ptr_vec_neighbors space filling codes of neighbors
*            of the center node.
*/
void GridManager3D::GridFindAllNeighborsVir(const DefSFBitset& bitset_in,
    std::vector<DefSFBitset>* const ptr_vec_neighbors) const {
    std::array<DefSFBitset, 27> array_neighbors;
    SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
    ptr_vec_neighbors->resize(27);
    memcpy(ptr_vec_neighbors->data(), array_neighbors.data(),
        27 * sizeof(DefSFBitset));
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
void GridManager3D::ResetExtendLayerBasedOnDomainSize(
    const DefAmrIndexUint i_level, const DefSFBitset& sfbitset_in,
    std::vector<DefAmrIndexLUint>* const ptr_vec_extend_neg,
    std::vector<DefAmrIndexLUint>* const ptr_vec_extend_pos) const {
    //  extended layer at the background refinement level
    DefAmrIndexLUint two_power_i_level = static_cast<DefAmrIndexLUint>(TwoPowerN(i_level));
    DefAmrIndexLUint index_xmin = k0MinIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrIndexLUint index_xmax = k0MaxIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrIndexLUint index_ymin = k0MinIndexOfBackgroundNode_[kYIndex] * two_power_i_level;
    DefAmrIndexLUint index_ymax = k0MaxIndexOfBackgroundNode_[kYIndex] * two_power_i_level;
    DefAmrIndexLUint index_zmin = k0MinIndexOfBackgroundNode_[kZIndex] * two_power_i_level;
    DefAmrIndexLUint index_zmax = k0MaxIndexOfBackgroundNode_[kZIndex] * two_power_i_level;

    std::array<DefAmrIndexLUint, 3> indices;
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
    if ((indices[kZIndex] - index_zmin) < ptr_vec_extend_neg->at(kZIndex)) {
        ptr_vec_extend_neg->at(kZIndex) = (indices[kZIndex] - index_zmin);
    }
    if ((index_zmax - indices[kZIndex])
        < ptr_vec_extend_pos->at(kZIndex)) {
        ptr_vec_extend_pos->at(kZIndex) = (index_zmax - indices[kZIndex]);
    }
}
/**
* @brief   function to find space filling code of nodes on a cell (3D)
* @param[in]  sfbitset_in   bitset of the node at the origin of a cell
* @param[in]  node_exist grid containing nodes at the same level
* @param[out] ptr_bitsets nodes of the cell
*/
bool GridManager3D::NodesBelongToOneCell(const DefSFBitset bitset_in,
    const DefMap<DefAmrIndexUint>& node_exist,
    std::vector<DefSFBitset>* const ptr_bitsets) {
    ptr_bitsets->clear();
    bool bool_cell;
    std::array<DefSFBitset, 8> bitset_cell;
    bool_cell = SFBitsetBelongToOneCell(bitset_in, node_exist, &bitset_cell);
    ptr_bitsets->assign(bitset_cell.begin(), bitset_cell.end());
    return bool_cell;
}
/**
* @brief   function to find conners of cells including the given node
* @param[in]  sfbitset_in   bitset of the given node
* @param[out] ptr_bitsets conners of neighboring cells
*/
void GridManager3D::FindCornersForNeighbourCells(const DefSFBitset bitset_in,
    std::vector<DefSFBitset>* const ptr_bitsets)  const {
    ptr_bitsets->resize(8);
    ptr_bitsets->at(0) = bitset_in;
    ptr_bitsets->at(1) = FindXNeg(bitset_in);
    ptr_bitsets->at(2) = FindYNeg(bitset_in);
    ptr_bitsets->at(3) = FindXNeg(ptr_bitsets->at(2));
    ptr_bitsets->at(4) = FindZNeg(bitset_in);
    ptr_bitsets->at(5) = FindXNeg(ptr_bitsets->at(4));
    ptr_bitsets->at(6) = FindYNeg(ptr_bitsets->at(4));
    ptr_bitsets->at(7) = FindXNeg(ptr_bitsets->at(6));
}
/**
* @brief   function to Identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node (lower level)
* @param[in] node_coarse_outmost nodes on the outmost layer of coarser grid
* @param[in] node_exist   existing nodes
* @param[out]  ptr_inner_layer map store nodes on the inner layer
* @param[out]  ptr_mid_layer map store nodes on the middle layer
* @param[out]  ptr_outer_layer map store nodes on the outer layer
*/
void GridManager3D::IdentifyInterfaceForACell(const DefSFBitset bitset_in,
    const DefMap<DefAmrUint>& node_coarse_outmost, const DefMap<DefAmrIndexUint>& node_exist,
    DefMap<DefAmrUint>* const ptr_inner_layer,
    DefMap<DefAmrUint>* const ptr_mid_layer,
    DefMap<DefAmrUint>* const ptr_outer_layer) {
    // bitset_neighbors[0]:(0, 0, 0);   bitset_neighbors[1]:(+x, 0, 0);
    // bitset_neighbors[2]:(0, +y, 0);  bitset_neighbors[3]:(+x, +y, 0);
    // bitset_neighbors[4]:(0, 0, +z);  bitset_neighbors[5]:(+x, 0, +z);
    // bitset_neighbors[6]:(0, +y, +z); bitset_neighbors[7]:(+x, +y, +z).
    std::array<DefSFBitset, 8> bitset_neighbors;
    DefSFBitset bitset_mid_higher;
    std::array<DefMap<DefAmrUint>* const, 3> arr_ptr_layer = {
    ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    if (SFBitsetBelongToOneCell<DefAmrIndexUint>(bitset_in, node_exist, &bitset_neighbors)) {
        // bottom surface
        // edge (0 + dx/2, 0, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[1] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (0, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[2] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (0 + dx/2, y, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[2], bitset_neighbors[3] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (x, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[1], bitset_neighbors[3] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (0 + x/2, 0 + dy/2, 0)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[0], bitset_neighbors[1],
              bitset_neighbors[2], bitset_neighbors[3] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);

        // top surface
        // mid edge (0 + dx/2, 0, z)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[4], bitset_neighbors[5] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (0, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[4], bitset_neighbors[6] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (0 + dx/2, y, z)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[6]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[6], bitset_neighbors[7] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (x, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[5]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[5], bitset_neighbors[7] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (0 + x/2, 0 + dy/2, z)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[4], bitset_neighbors[5],
              bitset_neighbors[6], bitset_neighbors[7] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);

        // middle surface
        // mid edge (0, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[4] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (0 + x/2, 0, 0 + dz/2)
        bitset_mid_higher = FindXPos(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[0], bitset_neighbors[1],
              bitset_neighbors[5], bitset_neighbors[4] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (x, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[1], bitset_neighbors[5] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (x, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYPos(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[1], bitset_neighbors[3],
              bitset_neighbors[7], bitset_neighbors[5] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (0, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[2], bitset_neighbors[6] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (0, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[0], bitset_neighbors[2],
              bitset_neighbors[6], bitset_neighbors[4] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // mid edge (x, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[3]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[3], bitset_neighbors[7] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);
        // diagonal (0 + dx/2, y, 0 + dz/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[2], bitset_neighbors[3],
              bitset_neighbors[7], bitset_neighbors[6] },
            bitset_mid_higher, node_coarse_outmost, arr_ptr_layer);

        // center (0 + dx/2, 0 + dy/2, 0 + dz/2)
        ptr_mid_layer->insert({ FindYNeg(bitset_mid_higher), kFlag0_ });
    }
}
/**
* @brief   function to identify types of interface node on the edge
* @param[in]  flag_interface   flag of marked interface node
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_outmost   existing nodes on the outermost layer of coarse grid
* @param[out]  arr_ptr_layer pointer to map store interface layers
*/
void GridManager3D::IdentifyInterfaceNodeOnEdge(
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
* @brief   function to Identify types of interface diagonal node
* @param[in]  flag_interface   flag of marked interface node
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_outmost   existing nodes on the outermost layer of coarse grid
* @param[out]  arr_ptr_layer pointer to map store interface layers
*/
void GridManager3D::IdentifyInterfaceNodeDiagonal(
    const std::array<DefSFBitset, 4>& arr_bitset_lower,
    const DefSFBitset bitset_center_higher,
    const DefMap<DefAmrUint>& node_outmost,
    const std::array<DefMap<DefAmrUint>* const, 3>& arr_ptr_layer) {
    DefAmrUint node0_flag = node_outmost.find(arr_bitset_lower[0]) != node_outmost.end(),
     node1_flag = node_outmost.find(arr_bitset_lower[1]) != node_outmost.end(),
     node2_flag = node_outmost.find(arr_bitset_lower[2]) != node_outmost.end(),
     node3_flag = node_outmost.find(arr_bitset_lower[3]) != node_outmost.end();
    if (node0_flag == node1_flag && node1_flag == node2_flag
        && node2_flag == node3_flag) {
        if (node0_flag) {
            arr_ptr_layer[2]->insert({ bitset_center_higher, kFlag0_ });
        } else {
            arr_ptr_layer[0]->insert({ bitset_center_higher, kFlag0_ });
        }
    } else {
        arr_ptr_layer[1]->insert({ bitset_center_higher, kFlag0_ });
    }
}
/**
* @brief   function to check if the node at a given level coincides
*          with the background node
* @param[in]  i_level   refinement level
* @param[in] bitset_higher   space filling code at the given refinement level
* @param[out] ptr_bitset   space filling code at the background level (level 0)
*/
bool GridManager3D::CheckCoincideBackground(const DefAmrIndexUint i_level,
    const DefSFBitset& bitset_higher, DefSFBitset* const ptr_bitset) const {
    DefSFBitset bitset_refine = SFBitsetBitsForRefinement(i_level);
    if ((bitset_higher & bitset_refine) == 0) {
        // a node at the current level whose coordiniates
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
void GridManager3D::FindAllNodesInACellAtLowerLevel(
    const std::vector<DefSFBitset> bitset_cell,
    std::vector<DefSFBitset>* const ptr_bitset_all) const {
    ptr_bitset_all->resize(27);
    // bitset_cell[0]:(0, 0, 0);   bitset_cell[1]:(+x, 0, 0);
    // bitset_cell[2]:(0, +y, 0);  bitset_cell[3]:(+x, +y, 0);
    // bitset_cell[4]:(0, 0, +z);  bitset_cell[5]:(+x, 0, +z);
    // bitset_cell[6]:(0, +y, +z); bitset_cell[7]:(+x, +y, +z).
    DefSFBitset bitset_temp;
    // bottom
    bitset_temp = SFBitsetToOneHigherLevel(bitset_cell.at(0));
    ptr_bitset_all->at(0) = bitset_temp;
    ptr_bitset_all->at(1) = FindXPos(bitset_temp);
    ptr_bitset_all->at(2) = SFBitsetToOneHigherLevel(bitset_cell.at(1));
    ptr_bitset_all->at(3) = FindYPos(ptr_bitset_all->at(0));
    ptr_bitset_all->at(4) = FindYPos(ptr_bitset_all->at(1));
    ptr_bitset_all->at(5) = FindYPos(ptr_bitset_all->at(2));
    bitset_temp = SFBitsetToOneHigherLevel(bitset_cell.at(2));
    ptr_bitset_all->at(6) = bitset_temp;
    ptr_bitset_all->at(7) = FindXPos(bitset_temp);
    ptr_bitset_all->at(8) = SFBitsetToOneHigherLevel(bitset_cell.at(3));

    // middle
    ptr_bitset_all->at(9) = FindZPos(ptr_bitset_all->at(0));
    ptr_bitset_all->at(10) = FindZPos(ptr_bitset_all->at(1));
    ptr_bitset_all->at(11) = FindZPos(ptr_bitset_all->at(2));
    ptr_bitset_all->at(12) = FindZPos(ptr_bitset_all->at(3));
    ptr_bitset_all->at(13) = FindZPos(ptr_bitset_all->at(4));
    ptr_bitset_all->at(14) = FindZPos(ptr_bitset_all->at(5));
    ptr_bitset_all->at(15) = FindZPos(ptr_bitset_all->at(6));
    ptr_bitset_all->at(16) = FindZPos(ptr_bitset_all->at(7));
    ptr_bitset_all->at(17) = FindZPos(ptr_bitset_all->at(8));
    // top
    bitset_temp = SFBitsetToOneHigherLevel(bitset_cell.at(4));
    ptr_bitset_all->at(18) = bitset_temp;
    ptr_bitset_all->at(19) = FindXPos(bitset_temp);
    ptr_bitset_all->at(20) = SFBitsetToOneHigherLevel(bitset_cell.at(5));
    ptr_bitset_all->at(21) = FindYPos(ptr_bitset_all->at(18));
    ptr_bitset_all->at(22) = FindYPos(ptr_bitset_all->at(19));
    ptr_bitset_all->at(23) = FindYPos(ptr_bitset_all->at(20));
    bitset_temp = SFBitsetToOneHigherLevel(bitset_cell.at(6));
    ptr_bitset_all->at(24) = bitset_temp;
    ptr_bitset_all->at(25) = FindXPos(bitset_temp);
    ptr_bitset_all->at(26) = SFBitsetToOneHigherLevel(bitset_cell.at(7));
}
/**
* @brief function to space filling code to n level lower
* @param[in]  n_level  number of levels need to be shrink.
* @param[in]  bitset_in   input space filling code.
* @return   space filling code at lower levels.
*/
DefSFBitset GridManager3D::NodeAtNLowerLevel(
    const DefAmrIndexUint n_level, const DefSFBitset& bitset_in) const {
    return SFBitsetToNLowerLevel(n_level, bitset_in);
}
/**
* @brief   function to find layer in the overlapping region based on the layer
*          at one level higher
* @param[in]  layer_high_level   layer in the overlapping region at high level
* @param[out] ptr_layer_low_level layer in the overlapping region at low level
*/
void GridManager3D::OverlapLayerFromHighToLow(
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
bool GridManager3D::CheckBackgroundOffset(const DefSFBitset& bitset_in) const {
    std::array<DefAmrIndexLUint, 3> indices;
    SFBitsetComputeIndices(bitset_in, &indices);
    if ((indices[kXIndex] < k0MinIndexOfBackgroundNode_[kXIndex])
        || (indices[kYIndex] < k0MinIndexOfBackgroundNode_[kYIndex])
        || (indices[kZIndex] < k0MinIndexOfBackgroundNode_[kZIndex])) {
        return true;
    } else {
        return false;
    }
}
/**
* @brief   function to instantiate background node.
* @param[in] bitset_min minimum space filling code of background nodes
* @param[in] bitset_max maximum space filling code of background nodes
* @param[in] map_occupied node exist in grid at high refinement level.
*/
void GridManager3D::InstantiateBackgroundGrid(const DefSFBitset bitset_min,
    const DefSFBitset bitset_max, const DefMap<DefAmrIndexUint>& map_occupied) {
    DefSFBitset bitset_temp;
    GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(0));
    DefSFCodeToUint i_code = bitset_min.to_ullong(), code_max = bitset_max.to_ullong();
    while (i_code <= code_max) {
        ResetIndicesExceedingDomain(k0MinIndexOfBackgroundNode_, k0MaxIndexOfBackgroundNode_, &i_code, &bitset_temp);
        if (map_occupied.find(bitset_temp) == map_occupied.end()) {
                grid_info.map_grid_node_.insert({ bitset_temp, grid_info.GridNodeCreator() });
        }
        ++i_code;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
