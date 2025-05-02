//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
* @brief function to setup size of the computational domain.
* @param[in] domain_size maximum coordinates of the computational domain.
*/
void GridManager2D::SetDomainSize(const std::vector<DefReal>& domain_size) {
    if (domain_size.size() == 2) {
        k0DomainSize_.at(kXIndex) = domain_size.at(kXIndex);
        k0DomainSize_.at(kYIndex) = domain_size.at(kYIndex);
    } else {
        LogManager::LogError("size of the input vector should be 2 in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief function to setup grid spacing the background level.
* @param[in] domain_grid_size grid size of the computational domain.
*/
void GridManager2D::SetDomainGridSize(const std::vector<DefReal>& domain_grid_size) {
    if (domain_grid_size.size() == 1) {
        k0DomainDx_.at(kXIndex) = domain_grid_size.at(kXIndex);
    } else if (domain_grid_size.size() == 2) {
        k0DomainDx_.at(kXIndex) = domain_grid_size.at(kXIndex);
        k0DomainDx_.at(kYIndex) = domain_grid_size.at(kYIndex);
    } else {
        LogManager::LogError("size of the input vector should be 1 or 2 in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief function to setup and check grid related parameters.
*/
void GridManager2D::SetupDependentGridParameters() {
    // check if length of computational domain is given
    if (k0DomainSize_.at(kXIndex) < kEps) {
        LogManager::LogError("Domain length in x direction (k0DomainSize_[0])"
            " should be a positive value in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else if (k0DomainSize_.at(kYIndex) < kEps) {
        LogManager::LogError("Domain length in x direction (k0DomainSize_[1])"
            " should be a positive value in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    // check if grid space is given
    if (k0DomainDx_.at(kXIndex) < kEps
        && k0DomainDx_.at(kYIndex) < kEps) {
        LogManager::LogError("Grid space of x or y(k0DomainDx_) in "
            " should be positive values"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    // set grid space if not all grid spaces are given
    if (k0DomainDx_.at(kXIndex) < kEps) {
        k0DomainDx_.at(kXIndex) = k0DomainDx_.at(kYIndex);
    }
    if (k0DomainDx_.at(kYIndex) < kEps) {
        k0DomainDx_.at(kYIndex) = k0DomainDx_.at(kXIndex);
    }
    k0SpaceBackground_ = { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex]};

    // calculate number of background nodes in each direction
    k0MaxIndexOfBackgroundNode_ = {
        static_cast<DefAmrLUint>(k0DomainSize_[kXIndex]
        / k0DomainDx_[kXIndex] + kEps) + k0MinIndexOfBackgroundNode_[kXIndex],
        static_cast<DefAmrLUint>(k0DomainSize_[kYIndex]
        / k0DomainDx_[kYIndex] + kEps) + k0MinIndexOfBackgroundNode_[kYIndex]};

    SFBitsetSetMinAndMaxBounds(k0MinIndexOfBackgroundNode_, k0MaxIndexOfBackgroundNode_);

    // set offsets
    k0RealMin_[kXIndex] = k0MinIndexOfBackgroundNode_[kXIndex] * k0DomainDx_[kXIndex];
    k0RealMin_[kYIndex] = k0MinIndexOfBackgroundNode_[kYIndex] * k0DomainDx_[kYIndex];

    k0SFBitsetDomainMin_ = SFBitsetEncoding({k0MinIndexOfBackgroundNode_[kXIndex],
        k0MinIndexOfBackgroundNode_[kYIndex]});
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
    DefInt bit_max = kSFBitsetBit / k0GridDims_ - k0MaxLevel_ - 1;
    DefAmrLUint index_max = TwoPowerN(bit_max);

    if (k0MaxIndexOfBackgroundNode_.at(kXIndex) > index_max) {
        LogManager::LogError("Domain size exceeds the limits of space filling code in"
            " x direction, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h. Error in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (k0MaxIndexOfBackgroundNode_.at(kYIndex) > index_max) {
        LogManager::LogError("Domain size exceeds the limits of space filling code in"
            " y direction, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h. Error in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
 * @brief function to calculates the domain bounds for a given refinement level.
 * @param[in] i_level level of refinement.
 * @param[out] ptr_domain_min pointer to the minimum domain bounds at a given refinement level.
 * @param[out] ptr_domain_max pointer to the maximum domain bounds at a given refinement level.
 */
void GridManager2D::CalDomainBoundsAtGivenLevel(const DefInt i_level,
    std::vector<DefSFBitset>* const ptr_domain_min, std::vector<DefSFBitset>* const ptr_domain_max) const {
    ptr_domain_min->resize(2);
    ptr_domain_max->resize(2);
    if (SFBitsetMin_.at(kXIndex) == ~DefSFCodeToUint(0)) {
         LogManager::LogError("SFBitsetMin_ should be set before calling CalDomainBoundsAtGivenLevel"
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        ptr_domain_min->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetMin_.at(kXIndex));
    }
    if (SFBitsetMin_.at(kYIndex) == ~DefSFCodeToUint(0)) {
         LogManager::LogError("SFBitsetMin_ should be set before calling CalDomainBoundsAtGivenLevel"
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        ptr_domain_min->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetMin_.at(kYIndex));
    }
    if (SFBitsetMax_.at(kXIndex) == 0) {
         LogManager::LogError("SFBitsetMax_ should be set before calling CalDomainBoundsAtGivenLevel"
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        ptr_domain_max->at(kXIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetMax_.at(kXIndex));
    }
    if (SFBitsetMax_.at(kYIndex) == 0) {
         LogManager::LogError("SFBitsetMax_ should be set before calling CalDomainBoundsAtGivenLevel"
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        ptr_domain_max->at(kYIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetMax_.at(kYIndex));
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
void GridManager2D::GridFindAllNeighborsVir(const DefSFBitset& bitset_in,
    std::vector<DefSFBitset>* const ptr_vec_neighbors) const {
    std::array<DefSFBitset, 9> array_neighbors;
    SFBitsetFindAllNeighbors(bitset_in, &array_neighbors);
    ptr_vec_neighbors->resize(9);
    memcpy(ptr_vec_neighbors->data(), array_neighbors.data(), 9 * sizeof(DefSFBitset));
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
    const DefInt i_level, const DefSFBitset& sfbitset_in,
    std::vector<DefAmrLUint>* const ptr_vec_extend_neg,
    std::vector<DefAmrLUint>* const ptr_vec_extend_pos) const {

    // extended layer at the background refinement level
    DefAmrLUint two_power_i_level = static_cast<DefAmrLUint>(TwoPowerN(i_level));
    DefAmrLUint index_xmin = k0MinIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrLUint index_xmax = k0MaxIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrLUint index_ymin = k0MinIndexOfBackgroundNode_[kYIndex] * two_power_i_level;
    DefAmrLUint index_ymax = k0MaxIndexOfBackgroundNode_[kYIndex] * two_power_i_level;

    std::array<DefAmrLUint, 2> indices;
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
* @brief   function to check node on which domain boundaries.
* @param[in] i_level refinement level.
* @param[in] sfbitset_in space filling code of a given node.
* @return flag indicate node on which domain boundaries, 1: x min, 2: x max, 4: y min, 8: y max
*/
int GridManager2D::CheckNodeOnDomainBoundary(
    const DefInt i_level, const DefSFBitset& sfbitset_in) const {
    int node_status = 0;
    std::array<DefSFBitset, 2> sfbitset_min, sfbitset_max;
    sfbitset_min[kXIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMin_[kXIndex]);
    sfbitset_min[kYIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMin_[kYIndex]);
    sfbitset_max[kXIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMax_[kXIndex]);
    sfbitset_max[kYIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMax_[kYIndex]);
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_min[kXIndex]) {
        node_status |= 1;
    }
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_max[kXIndex]) {
        node_status |= 2;  // 1<< 1
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_min[kYIndex]) {
        node_status |= 4;  // 1<< 2
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_max[kYIndex]) {
        node_status |= 8;  // 1<< 3
    }
    return node_status;
}
/**
* @brief   function to check if node is not outside a cubic domain boundary.
* @param[in] sfbitset_in space filling code of a given node.
* @param[in] sfbitset_min bitset corresponding to the minimum coordinate in each direction.
* @param[in] sfbitset_max bitset corresponding to the maximum coordinate in each direction.
* @return  if false, node is outside a cubic domain boundary
*/
bool GridManager2D::CheckNodeNotOutsideDomainBoundary(const DefSFBitset& sfbitset_in,
    const std::vector<DefSFCodeToUint>& sfbitset_min,
    const std::vector<DefSFCodeToUint>& sfbitset_max) const {
    if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeXRef_[kRefCurrent_])< sfbitset_min.at(kXIndex)) {
        return false;
    } else if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeXRef_[kRefCurrent_]) > sfbitset_max.at(kXIndex)) {
        return false;
    }
    if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeYRef_[kRefCurrent_]) < sfbitset_min.at(kYIndex)) {
        return false;
    } else if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeYRef_[kRefCurrent_]) > sfbitset_max.at(kYIndex)) {
        return false;
    }
    return true;
}
/**
* @brief   function to find space filling code of nodes on a cell (2D)
* @param[in]  sfbitset_in   space filling code of the node at the corner of a cell
* @param[in]  node_exist grid containing nodes at the same level
* @param[out] ptr_sfbitset all nodes in the cell
*/
bool GridManager2D::NodesBelongToOneCell(const DefSFBitset sfbitset_in,
    const DefMap<DefInt>& node_exist,
    std::vector<DefSFBitset>* const ptr_sfbitset) const {
    bool bool_cell;
    ptr_sfbitset->clear();
    std::array<DefSFBitset, 4> bitset_cell;
    bool_cell = SFBitsetBelongToOneCell(sfbitset_in, node_exist, &bitset_cell);
    ptr_sfbitset->assign(bitset_cell.begin(), bitset_cell.end());
    return bool_cell;
}
/**
* @brief   function to find space filling code of nodes on an edge (2D)
* @param[in]  sfbitset_in   space filling code of the node at the stating point of an edge
* @param[in]  dir_norm   normal to the edge: kXIndex or kYIndex
* @param[in]  map_node_exist grid containing nodes at the same level
* @param[out] ptr_sfbitset all nodes on the edge
*/
bool GridManager2D::NodesBelongToOneSurfAtHigherLevel(const DefSFBitset sfbitset_in,
    const DefInt dir_norm, const DefMap<DefInt>& map_node_exist,
    std::vector<DefSFBitset>* const ptr_sfbitset) const {
    ptr_sfbitset->clear();
    DefSFBitset sfbitset_tmp;
    switch (dir_norm) {
    case kXIndex:
        if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_in));
        // (+y, 0)
        sfbitset_tmp = FindYPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(FindYPos(ptr_sfbitset->at(0)));
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp));
        break;
    case kYIndex:
        if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_in));
        // (+x, 0)
        sfbitset_tmp = FindXPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(FindXPos(ptr_sfbitset->at(0)));
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp));
        break;
    default:
        return false;
    }
    return true;
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
* @brief   function to identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node (lower level)
* @param[in] node_coarse_interface nodes on the interface layer of coarser grid
* @param[in] node_exist_lower   existing nodes at lower level
* @param[out] ptr_inner_layer map storing nodes on the inner layer
* @param[out] ptr_mid_layer map storing nodes on the middle layer
* @param[out] ptr_outer_layer map storing nodes on the outer layer
* @note
*  o     o     o     o coarse grid \n
*
*  o  x  o  x  o  x  o outer layer \n
*  x  x  x  x  x  x  x mid layer \n
*  o  x  o  x  o  x  o inner layer (outmost coarse layer) \n
*  x  x  x  x  x  x  x fine grid \n
*  x is node at higher level and o is node at lower level
*/
// in node_exist_lower, only nodes on the refinement interface exist since grid generation hasn't been done
void GridManager2D::IdentifyInterfaceForACell(
    const DefSFBitset bitset_in,
    const DefMap<DefInt>& node_coarse_interface,
    const DefMap<DefInt>& node_exist_lower,
    DefMap<DefInt>* const ptr_inner_layer,
    DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer) {
    // bitset_neighbors[0]:(0, 0); bitset_neighbors[1]:(+x, 0);
    // bitset_neighbors[2]:(0, +y); bitset_neighbors[3]:(+x, +y);
    std::array<DefSFBitset, 4> bitset_neighbors;
    DefSFBitset bitset_mid_higher, bitset_node0_higher, bitset_node1_higher;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
        ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    bool belong_to_cell =  SFBitsetBelongToOneCell<DefInt>(
        bitset_in, node_exist_lower, &bitset_neighbors);
    if (belong_to_cell) {
        // mid (0 + dx/2, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[1]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid (0, 0 + dy/2)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid (0 + dx/2, y)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[2], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid (x, 0 + dy/2)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (0 + x/2, 0 + dy/2)
        ptr_mid_layer->insert({ FindXNeg(bitset_mid_higher), kFlag0_ });
    }
}
/**
* @brief   function to identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node (lower level)
* @param[in] node_coarse_interface_previous nodes on coarse to fine interface has been marked
* @param[in] node_coarse_interface_inner nodes on inner coarse to fine interface
* @param[in] node_exist_current existing nodes at current level
* @param[in] node_exist_lower   existing nodes at lower level
* @param[out] ptr_inner_layer map storing nodes on the inner layer
* @param[out] ptr_mid_layer map storing nodes on the middle layer
* @param[out] ptr_outer_layer map storing nodes on the outer layer
* @param[out] ptr_node_coarse_interface_outer pointer to map storing nodes on outer coarse to fine interface
*/
// in node_exist_current and node_exist_lower, all nodes exist since grid generation is done
// the aim is to add nodes at current refinement level to refinement interfaces
void GridManager2D::IdentifyInterfaceForACell(const DefSFBitset bitset_in,
    const DefMap<DefInt>& node_coarse_interface_previous,
    const DefMap<DefInt>& node_coarse_interface_inner,
    const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
    const DefMap<std::unique_ptr<GridNode>>& node_exist_lower,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_node_coarse_interface_outer) {
    // bitset_neighbors[0]:(0, 0); bitset_neighbors[1]:(+x, 0);
    // bitset_neighbors[2]:(0, +y); bitset_neighbors[3]:(+x, +y);
    std::array<DefSFBitset, 4> bitset_neighbors;
    DefSFBitset bitset_mid_higher, bitset_tmp;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
        ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    DefMap<DefInt> discard_layer;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_mid_layer = {
        &discard_layer, ptr_mid_layer, &discard_layer };
    bool belong_to_cell =  SFBitsetBelongToOneCell<std::unique_ptr<GridNode>>(
        bitset_in, node_exist_lower, &bitset_neighbors);
    if (belong_to_cell) {
        // mid (0 + dx/2, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[1]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid (0, 0 + dy/2)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid (0 + dx/2, y)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[2], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid (x, 0 + dy/2)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0 + dy/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_mid_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_mid_layer, ptr_node_coarse_interface_outer);
    }
}
/**
* @brief   function to identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node at lower level(level - 1)
* @param[in] node_coarse_innermost nodes on the innermost interface layer of coarser grid
* @param[in] map_node_exist   existing fine nodes at lower level (level - 1)
* @param[in] map_exist_coarse   existing coarse nodes at level - 1
* @param[out] ptr_inner_layer map storing nodes on the inner layer
* @param[out] ptr_mid_layer map storing nodes on the middle layer
* @param[out] ptr_outer_layer map storing nodes on the outer layer
*/
void GridManager2D::IdentifyInterfaceForACellAcrossTwoLevels(
    const DefSFBitset bitset_in, const DefMap<DefInt>& node_coarse_innermost,
    const DefMap<DefInt>& map_node_exist, const DefMap<DefInt>& map_exist_coarse,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* ptr_coarse_outer) {
    std::array<std::pair<DefSFBitset, DefInt>, 4> bitset_neighbors;
    DefSFBitset bitset_mid_higher, bitset_tmp;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
        ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    DefMap<DefInt> discard_layer;
    bool belong_to_cell = SFBitsetBelongToOneCellAcrossTwoLevels(
        bitset_in, map_node_exist, map_exist_coarse, &bitset_neighbors);
    if (belong_to_cell) {
        // mid (0 + dx/2, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels(
            { bitset_neighbors[0], bitset_neighbors[1]},
            bitset_mid_higher, *this, node_coarse_innermost, arr_ptr_layer, ptr_coarse_outer);

        // mid (0, 0 + dy/2)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels(
            { bitset_neighbors[0], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_innermost, arr_ptr_layer, ptr_coarse_outer);

        // mid (0 + dx/2, y)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels(
            { bitset_neighbors[2], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_innermost, arr_ptr_layer, ptr_coarse_outer);
        // mid (x, 0 + dy/2)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels(
            { bitset_neighbors[1], bitset_neighbors[3] },
            bitset_mid_higher, *this, node_coarse_innermost, arr_ptr_layer, ptr_coarse_outer);

        // diagonal (0 + x/2, 0 + dy/2)
        if (bitset_neighbors[0].second == 1 || bitset_neighbors[1].second == 1
            || bitset_neighbors[2].second == 1 || bitset_neighbors[3].second == 1) {
            ptr_mid_layer->insert({ FindXNeg(bitset_mid_higher), kFlag0_ });
        }
    }
}
/**
* @brief   function to identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node (lower level)
* @param[in] node_coarse_interface_innermost  nodes on the innermost coarse to fine interface
* @param[in] node_exist_current existing nodes at current level
* @param[in] node_exist_lower   existing nodes at lower level
* @param[out] ptr_inner_layer map storing nodes on the inner layer
* @param[out] ptr_mid_layer map storing nodes on the middle layer
* @param[out] ptr_outer_layer map storing nodes on the outer layer
* @param[out] ptr_node_coarse_interface_outer  pointer to map storing nodes on outer coarse to fine interface
*/
void GridManager2D::IdentifyInnermostInterfaceForACell(const DefSFBitset sfbitset_in,
    const DefMap<DefInt>& node_coarse_interface_innermost,
    const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
    const DefMap<DefInt>& node_exist_lower,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_node_coarse_interface_outer) {
    std::array<DefSFBitset, 4> bitset_neighbors;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
        ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    bool belong_to_cell =  SFBitsetBelongToOneCell<DefInt>(
        sfbitset_in, node_exist_lower, &bitset_neighbors);
    DefSFBitset bitset_mid_higher;
    if (belong_to_cell) {
        // mid (0 + dx/2, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdgeInnermost(
            { bitset_neighbors[0], bitset_neighbors[1] },
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid (0, 0 + dy/2)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdgeInnermost(
            { bitset_neighbors[0], bitset_neighbors[2] },
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid (0 + dx/2, y)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdgeInnermost(
            { bitset_neighbors[2], bitset_neighbors[3] },
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid (x, 0 + dy/2)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdgeInnermost(
            { bitset_neighbors[1], bitset_neighbors[3] },
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0 + dy/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost(
            { bitset_neighbors[0], bitset_neighbors[3] },
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
    }
}
/**
* @brief   function to check if the node at a given level coincides with the background node
* @param[in]  i_level   refinement level
* @param[in]  bitset_higher   space filling code at the given refinement level
* @param[out] ptr_bitset   space filling code at the background level (level 0)
*/
bool GridManager2D::CheckCoincideBackground(const DefInt i_level,
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
void GridManager2D::FindAllNodesInACellAtOneLevelLower(
    const std::vector<DefSFBitset> bitset_cell,
    std::vector<DefSFBitset>* const ptr_bitset_all) const {
    ptr_bitset_all->resize(9);
    // bitset_cell[0]:(0, 0); bitset_cell[1]:(+x, 0);
    // bitset_cell[2]:(0, +y); bitset_cell[3]:(+x, +y);
    DefSFBitset sfbitset_tmp;
    // bottom
    sfbitset_tmp = SFBitsetToOneHigherLevel(bitset_cell.at(0));
    ptr_bitset_all->at(0) = sfbitset_tmp;
    ptr_bitset_all->at(1) = FindXPos(sfbitset_tmp);
    ptr_bitset_all->at(2) = SFBitsetToOneHigherLevel(bitset_cell.at(1));
    // middle
    ptr_bitset_all->at(3) = FindYPos(ptr_bitset_all->at(0));
    ptr_bitset_all->at(4) = FindYPos(ptr_bitset_all->at(1));
    ptr_bitset_all->at(5) = FindYPos(ptr_bitset_all->at(2));
    // top
    sfbitset_tmp = SFBitsetToOneHigherLevel(bitset_cell.at(2));
    ptr_bitset_all->at(6) = sfbitset_tmp;
    ptr_bitset_all->at(7) = FindXPos(sfbitset_tmp);
    ptr_bitset_all->at(8) = SFBitsetToOneHigherLevel(bitset_cell.at(3));
}
/**
* @brief function to space filling code to n level lower
* @param[in]  n_level  number of levels need to be shrink.
* @param[in]  bitset_in   input space filling code.
* @return   space filling code at lower levels.
*/
DefSFBitset GridManager2D::NodeAtNLowerLevel(
    const DefInt n_level, const DefSFBitset& bitset_in) const {
    return SFBitsetToNLowerLevel(n_level, bitset_in);
}
/**
* @brief   function to find layer in the overlapping region based on the layer
*          at one level higher
* @param[in]  layer_high_level   layer in the overlapping region at high level
* @param[out] ptr_layer_low_level layer in the overlapping region at low level
*/
void GridManager2D::OverlapLayerFromHighToLow(
    const DefMap<DefInt>& layer_high_level,
    DefMap<DefInt>* const ptr_layer_low_level) {
    for (const auto& iter : layer_high_level) {
        if ((iter.first & k0SfBitsetCurrentLevelBits_) == 0) {
            ptr_layer_low_level->insert({
                SFBitsetToOneLowerLevel(iter.first), kFlag0_ });
        }
    }
}
/**
* @brief   function to find if a background node is in the offset region.
* @param[in] bitset_in space filling code of the node
* @return true if node is in the offset region, others false
*/
bool GridManager2D::CheckBackgroundOffset(const DefSFBitset& bitset_in) const {
    std::array<DefAmrLUint, 2> indices;
    SFBitsetComputeIndices(bitset_in, &indices);
    if ((indices[kXIndex] < k0MinIndexOfBackgroundNode_[kXIndex])
        || (indices[kYIndex] < k0MinIndexOfBackgroundNode_[kYIndex])) {
        return true;
    } else {
        return false;
    }
}
/**
* @brief   function to instantiate background node.
* @param[in] code_min the minimum space filling code at background level of current rank.
* @param[in] code_max the maximum space filling code at background level of current rank.
* @param[in] map_occupied node exist in grid at high refinement level.
*/
void GridManager2D::InstantiateBackgroundGrid(const DefSFCodeToUint code_min,
    const DefSFCodeToUint code_max, const DefMap<DefInt>& map_occupied) {
    DefSFBitset sfbitset_tmp;
    GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(0));
    DefSFCodeToUint i_code = code_min;
    DefInt flag_node;
    while (i_code <= code_max) {
        ResetIndicesExceedingDomain(k0MinIndexOfBackgroundNode_, k0MaxIndexOfBackgroundNode_, &i_code, &sfbitset_tmp);
        if (map_occupied.find(sfbitset_tmp) == map_occupied.end()) {
            InstantiateGridNode(sfbitset_tmp, &grid_info);
            flag_node = grid_info.CheckIfNodeOutsideCubicDomain(k0GridDims_, sfbitset_tmp, *this);
            (this->*ptr_func_insert_domain_boundary_)(flag_node, sfbitset_tmp, &grid_info);
        }
        ++i_code;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
