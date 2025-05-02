//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
#include "./auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to setup size of the computational domain.
* @param[in] domain_size maximum coordinates of the computational domain.
*/
void GridManager3D::SetDomainSize(const std::vector<DefReal>& domain_size) {
    if (domain_size.size() == 3) {
        k0DomainSize_.at(kXIndex) = domain_size.at(kXIndex);
        k0DomainSize_.at(kYIndex) = domain_size.at(kYIndex);
        k0DomainSize_.at(kZIndex) = domain_size.at(kZIndex);
    } else {
        LogManager::LogError("size of the input vector should be 3"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief function to setup grid spacing the background level.
* @param[in] domain_grid_size grid size of the computational domain.
*/
void GridManager3D::SetDomainGridSize(const std::vector<DefReal>& domain_grid_size) {
    if (domain_grid_size.size() == 1) {
        k0DomainDx_.at(kXIndex) = domain_grid_size.at(kXIndex);
    } else if (domain_grid_size.size() == 2) {
        k0DomainDx_.at(kXIndex) = domain_grid_size.at(kXIndex);
        k0DomainDx_.at(kYIndex) = domain_grid_size.at(kYIndex);
    } else if (domain_grid_size.size() == 3) {
        k0DomainDx_.at(kXIndex) = domain_grid_size.at(kXIndex);
        k0DomainDx_.at(kYIndex) = domain_grid_size.at(kYIndex);
        k0DomainDx_.at(kZIndex) = domain_grid_size.at(kZIndex);
    } else {
        LogManager::LogError("size of the input vector should be 1 2, or 3"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
* @brief function to setup and check grid related parameters.
*/
void GridManager3D::SetupDependentGridParameters() {
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
    k0SpaceBackground_ = { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex], k0DomainDx_[kZIndex] };

    // calculate number of background nodes in each direction
    k0MaxIndexOfBackgroundNode_ = {
        static_cast<DefAmrLUint>(k0DomainSize_[kXIndex]
        / k0DomainDx_[kXIndex] + kEps) + k0MinIndexOfBackgroundNode_[kXIndex],
        static_cast<DefAmrLUint>(k0DomainSize_[kYIndex]
        / k0DomainDx_[kYIndex] + kEps) + k0MinIndexOfBackgroundNode_[kYIndex],
        static_cast<DefAmrLUint>(k0DomainSize_[kZIndex]
        / k0DomainDx_[kZIndex] + kEps) + k0MinIndexOfBackgroundNode_[kZIndex]};

    SFBitsetSetMinAndMaxBounds(
        k0MinIndexOfBackgroundNode_, k0MaxIndexOfBackgroundNode_);

    // set offsets
    k0RealMin_[kXIndex] = k0MinIndexOfBackgroundNode_[kXIndex] * k0DomainDx_[kXIndex];
    k0RealMin_[kYIndex] = k0MinIndexOfBackgroundNode_[kYIndex] * k0DomainDx_[kYIndex];
    k0RealMin_[kZIndex] = k0MinIndexOfBackgroundNode_[kZIndex] * k0DomainDx_[kZIndex];

    k0SFBitsetDomainMin_ = SFBitsetEncoding({k0MinIndexOfBackgroundNode_[kXIndex],
        k0MinIndexOfBackgroundNode_[kYIndex], k0MinIndexOfBackgroundNode_[kZIndex]});
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
    DefInt bit_max = kSFBitsetBit / k0GridDims_ - k0MaxLevel_;
    DefSizet index_max = TwoPowerN(bit_max);
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
* @brief   function to check node on which domain boundaries.
* @param[in] sfbitset_in space filling code of a given node.
* @param[in] sfbitset_min bitset corresponding to the minimum coordinate in each direction.
* @param[in] sfbitset_max bitset corresponding to the maximum coordinate in each direction.
* @return flag indicate node on which domain boundaries, 1: x min, 2: x max, 4: y min, 8: y max
*/
void GridManager3D::CalDomainBoundsAtGivenLevel(const DefInt i_level,
    std::vector<DefSFBitset>* const ptr_domain_min, std::vector<DefSFBitset>* const ptr_domain_max) const {
    ptr_domain_min->resize(3);
    ptr_domain_max->resize(3);
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
    if (SFBitsetMin_.at(kZIndex) == ~DefSFCodeToUint(0)) {
         LogManager::LogError("SFBitsetMin_ should be set before calling CalDomainBoundsAtGivenLevel"
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        ptr_domain_min->at(kZIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetMin_.at(kZIndex));
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
    if (SFBitsetMax_.at(kZIndex) == 0) {
         LogManager::LogError("SFBitsetMax_ should be set before calling CalDomainBoundsAtGivenLevel"
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        ptr_domain_max->at(kZIndex) = SFBitsetToNHigherLevel(i_level, SFBitsetMax_.at(kZIndex));
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
    const DefInt i_level, const DefSFBitset& sfbitset_in,
    std::vector<DefAmrLUint>* const ptr_vec_extend_neg,
    std::vector<DefAmrLUint>* const ptr_vec_extend_pos) const {
    //  extended layer at the background refinement level
    DefAmrLUint two_power_i_level = static_cast<DefAmrLUint>(TwoPowerN(i_level));
    DefAmrLUint index_xmin = k0MinIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrLUint index_xmax = k0MaxIndexOfBackgroundNode_[kXIndex] * two_power_i_level;
    DefAmrLUint index_ymin = k0MinIndexOfBackgroundNode_[kYIndex] * two_power_i_level;
    DefAmrLUint index_ymax = k0MaxIndexOfBackgroundNode_[kYIndex] * two_power_i_level;
    DefAmrLUint index_zmin = k0MinIndexOfBackgroundNode_[kZIndex] * two_power_i_level;
    DefAmrLUint index_zmax = k0MaxIndexOfBackgroundNode_[kZIndex] * two_power_i_level;

    std::array<DefAmrLUint, 3> indices;
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
* @brief   function to check node on which domain boundaries.
* @param[in] i_level refinement level.
* @param[in] sfbitset_in space filling code of a given node.
* @return flag indicate node on which domain boundaries, 1: x min, 2: x max, 4: y min, 8: y max, 16: z min, 32: z max
*/
int GridManager3D::CheckNodeOnDomainBoundary(
    const DefInt i_level, const DefSFBitset& sfbitset_in) const {
    int node_status = 0;
    std::array<DefSFBitset, 3> sfbitset_min, sfbitset_max;
    sfbitset_min[kXIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMin_[kXIndex]);
    sfbitset_min[kYIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMin_[kYIndex]);
    sfbitset_min[kZIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMin_[kZIndex]);
    sfbitset_max[kXIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMax_[kXIndex]);
    sfbitset_max[kYIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMax_[kYIndex]);
    sfbitset_max[kZIndex] = SFBitsetToNHigherLevel(i_level, SFBitsetMax_[kZIndex]);
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_min[kXIndex]) {
        node_status |= 1;
    }
    if ((sfbitset_in & k0SFBitsetTakeXRef_[kRefCurrent_])
        == sfbitset_max[kXIndex]) {
        node_status |= 2;
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_min[kYIndex]) {
        node_status |= 4;
    }
    if ((sfbitset_in & k0SFBitsetTakeYRef_[kRefCurrent_])
        == sfbitset_max[kYIndex]) {
        node_status |= 8;
    }
    if ((sfbitset_in & k0SFBitsetTakeZRef_[kRefCurrent_])
        == sfbitset_min[kZIndex]) {
        node_status |= 16;
    }
    if ((sfbitset_in & k0SFBitsetTakeZRef_[kRefCurrent_])
        == sfbitset_max[kZIndex]) {
        node_status |= 32;
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
bool GridManager3D::CheckNodeNotOutsideDomainBoundary(const DefSFBitset& sfbitset_in,
    const std::vector<DefSFCodeToUint>& sfbitset_min,
    const std::vector<DefSFCodeToUint>& sfbitset_max) const {
    if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeXRef_[kRefCurrent_]) < sfbitset_min.at(kXIndex)) {
        return false;
    } else if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeXRef_[kRefCurrent_]) > sfbitset_max.at(kXIndex)) {
        return false;
    }
    if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeYRef_[kRefCurrent_])< sfbitset_min.at(kYIndex)) {
        return false;
    } else if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeYRef_[kRefCurrent_]) > sfbitset_max.at(kYIndex)) {
        return false;
    }
     if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeZRef_[kRefCurrent_]) < sfbitset_min.at(kZIndex)) {
        return false;
    } else if (SFBitsetoSFCode(sfbitset_in&k0SFBitsetTakeZRef_[kRefCurrent_]) > sfbitset_max.at(kZIndex)) {
        return false;
    }
    return true;
}
/**
* @brief   function to find space filling code of nodes on a cell (3D)
* @param[in]  sfbitset_in   space filling code of the node at the corner of a cell
* @param[in]  node_exist grid containing nodes at the same level
* @param[out] ptr_bitsets nodes of the cell
*/
bool GridManager3D::NodesBelongToOneCell(const DefSFBitset bitset_in,
    const DefMap<DefInt>& node_exist,
    std::vector<DefSFBitset>* const ptr_bitsets) const {
    bool bool_cell;
    ptr_bitsets->clear();
    std::array<DefSFBitset, 8> bitset_cell;
    bool_cell = SFBitsetBelongToOneCell(bitset_in, node_exist, &bitset_cell);
    ptr_bitsets->assign(bitset_cell.begin(), bitset_cell.end());
    return bool_cell;
}
/**
* @brief   function to find space filling code of nodes on a surface (3D)
* @param[in]  sfbitset_in   space filling code of the node at the corner of a surface
* @param[in]  dir_norm   normal to the surface: kXIndex, kYIndex or kZIndex
* @param[in]  map_node_exist existing nodes
* @param[out] ptr_sfbitset all nodes on the surface
*/
bool GridManager3D::NodesBelongToOneSurfAtHigherLevel(const DefSFBitset sfbitset_in,
    const DefInt dir_norm, const DefMap<DefInt>& map_node_exist,
    std::vector<DefSFBitset>* const ptr_sfbitset) const {
    ptr_sfbitset->clear();
    DefSFBitset sfbitset_tmp1, sfbitset_tmp2, sfbitset_tmp3, sfbitset_mid2;
    bool bool_vertex1, bool_vertex2, bool_vertex3;
    switch (dir_norm) {
    case kXIndex:
        if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_in));
        // (0, +y, 0)
        sfbitset_tmp1 = FindYPos(sfbitset_in);

        if (map_node_exist.find(sfbitset_tmp1) == map_node_exist.end()) {
            bool_vertex1 =  false;
        } else {
            bool_vertex1 = true;
            ptr_sfbitset->emplace_back(FindYPos(ptr_sfbitset->at(0)));
            ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp1));
        }
        // (0, 0, +z)
        sfbitset_tmp2 = FindZPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp2) == map_node_exist.end()) {
            bool_vertex2 =  false;
        } else {
            bool_vertex2 = true;
            ptr_sfbitset->emplace_back(FindZPos(ptr_sfbitset->at(0)));
            sfbitset_mid2 = SFBitsetToOneHigherLevel(sfbitset_tmp2);
            ptr_sfbitset->emplace_back(sfbitset_mid2);
        }
        // (0, +y, +z)
        sfbitset_tmp3 = FindZPos(sfbitset_tmp1);
        if (map_node_exist.find(sfbitset_tmp3) == map_node_exist.end()) {
            bool_vertex3 =  false;
        } else {
            bool_vertex3 = true;
            if (bool_vertex1) {
                ptr_sfbitset->emplace_back(FindZPos(ptr_sfbitset->at(2)));
            }
            if (bool_vertex2) {
                ptr_sfbitset->emplace_back(FindYPos(sfbitset_mid2));
            }
            if (bool_vertex1&&bool_vertex2) {
                ptr_sfbitset->emplace_back(FindZPos(ptr_sfbitset->at(1)));
            }
            ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp3));
        }
        if (bool_vertex1&&bool_vertex2&&bool_vertex3) {
            return true;
        } else {
            return false;
        }
        break;
    case kYIndex:
        if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_in));
        // (+x, 0, 0)
        sfbitset_tmp1 = FindXPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp1) == map_node_exist.end()) {
            bool_vertex1 =  false;
        } else {
            bool_vertex1 = true;
            ptr_sfbitset->emplace_back(FindXPos(ptr_sfbitset->at(0)));
            ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp1));
        }
        // (0, 0, +z)
        sfbitset_tmp2 = FindZPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp2) == map_node_exist.end()) {
            bool_vertex2 =  false;
        } else {
            bool_vertex2 = true;
            ptr_sfbitset->emplace_back(FindZPos(ptr_sfbitset->at(0)));
            sfbitset_mid2 = SFBitsetToOneHigherLevel(sfbitset_tmp2);
            ptr_sfbitset->emplace_back(sfbitset_mid2);
        }
        // (+x, 0, +z)
        sfbitset_tmp3 = FindZPos(sfbitset_tmp1);
        if (map_node_exist.find(sfbitset_tmp3) == map_node_exist.end()) {
            bool_vertex3 =  false;
        } else {
            bool_vertex3 = true;
            if (bool_vertex1) {
                ptr_sfbitset->emplace_back(FindZPos(ptr_sfbitset->at(2)));
            }
            if (bool_vertex2) {
                ptr_sfbitset->emplace_back(FindXPos(sfbitset_mid2));
            }
            if (bool_vertex1&&bool_vertex2) {
                ptr_sfbitset->emplace_back(FindZPos(ptr_sfbitset->at(1)));
            }
            ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp3));
        }
        if (bool_vertex1&&bool_vertex2&&bool_vertex3) {
            return true;
        } else {
            return false;
        }
        break;
    case kZIndex:
        if (map_node_exist.find(sfbitset_in) == map_node_exist.end()) {
            return false;
        }
        ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_in));
        // (+x, 0, 0)
        sfbitset_tmp1 = FindXPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp1) == map_node_exist.end()) {
            bool_vertex1 =  false;
        } else {
            bool_vertex1 = true;
            ptr_sfbitset->emplace_back(FindXPos(ptr_sfbitset->at(0)));
            ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp1));
        }
        // (0, +y, 0)
        sfbitset_tmp2 = FindYPos(sfbitset_in);
        if (map_node_exist.find(sfbitset_tmp2) == map_node_exist.end()) {
            bool_vertex2 =  false;
        } else {
            bool_vertex2 = true;
            ptr_sfbitset->emplace_back(FindYPos(ptr_sfbitset->at(0)));
            sfbitset_mid2 = SFBitsetToOneHigherLevel(sfbitset_tmp2);
            ptr_sfbitset->emplace_back(sfbitset_mid2);
        }
        // (+x, +y, 0)
        sfbitset_tmp3 = FindYPos(sfbitset_tmp1);
        if (map_node_exist.find(sfbitset_tmp3) == map_node_exist.end()) {
            bool_vertex3 =  false;
        } else {
            bool_vertex3 = true;
            if (bool_vertex1) {
                ptr_sfbitset->emplace_back(FindYPos(ptr_sfbitset->at(2)));
            }
            if (bool_vertex2) {
                ptr_sfbitset->emplace_back(FindXPos(sfbitset_mid2));
            }
            if (bool_vertex1&&bool_vertex2) {
                ptr_sfbitset->emplace_back(FindYPos(ptr_sfbitset->at(1)));
            }
            ptr_sfbitset->emplace_back(SFBitsetToOneHigherLevel(sfbitset_tmp3));
        }
        if (bool_vertex1&&bool_vertex2&&bool_vertex3) {
            return true;
        } else {
            return false;
        }
        break;
    default:
        return false;
    }
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
* @param[in] node_coarse_interface nodes on on interface layer of coarser grid
* @param[in] node_exist_lower   existing nodes at lower level
* @param[out]  ptr_inner_layer map store nodes on the inner layer
* @param[out]  ptr_mid_layer map store nodes on the middle layer
* @param[out]  ptr_outer_layer map store nodes on the outer layer
*/
// in node_exist_lower, only nodes on the refinement interface exist
void GridManager3D::IdentifyInterfaceForACell(const DefSFBitset bitset_in,
    const DefMap<DefInt>& node_coarse_interface,
    const DefMap<DefInt>& node_exist_lower,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer) {
    // bitset_neighbors[0]:(0, 0, 0);   bitset_neighbors[1]:(+x, 0, 0);
    // bitset_neighbors[2]:(0, +y, 0);  bitset_neighbors[3]:(+x, +y, 0);
    // bitset_neighbors[4]:(0, 0, +z);  bitset_neighbors[5]:(+x, 0, +z);
    // bitset_neighbors[6]:(0, +y, +z); bitset_neighbors[7]:(+x, +y, +z).
    std::array<DefSFBitset, 8> bitset_neighbors;
    DefSFBitset bitset_mid_higher;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
    ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    bool belong_to_cell =  SFBitsetBelongToOneCell<DefInt>(
                bitset_in, node_exist_lower, &bitset_neighbors);
    if (belong_to_cell) {
        // bottom surface
        // edge (0 + dx/2, 0, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[1] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid edge (0, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[2] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid edge (0 + dx/2, y, 0)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[2], bitset_neighbors[3] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid edge (x, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[1], bitset_neighbors[3] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (0 + x/2, 0 + dy/2, 0)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[0], bitset_neighbors[1],
              bitset_neighbors[2], bitset_neighbors[3] },
            bitset_mid_higher, node_coarse_interface, arr_ptr_layer);

        // top surface
        // mid edge (0 + dx/2, 0, z)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[4], bitset_neighbors[5] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid edge (0, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[4], bitset_neighbors[6] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid edge (0 + dx/2, y, z)
        bitset_mid_higher = FindXPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[6]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[6], bitset_neighbors[7] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // mid edge (x, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[5]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[5], bitset_neighbors[7] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (0 + x/2, 0 + dy/2, z)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[4], bitset_neighbors[5],
              bitset_neighbors[6], bitset_neighbors[7] },
            bitset_mid_higher, node_coarse_interface, arr_ptr_layer);

        // middle surface
        // mid edge (0, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[0], bitset_neighbors[4] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (0 + x/2, 0, 0 + dz/2)
        bitset_mid_higher = FindXPos(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[0], bitset_neighbors[1],
              bitset_neighbors[5], bitset_neighbors[4] },
            bitset_mid_higher, node_coarse_interface, arr_ptr_layer);
        // mid edge (x, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[1], bitset_neighbors[5] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (x, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYPos(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[1], bitset_neighbors[3],
              bitset_neighbors[7], bitset_neighbors[5] },
            bitset_mid_higher, node_coarse_interface, arr_ptr_layer);
        // mid edge (0, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[2], bitset_neighbors[6] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (0, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[0], bitset_neighbors[2],
              bitset_neighbors[6], bitset_neighbors[4] },
            bitset_mid_higher, node_coarse_interface, arr_ptr_layer);
        // mid edge (x, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(
            SFBitsetToOneHigherLevel(bitset_neighbors[3]));
        IdentifyInterfaceNodeOnEdge(
            { bitset_neighbors[3], bitset_neighbors[7] },
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer);
        // diagonal (0 + dx/2, y, 0 + dz/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonal(
            { bitset_neighbors[2], bitset_neighbors[3],
              bitset_neighbors[7], bitset_neighbors[6] },
            bitset_mid_higher, node_coarse_interface, arr_ptr_layer);

        // center (0 + dx/2, 0 + dy/2, 0 + dz/2)
        ptr_mid_layer->insert({ FindYNeg(bitset_mid_higher), kFlag0_ });
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
void GridManager3D::IdentifyInterfaceForACell(const DefSFBitset bitset_in,
    const DefMap<DefInt>& node_coarse_interface_previous,
    const DefMap<DefInt>& node_coarse_interface_inner,
    const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
    const DefMap<std::unique_ptr<GridNode>>& node_exist_lower,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_node_coarse_interface_outer) {
    // bitset_neighbors[0]:(0, 0, 0);   bitset_neighbors[1]:(+x, 0, 0);
    // bitset_neighbors[2]:(0, +y, 0);  bitset_neighbors[3]:(+x, +y, 0);
    // bitset_neighbors[4]:(0, 0, +z);  bitset_neighbors[5]:(+x, 0, +z);
    // bitset_neighbors[6]:(0, +y, +z); bitset_neighbors[7]:(+x, +y, +z).
    std::array<DefSFBitset, 8> bitset_neighbors;
    DefSFBitset bitset_mid_higher;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
    ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    DefMap<DefInt> discard_layer;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_mid_layer = {
        &discard_layer, ptr_mid_layer, &discard_layer };
    bool belong_to_cell =  SFBitsetBelongToOneCell<std::unique_ptr<GridNode>>(
        bitset_in, node_exist_lower, &bitset_neighbors);
    if (belong_to_cell) {
        // bottom surface
        // edge (0 + dx/2, 0, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[1]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0 + dx/2, y, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[2], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0 + dy/2, 0)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);

        // top surface
        // mid edge (0 + dx/2, 0, z)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[4], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[4], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0 + dx/2, y, z)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[6]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[6], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[5]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[5], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0 + dy/2, z)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[4], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[5], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);

        // middle surface
        // mid edge (0, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0, 0 + dz/2)
        bitset_mid_higher = FindXPos(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (x, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYPos(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[1], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[3], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[2], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[2], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[3]));
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[3], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + dx/2, y, 0 + dz/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[2], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[3], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_layer, ptr_node_coarse_interface_outer);

        // center (0 + dx/2, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[0], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_mid_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdge({bitset_neighbors[3], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_previous, node_coarse_interface_inner,
            node_exist_current, arr_ptr_mid_layer, ptr_node_coarse_interface_outer);
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
void GridManager3D::IdentifyInnermostInterfaceForACell(const DefSFBitset sfbitset_in,
    const DefMap<DefInt>& node_coarse_interface_innermost,
    const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
    const DefMap<DefInt>& node_exist_lower,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_node_coarse_interface_outer) {
    std::array<DefSFBitset, 8> bitset_neighbors;
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
        IdentifyInterfaceNodeOnEdgeInnermost(
            { bitset_neighbors[1], bitset_neighbors[2] },
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);

        // bottom surface
        // edge (0 + dx/2, 0, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[1]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0 + dx/2, y, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[2], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[1], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0 + dy/2, 0)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[1], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);

        // top surface
        // mid edge (0 + dx/2, 0, z)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[4], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[4]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[4], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0 + dx/2, y, z)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[6]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[6], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[5]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[5], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0 + dy/2, z)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[4], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[5], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);

        // middle surface
        // mid edge (0, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[0]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + x/2, 0, 0 + dz/2)
        bitset_mid_higher = FindXPos(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[1], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[1]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[1], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (x, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYPos(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[1], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[3], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (0, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[2]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[2], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[2], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // mid edge (x, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[3]));
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[3], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        // diagonal (0 + dx/2, y, 0 + dz/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[2], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[3], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);

        // center (0 + dx/2, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[0], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
        IdentifyInterfaceNodeOnEdgeInnermost({bitset_neighbors[3], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface_innermost, node_exist_current,
            arr_ptr_layer, ptr_node_coarse_interface_outer);
    }
}
/**
* @brief   function to identify interface for a given cell
* @param[in] sfbitset_in   space filling code of the given node at lower level(level - 1)
* @param[in] node_coarse_interface nodes on the interface layer of coarser grid
* @param[in] map_node_exist   existing fine nodes at lower level (level - 1)
* @param[in] map_exist_coarse   existing coarse nodes at level - 1
* @param[out] ptr_inner_layer map storing nodes on the inner layer
* @param[out] ptr_mid_layer map storing nodes on the middle layer
* @param[out] ptr_outer_layer map storing nodes on the outer layer
*/
void GridManager3D::IdentifyInterfaceForACellAcrossTwoLevels(
     const DefSFBitset bitset_in, const DefMap<DefInt>& node_coarse_interface,
    const DefMap<DefInt>& map_node_exist, const DefMap<DefInt>& map_exist_coarse,
    DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
    DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* ptr_coarse_outer) {
    std::array<std::pair<DefSFBitset, DefInt>, 8> bitset_neighbors;
    DefSFBitset bitset_mid_higher, bitset_tmp;
    std::array<DefMap<DefInt>* const, 3> arr_ptr_layer = {
        ptr_inner_layer, ptr_mid_layer, ptr_outer_layer };
    bool belong_to_cell =  SFBitsetBelongToOneCellAcrossTwoLevels(
        bitset_in, map_node_exist, map_exist_coarse, &bitset_neighbors);
    DefMap<DefInt> discard_layer;
    if (belong_to_cell) {
        // bottom surface
        // edge (0 + dx/2, 0, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[0].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[0], bitset_neighbors[1]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // mid edge (0, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[0].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[0], bitset_neighbors[2]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // mid edge (0 + dx/2, y, 0)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[2].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[2], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // mid edge (x, 0 + dy/2, 0)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[1].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[1], bitset_neighbors[3]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // diagonal (0 + x/2, 0 + dy/2, 0)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonalAcrossTwoLevels({bitset_neighbors[0], bitset_neighbors[1],
            bitset_neighbors[2], bitset_neighbors[3]}, bitset_mid_higher, node_coarse_interface,
            arr_ptr_layer, ptr_coarse_outer);

        // top surface
        // mid edge (0 + dx/2, 0, z)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[4].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[4], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // mid edge (0, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[4].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[4], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // mid edge (0 + dx/2, y, z)
        bitset_mid_higher = FindXPos(SFBitsetToOneHigherLevel(bitset_neighbors[6].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[6], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // mid edge (x, 0 + dy/2, z)
        bitset_mid_higher = FindYPos(SFBitsetToOneHigherLevel(bitset_neighbors[5].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[5], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // diagonal (0 + x/2, 0 + dy/2, z)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonalAcrossTwoLevels({bitset_neighbors[4], bitset_neighbors[5],
            bitset_neighbors[6], bitset_neighbors[7]}, bitset_mid_higher, node_coarse_interface,
            arr_ptr_layer, ptr_coarse_outer);

        // middle surface
        // mid edge (0, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[0].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[0], bitset_neighbors[4]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // diagonal (0 + x/2, 0, 0 + dz/2)
        bitset_mid_higher = FindXPos(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonalAcrossTwoLevels({bitset_neighbors[0], bitset_neighbors[1],
            bitset_neighbors[4], bitset_neighbors[5]}, bitset_mid_higher, node_coarse_interface,
            arr_ptr_layer, ptr_coarse_outer);
        // mid edge (x, 0, 0 + dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[1].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[1], bitset_neighbors[5]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // diagonal (x, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYPos(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonalAcrossTwoLevels({bitset_neighbors[1], bitset_neighbors[3],
            bitset_neighbors[5], bitset_neighbors[7]}, bitset_mid_higher, node_coarse_interface,
            arr_ptr_layer, ptr_coarse_outer);
        // mid edge (0, y, 0+ dz/2)
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[2].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[2], bitset_neighbors[6]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // diagonal (0, 0 + dy/2, 0 + dz/2)
        bitset_mid_higher = FindYNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonalAcrossTwoLevels({bitset_neighbors[0], bitset_neighbors[2],
            bitset_neighbors[4], bitset_neighbors[6]}, bitset_mid_higher, node_coarse_interface,
            arr_ptr_layer, ptr_coarse_outer);
        bitset_mid_higher = FindZPos(SFBitsetToOneHigherLevel(bitset_neighbors[3].first));
        IdentifyInterfaceNodeOnEdgeAcrossTwoLevels({bitset_neighbors[3], bitset_neighbors[7]},
            bitset_mid_higher, *this, node_coarse_interface, arr_ptr_layer, ptr_coarse_outer);
        // diagonal (0 + dx/2, y, 0 + dz/2)
        bitset_mid_higher = FindXNeg(bitset_mid_higher);
        IdentifyInterfaceNodeDiagonalAcrossTwoLevels({bitset_neighbors[2], bitset_neighbors[3],
            bitset_neighbors[6], bitset_neighbors[7]}, bitset_mid_higher, node_coarse_interface,
            arr_ptr_layer, ptr_coarse_outer);

        // center (0 + dx/2, 0 + dy/2, 0 + dz/2)
        if (bitset_neighbors[0].second == 1 || bitset_neighbors[3].second == 1
            || bitset_neighbors[4].second == 1 || bitset_neighbors[7].second == 1) {
            ptr_mid_layer->insert({FindYNeg(bitset_mid_higher), kFlag0_});
            if (node_coarse_interface.find(bitset_neighbors[0].first) == node_coarse_interface.end()) {
                ptr_coarse_outer->insert({bitset_neighbors[0].first, kFlag0_});
            }
            if (node_coarse_interface.find(bitset_neighbors[3].first) == node_coarse_interface.end()) {
                ptr_coarse_outer->insert({bitset_neighbors[3].first, kFlag0_});
            }
            if (node_coarse_interface.find(bitset_neighbors[4].first) == node_coarse_interface.end()) {
                ptr_coarse_outer->insert({bitset_neighbors[4].first, kFlag0_});
            }
            if (node_coarse_interface.find(bitset_neighbors[7].first) == node_coarse_interface.end()) {
                ptr_coarse_outer->insert({bitset_neighbors[7].first, kFlag0_});
            }
        }
    }
}
/**
* @brief   function to Identify types of interface diagonal node
* @param[in]  flag_interface   flag of marked interface node
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_coarse_interface   existing nodes on the outermost layer of coarse grid
* @param[out]  arr_ptr_layer pointer to map store interface layers
*/
void GridManager3D::IdentifyInterfaceNodeDiagonal(
    const std::array<DefSFBitset, 4>& arr_bitset_lower,
    const DefSFBitset bitset_center_higher,
    const DefMap<DefInt>& node_coarse_interface,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer) {
    DefInt node0_flag = node_coarse_interface.find(arr_bitset_lower[0]) != node_coarse_interface.end(),
        node1_flag = node_coarse_interface.find(arr_bitset_lower[1]) != node_coarse_interface.end(),
        node2_flag = node_coarse_interface.find(arr_bitset_lower[2]) != node_coarse_interface.end(),
        node3_flag = node_coarse_interface.find(arr_bitset_lower[3]) != node_coarse_interface.end();
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
* @brief   function to Identify types of interface diagonal node
* @param[in]  flag_interface   flag of marked interface node
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_coarse_interface   existing nodes on the outermost layer of coarse grid
* @param[out]  arr_ptr_layer pointer to map store interface layers
*/
void GridManager3D::IdentifyInterfaceNodeDiagonalAcrossTwoLevels(
    const std::array<std::pair<DefSFBitset, DefInt>, 4>& arr_bitset_lower,
    const DefSFBitset bitset_center_higher,
    const DefMap<DefInt>& node_coarse_interface,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer, DefMap<DefInt>* const ptr_coarse_outer) {
    const DefInt& flag_node_0 = arr_bitset_lower.at(0).second, flag_node_1 = arr_bitset_lower.at(1).second,
        flag_node_2 = arr_bitset_lower.at(2).second, flag_node_3 = arr_bitset_lower.at(3).second;
    if (((flag_node_0 == 1) || (flag_node_1 == 1)|| (flag_node_2 == 1)|| (flag_node_3 == 1))
        && ((flag_node_0 > 0) && (flag_node_1 > 0) && (flag_node_2 > 0) && (flag_node_3 > 0))) {
        bool node0_interface = node_coarse_interface.find(arr_bitset_lower[0].first) != node_coarse_interface.end(),
            node1_interface = node_coarse_interface.find(arr_bitset_lower[1].first) != node_coarse_interface.end(),
            node2_interface = node_coarse_interface.find(arr_bitset_lower[2].first) != node_coarse_interface.end(),
            node3_interface = node_coarse_interface.find(arr_bitset_lower[3].first) != node_coarse_interface.end();
        if (node0_interface == node1_interface && node1_interface == node2_interface
            && node2_interface == node3_interface) {
            if (node0_interface) {
                arr_ptr_layer[2]->insert({ bitset_center_higher, kFlag0_ });
            } else {
                arr_ptr_layer[0]->insert({ bitset_center_higher, kFlag0_ });
                ptr_coarse_outer->insert({ arr_bitset_lower[0].first, kFlag0_ });
                ptr_coarse_outer->insert({ arr_bitset_lower[1].first, kFlag0_ });
                ptr_coarse_outer->insert({ arr_bitset_lower[2].first, kFlag0_ });
                ptr_coarse_outer->insert({ arr_bitset_lower[3].first, kFlag0_ });
            }
        } else {
            arr_ptr_layer[1]->insert({ bitset_center_higher, kFlag0_ });
            if (!node0_interface) {
                ptr_coarse_outer->insert({ arr_bitset_lower[0].first, kFlag0_ });
            }
            if (!node1_interface) {
                ptr_coarse_outer->insert({ arr_bitset_lower[1].first, kFlag0_ });
            }
            if (!node2_interface) {
                ptr_coarse_outer->insert({ arr_bitset_lower[2].first, kFlag0_ });
            }
            if (!node3_interface) {
                ptr_coarse_outer->insert({ arr_bitset_lower[3].first, kFlag0_ });
            }
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
bool GridManager3D::CheckCoincideBackground(const DefInt i_level,
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
void GridManager3D::FindAllNodesInACellAtOneLevelLower(
    const std::vector<DefSFBitset> bitset_cell,
    std::vector<DefSFBitset>* const ptr_bitset_all) const {
    ptr_bitset_all->resize(27);
    // bitset_cell[0]:(0, 0, 0);   bitset_cell[1]:(+x, 0, 0);
    // bitset_cell[2]:(0, +y, 0);  bitset_cell[3]:(+x, +y, 0);
    // bitset_cell[4]:(0, 0, +z);  bitset_cell[5]:(+x, 0, +z);
    // bitset_cell[6]:(0, +y, +z); bitset_cell[7]:(+x, +y, +z).
    DefSFBitset sfbitset_tmp;
    // bottom
    sfbitset_tmp = SFBitsetToOneHigherLevel(bitset_cell.at(0));
    ptr_bitset_all->at(0) = sfbitset_tmp;
    ptr_bitset_all->at(1) = FindXPos(sfbitset_tmp);
    ptr_bitset_all->at(2) = SFBitsetToOneHigherLevel(bitset_cell.at(1));
    ptr_bitset_all->at(3) = FindYPos(ptr_bitset_all->at(0));
    ptr_bitset_all->at(4) = FindYPos(ptr_bitset_all->at(1));
    ptr_bitset_all->at(5) = FindYPos(ptr_bitset_all->at(2));
    sfbitset_tmp = SFBitsetToOneHigherLevel(bitset_cell.at(2));
    ptr_bitset_all->at(6) = sfbitset_tmp;
    ptr_bitset_all->at(7) = FindXPos(sfbitset_tmp);
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
    sfbitset_tmp = SFBitsetToOneHigherLevel(bitset_cell.at(4));
    ptr_bitset_all->at(18) = sfbitset_tmp;
    ptr_bitset_all->at(19) = FindXPos(sfbitset_tmp);
    ptr_bitset_all->at(20) = SFBitsetToOneHigherLevel(bitset_cell.at(5));
    ptr_bitset_all->at(21) = FindYPos(ptr_bitset_all->at(18));
    ptr_bitset_all->at(22) = FindYPos(ptr_bitset_all->at(19));
    ptr_bitset_all->at(23) = FindYPos(ptr_bitset_all->at(20));
    sfbitset_tmp = SFBitsetToOneHigherLevel(bitset_cell.at(6));
    ptr_bitset_all->at(24) = sfbitset_tmp;
    ptr_bitset_all->at(25) = FindXPos(sfbitset_tmp);
    ptr_bitset_all->at(26) = SFBitsetToOneHigherLevel(bitset_cell.at(7));
}
/**
* @brief function to space filling code to n level lower
* @param[in]  n_level  number of levels need to be shrink.
* @param[in]  bitset_in   input space filling code.
* @return   space filling code at lower levels.
*/
DefSFBitset GridManager3D::NodeAtNLowerLevel(
    const DefInt n_level, const DefSFBitset& bitset_in) const {
    return SFBitsetToNLowerLevel(n_level, bitset_in);
}
/**
* @brief   function to find layer in the overlapping region based on the layer
*          at one level higher
* @param[in]  layer_high_level   layer in the overlapping region at high level
* @param[out] ptr_layer_low_level layer in the overlapping region at low level
*/
void GridManager3D::OverlapLayerFromHighToLow(
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
* @brief   function to find background node in the offset region.
* @param[in] bitset_in space filling code of the node
* @return if node is in the offset region
*/
bool GridManager3D::CheckBackgroundOffset(const DefSFBitset& bitset_in) const {
    std::array<DefAmrLUint, 3> indices;
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
void GridManager3D::InstantiateBackgroundGrid(const DefSFCodeToUint code_min,
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
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
