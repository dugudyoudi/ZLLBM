//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
* @note  functions from other managers will be called.
*/
#include <string>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
namespace grid {
/**
* @brief   function to initialize grid noes.
* @param[in] sfbitset_in   space filling code of the grid node.
* @note
*/
void GridInfoInterface::InitialGridNode(const DefSFBitset& sfbitset_in) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::GetInstance()->rank_id_;
#endif  // ENABLE_MPI
    io::LogWarning("Function for initializing grid node "
        "(GridInfoInterface::InitialGridNode) does not given");
}
/**
* @brief   function to create instance based on the type of solver.
* @note
*/
void GridManagerInterface::CreateGridInstance(const DefUint i_level,
    std::shared_ptr<SolverInterface> ptr_sovler) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::GetInstance()->rank_id_;
#endif  // ENABLE_MPI
}

/**
* @brief function to setup default grid related parameters.
*/
void GridManagerInterface::DefaultInitialization(const DefSizet max_level) {
    k0MaxLevel_ = max_level;
    k0IntExtend_ = std::vector<DefLUint>(max_level + 1, 8);
    k0XIntExtendPositive_ = std::vector<DefLUint>(max_level + 1, 0);
    k0XIntExtendNegative_ = std::vector<DefLUint>(max_level + 1, 0);
    k0YIntExtendPositive_ = std::vector<DefLUint>(max_level + 1, 0);
    k0YIntExtendNegative_ = std::vector<DefLUint>(max_level + 1, 0);
    k0IntInnerExtend_ = std::vector<DefLUint>(max_level + 1, 2);
}
/**
* @brief   function to check if the geometry is out of domain.
* @param[in]  coordinate_min  minimum coordinates.
* @param[in]  coordinate_min  maximum coordinates.
* @param[in]  domain_offset  offset of the background grid.
* @param[in]  domain_size  size of the computational domain.
* @return  0: coordinates does not exceed the computational domain; 1, 3: 
*          coordinates are less than the lower bound; 2, 4: coordinates
*          greater than the upper bound;
*/
int GridManagerInterface::CheckIfPointOutsideDomain(
    const std::array<DefReal, 2>& coordinate_min,
    const std::array<DefReal, 2>& coordinate_max,
    const std::array<DefReal, 2>& domain_offset,
    const std::array<DefReal, 2>& domain_size) {
    if (coordinate_min[kXIndex] < domain_offset[kXIndex]) {
        return 1;
    } else if (coordinate_min[kXIndex] > domain_offset[kXIndex]) {
        return 2;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        return 3;
    } else if (coordinate_min[kYIndex] > domain_offset[kYIndex]) {
        return 4;
    }  
    return 0;
}
/**
* @brief   function to check if the geometry is out of domain.
* @param[in]  coordinate_min  minimum coordinates.
* @param[in]  coordinate_min  maximum coordinates.
* @param[in]  domain_offset  offset of the background grid.
* @param[in]  domain_size  size of the computational domain.
* @return  0: coordinates does not exceed the computational domain; 1, 3, 5:
*          coordinates are less than the lower bound; 2, 4, 6: coordinates
*          greater than the upper bound;
*/
int GridManagerInterface::CheckIfPointOutsideDomain(
    const std::array<DefReal, 3>& coordinate_min,
    const std::array<DefReal, 3>& coordinate_max,
    const std::array<DefReal, 3>& domain_offset,
    const std::array<DefReal, 3>& domain_size) {
    if (coordinate_min[kXIndex] < domain_offset[kXIndex]) {
        return 1;
    } else if (coordinate_min[kXIndex] > domain_offset[kXIndex]) {
        return 2;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        return 3;
    } else if (coordinate_min[kYIndex] > domain_offset[kYIndex]) {
        return 4;
    }
    if (coordinate_min[kZIndex] < domain_offset[kZIndex]) {
        return 5;
    } else if (coordinate_min[kZIndex] > domain_offset[kZIndex]) {
        return 6;
    }
    return 0;
}
/**
* @brief   function to check if the geometry is out of domain (3D).
* @param[in]  sfbitset_lower_level  sfbitset at lower level.
* @param[in]  map_sfbitset_at_lower_level  existing nodes at lower level.
* @param[in]  ptr_grid_info pointer to grids at the current level.
*/
void GridManagerInterface::AddNodesInstanceBasedOnLowerLevel(
    const DefSFBitset& sfbitset_lower_level,
    const DefMap<DefUint>& map_sfbitset_at_lower_level,
    std::shared_ptr<GridInfoInterface> ptr_grid_info) {

//    DefSFBitset sfbitset_temp0, sfbitset_temp1, sfbitset_temp2;
//    GridNode node_temp = ptr_grid_info->grid_node_instance_;
//    DefMap<GridNode>* ptr_map_nodes = &ptr_grid_info->map_grid_node_;
//    DefSFBitset sfbitset_current =
//        SFBitsetToNHigherLevel(1, sfbitset_lower_level);
//    DefSFBitset sfbitset_center;
//    std::vector<DefSFBitset> vec_neighbour_sfbitset(k0NumNeighbours_, 0);
//    DefSizet size_vec = TwoPowerN(static_cast<DefSizet>(k0GridDims_));
//    std::vector<DefSFBitset> vec_sfbitset(size_vec, 0);
//
//    // pointer to dimesional related function
//    bool(*ptrfunc_bitset_belong_to_one_cell)(
//        const DefSFBitset&, const DefMap<DefUint>&,
//        std::vector<DefSFBitset>*) = nullptr;
//    void (*bitset_find_all_neighbours_ptr)(const DefSFBitset&,
//        std::vector<DefSFBitset>*) = nullptr;
//    if (k0GridDims_ == 2) {
//#ifndef DEBUG_DISABLE_2D_FUNCTIONS
//        ptrfunc_bitset_belong_to_one_cell =
//            &SFBitsetBelongToOneCell2D<DefUint>;
//        bitset_find_all_neighbours_ptr =
//            &SFBitsetFindAllNeighbours2D;
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//    } else {
//#ifndef DEBUG_DISABLE_3D_FUNCTIONS
//        ptrfunc_bitset_belong_to_one_cell =
//            SFBitsetBelongToOneCell3D<DefUint>;
//        bitset_find_all_neighbours_ptr =
//            &SFBitsetFindAllNeighbours3D;
//    }
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
//
//    // check if node given node is at the lower corner of a cell
//    if ((*ptrfunc_bitset_belong_to_one_cell)(sfbitset_lower_level,
//        map_sfbitset_at_lower_level, &vec_sfbitset)) {
//        // sfbitset of the node located at the cell center
//        sfbitset_center = sfbitset_current
//            | sfbitset_aux_.k0SfBItsetTakeCenter_;
//        (*bitset_find_all_neighbours_ptr)(
//            sfbitset_center, &vec_neighbour_sfbitset);
//        for (const auto& iter : vec_neighbour_sfbitset) {
//            if (ptr_map_nodes->find(iter) == ptr_map_nodes->end()) {
//                ptr_map_nodes->insert({ iter, node_temp });
//                ptr_grid_info->InitialGridNode(iter);
//            }
//        }
//    }
}
//#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
///**
//* @brief   function to find tracking node based on a given 2D geo point.
//* @param[in]  grid_space   grid spacing in x, and y direction.
//* @param[in]  geo_coordi  coordinate of geometry point.
//* @param[out] space filling codes of nodes containing the input geometry point.
//*/
//void GridManagerInterface::FindTrackingNodesContainingAGeometryCoordinate2D(
//    const std::vector<DefReal>& grid_space,
//    const std::vector<DefReal>& geo_coordi,
//    std::vector<DefSFBitset>* ptr_vec_sfbitset) {
//    DefLUint x_index = static_cast<DefLUint>
//        (geo_coordi[kXIndex] / grid_space[kXIndex] + kEps);
//    DefLUint y_index = static_cast<DefLUint>
//        (geo_coordi[kYIndex] / grid_space[kYIndex] + kEps);
//    ptr_vec_sfbitset->at(0) = SFBitsetEncoding2D({ x_index, y_index });
//    ptr_vec_sfbitset->at(1) = FindXPos2D(ptr_vec_sfbitset->at(0));
//    ptr_vec_sfbitset->at(2) = FindYPos2D(ptr_vec_sfbitset->at(1));
//    ptr_vec_sfbitset->at(3) = FindYPos2D(ptr_vec_sfbitset->at(0));
//}
///**
//* @brief   function to reset extened layers if the searching region
//           if ouside the computational domain (2D).
//* @param[in]  i_level   level of refinement.
//* @param[in]  sfbitset_in bitset of node at the origin of the searching region.
//* @param[out]  extend_pos  extended layer in positve direction.
//* @param[out]  extend_neg  extended layer in negative direction.
//*/
//void GridManagerInterface::ResetExtendLayerBasedOnDomainSize2D(
//    const DefSizet i_level, const DefSFBitset& sfbitset_in,
//    std::vector<DefLUint>* ptr_vec_extend_neg,
//    std::vector<DefLUint>* ptr_vec_extend_pos) {
//
//    // extended layer at the background refinement level
//    DefLUint two_power_i_level = TwoPowerN(static_cast<DefLUint>(i_level));
//    DefLUint index_xmin = k0IntOffest_[kXIndex] * two_power_i_level;
//    DefLUint index_xmax = k0MaxIndexOfBackgroundNode_[kXIndex]
//        * two_power_i_level;
//    DefLUint index_ymin = k0IntOffest_[kYIndex] * two_power_i_level;
//    DefLUint index_ymax = k0MaxIndexOfBackgroundNode_[kYIndex]
//        * two_power_i_level;
//
//    std::array<DefLUint, 2> indices;
//    SFBitsetComputeIndices2D(sfbitset_in, &indices);
//
//    // reset extended layer
//    if ((indices[kXIndex] - index_xmin) < ptr_vec_extend_neg->at(kXIndex)) {
//        ptr_vec_extend_neg->at(kXIndex) = (indices[kXIndex] - index_xmin);
//    }
//    if ((index_xmax - indices[kXIndex])
//        < ptr_vec_extend_pos->at(kXIndex)) {
//        ptr_vec_extend_pos->at(kXIndex) = (index_xmax - indices[kXIndex]);
//    }
//    if ((indices[kYIndex] - index_ymin) < ptr_vec_extend_neg->at(kYIndex)) {
//        ptr_vec_extend_neg->at(kYIndex) = (indices[kYIndex] - index_ymin);
//    }
//    if ((index_ymax - indices[kYIndex])
//        < ptr_vec_extend_pos->at(kYIndex)) {
//        ptr_vec_extend_pos->at(kYIndex) = (index_ymax - indices[kYIndex]);
//    }
//}
//#endif  // DEBUG_DISABLE_2D_FUNCTIONS
//#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
///**
//* @brief   function to find tracking node based on a given 3D geo point.
//* @param[in]  grid_space   grid spacing in x, y and z direction.
//* @param[in]  geo_coordi  coordinate of geometry point.
//* @param[out] space filling codes of nodes containing the input geometry point.
//*/
//void GridManagerInterface::FindTrackingNodesContainingAGeometryCoordinate3D(
//    const std::vector<DefReal>& grid_space,
//    const std::vector<DefReal>& geo_coordi,
//    std::vector<DefSFBitset>* ptr_vec_sfbitset) {
//    DefLUint x_index = static_cast<DefLUint>
//        (geo_coordi[kXIndex] / grid_space[kXIndex] + kEps);
//    DefLUint y_index = static_cast<DefLUint>
//        (geo_coordi[kYIndex] / grid_space[kYIndex] + kEps);
//    DefLUint z_index = static_cast<DefLUint>
//        (geo_coordi[kZIndex] / grid_space[kZIndex] + kEps);
//    ptr_vec_sfbitset->at(0) = SFBitsetEncoding3D(
//        { x_index, y_index, z_index });
//    ptr_vec_sfbitset->at(1) = FindXPos3D(ptr_vec_sfbitset->at(0));
//    ptr_vec_sfbitset->at(2) = FindYPos3D(ptr_vec_sfbitset->at(1));
//    ptr_vec_sfbitset->at(3) = FindYPos3D(ptr_vec_sfbitset->at(0));
//    ptr_vec_sfbitset->at(4) = FindZPos3D(ptr_vec_sfbitset->at(0));
//    ptr_vec_sfbitset->at(5) = FindZPos3D(ptr_vec_sfbitset->at(1));
//    ptr_vec_sfbitset->at(6) = FindZPos3D(ptr_vec_sfbitset->at(2));
//    ptr_vec_sfbitset->at(7) = FindZPos3D(ptr_vec_sfbitset->at(3));
//}
///**
//* @brief   function to reset extened layers if the searching region
//           if ouside the computational domain (3D).
//* @param[in]  i_level   level of refinement.
//* @param[in]   sfbitset_in   space filling code at the origin
//*                          of the searching region.
//* @param[out]  extend_neg  number of extended layer in negative direction.
//* @param[out]  extend_pos  number of extended layer in positve direction.
//*/
//void GridManagerInterface::ResetExtendLayerBasedOnDomainSize3D(
//    const DefSizet i_level, const DefSFBitset& sfbitset_in,
//    std::vector<DefLUint>* ptr_vec_extend_neg,
//    std::vector<DefLUint>* ptr_vec_extend_pos) {
//
//    // extended layer at the background refinement level
//    DefLUint two_power_i_level = TwoPowerN(static_cast<DefLUint>(i_level));
//    DefLUint index_xmin = k0IntOffest_[kXIndex] * two_power_i_level;
//    DefLUint index_xmax = k0MaxIndexOfBackgroundNode_[kXIndex]
//        * two_power_i_level;
//    DefLUint index_ymin = k0IntOffest_[kYIndex] * two_power_i_level;
//    DefLUint index_ymax = k0MaxIndexOfBackgroundNode_[kYIndex]
//        * two_power_i_level;
//    DefLUint index_zmin = k0IntOffest_[kZIndex] * two_power_i_level;
//    DefLUint index_zmax = k0MaxIndexOfBackgroundNode_[kZIndex]
//        * two_power_i_level;
//
//    std::array<DefLUint, 3> indices;
//    SFBitsetComputeIndices3D(sfbitset_in, &indices);
//
//    // reset extended layer
//    if ((indices[kXIndex] - index_xmin) < ptr_vec_extend_neg->at(kXIndex)) {
//        ptr_vec_extend_neg->at(kXIndex) = (indices[kXIndex] - index_xmin);
//    }
//    if ((index_xmax - indices[kXIndex])
//        < ptr_vec_extend_pos->at(kXIndex)) {
//        ptr_vec_extend_pos->at(kXIndex) = (index_xmax - indices[kXIndex]);
//    }
//    if ((indices[kYIndex] - index_ymin) < ptr_vec_extend_neg->at(kYIndex)) {
//        ptr_vec_extend_neg->at(kYIndex) = (indices[kYIndex] - index_ymin);
//    }
//    if ((index_ymax - indices[kYIndex])
//        < ptr_vec_extend_pos->at(kYIndex)) {
//        ptr_vec_extend_pos->at(kYIndex) = (index_ymax - indices[kYIndex]);
//    }
//    if ((indices[kZIndex] - index_zmin) < ptr_vec_extend_neg->at(kZIndex)) {
//        ptr_vec_extend_neg->at(kZIndex) = (indices[kZIndex] - index_zmin);
//    }
//    if ((index_zmax - indices[kZIndex])
//        < ptr_vec_extend_pos->at(kZIndex)) {
//        ptr_vec_extend_pos->at(kZIndex) = (index_zmax - indices[kZIndex]);
//    }
//}
//#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namsapce grid
}  // end namespace amrproject
}  // end namespace rootproject
