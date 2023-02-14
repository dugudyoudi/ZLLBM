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
///**
//* @brief   function to initialize grid noes.
//* @param[in] sfbitset_in   space filling code of the grid node.
//* @note
//*/
//void GridInfoInterface::InitialGridNode(const DefSFBitset& sfbitset_in) {
//    int rank_id = 0;
//#ifdef ENABLE_MPI
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
//#endif  // ENABLE_MPI
//    if (rank_id == 0) {
//        LogWarning("Function for initializing grid node "
//            "(GridInfoInterface::InitialGridNode) does not given");
//    }
//}
/**
* @brief   function to create instance based on the type of solver.
*/
void GridManagerInterface::CreateSameGridInstanceForAllLevel(
    GridInfoCreatorInterface* ptr_grid_creator) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
    for (DefSizet i_level = 0; i_level < k0MaxLevel_ + 1; ++i_level) {
        vec_ptr_grid_info_.emplace_back(ptr_grid_creator->CreateGridInfo());
        GridInfoInterface& grid_ref = *(vec_ptr_grid_info_).back();
        grid_ref.i_level_ = i_level;
        grid_ref.grid_space_ = std::vector<DefReal>(k0GridDims_, 0.);
        // set computational cost for each node 2^i_level
        grid_ref.computational_cost_ =
            static_cast<DefReal>(TwoPowerN(i_level));
        std::vector<DefReal> domain_dx_ = GetDomainDxArrAsVec();
        for (DefUint idim = 0; idim < k0GridDims_; ++idim) {
            grid_ref.grid_space_.at(idim) =
                domain_dx_.at(idim) /
                static_cast<DefReal>(TwoPowerN(i_level));
        }
        grid_ref.k0GridNodeInstance_.flag_status_ = kFlagExist_;
        grid_ref.set_number_of_vec_elements();
    }
}
/**
* @brief function to setup default grid related parameters.
*/
void GridManagerInterface::DefaultInitialization(const DefSizet max_level) {
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
* @brief   function to set extened layer of grid based on geomtery information.
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
void GridManagerInterface::SetNumberOfExtendLayerForGrid(const DefSizet i_level,
    const GeometryInfoInterface& geo_info,
    std::vector<DefLUint>* const ptr_inner_layer_neg,
    std::vector<DefLUint>* const ptr_inner_layer_pos,
    std::vector<DefLUint>* const ptr_outer_layer_neg,
    std::vector<DefLUint>* const ptr_outer_layer_pos){
    std::vector<DefLUint> layer_min(k0GridDims_, k0IntExtendMin_);
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
        && (geo_info.k0XIntExtendNegative_.at(i_level) > k0IntExtendMin_)){
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
* @brief   function to check if the geometry is out of domain (3D).
* @param[in]  sfbitset_lower_level  sfbitset at lower level.
* @param[in]  map_sfbitset_at_lower_level  existing nodes at lower level.
* @param[in]  ptr_grid_info pointer to grids at the current level.
*/
//void GridManagerInterface::AddNodesInstanceBasedOnLowerLevel(
//    const DefSFBitset& sfbitset_lower_level,
//    const DefMap<DefUint>& map_sfbitset_at_lower_level,
//    std::shared_ptr<GridInfoInterface> ptr_grid_info) {
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
//            | sfbitset_aux_.k0SfBitsetCurrentLevelBits_;
//        (*bitset_find_all_neighbours_ptr)(
//            sfbitset_center, &vec_neighbour_sfbitset);
//        for (const auto& iter : vec_neighbour_sfbitset) {
//            if (ptr_map_nodes->find(iter) == ptr_map_nodes->end()) {
//                ptr_map_nodes->insert({ iter, node_temp });
//                ptr_grid_info->InitialGridNode(iter);
//            }
//        }
//    }
//}
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
//    std::vector<DefLUint>* const ptr_vec_extend_neg,
//    std::vector<DefLUint>* const ptr_vec_extend_pos) {
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
//    std::vector<DefLUint>* const ptr_vec_extend_neg,
//    std::vector<DefLUint>* const ptr_vec_extend_pos) {
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
}  // end namespace amrproject
}  // end namespace rootproject
