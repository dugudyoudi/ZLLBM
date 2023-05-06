//  Copyright (c) 2022, Zhengliang Liu
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
    for (DefSizet i_level = 0; i_level < k0MaxLevel_ + 1; ++i_level) {
        vec_ptr_grid_info_.emplace_back(ptr_grid_creator->CreateGridInfo());
        GridInfoInterface& grid_ref = *(vec_ptr_grid_info_).back();
        grid_ref.i_level_ = i_level;
        grid_ref.grid_space_ = std::vector<DefReal>(k0GridDims_, 0.);
        // set computational cost for each node 2^i_level
        grid_ref.computational_cost_ = static_cast<DefUint>(TwoPowerN(i_level));
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
void GridManagerInterface::SetNumberOfExtendLayerForGrid(const DefSizet i_level,
    const GeometryInfoInterface& geo_info,
    std::vector<DefLUint>* const ptr_inner_layer_neg,
    std::vector<DefLUint>* const ptr_inner_layer_pos,
    std::vector<DefLUint>* const ptr_outer_layer_neg,
    std::vector<DefLUint>* const ptr_outer_layer_pos) {
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
}  // end namespace amrproject
}  // end namespace rootproject
