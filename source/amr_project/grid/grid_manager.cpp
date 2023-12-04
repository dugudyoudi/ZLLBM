//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
*/
#include <string>
#include <set>
#include <mpi.h>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
#include "io/debug_write.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to get background time step of current level.
*/
DefReal MultiTimeSteppingC2F::GetCurrentTimeStep(const DefAmrIndexUint i_level,
    const DefAmrIndexLUint time_step_background, const DefAmrIndexUint time_step_level) const {
    return time_step_background
        + static_cast<DefReal>(time_step_level) / static_cast<DefReal>(TwoPowerN(i_level));
}
/**
* @brief function to setup default grid related parameters.
* @param[in]  max_level  maximum refinement level.
*/
void GridManagerInterface::DefaultInitialization(const DefAmrIndexUint max_level) {
    k0MaxLevel_ = max_level;
}
/**
* @brief   function tp create a tracking grid instance for a given geometry.
* @param[in] i_geo  index of the geometry.
* @param[in] geo_info object containing the geometry information.
*/ 
void GridManagerInterface::CreateTrackingGridInstanceForAGeo(const DefAmrIndexUint i_geo,
    const GeometryInfoInterface& geo_info) {
    std::pair<ECriterionType, DefAmrIndexUint> key_tracking_grid = { ECriterionType::kGeometry, i_geo };
    if (vec_ptr_grid_info_.at(geo_info.i_level_)->map_ptr_tracking_grid_info_.find(key_tracking_grid)
        == vec_ptr_grid_info_.at(geo_info.i_level_)->map_ptr_tracking_grid_info_.end()) {
        vec_ptr_grid_info_.at(geo_info.i_level_)->map_ptr_tracking_grid_info_.insert(
        { key_tracking_grid, geo_info.ptr_tracking_grid_info_creator_->CreateTrackingGridInfo() });
    }
}
/**
* @brief   function to create instance based on the type of solver.
* @param[in] ptr_solver pointer to the solver.
* @param[in] grid_creator instance for creating grid infomation.
*/
void GridManagerInterface::CreateSameGridInstanceForAllLevel(const std::shared_ptr<SolverInterface>& ptr_solver,
    const GridInfoCreatorInterface& grid_creator) {
    ptr_solver->ptr_grid_manager_ = this;
    ptr_solver->SetNodeFlagForSolver();
    for (DefAmrIndexUint i_level = 0; i_level < k0MaxLevel_ + 1; ++i_level) {
        vec_ptr_grid_info_.emplace_back(grid_creator.CreateGridInfo());
        GridInfoInterface& grid_ref = *(vec_ptr_grid_info_).back();
        grid_ref.i_level_ = i_level;
        grid_ref.ptr_solver_ = ptr_solver;
        // set computational cost for each node 2^i_level
        grid_ref.computational_cost_ = static_cast<DefAmrUint>(TwoPowerN(i_level));

        // grid spacing
        grid_ref.grid_space_ = std::vector<DefReal>(k0GridDims_, 0.);
        std::vector<DefReal> domain_dx_ = GetDomainDxArrAsVec();
        for (DefAmrIndexUint idim = 0; idim < k0GridDims_; ++idim) {
            grid_ref.grid_space_.at(idim) =
                domain_dx_.at(idim) / static_cast<DefReal>(TwoPowerN(i_level));
        }

        // domain boundary related
        CalDomainBoundsAtGivenLevel(i_level,
            &grid_ref.k0VecBitsetDomainMin_, &grid_ref.k0VecBitsetDomainMax_);
        if (this->ptr_func_insert_domain_boundary_ == &GridManagerInterface::InsertCubicDomainBoundary) {
            grid_ref.domain_boundary_min_.resize(k0GridDims_);
            grid_ref.domain_boundary_max_.resize(k0GridDims_);
        }
    }
}
/**
 * @brief function to insert boundary node on a cubic computational domain.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[out] ptr_grid_info pointer to grid information including records of computational domain.
 */
void GridManagerInterface::InstantiateGridNode(const DefSFBitset& bitset_in,
    GridInfoInterface* const ptr_grid_info) {
    ptr_grid_info->map_grid_node_.insert({bitset_in, ptr_grid_info->GridNodeCreator()});
    (this->*ptr_func_insert_domain_boundary_)(bitset_in, ptr_grid_info);
}
/**
 * @brief function to insert boundary node on a cubic computational domain.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[out] ptr_grid_info pointer to grid information including records of computational domain.
 */
void GridManagerInterface::InsertCubicDomainBoundary(const DefSFBitset& bitset_in,
    GridInfoInterface* const ptr_grid_info) const {
    DefAmrIndexUint flag_node = CheckNodeOnDomainBoundary(bitset_in,
        ptr_grid_info->k0VecBitsetDomainMin_, ptr_grid_info->k0VecBitsetDomainMax_);
    // flag_node 1: x min, 2: x max, 4: y min, 8: y max, 16: z min, 32: z max
    if (flag_node > 0) {  // flag_node == 0 indicates not on domain boundary
        if ((flag_node & 1) == 1) {  // on x min boundary
        ptr_grid_info->domain_boundary_min_[kXIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 2) == 2) {  // on x max boundary
        ptr_grid_info->domain_boundary_max_[kXIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 4) == 4) {  // on y min boundary
        ptr_grid_info->domain_boundary_min_[kYIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 8) == 8) {  // on y max boundary
        ptr_grid_info->domain_boundary_max_[kYIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 16) == 16) {  // on z min boundary
            ptr_grid_info->domain_boundary_min_[kZIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 32) == 32) {  // on z max boundary
            ptr_grid_info->domain_boundary_max_[kZIndex].insert({bitset_in, 0});
        }
    }
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
    if (coordinate_max[kXIndex] > domain_size[kXIndex] + domain_offset[kXIndex]) {
        status |= 2;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        status |= 4;
    }
    if (coordinate_max[kYIndex] > domain_size[kYIndex] + domain_offset[kYIndex]) {
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
    if (coordinate_max[kXIndex] > domain_size[kXIndex] + domain_offset[kXIndex]) {
        status |= 2;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        status |= 4;
    }
    if (coordinate_max[kYIndex] > domain_size[kYIndex] + domain_offset[kYIndex]) {
        status |= 8;
    }
    if (coordinate_min[kZIndex] < domain_offset[kZIndex]) {
        status |= 16;
    }
    if (coordinate_max[kZIndex] > domain_size[kZIndex] + domain_offset[kZIndex]) {
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
void GridManagerInterface::SetNumberOfExtendLayerForGrid(const DefAmrIndexUint i_level,
    const GeometryInfoInterface& geo_info,
    std::vector<DefAmrIndexLUint>* const ptr_inner_layer_neg,
    std::vector<DefAmrIndexLUint>* const ptr_inner_layer_pos,
    std::vector<DefAmrIndexLUint>* const ptr_outer_layer_neg,
    std::vector<DefAmrIndexLUint>* const ptr_outer_layer_pos) {
    std::vector<DefAmrIndexLUint> layer_min(k0GridDims_, k0IntExtendMin_);
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
/**
 * @brief function to find overlapping layers between grid of adjacent refinement levels based on the outermost coarse layer.
 * @param[in]  layer_coarse_m1  nodes in the second outermost coarse layer.
 * @param[in]  sfbitset_exist   existent nodes.
 * @param[out] ptr_layer_coarse_0 pointer to nodes in the outermost coarse layer.
 * @param[out] ptr_layer_fine_0 pointer to nodes in the outermost fine layer.
 * @param[out] ptr_layer_fine_m1 pointer to nodes in the second outermost fine layer.
 * @param[out] ptr_layer_fine_m2 pointer to nodes in the third outermost fine layer.
 */
void GridManagerInterface::FindOverlappingLayersBasedOnOutermostCoarse(
    const DefMap<DefAmrUint>& layer_coarse_m1, const DefMap<DefAmrIndexUint>& sfbitset_exist,
    DefMap<DefAmrUint>* const ptr_layer_coarse_0, DefMap<DefAmrUint>* const ptr_layer_fine_0,
    DefMap<DefAmrUint>* const ptr_layer_fine_m1, DefMap<DefAmrUint>* const ptr_layer_fine_m2) {
#ifdef DEBUG_CHECK_GRID
    if (&layer_coarse_m1 == ptr_layer_coarse_0) {
        LogManager::LogError("input (layer_coarse_0)"
         " should not be the same as output (ptr_layer_coarse_m1) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID

    std::vector<DefSFBitset> corner_bitsets;
    for (const auto& iter_node : layer_coarse_m1) {
        FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
        for (auto& iter_conner : corner_bitsets) {
            IdentifyInterfaceForACell(iter_conner, layer_coarse_m1, sfbitset_exist,
             ptr_layer_fine_m2, ptr_layer_fine_m1, ptr_layer_fine_0);
        }
    }
    // coarse node on the layer_coarse_m1 layer, which overlap with
    // layer0 layer at one level higher
#ifdef DEBUG_CHECK_GRID
    if (ptr_layer_fine_0 == ptr_layer_coarse_0) {
        LogManager::LogError("input (ptr_layer_fine_0)"
         " should not be the same as output (ptr_layer_coarse_0) in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID
    OverlapLayerFromHighToLow(*ptr_layer_fine_m2, ptr_layer_coarse_0);
}
/**
 * @brief function to instantiates nodes on the overlap layers of refinement interface for all levels.
 * @param sfbitset_one_lower_level nodes indices encoded by space filling code at one lower level.
 */
void GridManagerInterface::InstantiateOverlapLayerOfRefinementInterface(
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level) {
    InterfaceLayerInfo* ptr_interface_info = nullptr;
    InterfaceLayerInfo* ptr_interface_info_lower = nullptr;
    DefAmrIndexUint layer_coarse0, layer_coarse_m1, layer0, layer_m1, layer_m2;
    DefAmrUint flag_temp, flag_refined = kNodeStatus0_;
    int maxlayer;
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        GridInfoInterface& grid_info_lower =
            *(vec_ptr_grid_info_.at(i_level - 1));
        DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
        DefMap<std::unique_ptr<GridNode>>& map_grid_lower = grid_info_lower.map_grid_node_;
#ifdef DEBUG_CHECK_GRID
        if (grid_info_lower.k0NumCoarse2FineLayer_ < 2) {
            LogManager::LogError("number of coarse to fine layers at level "
                + std::to_string(i_level - 1) + " is at least 1"
                + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        if (grid_info.k0NumFine2CoarseLayer_ < 3) {
            LogManager::LogError("number of fine to coarse layers at level "
                + std::to_string(i_level) + " is at least 3"
                + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
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
                    ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_m1),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse0),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
                // insert node instance
                maxlayer = static_cast<int>(ptr_interface_info->vec_inner_fine2coarse_.size());
                for (int ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_inner_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                            InstantiateGridNode(iter_layer_node.first, &grid_info);
                            map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        } else {
                            map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                        }
                    }
                }
                maxlayer = static_cast<int>(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_inner_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first)
                            == map_grid_lower.end()) {
                            InstantiateGridNode(iter_layer_node.first, &grid_info_lower);
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ =
                                flag_temp;
                        } else {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ |=
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
                    ->vec_outer_coarse2fine_.at(layer_coarse_m1),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(layer_coarse0),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m2));
                // insert node instance
                DefAmrIndexUint maxlayer = DefAmrIndexUint(ptr_interface_info->vec_outer_fine2coarse_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_outer_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                            InstantiateGridNode(iter_layer_node.first, &grid_info);
                            map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        } else {
                            map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                        }
                    }
                }
                maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_outer_coarse2fine_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | kNodeStatusCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | kNodeStatusCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_outer_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first)
                            == map_grid_lower.end()) {
                            InstantiateGridNode(iter_layer_node.first, &grid_info_lower);
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_temp;
                        } else {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                        }
                    }
                }
            }
        }
    }
}
/**
* @brief function to instantiate interface layers between grids of adjacent refinement level.
* @param[in] i_level level of refinement.
* @param[in] num_partition_outer_layer  number of outer layers for mpi communication.
* @param[in] flag_refinement flag indicates node is on the refinement interface.
* @param[in] code_min the minimum space filling code assigned to current rank.
* @param[in] code_max the maximum space filling code assigned to current rank.
* @param[in] sfbitset_aux class manages space filling curves.
* @param[in] map_sfbitset_one_lower_level  nodes indices encoded by space filling code at one lower level.
* @param[in] sfbitset_partition_interface_background nodes on mpi partitioned interface at background level. 
* @param[out] ptr_outer_layer_current_level pointer to nodes on outer layers for mpi communication at current level. 
* @param[out] ptr_outer_layer_lower_level pointer to nodes on outer layers for mpi communication at current level. 
*/

void GridManagerInterface::InstantiateOverlapLayerOfRefinementInterfaceMpi(
    const DefAmrIndexUint i_level, const DefAmrIndexUint num_partition_outer_layer,
    const DefAmrUint flag_refinement,
    const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
    const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefAmrIndexUint>& map_sfbitset_one_lower_level,
    const DefMap<DefAmrIndexUint>& sfbitset_partition_interface_background,
    DefMap<DefAmrUint>* const ptr_outer_layer_current_level,
    DefMap<DefAmrUint>* const ptr_outer_layer_lower_level) {
    InterfaceLayerInfo* ptr_interface_info = nullptr;
    InterfaceLayerInfo* ptr_interface_info_lower = nullptr;
    DefAmrUint flag_temp, flag_refined = 0;
    int maxlayer;
    if (i_level == 0) {
        LogManager::LogError("input refinement level should be greater than 0 in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    DefAmrIndexUint i_lower_level = i_level - 1;
    GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
    GridInfoInterface& grid_info_lower = *(vec_ptr_grid_info_.at(i_lower_level));
    DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
    DefMap<std::unique_ptr<GridNode>>& map_grid_lower = grid_info_lower.map_grid_node_;
    DefSFCodeToUint code_background;
#ifdef DEBUG_CHECK_GRID
    if (grid_info_lower.k0NumCoarse2FineLayer_ < 2) {
        LogManager::LogError("number of coarse to fine layers at level "
            + std::to_string(i_level - 1) + " is at least 1"
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (grid_info.k0NumFine2CoarseLayer_ < 3) {
        LogManager::LogError("number of fine to coarse layers at level "
            + std::to_string(i_level) + " is at least 3"
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID
    // find interface node at current level
    DefAmrIndexUint layer_coarse0 = grid_info_lower.k0NumCoarse2FineLayer_ - 1,
        layer_coarse_m1 = layer_coarse0 - 1,
        layer0 = grid_info.k0NumFine2CoarseLayer_ - 1,
        layer_m1 = layer0 - 1,
        layer_m2 = layer_m1 - 1;
    std::vector<DefSFBitset> sfbitset_neighbors;
    for (auto& iter_interface : grid_info_lower.map_ptr_interface_layer_info_) {
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
                ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_m1),
                map_sfbitset_one_lower_level,
                &ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse0),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
            // insert node instance
            maxlayer = static_cast<int>(ptr_interface_info->vec_inner_fine2coarse_.size());
            for (int ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_refined | kNodeStatusFine2Coarse0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_refined | kNodeStatusFine2CoarseM1_;
                } else {
                    flag_temp = flag_refined;
                }
                for (const auto& iter_layer_node : ptr_interface_info
                    ->vec_inner_fine2coarse_.at(ilayer)) {
                    if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                        InstantiateGridNode(iter_layer_node.first, &grid_info);
                        map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                    } else {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_level, iter_layer_node.first).to_ullong();
                    if (ilayer == maxlayer - 1 && (code_background < code_min || code_background > code_max)) {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                        if (ptr_outer_layer_current_level->find(iter_layer_node.first)
                            == ptr_outer_layer_current_level->end()) {
                            ptr_outer_layer_current_level->insert({iter_layer_node.first, flag_refinement});
                        } else {
                            ptr_outer_layer_current_level->at(iter_layer_node.first) = flag_refinement;
                        }
                    }
                }
            }
            maxlayer = static_cast<int>(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_refined | kNodeStatusCoarse2Fine0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_refined | kNodeStatusCoarse2FineM1_;
                } else {
                    flag_temp = flag_refined;
                }
                for (const auto& iter_layer_node : ptr_interface_info_lower
                    ->vec_inner_coarse2fine_.at(ilayer)) {
                    if (map_grid_lower.find(iter_layer_node.first)
                        == map_grid_lower.end()) {
                        InstantiateGridNode(iter_layer_node.first, &grid_info_lower);
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ =  flag_temp;
                    } else {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_lower_level, iter_layer_node.first).to_ullong();
                    if (ilayer == maxlayer - 1 && (code_background < code_min || code_background > code_max)) {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                        if (ptr_outer_layer_lower_level->find(iter_layer_node.first)
                            == ptr_outer_layer_lower_level->end()) {
                            ptr_outer_layer_lower_level->insert({iter_layer_node.first, flag_refinement});
                        } else {
                            ptr_outer_layer_lower_level->at(iter_layer_node.first) = flag_refinement;
                        }
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
                ->vec_outer_coarse2fine_.at(layer_coarse_m1),
                map_sfbitset_one_lower_level,
                &ptr_interface_info_lower
                ->vec_outer_coarse2fine_.at(layer_coarse0),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer0),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m1),
                &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m2));
            // insert node instance
            DefAmrIndexUint maxlayer = DefAmrIndexUint(ptr_interface_info->vec_outer_fine2coarse_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_refined | kNodeStatusFine2Coarse0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_refined | kNodeStatusFine2CoarseM1_;
                } else {
                    flag_temp = flag_refined;
                }
                for (const auto& iter_layer_node : ptr_interface_info
                    ->vec_outer_fine2coarse_.at(ilayer)) {
                    if (map_grid.find(iter_layer_node.first)
                        == map_grid.end()) {
                        InstantiateGridNode(iter_layer_node.first, &grid_info);
                        map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                    } else {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_level, iter_layer_node.first).to_ullong();
                    if (ilayer == maxlayer - 1 && (code_background < code_min || code_background > code_max)) {
                        map_grid.at(iter_layer_node.first)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                        if (ptr_outer_layer_current_level->find(iter_layer_node.first)
                            == ptr_outer_layer_current_level->end()) {
                            ptr_outer_layer_current_level->insert({iter_layer_node.first, flag_refinement});
                        } else {
                            ptr_outer_layer_current_level->at(iter_layer_node.first) = flag_refinement;
                        }
                    }
                }
            }
            maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_outer_coarse2fine_.size());
            for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                if (ilayer == maxlayer - 1) {
                    flag_temp = flag_refined | kNodeStatusCoarse2Fine0_;
                } else if (ilayer == maxlayer - 2) {
                    flag_temp = flag_refined | kNodeStatusCoarse2FineM1_;
                } else {
                    flag_temp = flag_refined;
                }
                for (const auto& iter_layer_node : ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(ilayer)) {
                    if (map_grid_lower.find(iter_layer_node.first)
                        == map_grid_lower.end()) {
                        InstantiateGridNode(iter_layer_node.first, &grid_info_lower);
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_temp;
                    } else {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                    }
                    // check if node on outer mpi communication layer
                    code_background = sfbitset_aux.SFBitsetToNLowerLevelVir(
                        i_lower_level, iter_layer_node.first).to_ullong();
                    if (ilayer == maxlayer - 1 && (code_background < code_min || code_background > code_max)) {
                        map_grid_lower.at(iter_layer_node.first)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                        if (ptr_outer_layer_lower_level->find(iter_layer_node.first)
                            == ptr_outer_layer_lower_level->end()) {
                            ptr_outer_layer_lower_level->insert({iter_layer_node.first, flag_refinement});
                        } else {
                            ptr_outer_layer_lower_level->at(iter_layer_node.first) = flag_refinement;
                        }
                    }
                }
            }
        }
    }
}
/**
 * @brief function to instantiate grid noes for all refinement levels (only works for serial).
 * @param[in] sfbitset_min  minimum space filling code of the current rank.
 * @param[in] sfbitset_max  maximum space filling code of the current rank.
 * @param[in] sfbitset_one_lower_level space filling codes at one lower refinement level.
 */
void GridManagerInterface::InstantiateGridNodeAllLevel(const DefSFBitset sfbitset_min,
    const DefSFBitset sfbitset_max, const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level) {
    DefMap<DefAmrIndexUint> background_occupied;
    InstantiateOverlapLayerOfRefinementInterface(sfbitset_one_lower_level);
    DefSFCodeToUint code_min = sfbitset_min.to_ullong(), code_max = sfbitset_max.to_ullong();
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
        DefAmrIndexUint i_level_lower = i_level - 1;
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        // initialize grid information
        grid_info.InitialGridInfo();
        DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
        std::vector<DefSFBitset> bitset_cell_lower, bitset_all;
        DefSFBitset bitset_background;
        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            if (CheckCoincideBackground(i_level_lower, iter_low.first, &bitset_background)) {
                if (background_occupied.find(bitset_background)
                 == background_occupied.end()) {
                    background_occupied.insert({ bitset_background, kFlagSize0_ });
                }
            }
            if (NodesBelongToOneCell(iter_low.first,
                sfbitset_one_lower_level.at(i_level), &bitset_cell_lower)) {
                FindAllNodesInACellAtLowerLevel(bitset_cell_lower, &bitset_all);
                for (const auto& iter_node : bitset_all) {
                    if (map_grid.find(iter_node) == map_grid.end()) {
                        InstantiateGridNode(iter_node, &grid_info);
                    }
                }
            }
        }
    }
    // find overlapping node for refinement levels of 0 and 1
    for (const auto & iter_interfaces : vec_ptr_grid_info_.at(0)->map_ptr_interface_layer_info_) {
        for (const auto & iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
            for (const auto & iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
        }
    }
    vec_ptr_grid_info_.at(0)->InitialGridInfo();
    InstantiateBackgroundGrid(code_min, code_max, background_occupied);
}
/**
 * @brief function to instantiate grid noes for all refinement levels and find layers for mpi communication.
 * @param[in] i_rank current rank id
 * @param[in] num_partition_inner_layer number of inner layers for mpi communication
 * @param[in] num_partition_outer_layer number of outer layers for mpi communication
 * @param[in] vec_sfbitset_min  minimum space filling code of all ranks.
 * @param[in] vec_sfbitset_max  maximum space filling code of all ranks.
 * @param[in] sfbitset_aux  class manage space filling curves.
 * @param[in] sfbitset_one_lower_level space filling codes at one lower refinement level.
 * @param[in] sfbitset_ghost_one_lower_level space filling codes of mpi communication node near coarse to fine refinement interface.
 * @param[in] sfbitset_partition_interface_0  nodes on the background partitioned interface of background level on current rank.
 * @param[out] ptr_mpi_inner_layer pointer to nodes on the inner layer for mpi communication (sending).
 * @param[out] ptr_mpi_outer_layer pointer to nodes on the outer layer for mpi communication (sending). 
 */
void GridManagerInterface::InstantiateGridNodeAllLevelMpi(const int i_rank,
    const DefAmrIndexUint num_partition_inner_layer, const DefAmrIndexUint num_partition_outer_layer,
    const std::vector<DefSFBitset>& vec_sfbitset_min,
    const std::vector<DefSFBitset>& vec_sfbitset_max, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level,
    const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_ghost_one_lower_level,
    const DefMap<DefAmrIndexUint>& sfbitset_partition_interface_0,
    std::vector<DefMap<std::set<int>>>* const ptr_mpi_inner_layer,
    std::vector<DefMap<DefAmrUint>>* const ptr_mpi_outer_layer) {
    DefAmrUint flag_0 = 0, flag_near_refinement = 1, flag_refinement = 2, flag_background_tmp = 3;
    DefMap<DefAmrIndexUint> background_occupied;  // background nodes coincide with those at higher levels
    ptr_mpi_inner_layer->resize(k0MaxLevel_ + 1);
    ptr_mpi_outer_layer->resize(k0MaxLevel_ + 1);
    DefSFCodeToUint code_min = vec_sfbitset_min.at(i_rank).to_ullong(),
     code_max = vec_sfbitset_max.at(i_rank).to_ullong();
    DefSFCodeToUint code_tmp, code_neighbor_tmp;
    bool bool_partition_outside, bool_partition_inside;
    std::vector<DefSFBitset> domain_min_m1_n_level(k0GridDims_), domain_max_p1_n_level(k0GridDims_);
    std::vector<DefSFBitset> vec_in_region;
    std::vector<DefSFCodeToUint>::iterator iter_index;
    int node_rank;
    std::vector<DefSFCodeToUint> ull_max(vec_sfbitset_max.size());
    for (DefSizet i = 0; i < vec_sfbitset_max.size(); ++i) {
        ull_max.at(i) = vec_sfbitset_max.at(i).to_ullong();
    }
    std::vector<DefAmrIndexLUint> indices_min = GetMinIndexOfBackgroundNodeArrAsVec(),
        indices_max = GetMaxIndexOfBackgroundNodeArrAsVec();
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
        DefAmrIndexUint i_level_lower = i_level - 1;
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        // initialize grid information
        grid_info.InitialGridInfo();
        DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
        std::vector<DefSFBitset> bitset_cell_lower, bitset_all, bitset_neighbors;
        DefSFBitset bitset_background;
        std::vector<DefSFBitset> domain_min_n_level(k0GridDims_), domain_max_n_level(k0GridDims_);
        sfbitset_aux.GetMinAtGivenLevel(i_level, indices_min, &domain_min_n_level);
        sfbitset_aux.GetMaxAtGivenLevel(i_level, indices_max, &domain_max_n_level);

        if (grid_info.ptr_solver_ == nullptr) {
            LogManager::LogError("solver has not been assigned to grid info"
            " of refinement level " + std::to_string(i_level)
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }

        // instantiate nodes on outer mpi communication layer and refinement interface
        // The spatial step of node on the coarse to fine interface at i_level
        // is half of those in sfbitset_one_lower_level used for instantiation.
        // Thus nodes on coarse to fine interfaces are sent from rank0 individually
        for (const auto& iter_ghost : sfbitset_ghost_one_lower_level.at(i_level)) {
            InstantiateGridNode(iter_ghost.first, vec_ptr_grid_info_.at(i_level_lower).get());
            vec_ptr_grid_info_.at(i_level_lower)->map_grid_node_.at(iter_ghost.first)->flag_status_
                |= kNodeStatusMpiPartitionOutside_;
            ptr_mpi_outer_layer->at(i_level_lower).insert({iter_ghost.first, flag_near_refinement});
        }
        // instantiate refinement interfaces
        InstantiateOverlapLayerOfRefinementInterfaceMpi(i_level, num_partition_outer_layer, flag_refinement,
            code_min, code_max, sfbitset_aux, sfbitset_one_lower_level.at(i_level),
            sfbitset_partition_interface_0,
            &ptr_mpi_outer_layer->at(i_level), &ptr_mpi_outer_layer->at(i_level_lower));
        // instantiate nodes stored in sfbitset_one_lower_level
        DefMap<DefAmrUint> partition_interface_level;   // nodes on partition interface at current level
        for (const auto& iter_low : sfbitset_one_lower_level.at(i_level)) {
            if (CheckCoincideBackground(i_level_lower, iter_low.first, &bitset_background)) {
                if (background_occupied.find(bitset_background)
                 == background_occupied.end()) {
                    background_occupied.insert({ bitset_background, kFlagSize0_ });
                }
            }
            if (NodesBelongToOneCell(iter_low.first,
                sfbitset_one_lower_level.at(i_level), &bitset_cell_lower)) {
                bool_partition_outside = false;
                bool_partition_inside = false;
                FindAllNodesInACellAtLowerLevel(bitset_cell_lower, &bitset_all);
                for (const auto& iter_node : bitset_all) {
                    if (map_grid.find(iter_node) == map_grid.end()) {
                        InstantiateGridNode(iter_node, &grid_info);
                    }
                    code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                    if (code_tmp < code_min || code_tmp > code_max) {
                        map_grid.at(iter_node)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                        if (ptr_mpi_outer_layer->at(i_level).find(iter_node)
                         == ptr_mpi_outer_layer->at(i_level).end()) {
                            ptr_mpi_outer_layer->at(i_level).insert({iter_node, flag_0});
                        }
                        bool_partition_outside = true;
                    } else {
                        bool_partition_inside = true;
                    }
                }
                // one or more nodes in the cell should be on the partition interface.
                // check if one node is on the partition interface by finding if it is inside
                // and at least one of its neighbors is outside the partition block
                if (bool_partition_outside && bool_partition_inside) {
                    for (const auto& iter_node : bitset_all) {
                        code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                        if (code_tmp >= code_min && code_tmp <= code_max
                            && partition_interface_level.find(iter_node) == partition_interface_level.end()) {
                            sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_node,
                                domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                            for (const auto& iter_neighbor : bitset_neighbors) {
                                code_neighbor_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(
                                    i_level, iter_neighbor).to_ullong();
                                if (code_neighbor_tmp < code_min || code_neighbor_tmp > code_max) {
                                    partition_interface_level.insert({iter_node, 0});
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (const auto& iter_interface : partition_interface_level) {
            // extend mpi interface nodes for a given number of inner communication layers
            sfbitset_aux.FindNodesInReginOfGivenLength(iter_interface.first, num_partition_inner_layer,
                domain_min_m1_n_level, domain_max_p1_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                if ((code_tmp < code_min || code_tmp > code_max)
                    && ptr_mpi_inner_layer->at(i_level).find(iter_node) == ptr_mpi_inner_layer->at(i_level).end()
                    && map_grid.find(iter_node) != map_grid.end()) {
                    ptr_mpi_inner_layer->at(i_level).insert({iter_node, {}});
                }
            }
        }
        DefMap<DefAmrIndexUint> map_outmost_tmp;
        bool bool_interface_upper_extra = ((num_partition_outer_layer%2) == 0);
        for (const auto& iter_interface : ptr_mpi_outer_layer->at(i_level)) {
            // set ranks that should be sent by mpi inner layer nodes
            sfbitset_aux.FindNodesInReginOfGivenLength(iter_interface.first, num_partition_outer_layer,
                domain_min_m1_n_level, domain_max_p1_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = sfbitset_aux.SFBitsetToNLowerLevelVir(i_level, iter_node).to_ullong();
                if ((code_tmp >= code_min && code_tmp <= code_max)
                    && ptr_mpi_inner_layer->at(i_level).find(iter_node) != ptr_mpi_inner_layer->at(i_level).end()) {
                    // set indices for nodes on inner communication layers
                    iter_index = std::lower_bound(ull_max.begin(), ull_max.end(),
                        iter_interface.first.to_ullong());
                    node_rank = static_cast<int>(iter_index - ull_max.begin());
                    if (node_rank != i_rank
                        && ptr_mpi_inner_layer->at(i_level).at(iter_node).find(node_rank)
                        == ptr_mpi_inner_layer->at(i_level).at(iter_node).end()) {
                        ptr_mpi_inner_layer->at(i_level).at(iter_node).insert(node_rank);
                    }
                }
            }
            if ((bool_interface_upper_extra && code_tmp > code_max)
                && iter_interface.second == flag_0) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                    domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                for (const auto& iter_neighbor : bitset_neighbors) {
                    if (map_outmost_tmp.find(iter_neighbor) == map_outmost_tmp.end()
                        && map_grid.find(iter_neighbor) == map_grid.end()) {
                        InstantiateGridNode(iter_neighbor, &grid_info);
                        map_grid.at(iter_neighbor)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                        map_outmost_tmp.insert({iter_neighbor, 0});
                    }
                }
            } else if ((!bool_interface_upper_extra && code_tmp < code_min)
                && iter_interface.second == flag_0) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                    domain_min_n_level, domain_max_n_level, &bitset_neighbors);
                for (const auto& iter_neighbor : bitset_neighbors) {
                    if (map_outmost_tmp.find(iter_neighbor) == map_outmost_tmp.end()
                        && map_grid.find(iter_neighbor) == map_grid.end()) {
                            InstantiateGridNode(iter_neighbor, &grid_info);
                            map_grid.at(iter_neighbor)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                            map_outmost_tmp.insert({iter_neighbor, 0});
                    }
                }
            }
        }
    }
    // find overlapping node for refinement levels of 0 and 1
    for (const auto& iter_interfaces : vec_ptr_grid_info_.at(0)->map_ptr_interface_layer_info_) {
        for (const auto& iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
            for (const auto& iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
        }
    }
    vec_ptr_grid_info_.at(0)->InitialGridInfo();

    InstantiateBackgroundGrid(code_min, code_max, background_occupied);

    // find background nodes on mpi communication interface
    DefMap<std::unique_ptr<GridNode>>& map_grid_level_0 = vec_ptr_grid_info_.at(0)->map_grid_node_;
    sfbitset_aux.GetMinM1AtGivenLevel(0, indices_min, &domain_min_m1_n_level);
    sfbitset_aux.GetMaxP1AtGivenLevel(0, indices_max, &domain_max_p1_n_level);
    DefMap<DefAmrIndexUint> map_outmost_current, map_outmost_pre;
    std::vector<DefSFBitset> vec_neighbors;
    for (const auto& iter_interface : sfbitset_partition_interface_0) {
        if (background_occupied.find(iter_interface.first) == background_occupied.end()) {
            map_outmost_current.insert({iter_interface.first, 0});
        }
    }
    // find nodes on inner and outer communication layers
    for (const auto& iter_interface : map_outmost_current) {
        if (map_grid_level_0.find(iter_interface.first) != map_grid_level_0.end()) {
            sfbitset_aux.FindNodesInReginOfGivenLength(iter_interface.first, num_partition_inner_layer,
                domain_min_m1_n_level, domain_max_p1_n_level, &vec_in_region);
            for (const auto& iter_node : vec_in_region) {
                code_tmp = iter_node.to_ullong();
                if ((code_tmp >= code_min && code_tmp <= code_max)
                    && ptr_mpi_inner_layer->at(0).find(iter_node) == ptr_mpi_inner_layer->at(0).end()
                    && background_occupied.find(iter_node) == background_occupied.end()) {
                    ptr_mpi_inner_layer->at(0).insert({iter_node, {}});
                }
            }
        }
    }

    std::vector<DefSFBitset> domain_min_n_level(k0GridDims_), domain_max_n_level(k0GridDims_);
        sfbitset_aux.GetMinAtGivenLevel(0, indices_min, &domain_min_n_level);
        sfbitset_aux.GetMaxAtGivenLevel(0, indices_max, &domain_max_n_level);
    // find nodes on outer communication layers layer by layer
    for (auto i_layer = 0; i_layer < num_partition_outer_layer; ++i_layer) {
        for (const auto& iter_interface : map_outmost_current) {
            if (map_grid_level_0.find(iter_interface.first) != map_grid_level_0.end()) {
                sfbitset_aux.SFBitsetFindAllBondedNeighborsVir(iter_interface.first,
                domain_min_n_level, domain_max_n_level, &vec_neighbors);
                for (const auto& iter_node : vec_neighbors) {
                    code_tmp = iter_node.to_ullong();
                    if ((code_tmp < code_min || code_tmp > code_max)
                        && background_occupied.find(iter_node) == background_occupied.end()) {
                        if (ptr_mpi_outer_layer->at(0).find(iter_node) == ptr_mpi_outer_layer->at(0).end()) {
                            ptr_mpi_outer_layer->at(0).insert({iter_node, flag_background_tmp});
                            if (map_grid_level_0.find(iter_node) == map_grid_level_0.end()) {
                                InstantiateGridNode(iter_node, vec_ptr_grid_info_.at(0).get());
                                map_grid_level_0.at(iter_node)->flag_status_ |= kNodeStatusMpiPartitionOutside_;
                            }
                            map_outmost_pre.insert({iter_node, 0});
                        } else if (ptr_mpi_outer_layer->at(0).at(iter_node) != flag_refinement) {
                            ptr_mpi_outer_layer->at(0).at(iter_node) = flag_background_tmp;
                            map_outmost_pre.insert({iter_node, 0});
                        }
                    }
                }
            }
        }
        map_outmost_current = map_outmost_pre;
        map_outmost_pre.clear();
    }

    // set ranks that should be sent by mpi inner layer nodes (level 0)
    for (const auto& iter_interface : ptr_mpi_outer_layer->at(0)) {
        sfbitset_aux.FindNodesInReginOfGivenLength(iter_interface.first, num_partition_outer_layer,
            domain_min_m1_n_level, domain_max_p1_n_level, &vec_in_region);
        for (const auto& iter_node : vec_in_region) {
            code_tmp = iter_node.to_ullong();
            if ((code_tmp >= code_min && code_tmp <= code_max)
                && ptr_mpi_inner_layer->at(0).find(iter_node) != ptr_mpi_inner_layer->at(0).end()) {
                iter_index = std::lower_bound(ull_max.begin(), ull_max.end(), iter_interface.first.to_ullong());
                node_rank = static_cast<int>(iter_index - ull_max.begin());
                if (node_rank != i_rank
                    && ptr_mpi_inner_layer->at(0).at(iter_node).find(node_rank)
                    != ptr_mpi_inner_layer->at(0).at(iter_node).end()) {
                    ptr_mpi_inner_layer->at(0).at(iter_node).insert(node_rank);
                }
            }
        }
    }
}
/**
* @brief function to set k0TimeSteppingOrder_ as multi-stepping scheme
*        from the background level to the finest (ratio of 2).
*/
MultiTimeSteppingC2F::MultiTimeSteppingC2F(const DefAmrIndexUint max_level) {
    if (max_level == 0) {
        k0TimeSteppingOrder_ = {0};
    } else if (max_level == 1) {
        k0TimeSteppingOrder_ = {0, 1, 1};
    } else {
        std::vector<DefSizet> accumulate_t(max_level + 1, 0);
        DefAmrIndexUint i_level = 1;
        accumulate_t[0] = TwoPowerN(max_level);
        k0TimeSteppingOrder_ = {0};
        while (i_level > 0) {
            if (accumulate_t[i_level- 1] == accumulate_t[i_level]) {
                --i_level;
            } else {
                k0TimeSteppingOrder_.push_back(i_level);
                accumulate_t[i_level] += TwoPowerN(max_level - i_level);
                if (i_level + 1 < max_level) {
                    ++i_level;
                } else {
                    ++i_level;
                    k0TimeSteppingOrder_.push_back(i_level);
                    k0TimeSteppingOrder_.push_back(i_level);
                    accumulate_t[i_level] += TwoPowerN(max_level - i_level + 1);
                }
            }
        }
    }
}
/**
* @brief function to setup default grid related parameters.
* @param[in]  max_level  maximum refinement level.
*/
void GridManagerInterface::InstantiateTimeSteppingScheme() {
    switch (k0TimeSteppingType_) {
    case ETimeSteppingScheme::kMultiSteppingC2F:
        ptr_time_stepping_scheme_ = std::make_unique<MultiTimeSteppingC2F>(k0MaxLevel_);
        break;
    default:
        LogManager::LogError("undefined time stepping type in"
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        break;
    }
}
/**
* @brief      function to manage grid information during one time step at the background refinement level.
* @param[in]  time_step_background   current background time step.
* @param[in]  sfbitset_aux   class manages space filling curves.
*/
void GridManagerInterface::TimeMarching(const DefAmrIndexLUint time_step_background) {
    // record number of time step at i_level
    std::vector<DefAmrIndexUint> time_step_level(k0MaxLevel_ + 1, 0);
    DefAmrIndexUint i_level;
    DefReal time_step_current;
    for (auto iter_level = ptr_time_stepping_scheme_->k0TimeSteppingOrder_.begin();
        iter_level != ptr_time_stepping_scheme_->k0TimeSteppingOrder_.end(); ++iter_level) {
        i_level = *iter_level;
        ++time_step_level[i_level];
        time_step_current = ptr_time_stepping_scheme_->GetCurrentTimeStep(
            i_level, time_step_background, time_step_level[i_level]);
        GridInfoInterface& grid_ref = *vec_ptr_grid_info_.at(*iter_level);
        grid_ref.ptr_solver_->RunSolver(time_step_current, *this->GetSFBitsetAuxPtr(), &grid_ref);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
