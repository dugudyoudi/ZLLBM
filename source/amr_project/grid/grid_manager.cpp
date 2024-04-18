//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
    for (DefAmrIndexUint i_level = 0; i_level < k0MaxLevel_ + 1; ++i_level) {
        vec_ptr_grid_info_.emplace_back(grid_creator.CreateGridInfo());
        GridInfoInterface& grid_ref = *(vec_ptr_grid_info_).back();
        grid_ref.i_level_ = i_level;
        grid_ref.ptr_solver_ = ptr_solver;
        if (k0GridDims_ == 2) {
            grid_ref.ptr_sfbitset_aux_ = dynamic_cast<GridManager2D*>(this);
        } else {
            grid_ref.ptr_sfbitset_aux_ = dynamic_cast<GridManager3D*>(this);
        }

        // set computational cost for each node 2^i_level
        grid_ref.computational_cost_ = static_cast<DefAmrUint>(TwoPowerN(i_level));
        grid_ref.InitialNotComputeNodeFlag();

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
 * @brief function to create instance of a grid node.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[out] ptr_grid_info pointer to grid information including records of computational domain.
 * @return if true, node is instantiated
 */
bool GridManagerInterface::InstantiateGridNode(const DefSFBitset& bitset_in,
    GridInfoInterface* const ptr_grid_info) {
    int flag_node = ptr_grid_info->CheckIfNodeOutsideCubicDomain(k0GridDims_, bitset_in, *GetSFBitsetAuxPtr());
    if (flag_node >= 0) {
        ptr_grid_info->map_grid_node_.insert({bitset_in, ptr_grid_info->GridNodeCreator()});
        return true;
    }
    return false;
}
/**
 * @brief function to insert boundary node on a cubic computational domain.
 * @param[in] flag_node flag indicates node on which boundary of a cubic domain.
 * @param[in] bitset_in spacing filling code of a grid node.
 * @param[out] ptr_grid_info pointer to grid information including records of computational domain.
 */
void GridManagerInterface::InsertCubicDomainBoundary(const int flag_node,
    const DefSFBitset& bitset_in, GridInfoInterface* const ptr_grid_info) const {
    // flag_node 1: x min, 8: x max, 2: y min, 16: y max, 4: z min, 32: z max
    if (flag_node > 0) {  // flag_node == 0 indicates not on domain boundary
        if ((flag_node & 1) == 1) {  // on x min boundary
        ptr_grid_info->domain_boundary_min_[kXIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 8) == 8) {  // on x max boundary
        ptr_grid_info->domain_boundary_max_[kXIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 2) == 2) {  // on y min boundary
        ptr_grid_info->domain_boundary_min_[kYIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 16) == 16) {  // on y max boundary
        ptr_grid_info->domain_boundary_max_[kYIndex].insert({bitset_in, 0});
        }
        if ((flag_node & 4) == 4) {  // on z min boundary
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
* @return  0: coordinates does not exceed the computational domain; 1, 2:
*          coordinates are less than the lower bound; 8, 16: coordinates
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
        status |= 8;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        status |= 2;
    }
    if (coordinate_max[kYIndex] > domain_size[kYIndex] + domain_offset[kYIndex]) {
        status |= 16;
    }
    return status;
}
/**
* @brief   function to check if the geometry is out of domain.
* @param[in]  coordinate_min  minimum coordinates.
* @param[in]  coordinate_min  maximum coordinates.
* @param[in]  domain_offset  offset of the background grid.
* @param[in]  domain_size  size of the computational domain.
* @return  0: coordinates does not exceed the computational domain; 1, 2, 4:
*          coordinates are less than the lower bound; 8 16 32: coordinates
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
        status |= 8;
    }
    if (coordinate_min[kYIndex] < domain_offset[kYIndex]) {
        status |= 2;
    }
    if (coordinate_max[kYIndex] > domain_size[kYIndex] + domain_offset[kYIndex]) {
        status |= 16;
    }
    if (coordinate_min[kZIndex] < domain_offset[kZIndex]) {
        status |= 4;
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
    DefAmrUint flag_temp, flag_refined = NodeBitStatus::kNodeStatus0_;
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
            if (grid_info.map_ptr_interface_layer_info_.find(iter_interface.first)
                == grid_info.map_ptr_interface_layer_info_.end()) {
                grid_info.map_ptr_interface_layer_info_.insert(
                    { iter_interface.first, std::make_shared<InterfaceLayerInfo>() });
            }
            ptr_interface_info = grid_info.map_ptr_interface_layer_info_.at(iter_interface.first).get();
            // interface inside the geometry
            if (ptr_interface_info_lower->vec_inner_coarse2fine_.size() > 0) {
                ptr_interface_info->vec_inner_fine2coarse_.resize(grid_info.k0NumFine2CoarseLayer_);
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
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_inner_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                                map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                            }
                        } else {
                            map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                        }
                    }
                }
                maxlayer = static_cast<int>(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_inner_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                                map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_temp;
                            }
                        } else {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_temp;
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
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusFine2Coarse0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusFine2CoarseM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_outer_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                                map_grid.at(iter_layer_node.first)->flag_status_ = flag_temp;
                            }
                        } else {
                            map_grid.at(iter_layer_node.first)->flag_status_ |= flag_temp;
                        }
                    }
                }
                maxlayer = DefAmrIndexUint(ptr_interface_info_lower->vec_outer_coarse2fine_.size());
                for (DefAmrIndexUint ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    if (ilayer == maxlayer - 1) {
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusCoarse2Fine0_;
                    } else if (ilayer == maxlayer - 2) {
                        flag_temp = flag_refined | NodeBitStatus::kNodeStatusCoarse2FineM1_;
                    } else {
                        flag_temp = flag_refined;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_outer_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                                map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_temp;
                            }
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
 * @brief function to instantiate grid noes for all refinement levels (only works for serial).
 * @param[in] sfbitset_min  minimum space filling code of the current rank.
 * @param[in] sfbitset_max  maximum space filling code of the current rank.
 * @param[in] sfbitset_one_lower_level space filling codes at one lower refinement level.
 */
void GridManagerInterface::InstantiateGridNodeAllLevel(const DefSFBitset sfbitset_min,
    const DefSFBitset sfbitset_max, const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level) {
    DefMap<DefAmrIndexUint> background_occupied;
    // initialize grid information, will be used for instantiate grid nodes
    for (DefAmrIndexUint i_level = 0; i_level <= k0MaxLevel_; ++i_level) {
        vec_ptr_grid_info_.at(i_level)->InitialGridInfo();
    }
    InstantiateOverlapLayerOfRefinementInterface(sfbitset_one_lower_level);
    DefSFCodeToUint code_min = sfbitset_min.to_ullong(), code_max = sfbitset_max.to_ullong();
    for (DefAmrIndexUint i_level = k0MaxLevel_; i_level > 0; --i_level) {
        DefAmrIndexUint i_level_lower = i_level - 1;
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        // initialize grid information
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
                FindAllNodesInACellAtOneLevelLower(bitset_cell_lower, &bitset_all);
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
    InstantiateBackgroundGrid(code_min, code_max, background_occupied);
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
}  // end namespace amrproject
}  // end namespace rootproject
