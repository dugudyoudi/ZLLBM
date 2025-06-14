//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
*/
#include <string>
#include <set>
#include "./auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "grid/sfbitset_aux.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
#include "io/input_parser.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to get background time step of current level.
*/
DefReal MultiTimeSteppingC2F::GetCurrentTimeStep(const DefInt i_level,
    const DefAmrLUint time_step_background, const DefInt time_step_level) const {
    return time_step_background
        + static_cast<DefReal>(time_step_level) / static_cast<DefReal>(TwoPowerN(i_level));
}
/**
* @brief function to get computational cost of a grid node for all levels.
* @return   computational cost of a grid node at each level.
*/
std::vector<DefInt> GridManagerInterface::GetNodeCostAtEachLevel() const {
    std::vector<DefInt> vec_cost;
    for (auto iter_grid : vec_ptr_grid_info_) {
        if (k0MaxLevel_== 0) {
            vec_cost.emplace_back(iter_grid->GetComputationalCost());
        } else {
            vec_cost.emplace_back(iter_grid->GetComputationalCost() * 2 * k0GridDims_);
        }
    }
    return vec_cost;
}
/**
* @brief      function to read grid related parameters.
* @param[in, out] ptr_input_parser pointer to class for input parsing. 
*/
void GridManagerInterface::ReadAndSetupGridParameters(InputParser* const ptr_input_parser) {
    // parameters could be read
    std::vector<DefReal> grid_size;   // compulsory
    std::vector<DefReal> domain_size;   // compulsory
    DefInt max_level = 0;  // optional

    if (ptr_input_parser->GetValue<DefInt>("grid.max_level", &max_level)
        && ptr_input_parser->print_values_when_read_) {
        LogManager::LogInfo("Read and set max refinement level: " + std::to_string(max_level));
    }
    SetMaxLevel(max_level);

    std::string values_str;
    if (ptr_input_parser->GetValue<DefReal>("grid.domain_size", &domain_size)) {
        SetDomainSize(domain_size);
        if (ptr_input_parser->print_values_when_read_) {
            values_str = ptr_input_parser->ValuesToOutputStr<DefReal>(domain_size);
            LogManager::LogInfo("Read and set computational domain size: " + values_str);
        }
    } else {
        LogManager::LogError("Domain size is not given.");
    }

    if (ptr_input_parser->GetValue<DefReal>("grid.grid_size", &grid_size)) {
        SetDomainGridSize(grid_size);
        if (ptr_input_parser->print_values_when_read_) {
            values_str = ptr_input_parser->ValuesToOutputStr<DefReal>(grid_size);
            LogManager::LogInfo("Read and set grid size of background level: " + values_str);
        }
    } else {
        LogManager::LogError("Grid size (grid.grid_size) is not given");
    }
}
/**
* @brief      function to print grid related parameters.
*/
void GridManagerInterface::PrintGridSetting() const {
    printf("Max refinement level: %d \n", GetMaxLevel());

    std::vector<DefReal> vec_tmp;

    vec_tmp.clear();
    vec_tmp = GetDomainSize();
    printf("Size of computational domain:");
    for (const auto& iter : vec_tmp) {
        printf(" %f", iter);
    }
    printf("\n");

    vec_tmp.clear();
    vec_tmp = GetDomainGridSize();
    printf("Grid size at background level:");
    for (const auto& iter : vec_tmp) {
        printf(" %f", iter);
    }
    printf("\n");
}
/**
* @brief   function tp create a tracking grid instance for a given geometry.
* @param[in] i_geo  index of the geometry.
* @param[in] geo_info object containing the geometry information.
*/ 
void GridManagerInterface::CreateTrackingGridInstanceForAGeo(const DefInt i_geo,
    const GeometryInfoInterface& geo_info) {
    std::pair<ECriterionType, DefInt> key_tracking_grid = { ECriterionType::kGeometry, i_geo };
    DefInt geo_level = geo_info.GetLevel();
    TrackingGridInfoCreatorInterface* ptr_creator = geo_info.GetPtrTrackingGridInfoCreatorInterface();
    if (vec_ptr_grid_info_.at(geo_level)->map_ptr_tracking_grid_info_.find(key_tracking_grid)
        == vec_ptr_grid_info_.at(geo_level)->map_ptr_tracking_grid_info_.end()) {
        vec_ptr_grid_info_.at(geo_level)->map_ptr_tracking_grid_info_.insert(
            {key_tracking_grid, ptr_creator->CreateTrackingGridInfo()});
    }
}
/**
* @brief   function to create instance based on the type of solver.
* @param[in] ptr_solver pointer to the solver.
* @param[in] grid_creator instance for creating grid infomation.
*/
void GridManagerInterface::CreateSameGridInstanceForAllLevel(const std::weak_ptr<SolverInterface>& ptr_solver,
    const GridInfoCreatorInterface& grid_creator) {
    if (auto ptr_tmp = ptr_solver.lock()) {
        ptr_tmp->SetPtrToGridManager(this);
    } else {
        LogManager::LogError("Solver is not created.");
        return;
    }
    for (DefInt i_level = 0; i_level < k0MaxLevel_ + 1; ++i_level) {
        vec_ptr_grid_info_.emplace_back(grid_creator.CreateGridInfo());
        GridInfoInterface& grid_ref = *(vec_ptr_grid_info_).back();
        grid_ref.SetPtrToParentGridManager(this);
        grid_ref.SetGridLevel(i_level);
        grid_ref.SetPtrSolver(ptr_solver);
        if (k0GridDims_ == 2) {
            grid_ref.SetPtrSFBitsetAux(dynamic_cast<GridManager2D*>(this));
        } else {
            grid_ref.SetPtrSFBitsetAux(dynamic_cast<GridManager3D*>(this));
        }

        // set computational cost for each node 2^i_level
        grid_ref.SetComputationalCost(TwoPowerN(i_level));
        grid_ref.InitialNotComputeNodeFlag();

        // grid spacing
        std::vector<DefReal> grid_space = std::vector<DefReal>(k0GridDims_, 0.);
        std::vector<DefReal> domain_dx_ = GetDomainDxArrAsVec();
        for (DefInt idim = 0; idim < k0GridDims_; ++idim) {
            grid_space.at(idim) =
                domain_dx_.at(idim) / static_cast<DefReal>(TwoPowerN(i_level));
        }
        grid_ref.SetGridSpace(grid_space);

        // domain boundary related
        CalDomainBoundsAtGivenLevel(i_level,
            &grid_ref.k0VecBitsetDomainMin_, &grid_ref.k0VecBitsetDomainMax_);
        if (this->ptr_func_insert_domain_boundary_ == &GridManagerInterface::InsertCubicDomainBoundary) {
            domain_boundary_type_ = EDomainBoundaryType::kCubic;
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
    DefInt flag_node = ptr_grid_info->CheckIfNodeOutsideCubicDomain(k0GridDims_, bitset_in, *GetPtrToSFBitsetAux());
    if (flag_node >= GridInfoInterface::kFlagInsideDomain_) {
        ptr_grid_info->map_grid_node_.insert({bitset_in, ptr_grid_info->GridNodeCreator()});
        (this->*ptr_func_insert_domain_boundary_)(flag_node, bitset_in, ptr_grid_info);
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
    if (flag_node > GridInfoInterface::kFlagInsideDomain_) {  // flag_node == 0 indicates not on domain boundary
        if ((flag_node & GridInfoInterface::kFlagXMinBoundary_)
            == GridInfoInterface::kFlagXMinBoundary_) {  // on x min boundary
            ptr_grid_info->domain_boundary_min_[kXIndex].insert({bitset_in, 0});
        }
        if ((flag_node & GridInfoInterface::kFlagXMaxBoundary_)
            == GridInfoInterface::kFlagXMaxBoundary_) {  // on x max boundary
            ptr_grid_info->domain_boundary_max_[kXIndex].insert({bitset_in, 0});
        }
        if ((flag_node & GridInfoInterface::kFlagYMinBoundary_)
            == GridInfoInterface::kFlagYMinBoundary_) {  // on y min boundary
            ptr_grid_info->domain_boundary_min_[kYIndex].insert({bitset_in, 0});
        }
        if ((flag_node & GridInfoInterface::kFlagYMaxBoundary_)
            == GridInfoInterface::kFlagYMaxBoundary_) {  // on y max boundary
            ptr_grid_info->domain_boundary_max_[kYIndex].insert({bitset_in, 0});
        }
        if ((flag_node & GridInfoInterface::kFlagZMinBoundary_)
            == GridInfoInterface::kFlagZMinBoundary_) {  // on z min boundary
            ptr_grid_info->domain_boundary_min_[kZIndex].insert({bitset_in, 0});
        }
        if ((flag_node & GridInfoInterface::kFlagZMaxBoundary_)
            == GridInfoInterface::kFlagZMaxBoundary_) {  // on z max boundary
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
void GridManagerInterface::SetNumberOfExtendLayerForGrid(const DefInt i_level,
    const GeometryInfoInterface& geo_info,
    std::vector<DefAmrLUint>* const ptr_inner_layer_neg,
    std::vector<DefAmrLUint>* const ptr_inner_layer_pos,
    std::vector<DefAmrLUint>* const ptr_outer_layer_neg,
    std::vector<DefAmrLUint>* const ptr_outer_layer_pos) {
    std::vector<DefAmrLUint> layer_min(k0GridDims_, k0IntExtendMin_);
    ptr_inner_layer_neg->assign(layer_min.begin(), layer_min.end());
    ptr_inner_layer_pos->assign(layer_min.begin(), layer_min.end());
    ptr_outer_layer_neg->assign(layer_min.begin(), layer_min.end());
    ptr_outer_layer_pos->assign(layer_min.begin(), layer_min.end());
    std::vector<DefInt> geo_extend = geo_info.GetInnerExtend();
    if ((static_cast<DefInt>(geo_extend.size()) > i_level)
        && (geo_extend.at(i_level) > k0IntExtendMin_)) {
        ptr_inner_layer_neg->at(kXIndex) = geo_extend.at(i_level);
        ptr_inner_layer_neg->at(kYIndex) = geo_extend.at(i_level);
        ptr_inner_layer_pos->at(kXIndex) = geo_extend.at(i_level) + 1;
        ptr_inner_layer_pos->at(kYIndex) = geo_extend.at(i_level) + 1;
        if (k0GridDims_ == 3) {
            ptr_inner_layer_neg->at(kZIndex) = geo_extend.at(i_level);
            ptr_inner_layer_pos->at(kZIndex) = geo_extend.at(i_level) + 1;
        }
    } else {
        ptr_inner_layer_neg->at(kXIndex) = k0IntExtendMin_;
        ptr_inner_layer_neg->at(kYIndex) = k0IntExtendMin_;
        ptr_inner_layer_pos->at(kXIndex) = k0IntExtendMin_ + 1;
        ptr_inner_layer_pos->at(kYIndex) = k0IntExtendMin_ + 1;
        if (k0GridDims_ == 3) {
            ptr_inner_layer_neg->at(kZIndex) = k0IntExtendMin_;
            ptr_inner_layer_pos->at(kZIndex) = k0IntExtendMin_ + 1;
        }
    }
    geo_extend = geo_info.GetXExtendNegative();
    if ((static_cast<DefInt>(geo_extend.size()) > i_level) && (geo_extend.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_neg->at(kXIndex) = geo_extend.at(i_level);
    } else {
        ptr_outer_layer_neg->at(kXIndex) = k0IntExtendMin_;
    }
    geo_extend = geo_info.GetXExtendPositive();
    if ((static_cast<DefInt>(geo_extend.size()) > i_level) && (geo_extend.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_pos->at(kXIndex) = geo_extend.at(i_level) + 1;
    } else {
        ptr_outer_layer_pos->at(kXIndex) = k0IntExtendMin_ + 1;
    }
    geo_extend = geo_info.GetYExtendNegative();
    if ((static_cast<DefInt>(geo_extend.size()) > i_level) && (geo_extend.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_neg->at(kYIndex) = geo_extend.at(i_level);
    } else {
        ptr_outer_layer_neg->at(kYIndex) = k0IntExtendMin_;
    }
    geo_extend = geo_info.GetYExtendPositive();
    if ((static_cast<DefInt>(geo_extend.size()) > i_level) && (geo_extend.at(i_level) > k0IntExtendMin_)) {
        ptr_outer_layer_pos->at(kYIndex) = geo_extend.at(i_level) + 1;
    } else {
        ptr_outer_layer_pos->at(kYIndex) = k0IntExtendMin_ + 1;
    }
    if (k0GridDims_ == 3) {
        geo_extend = geo_info.GetZExtendNegative();
        if ((static_cast<DefInt>(geo_extend.size()) > i_level) && (geo_extend.at(i_level) > k0IntExtendMin_)) {
            ptr_outer_layer_neg->at(kZIndex) = geo_extend.at(i_level);
        } else {
            ptr_outer_layer_neg->at(kZIndex) = k0IntExtendMin_;
        }
        geo_extend = geo_info.GetZExtendPositive();
        if ((static_cast<DefInt>(geo_extend.size()) > i_level) && (geo_extend.at(i_level) >
            k0IntExtendMin_)) {
            ptr_outer_layer_pos->at(kZIndex) = geo_extend.at(i_level) + 1;
        } else {
            ptr_outer_layer_pos->at(kZIndex) = k0IntExtendMin_ + 1;
        }
    }
}
/**
 * @brief function to find overlapping layers between grid of adjacent refinement levels based on the outermost coarse layer.
 * @param[in]  layer_coarse_innermost  nodes in the coarse layer overlapping with ptr_layer_fine_outmost.
 * @param[in]  sfbitset_exist   existent nodes.
 * @param[out] ptr_layer_coarse_outer pointer to nodes in the coarse layer overlapping with ptr_layer_fine_inner.
 * @param[out] ptr_layer_fine_outmost pointer to nodes in the outermost fine layer.
 * @param[out] ptr_layer_fine_mid pointer to nodes in the second outermost fine layer.
 * @param[out] ptr_layer_fine_inner pointer to nodes in the third outermost fine layer.
 */
void GridManagerInterface::FindOverlappingLayersBasedOnOutermostCoarse(
    const DefMap<DefInt>& layer_coarse_innermost, const DefMap<DefInt>& sfbitset_exist,
    DefMap<DefInt>* const ptr_layer_coarse_outer, DefMap<DefInt>* const ptr_layer_fine_outmost,
    DefMap<DefInt>* const ptr_layer_fine_mid, DefMap<DefInt>* const ptr_layer_fine_inner) {
#ifdef DEBUG_CHECK_GRID
    if (&layer_coarse_innermost == ptr_layer_coarse_outer) {
        LogManager::LogError("input (layer_coarse_innermost)"
         " should not be the same as output (ptr_layer_coarse_outer) in "
        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID

    std::vector<DefSFBitset> corner_bitsets;
    for (const auto& iter_node : layer_coarse_innermost) {
        FindCornersForNeighbourCells(iter_node.first, &corner_bitsets);
        for (auto& iter_conner : corner_bitsets) {
            IdentifyInterfaceForACell(iter_conner, layer_coarse_innermost, sfbitset_exist,
                ptr_layer_fine_inner, ptr_layer_fine_mid, ptr_layer_fine_outmost);
        }
    }
#ifdef DEBUG_CHECK_GRID
    if (ptr_layer_fine_outmost == ptr_layer_coarse_outer) {
        LogManager::LogError("input (ptr_layer_fine_outmost)"
         " should not be the same as output (ptr_layer_coarse_0) in "
         + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
#endif  // DEBUG_CHECK_GRID
    OverlapLayerFromHighToLow(*ptr_layer_fine_inner, ptr_layer_coarse_outer);
}
/**
 * @brief function to instantiates nodes on the overlap layers of refinement interface for all levels.
 * @param sfbitset_one_lower_level nodes indices encoded by space filling code at one lower level.
 */
void GridManagerInterface::InstantiateOverlapLayerOfRefinementInterface(
    const std::vector<DefMap<DefInt>>& sfbitset_one_lower_level) {
    InterfaceLayerInfo* ptr_interface_info = nullptr;
    InterfaceLayerInfo* ptr_interface_info_lower = nullptr;
    DefInt layer_coarse_outer, layer_coarse_innermost, layer0, layer_m1, layer_m2;
    DefInt flag_tmp, flag_refined = NodeBitStatus::kNodeStatus0_;
    int maxlayer;
    for (DefInt i_level = k0MaxLevel_; i_level > 0; --i_level) {
        GridInfoInterface& grid_info = *(vec_ptr_grid_info_.at(i_level));
        GridInfoInterface& grid_info_lower =
            *(vec_ptr_grid_info_.at(i_level - 1));
        DefMap<std::unique_ptr<GridNode>>& map_grid = grid_info.map_grid_node_;
        DefMap<std::unique_ptr<GridNode>>& map_grid_lower = grid_info_lower.map_grid_node_;
        DefInt num_f2c_layer = grid_info.GetNumFine2CoarseLayer();
#ifdef DEBUG_CHECK_GRID
        if (grid_info_lower.GetNumCoarse2FineLayer() < 2) {
            LogManager::LogError("number of coarse to fine layers at level "
                + std::to_string(i_level - 1) + " is at least 1"
                + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
        if (num_f2c_layer < 3) {
            LogManager::LogError("number of fine to coarse layers at level "
                + std::to_string(i_level) + " is at least 3"
                + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        }
#endif  // DEBUG_CHECK_GRID
        // find interface node at current level
        DefInt num_f2c_ghost_layer = grid_info.GetNumFine2CoarseGhostLayer();
        layer_coarse_outer = 1;
        layer_coarse_innermost = 0;
        layer0 = num_f2c_layer - 1;
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
                ptr_interface_info->vec_inner_fine2coarse_.resize(num_f2c_layer);
                FindOverlappingLayersBasedOnOutermostCoarse(
                    ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_innermost),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower->vec_inner_coarse2fine_.at(layer_coarse_outer),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_inner_fine2coarse_.at(layer_m2));
                // insert node instance
                maxlayer = static_cast<int>(ptr_interface_info->vec_inner_fine2coarse_.size());
                for (int ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    flag_tmp = flag_refined;
                    if (ilayer == maxlayer - 1) {
                        flag_tmp |= NodeBitStatus::kNodeStatusFine2Coarse0_;
                    }
                    if (ilayer >= maxlayer - num_f2c_ghost_layer) {
                        flag_tmp |=  NodeBitStatus::kNodeStatusFine2CoarseGhost_;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_inner_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                                map_grid.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                            }
                        } else {
                            map_grid.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                        }
                    }
                }
                maxlayer = DefInt(ptr_interface_info_lower->vec_inner_coarse2fine_.size());
                for (DefInt ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    flag_tmp = flag_refined;
                    if (ilayer == maxlayer - 1) {
                        flag_tmp |= NodeBitStatus::kNodeStatusCoarse2Fine0_;
                    }
                    if (ilayer >= maxlayer - num_f2c_ghost_layer) {
                        flag_tmp |= NodeBitStatus::kNodeStatusCoarse2FineGhost_;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_inner_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                                map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                            }
                        } else {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                        }
                    }
                }
            }
            // interface outside the geometry
            if (ptr_interface_info_lower->vec_outer_coarse2fine_.size() > 0) {
                ptr_interface_info->vec_outer_fine2coarse_.resize(num_f2c_layer);
                FindOverlappingLayersBasedOnOutermostCoarse(
                    ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(layer_coarse_innermost),
                    sfbitset_one_lower_level.at(i_level),
                    &ptr_interface_info_lower
                    ->vec_outer_coarse2fine_.at(layer_coarse_outer),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer0),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m1),
                    &ptr_interface_info->vec_outer_fine2coarse_.at(layer_m2));
                // insert node instance
                DefInt maxlayer = DefInt(ptr_interface_info->vec_outer_fine2coarse_.size());
                for (DefInt ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    flag_tmp = flag_refined;
                    if (ilayer == maxlayer - 1) {
                        flag_tmp |= NodeBitStatus::kNodeStatusFine2Coarse0_;
                    }
                    if (ilayer >= maxlayer - num_f2c_ghost_layer) {
                        flag_tmp |=  NodeBitStatus::kNodeStatusFine2CoarseGhost_;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info
                        ->vec_outer_fine2coarse_.at(ilayer)) {
                        if (map_grid.find(iter_layer_node.first) == map_grid.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info)) {
                                map_grid.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                            }
                        } else {
                            map_grid.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
                        }
                    }
                }
                maxlayer = DefInt(ptr_interface_info_lower->vec_outer_coarse2fine_.size());
                for (DefInt ilayer = 0; ilayer < maxlayer; ++ilayer) {
                    flag_tmp = flag_refined;
                    if (ilayer == maxlayer - 1) {
                        flag_tmp |= NodeBitStatus::kNodeStatusCoarse2Fine0_;
                    }
                    if (ilayer >= maxlayer - grid_info_lower.GetNumCoarse2FineGhostLayer()) {
                        flag_tmp |= NodeBitStatus::kNodeStatusCoarse2FineGhost_;
                    }
                    for (const auto& iter_layer_node : ptr_interface_info_lower
                        ->vec_outer_coarse2fine_.at(ilayer)) {
                        if (map_grid_lower.find(iter_layer_node.first) == map_grid_lower.end()) {
                            if (InstantiateGridNode(iter_layer_node.first, &grid_info_lower)) {
                                map_grid_lower.at(iter_layer_node.first)->flag_status_ = flag_tmp;
                            }
                        } else {
                            map_grid_lower.at(iter_layer_node.first)->flag_status_ |= flag_tmp;
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
    const DefSFBitset sfbitset_max, const std::vector<DefMap<DefInt>>& sfbitset_one_lower_level) {
    DefMap<DefInt> background_occupied;
    // initialize grid information, will be used for instantiate grid nodes
    for (DefInt i_level = 0; i_level <= k0MaxLevel_; ++i_level) {
        vec_ptr_grid_info_.at(i_level)->InitialGridInfo(k0GridDims_);
    }
    InstantiateOverlapLayerOfRefinementInterface(sfbitset_one_lower_level);
    DefSFCodeToUint code_min, code_max;
    if (k0GridDims_ == 2) {
        SFBitsetAux2D sfbitset_aux;
        code_min = sfbitset_aux.SFBitsetToSFCode(sfbitset_min);
        code_max = sfbitset_aux.SFBitsetToSFCode(sfbitset_max);
    } else {
        SFBitsetAux3D sfbitset_aux;
        code_min = sfbitset_aux.SFBitsetToSFCode(sfbitset_min);
        code_max = sfbitset_aux.SFBitsetToSFCode(sfbitset_max);
    }
    for (DefInt i_level = k0MaxLevel_; i_level > 0; --i_level) {
        DefInt i_level_lower = i_level - 1;
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
        for (const auto & iter_coarse2fine : iter_interfaces.second->vec_inner_coarse2fine_) {
            for (const auto & iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
        }
        for (const auto & iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
            for (const auto & iter_node : iter_coarse2fine) {
                background_occupied.erase(iter_node.first);
            }
        }
    }
    InstantiateBackgroundGrid(code_min, code_max, background_occupied);
}
/**
* @brief   function to identify types of interface nodes
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_coarse_interface   nodes on on interface layer of coarser grid
* @param[out] arr_ptr_layer pointer to maps storing interface layers
*/
void GridManagerInterface::IdentifyInterfaceNodeOnEdge(
    const std::array<DefSFBitset, 2>& arr_bitset_lower,
    const DefSFBitset bitset_mid_higher,
    const SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<DefInt>& node_coarse_interface,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer) {
    bool node0_flag = node_coarse_interface.find(arr_bitset_lower[0]) != node_coarse_interface.end(),
        node1_flag = node_coarse_interface.find(arr_bitset_lower[1]) != node_coarse_interface.end();
    if (node0_flag == node1_flag) {
        if (node0_flag) {
            arr_ptr_layer[2]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[2]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[2]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]), kFlag0_ });
        } else {
            arr_ptr_layer[0]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[0]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[0]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]), kFlag0_ });
        }
    } else {
        if (node0_flag) {
            arr_ptr_layer[2]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[1]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[0]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]), kFlag0_ });
        } else {
            arr_ptr_layer[0]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]), kFlag0_ });
            arr_ptr_layer[1]->insert({ bitset_mid_higher, kFlag0_ });
            arr_ptr_layer[2]->insert({
                sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]), kFlag0_ });
        }
    }
}

/**
* @brief   function to identify types of interface nodes which may cross two refinement levels
* @param[in] arr_bitset_lower   two nodes of an edge at lower level
* @param[in] arr_bitset_current   three nodes on the edge at current level
* @param[out] arr_ptr_layer pointer to maps storing interface layers
* @note noting that it assumes input nodes at lower level should exist on current or lower level
*/
// arr_bitset_lower[0] overlaps arr_bitset_current[0] and arr_bitset_lower[1] overlaps arr_bitset_current[2]
void GridManagerInterface::IdentifyInterfaceNodeOnEdge(
    const std::array<std::pair<DefSFBitset, std::array<bool, 2>>, 2>& arr_bitset_lower,
    const std::array<DefSFBitset, 3>& arr_bitset_current,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer) {
    // first bit: if on current rank and refinement level, second bit: if on given interface
    // at least one node exist in desired grid
    if (arr_bitset_lower.at(0).second.at(0) || arr_bitset_lower.at(1).second.at(0)) {
        if (arr_bitset_lower.at(0).second.at(1) && arr_bitset_lower.at(1).second.at(1)) {
            // both nodes on given interface
            if (arr_bitset_lower.at(0).second.at(0)) {
                arr_ptr_layer[2]->insert({arr_bitset_current.at(0), kFlag0_ });
            }
            if (arr_bitset_lower.at(0).second.at(0) || arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[2]->insert({arr_bitset_current.at(1), kFlag0_});
            }
            if (arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[2]->insert({arr_bitset_current.at(2), kFlag0_ });
            }
            return;
        }
        if ((!arr_bitset_lower.at(0).second.at(1)) && (!arr_bitset_lower.at(1).second.at(1))) {
            // neither node on given interface
            if (arr_bitset_lower.at(0).second.at(0)) {
                arr_ptr_layer[0]->insert({arr_bitset_current.at(0), kFlag0_ });
            }
            if (arr_bitset_lower.at(0).second.at(0) || arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[0]->insert({arr_bitset_current.at(1), kFlag0_});
            }
            if (arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[0]->insert({arr_bitset_current.at(2), kFlag0_ });
            }
            return;
        }
        if ((arr_bitset_lower.at(0).second.at(1)) && (!arr_bitset_lower.at(1).second.at(1))) {
            // node 0 on given interface
            if (arr_bitset_lower.at(0).second.at(0)) {
                arr_ptr_layer[2]->insert({arr_bitset_current.at(0), kFlag0_ });
            }
            if (arr_bitset_lower.at(0).second.at(0) || arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[1]->insert({arr_bitset_current.at(1), kFlag0_});
            }
            if (arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[0]->insert({arr_bitset_current.at(2), kFlag0_ });
            }
            return;
        }
        if ((!arr_bitset_lower.at(0).second.at(1)) && (arr_bitset_lower.at(1).second.at(1))) {
            // node 1 on given interface
            if (arr_bitset_lower.at(0).second.at(0)) {
                arr_ptr_layer[0]->insert({arr_bitset_current.at(0), kFlag0_ });
            }
            if (arr_bitset_lower.at(0).second.at(0) || arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[1]->insert({arr_bitset_current.at(1), kFlag0_});
            }
            if (arr_bitset_lower.at(1).second.at(0)) {
                arr_ptr_layer[2]->insert({arr_bitset_current.at(2), kFlag0_ });
            }
            return;
        }
    }
}
/**
* @brief   function to identify types of interface nodes based on the innermost coarse to fine interface
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_coarse_interface_innermost nodes on the innermost coarse to fine interface
* @param[in] node_current nodes at current level
* @param[out] arr_ptr_layer pointer to maps storing interface layers
* @param[out] ptr_node_coarse_interface_outer pointer to map storing nodes on outer coarse to fine interface
*/
void GridManagerInterface::IdentifyInterfaceNodeOnEdgeInnermost(
    const std::array<DefSFBitset, 2>& arr_bitset_lower,
    const DefSFBitset bitset_mid_higher,
    const SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<DefInt>& node_coarse_interface_innermost,
    const DefMap<std::unique_ptr<GridNode>>& node_current,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer,
    DefMap<DefInt>* const ptr_node_coarse_interface_outer) {
    bool node0_inner =
        node_coarse_interface_innermost.find(arr_bitset_lower[0]) != node_coarse_interface_innermost.end();
    bool node1_inner =
        node_coarse_interface_innermost.find(arr_bitset_lower[1]) != node_coarse_interface_innermost.end();
    DefSFBitset sfbitset_tmp;
    if ((node0_inner) && (node1_inner)) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[2]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        return;
    }
    if ((!node0_inner) && (!node1_inner)) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
        }
        // node may exist on coarse grid but not on fine grid
        ptr_node_coarse_interface_outer->insert({arr_bitset_lower[0], kFlag0_ });
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[0]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
        }
        ptr_node_coarse_interface_outer->insert({arr_bitset_lower[1], kFlag0_ });
        return;
    }
    if (node0_inner && (!node1_inner)) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[1]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
        }
        ptr_node_coarse_interface_outer->insert({arr_bitset_lower[1], kFlag0_ });
        return;
    }
    if ((!node0_inner) && node1_inner) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
        }
        ptr_node_coarse_interface_outer->insert({arr_bitset_lower[0], kFlag0_ });
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[1]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        return;
    }
}
/**
* @brief   function to identify types of interface nodes
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_coarse_interface_previous nodes on coarse to fine interface has been marked
* @param[in] node_coarse_interface_inner nodes on inner coarse to fine interface
* @param[in] node_current nodes at current level
* @param[out] arr_ptr_layer pointer to maps storing interface layers
* @param[out] ptr_node_coarse_interface_outer pointer to map storing nodes on outer coarse to fine interface
*/
void GridManagerInterface::IdentifyInterfaceNodeOnEdge(
    const std::array<DefSFBitset, 2>& arr_bitset_lower,
    const DefSFBitset bitset_mid_higher,
    const SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<DefInt>& node_coarse_interface_previous,
    const DefMap<DefInt>& node_coarse_interface_inner,
    const DefMap<std::unique_ptr<GridNode>>& node_current,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer,
    DefMap<DefInt>* const ptr_node_coarse_interface_outer) {
    DefSFBitset sfbitset_tmp;
    bool node0_previous =
        node_coarse_interface_previous.find(arr_bitset_lower[0]) != node_coarse_interface_previous.end(),
        node1_previous =
        node_coarse_interface_previous.find(arr_bitset_lower[1]) != node_coarse_interface_previous.end(),
        node0_inner = false, node1_inner = false;
    node0_inner = node_coarse_interface_inner.find(arr_bitset_lower[0]) != node_coarse_interface_inner.end();
    node1_inner = node_coarse_interface_inner.find(arr_bitset_lower[1]) != node_coarse_interface_inner.end();

    if ((!node0_previous) && (!node1_previous)) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
            ptr_node_coarse_interface_outer->insert({arr_bitset_lower[0], kFlag0_ });
        }
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[0]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
            ptr_node_coarse_interface_outer->insert({arr_bitset_lower[1], kFlag0_ });
        }
        return;
    }
    if (node0_inner && node1_inner) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[2]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        return;
    }
    if (!node0_previous && node1_inner) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
            ptr_node_coarse_interface_outer->insert({arr_bitset_lower[0], kFlag0_ });
        }
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[1]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        return;
    }
    if (node0_inner && !node1_previous) {
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[2]->insert({sfbitset_tmp, kFlag0_ });
        }
        if (node_current.find(bitset_mid_higher) != node_current.end()) {
            arr_ptr_layer[1]->insert({bitset_mid_higher, kFlag0_});
        }
        sfbitset_tmp = sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1]);
        if (node_current.find(sfbitset_tmp) != node_current.end()) {
            arr_ptr_layer[0]->insert({sfbitset_tmp, kFlag0_ });
            ptr_node_coarse_interface_outer->insert({arr_bitset_lower[1], kFlag0_ });
        }
        return;
    }
}
/**
* @brief   function to identify types of interface nodes which may cross two refinement levels
* @param[in] arr_bitset_lower   two nodes of an edge
* @param[in] bitset_mid_higher   node at the mid point of the edge
* @param[in] node_coarse_interface   nodes on interface layer of coarser grid
* @param[out] arr_ptr_layer pointer to maps storing interface layers
*/
void GridManagerInterface::IdentifyInterfaceNodeOnEdgeAcrossTwoLevels(
    const std::array<std::pair<DefSFBitset, DefInt>, 2>& arr_bitset_lower,
    const DefSFBitset bitset_mid_higher, const SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<DefInt>& node_coarse_interface,
    const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer, DefMap<DefInt>* const ptr_coarse_outer) {
    const DefInt& flag_node_0 = arr_bitset_lower.at(0).second, flag_node_1 = arr_bitset_lower.at(1).second;
    // flag: 1, nodes exist at current level; flag: 2, nodes exist at coarse level;
    // at least one node exist in desired grid
    if (((flag_node_0 == 1) || (flag_node_1 == 1)) && ((flag_node_0 > 0) && (flag_node_1 > 0))) {
        bool node0_flag = node_coarse_interface.find(arr_bitset_lower[0].first) != node_coarse_interface.end(),
            node1_flag = node_coarse_interface.find(arr_bitset_lower[1].first) != node_coarse_interface.end();
        if (node0_flag == node1_flag) {
            if (node0_flag) {
                if (flag_node_0 == 1) {
                    arr_ptr_layer[2]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0].first), kFlag0_ });
                }
                if (flag_node_0 == 1 || flag_node_1 == 1) {
                     arr_ptr_layer[2]->insert({ bitset_mid_higher, kFlag0_ });
                }
                if (flag_node_1 == 1) {
                    arr_ptr_layer[2]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1].first), kFlag0_ });
                }
            } else {
                if (flag_node_0 == 1) {
                    arr_ptr_layer[0]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0].first), kFlag0_ });
                }
                if (flag_node_0 == 1 || flag_node_1 == 1) {
                     arr_ptr_layer[0]->insert({ bitset_mid_higher, kFlag0_ });
                }
                if (flag_node_1 == 1) {
                    arr_ptr_layer[0]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1].first), kFlag0_ });
                }
                ptr_coarse_outer->insert({ arr_bitset_lower[0].first, kFlag0_ });
                ptr_coarse_outer->insert({ arr_bitset_lower[1].first, kFlag0_ });
            }
        } else {
            if (node0_flag) {
                if (flag_node_0 == 1) {
                    arr_ptr_layer[2]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0].first), kFlag0_ });
                }
                if (flag_node_0 == 1 || flag_node_1 == 1) {
                     arr_ptr_layer[1]->insert({ bitset_mid_higher, kFlag0_ });
                }
                if (flag_node_1 == 1) {
                    arr_ptr_layer[0]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1].first), kFlag0_ });
                }
                ptr_coarse_outer->insert({ arr_bitset_lower[1].first, kFlag0_ });
            } else {
                if (flag_node_0 == 1) {
                    arr_ptr_layer[0]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[0].first), kFlag0_ });
                }
                if (flag_node_0 == 1 || flag_node_1 == 1) {
                     arr_ptr_layer[1]->insert({ bitset_mid_higher, kFlag0_ });
                }
                if (flag_node_1 == 1) {
                    arr_ptr_layer[2]->insert({
                        sfbitset_aux.SFBitsetToNHigherLevelVir(1, arr_bitset_lower[1].first), kFlag0_ });
                }
                ptr_coarse_outer->insert({ arr_bitset_lower[0].first, kFlag0_ });
            }
        }
    }
}
/**
* @brief function to set k0TimeSteppingOrder_ as multi-stepping scheme
*        from the background level to the finest (ratio of 2).
*/
MultiTimeSteppingC2F::MultiTimeSteppingC2F(const DefInt max_level) {
    if (max_level == 0) {
        k0TimeSteppingOrder_ = {0};
    } else if (max_level == 1) {
        k0TimeSteppingOrder_ = {0, 1, 1};
    } else {
        std::vector<DefSizet> accumulate_t(max_level + 1, 0);
        DefInt i_level = 1;
        accumulate_t[0] = TwoPowerN(max_level);
        k0TimeSteppingOrder_ = {0};
        while (i_level > 0) {
            if (accumulate_t[i_level- 1] == accumulate_t[i_level]) {
                --i_level;
            } else {
                k0TimeSteppingOrder_.emplace_back(i_level);
                accumulate_t[i_level] += TwoPowerN(max_level - i_level);
                if (i_level + 1 < max_level) {
                    ++i_level;
                } else {
                    ++i_level;
                    k0TimeSteppingOrder_.emplace_back(i_level);
                    k0TimeSteppingOrder_.emplace_back(i_level);
                    accumulate_t[i_level] += TwoPowerN(max_level - i_level + 1);
                }
            }
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
