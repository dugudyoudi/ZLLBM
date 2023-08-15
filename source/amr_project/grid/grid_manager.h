//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file def_libs.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief contain definations and libraries will be used in all modules.
*/

#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
#include <concepts>
#include <array>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include "../defs_libs.h"
#include "./solver_info_interface.h"
#include "grid/grid_info_interface.h"
#include "grid/sfbitset_aux.h"
#include "criterion/geometry_info_interface.h"
#ifdef DEBUG_UNIT_TEST
#include "../../googletest-main/googletest/include/gtest/gtest_prod.h"
#endif  // DEBUG_UNIT_TEST
namespace rootproject {
namespace amrproject {
template<typename InterfaceInfo>
concept InterfaceInfoHasType = requires(
    InterfaceInfo type_interface) {
    type_interface.node_type_;
};
/**
* @class GridManager
* @brief class used to manage adaptive mesh refinement.
* @note parameters with prefix 'k0'will be considered
*       as constants after project initialization.
*/
class GridManagerInterface{
 public:
    /**<Constant: index corresponding to x, y, and z in vectors */
    // 2 or 3
    DefAmrIndexUint  k0GridDims_ = 0;  ///< dimension
    // 9 for 2D, 27 for 3D
    DefAmrIndexUint  k0NumNeighbors_ = 0;  ///< number of neighboring nodes
    // 0 refers to uniform mesh
    DefAmrIndexUint  k0MaxLevel_ = 0;  ///< maximum levels of refinement

    /* this setup is used to avoid issues caused by rapid change of
    refinement level in a small region. Number of layers is counted
    at the lower refinement level.*/
    DefAmrIndexUint k0IntExtendMin_ = 3;  ///< a minimum number of extending layers
    DefAmrIndexUint k0IntExtendRemoved_ = 1;
    ///< a number of extending layers when inside nodes will be removed
    /* number of nodes extended from position of given criterion,
    e.g. solid boundary and free surface*/

    // overlapping region between different refinement levels
    const DefAmrIndexUint kOverlappingOutmostM1_ = 1;
    ///< index of the outmost layer of the overlapping region
    const DefAmrIndexUint kOverlappingOutmostM2_ = 2;
    ///< index of the second outmost layer of the overlapping region

    // grid status
    const DefAmrUint kNodeStatusExist_ = 1;
    const DefAmrUint kNodeStatusCoarse2Fine0_ = 1 << 1;
    const DefAmrUint kNodeStatusCoarse2FineM1_ = 1 << 2;
    const DefAmrUint kNodeStatusFine2Coarse0_ = 1 << 3;
    const DefAmrUint kNodeStatusFine2CoarseM1_ = 1 << 4;
    const DefAmrUint kNodeStatusGhostCommunication_ = 1 << 5;
    const DefAmrIndexUint kFlagSize0_ = 0;  // flag initialize size as 0

    // computational domain related information
    DefSFBitset k0SFBitsetDomainMin_ = 0;
    DefSFBitset k0SFBitsetDomainMax_ = ~0;

    // list of creators to construct classes
    std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>> vec_ptr_tracking_info_creator_;

    /*the method is used to create grid instance for all refinement levels.
      Give different methods for each grid is possible by using other functions
      instead of SetGridInfoInterfaceAllLevels(EGridNodeType node_type) during
      initialization in function SetGridParameters(). */
    std::vector<std::shared_ptr<GridInfoInterface>> vec_ptr_grid_info_;
    void CreateSameGridInstanceForAllLevel(
        GridInfoCreatorInterface* ptr_grid_creator);

    void DefaultInitialization(const DefAmrIndexUint max_level);
    int CheckIfPointOutsideDomain(
        const std::array<DefReal, 2>& coordinate_min,
        const std::array<DefReal, 2>& coordinate_max,
        const std::array<DefReal, 2>& domain_min,
        const std::array<DefReal, 2>& domain_max) const;
    int CheckIfPointOutsideDomain(
        const std::array<DefReal, 3>& coordinate_min,
        const std::array<DefReal, 3>& coordinate_max,
        const std::array<DefReal, 3>& domain_min,
        const std::array<DefReal, 3>& domain_max) const;
    void SetNumberOfExtendLayerForGrid(const DefAmrIndexUint i_level,
        const GeometryInfoInterface& geo_info,
        std::vector<DefAmrIndexLUint>* const num_extend_inner_layer_neg,
        std::vector<DefAmrIndexLUint>* const num_extend_inner_layer_pos,
        std::vector<DefAmrIndexLUint>* const num_extend_outer_layer_neg,
        std::vector<DefAmrIndexLUint>* const num_extend_outer_layer_pos);

    virtual void PrintGridInfo(void) const = 0;
    virtual void SetGridParameters(void) = 0;

    virtual SFBitsetAuxInterface* GetSFBitsetAuxPtr() = 0;
    virtual std::vector<DefReal> GetDomainDxArrAsVec() const = 0;

    virtual bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefAmrIndexUint>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) = 0;

    // node related functions
    virtual void FindAllNeighborsSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const = 0;
    virtual void FindAllNodesInACellAtLowerLevel(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const = 0;
    virtual DefSFBitset NodeAtNLowerLevel(
        const DefAmrIndexUint n_level, const DefSFBitset& bitset_in) const = 0;
    virtual DefAmrIndexLUint CalNumOfBackgroundNode() const = 0;
    virtual bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const = 0;

    // generate grid
    void GenerateGridFromHighToLowLevelSerial(const std::vector<std::shared_ptr
        <GeometryInfoInterface>>&vec_geo_info,
        std::vector<DefMap<DefAmrIndexUint>>* const ptr_sfbitset_one_lower_level);
    void InstantiateGridNodeAllLevel(const DefSFBitset sfbitset_min, const DefSFBitset sfbitset_max,
        const DefMap<DefAmrIndexUint>& partitioned_interface_background,
        const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level);

    template<InterfaceInfoHasType InterfaceInfo>
    DefSizet CheckExistenceOfTypeAtGivenLevel(
        const DefAmrIndexUint i_level, const std::string& type,
        const std::vector<std::shared_ptr<InterfaceInfo>>& vec_ptr_interface) const;

    virtual ~GridManagerInterface() {}
    GridManagerInterface() {
        static_assert(sizeof(DefSFBitset) == sizeof(DefSFCodeToUint),
            "size of DefSFBitset and DefSFCodeToU must be the same");
    }

 protected:
    //// Generate grid
    const DefAmrUint kFlag0_ = 0;  // flag to initialize status of a node

    // functions to generate grid
    void InstantiateOverlapLayerOfRefinementInterface(
        const std::vector<DefMap<DefAmrIndexUint>>& sfbitset_one_lower_level);

    virtual void OverlapLayerFromHighToLow(
        const DefMap<DefAmrUint>& layer_high_level,
        DefMap<DefAmrUint>* const ptr_layer_low_level) = 0;
    virtual void InstantiateBackgroundGrid(const DefSFBitset bitset_min,
        const DefSFBitset bitset_max, const DefMap<DefAmrIndexUint>& map_occupied) = 0;
    virtual bool CheckCoincideBackground(const DefAmrIndexUint i_level,
        const DefSFBitset& bitset_higher,
        DefSFBitset* const ptr_bitset) const = 0;
    // only run on rank 0
    int FloodFillForInAndOut(const DefSFBitset& sfbitset_start,
        const DefMap<DefAmrUint>& map_nodes_exist,
        DefMap<DefAmrUint>* const ptr_map_nodes_outside,
        DefMap<DefAmrUint>* const ptr_map_nodes_inside) const;

 private:
    const DefAmrUint K0IMaxFloodFill_ = 90000;
    ///<  maximum iteration for flood fill

    void ExtendGivenNumbOfLayer(DefAmrIndexUint i_level,
        const std::vector<DefAmrIndexLUint> num_extend_neg,
        const std::vector<DefAmrIndexLUint> num_extend_pos,
        const DefMap<DefAmrUint>& map_start_layer,
        DefMap<DefAmrIndexUint>* const ptr_map_exist,
        DefMap<DefAmrUint>* const ptr_map_outmost) const;
    void FindOverlappingLayersBasedOnOutermostCoarse(
        const DefMap<DefAmrUint>& layer_coarse_0,
        const DefMap<DefAmrIndexUint>& sfbitset_exist,
        DefMap<DefAmrUint>* const ptr_layer_coarse_m1,
        DefMap<DefAmrUint>* const ptr_layer_fine_0,
        DefMap<DefAmrUint>* const ptr_layer_fine_m1,
        DefMap<DefAmrUint>* const ptr_layer_fine_m2);

    virtual void ResetExtendLayerBasedOnDomainSize(
        const DefAmrIndexUint i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefAmrIndexLUint>* const ptr_vec_extend_neg,
        std::vector<DefAmrIndexLUint>* const ptr_vec_extend_pos) const = 0;
    virtual void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const = 0;
    virtual void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets) const = 0;
    virtual void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefAmrUint>& node_coarse_outmost, const DefMap<DefAmrIndexUint>& node_exist,
        DefMap<DefAmrUint>* const ptr_inner_layer,
        DefMap<DefAmrUint>* const ptr_mid_layer,
        DefMap<DefAmrUint>* const ptr_outer_layer) = 0;

    // functions to generate grid (Serial)
    virtual void GenerateGridNodeNearTrackingNode(const DefAmrIndexUint i_level,
        const std::pair<ECriterionType, DefAmrIndexUint>& tracking_grid_key,
        DefMap<DefAmrUint>* const ptr_map_node_temp) const = 0;
    virtual void IdentifyTypeOfLayerByFloodFill(
        const DefAmrIndexUint i_level, const DefAmrIndexUint i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefAmrUint>& map_nodes_exist,
        DefMap<DefAmrUint>* const ptr_map_nodes_outside,
        DefMap<DefAmrUint>* const ptr_map_nodes_inside) const = 0;
    virtual void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk) const = 0;
    virtual void ExtendOneLayerGrid(
        const DefMap<DefAmrUint>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefAmrUint>* const ptr_map_output_layer,
        DefMap<DefAmrIndexUint>* const ptr_map_exist,
        DefMap<DefAmrUint>* const ptr_outmost_layer,
        std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_max) const = 0;
    virtual void FindInterfaceBetweenGrid(
        const DefAmrIndexUint i_level, const DefMap<DefAmrUint>& map_outmost_layer,
        DefMap<DefAmrIndexUint>* const map_exist,
        DefMap<DefAmrUint>* const ptr_interface_outmost,
        DefMap<DefAmrIndexUint>* const ptr_layer_lower_level,
        DefMap<DefAmrUint>* const ptr_layer_lower_level_outer) = 0;

#ifdef ENABLE_MPI

 public:
    virtual void GetNLevelCorrespondingOnes(
        const DefAmrIndexUint i_level, std::vector<DefSFBitset>* const ptr_last_ones) const = 0;
    virtual bool CheckNodeOnOuterBoundaryOfBackgroundCell(DefAmrIndexUint i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max, const DefSFBitset bitset_in,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefAmrIndexUint>& partitioned_interface_background) const = 0;
    virtual void GetMinM1AtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const = 0;
    virtual void GetMaxP1AtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const = 0;
    virtual void SearchForGhostLayerForMinNMax(const DefSFBitset bitset_in,
        const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        std::vector<DefSFBitset>* const ptr_vec_ghost_layer) const = 0;
#endif  // ENABLE_MPI

#ifdef DEBUG_UNIT_TEST

 private:
     // gtest to access private member functions
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
     FRIEND_TEST(GridManagerGeneration2D, ExtendNumbOfLayerAndFindInterface);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
     FRIEND_TEST(GridManagerGeneration3D, ExtendNumbOfLayerAndFindInterface);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
#endif  // DEBUG_UNIT_TEST
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @class GridInfo2DInterface
* @brief class used to store 2D grid information
* @note
* @date  2022-11-15
*/
class GridManager2D :public  GridManagerInterface, public SFBitsetAux2D {
 public:
    std::array<DefAmrIndexLUint, 2> k0IntOffset_ = {1, 1};
    ///< number of offset background node
    std::array<DefReal, 2> k0DomainSize_{};  ///< domain size
    std::array<DefReal, 2> k0DomainDx_{};  ///< grid space
    /* offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (k0IntOffset_ * kDomainDx),
    and the default value is 1. */
    std::array<DefReal, 2> k0RealOffset_{};  ///< k0IntOffset_ * kDomainDx
    std::array<DefAmrIndexLUint, 2> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0IntOffset_ is included in k0MaxIndexOfBackgroundNode_

    void PrintGridInfo(void) const override;
    void SetGridParameters(void) override;
    SFBitsetAuxInterface* GetSFBitsetAuxPtr() override {
        return dynamic_cast<SFBitsetAuxInterface*>(this);
    }
    std::vector<DefReal> GetDomainDxArrAsVec() const final {
        return { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex] };
    };
    // node related function
    void FindAllNeighborsSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const override;
    void FindAllNodesInACellAtLowerLevel(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const override;
    DefSFBitset NodeAtNLowerLevel(
        const DefAmrIndexUint n_level, const DefSFBitset& biset_in) const override;
    DefAmrIndexLUint CalNumOfBackgroundNode() const override {
        return (k0MaxIndexOfBackgroundNode_[0] - k0IntOffset_[0] + 1) *
            (k0MaxIndexOfBackgroundNode_[1] - k0IntOffset_[1] + 1);
    }
    bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const override;

    GridManager2D() {
        k0GridDims_ = 2;
        k0NumNeighbors_ = 9;
    }

 protected:
    void InstantiateBackgroundGrid(const DefSFBitset bitset_min,
        const DefSFBitset bitset_max, const DefMap<DefAmrIndexUint>& map_occupied) override;
    void OverlapLayerFromHighToLow(
        const DefMap<DefAmrUint>& layer_high_level,
        DefMap<DefAmrUint>* const ptr_layer_low_level) override;

    bool CheckCoincideBackground(const DefAmrIndexUint i_level,
        const DefSFBitset& bitset_higher,
        DefSFBitset* const ptr_bitset) const override;
    void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const override;
    void ResetExtendLayerBasedOnDomainSize(
        const DefAmrIndexUint i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefAmrIndexLUint>* const ptr_vec_extend_neg,
        std::vector<DefAmrIndexLUint>* const ptr_vec_extend_pos) const override;
    void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
       std::vector<DefSFBitset>* const ptr_bitsets)  const override;
    void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefAmrUint>& node_coarse_outmost, const DefMap<DefAmrIndexUint>& node_exist,
        DefMap<DefAmrUint>* const ptr_inner_layer,
        DefMap<DefAmrUint>* const ptr_mid_layer,
        DefMap<DefAmrUint>* const ptr_outer_layer) override;
    bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefAmrIndexUint>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets)override;
    // functions to generate grid (Serial)
    void GenerateGridNodeNearTrackingNode(const DefAmrIndexUint i_level,
        const std::pair<ECriterionType, DefAmrIndexUint>& tracking_grid_key,
        DefMap<DefAmrUint>* const ptr_map_node_temp) const override;
    void IdentifyTypeOfLayerByFloodFill(
        const DefAmrIndexUint i_level, const DefAmrIndexUint i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefAmrUint>& map_nodes_exist,
        DefMap<DefAmrUint>* const ptr_map_nodes_outside,
        DefMap<DefAmrUint>* const ptr_map_nodes_inside) const override;
    void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk) const override;
    void ExtendOneLayerGrid(
        const DefMap<DefAmrUint>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefAmrUint>* const ptr_map_output_layer,
        DefMap<DefAmrIndexUint>* const ptr_map_exist,
        DefMap<DefAmrUint>* const ptr_outmost_layer,
        std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_max)
        const override;
    void FindInterfaceBetweenGrid(
        const DefAmrIndexUint i_level, const DefMap<DefAmrUint>& map_outmost_layer,
        DefMap<DefAmrIndexUint>* const map_exist,
        DefMap<DefAmrUint>* const ptr_interface_outmost,
        DefMap<DefAmrIndexUint>* const ptr_layer_lower_level,
        DefMap<DefAmrUint>* const ptr_layer_lower_level_outer) override;

 private:
    DefAmrUint kFlagCurrentNodeXNeg_ = 1, kFlagCurrentNodeXPos_ = 1 << 1,
        kFlagCurrentNodeYNeg_ = 1 << 2, kFlagCurrentNodeYPos_ = 1 << 3;
    DefAmrUint FindAllNeighborsWithSpecifiedDirection(const DefSFBitset bitset_in,
        const std::array<bool, 2>& bool_neg, const std::array<bool, 2>& bool_pos,
        std::vector <DefSFBitset>* const ptr_vec_neighbors) const;
    void IdentifyInterfaceNodeOnEdge(
        const std::array<DefSFBitset, 2>& arr_bitset_lower,
        const DefSFBitset bitset_mid_higher,
        const DefMap<DefAmrUint>& node_outmost,
        const std::array<DefMap<DefAmrUint>* const, 3>& arr_ptr_layer);

#ifdef ENABLE_MPI

 public:
    void GetNLevelCorrespondingOnes(
        const DefAmrIndexUint i_level, std::vector<DefSFBitset>* const ptr_last_ones) const final;
    bool CheckNodeOnOuterBoundaryOfBackgroundCell(DefAmrIndexUint i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max, const DefSFBitset bitset_in,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefAmrIndexUint>& partitioned_interface_background) const final;
    void GetMinM1AtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const final;
    void GetMaxP1AtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const final;
    void SearchForGhostLayerForMinNMax(const DefSFBitset bitset_in,
        const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        std::vector<DefSFBitset>* const ptr_vec_ghost_layer) const final;
#endif  // ENABLE_MPI

#ifdef DEBUG_UNIT_TEST

 private:
        // gtest to access private member functions
        FRIEND_TEST(GridManagerGeneration2D, ResetExtendLayer);
        FRIEND_TEST(GridManagerGeneration2D, IdentifyLayerThroughFloodFill);
        FRIEND_TEST(GridManagerGeneration2D, ExtendGridOneLayer);
        FRIEND_TEST(GridManagerGeneration2D, ExtendNumbOfLayerAndFindInterface);
        FRIEND_TEST(GridManagerGeneration2D, FindInterfaceCurrentLevel);
#endif  // DEBUG_UNIT_TEST
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GridManager3D :public  GridManagerInterface, public SFBitsetAux3D {
 public:
    std::array<DefAmrIndexLUint, 3> k0IntOffset_ = { 1, 1, 1};
    ///< number of offset background node
    std::array<DefReal, 3> k0DomainSize_{};  ///< domain size
    std::array<DefReal, 3> k0DomainDx_{};  ///< grid space
    /* offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (k0IntOffset_ * kDomainDx),
    and the default value is 1. */
    std::array<DefReal, 3> k0RealOffset_{};  ///< k0IntOffset_ * kDomainDx
    std::array<DefAmrIndexLUint, 3> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0IntOffset_ is included in k0MaxIndexOfBackgroundNode_

    void PrintGridInfo(void) const override;
    void SetGridParameters(void) override;

    SFBitsetAuxInterface* GetSFBitsetAuxPtr() override {
        return dynamic_cast<SFBitsetAuxInterface*>(this);
    }
    std::vector<DefReal> GetDomainDxArrAsVec() const final {
        return { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex], k0DomainDx_[kZIndex] };
    };

    // node related function
    void FindAllNeighborsSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const override;
    void FindAllNodesInACellAtLowerLevel(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const override;
    DefSFBitset NodeAtNLowerLevel(
        const DefAmrIndexUint n_level, const DefSFBitset& biset_in) const override;
    DefAmrIndexLUint CalNumOfBackgroundNode() const override {
        return (k0MaxIndexOfBackgroundNode_[0] - k0IntOffset_[0] + 1)
         * (k0MaxIndexOfBackgroundNode_[1] - k0IntOffset_[1] + 1)
         * (k0MaxIndexOfBackgroundNode_[2] - k0IntOffset_[2] + 1);
    }
    bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const override;

    GridManager3D() {
        k0GridDims_ = 3;
        k0NumNeighbors_ = 27;
    }

 protected:
    void InstantiateBackgroundGrid(const DefSFBitset bitset_min,
        const DefSFBitset bitset_max, const DefMap<DefAmrIndexUint>& map_occupied) override;
    void OverlapLayerFromHighToLow(
        const DefMap<DefAmrUint>& layer_high_level,
        DefMap<DefAmrUint>* const ptr_layer_low_level) override;

    bool CheckCoincideBackground(const DefAmrIndexUint i_level,
         const DefSFBitset& bitset_higher,
         DefSFBitset* const ptr_bitset) const override;
    void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const override;
    void ResetExtendLayerBasedOnDomainSize(
        const DefAmrIndexUint i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefAmrIndexLUint>* const ptr_vec_extend_neg,
        std::vector<DefAmrIndexLUint>* const ptr_vec_extend_pos) const override;
    void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets)  const override;
    void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefAmrUint>& node_coarse_outmost, const DefMap<DefAmrIndexUint>& node_exist,
        DefMap<DefAmrUint>* const ptr_inner_layer,
        DefMap<DefAmrUint>* const ptr_mid_layer,
        DefMap<DefAmrUint>* const ptr_outer_layer) override;
    bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefAmrIndexUint>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) override;
    // functions to generate grid (Serial)
    void GenerateGridNodeNearTrackingNode(const DefAmrIndexUint i_level,
        const std::pair<ECriterionType, DefAmrIndexUint>& tracking_grid_key,
        DefMap<DefAmrUint>* const ptr_map_node_temp) const override;
    void IdentifyTypeOfLayerByFloodFill(
        const DefAmrIndexUint i_level, const DefAmrIndexUint i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefAmrUint>& map_nodes_exist,
        DefMap<DefAmrUint>* const ptr_map_nodes_outside,
        DefMap<DefAmrUint>* const ptr_map_nodes_inside) const override;
    void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk)const  override;
    void ExtendOneLayerGrid(
        const DefMap<DefAmrUint>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefAmrUint>* const ptr_map_output_layer,
        DefMap<DefAmrIndexUint>* const ptr_map_exist,
        DefMap<DefAmrUint>* const ptr_outmost_layer,
        std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefAmrUint>>* const ptr_vector_boundary_max)
        const override;
    void FindInterfaceBetweenGrid(
        const DefAmrIndexUint i_level, const DefMap<DefAmrUint>& map_outmost_layer,
        DefMap<DefAmrIndexUint>* const map_exist,
        DefMap<DefAmrUint>* const ptr_interface_outmost,
        DefMap<DefAmrIndexUint>* const ptr_layer_lower_level,
        DefMap<DefAmrUint>* const ptr_layer_lower_level_outer) override;

 private:
    DefAmrUint kFlagCurrentNodeXNeg_ = 1, kFlagCurrentNodeXPos_ = 1 << 1,
     kFlagCurrentNodeYNeg_ = 1 << 2, kFlagCurrentNodeYPos_ = 1 << 3,
     kFlagCurrentNodeZNeg_ = 1 << 4, kFlagCurrentNodeZPos_ = 1 << 5;
    DefAmrUint FindAllNeighborsWithSpecifiedDirection(const DefSFBitset bitset_in,
        const std::array<bool, 3>& bool_neg, const std::array<bool, 3>& bool_pos,
        std::vector <DefSFBitset>* const ptr_vec_neigbours) const;
    void IdentifyInterfaceNodeOnEdge(
        const std::array<DefSFBitset, 2>& arr_bitset_lower,
        const DefSFBitset bitset_mid_higher,
        const DefMap<DefAmrUint>& node_outmost,
        const std::array<DefMap<DefAmrUint>* const, 3>& arr_ptr_layer);
    void IdentifyInterfaceNodeDiagonal(
        const std::array<DefSFBitset, 4>& arr_bitset_lower,
        const DefSFBitset bitset_center_higher,
        const DefMap<DefAmrUint>& node_outmost,
        const std::array<DefMap<DefAmrUint>* const, 3>& arr_ptr_layer);

#ifdef ENABLE_MPI

 public:
    void GetNLevelCorrespondingOnes(
        const DefAmrIndexUint i_level, std::vector<DefSFBitset>* const ptr_last_ones) const final;
    bool CheckNodeOnOuterBoundaryOfBackgroundCell(DefAmrIndexUint i_level,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max, const DefSFBitset bitset_in,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        const std::vector<DefSFBitset>& bitset_level_ones,
        const DefMap<DefAmrIndexUint>& partitioned_interface_background) const final;
    void GetMinM1AtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_min_m1_bitsets) const final;
    void GetMaxP1AtGivenLevel(const DefAmrIndexUint i_level,
        std::vector<DefSFBitset>* const ptr_max_p1_bitsets) const final;
    void SearchForGhostLayerForMinNMax(const DefSFBitset bitset_in,
        const DefAmrIndexUint num_of_ghost_layers, const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const std::vector<DefSFBitset>& domain_min_m1_n_level,
        const std::vector<DefSFBitset>& domain_max_p1_n_level,
        std::vector<DefSFBitset>* const ptr_vec_ghost_layer) const final;
#endif  // ENABLE_MPI
#ifdef DEBUG_UNIT_TEST

 private:
     // gtest to access private member functions
     FRIEND_TEST(GridManagerGeneration3D, ResetExtendLayer);
     FRIEND_TEST(GridManagerGeneration3D, IdentifyLayerThroughFloodFill);
     FRIEND_TEST(GridManagerGeneration3D, ExtendGridOneLayer);
     FRIEND_TEST(GridManagerGeneration3D, ExtendNumbOfLayerAndFindInterface);
     FRIEND_TEST(GridManagerGeneration3D, FindInterfaceCurrentLevel);
#endif  // DEBUG_UNIT_TEST
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to check if node type exists at a given level.
* @param[in]  i_level level of refinement.
* @param[in]  type the given type.
* @param[in]  vec_ptr_interface pointer to vector of information
*             with a specified node type.
* @return  if return 0, then the node type does not exist in the vector 
*          elements; else return the (indices + 1) of the element with
*          the given type.
*/
template<InterfaceInfoHasType InterfaceInfo>
DefSizet GridManagerInterface::CheckExistenceOfTypeAtGivenLevel(
    const DefAmrIndexUint i_level, const std::string& type,
    const std::vector<std::shared_ptr<InterfaceInfo>>&
    vec_ptr_interface) const {
    DefSizet iter_count = 0;
    if (vec_ptr_interface.empty()) {
        return 0;
    } else {
        for (const auto& iter : vec_ptr_interface) {
            ++iter_count;
            if ((iter->i_level_ == i_level) && (iter->node_type_ == type)) {
                return iter_count;
            }
        }
    }
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
