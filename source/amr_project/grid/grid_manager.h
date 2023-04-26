//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file def_libs.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief contain definations and libaries will be used in all modules.
*/

#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
#include <concepts>
#include <array>
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
concept InterfaceInfoHasType =requires (
    InterfaceInfo type_inferface){
    type_inferface.node_type_;
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
    DefUint  k0GridDims_ = 0;  ///< dimension
    // 9 for 2D, 27 for 3D
    DefUint  k0NumNeighbours_ = 0;  ///< number of neighbouring nodes
    // 0 refers to uniform mesh
    DefSizet  k0MaxLevel_ = 0;  ///< maximum levels of refinement
    
    /* this setup is used to avoid issues caused by rapid change of
    refinement level in a small region. Number of layers is counted
    at the lower refinement level.*/
    DefUint k0IntExtendMin_ = 3;  ///< a minimum number of extending layers
    DefUint k0IntExtendRemoved_ = 1;
    ///< a number of extending layers when inside nodes will be removed
    /* number of nodes extended from position of given criterion,
    e.g. solid boundary and free surface*/

    // overlapping region between differen refinement level
    const DefUint kOverlappingOutmostM1_ = 1;
    ///< index of the outmost layer of the overlapping region 
    const DefUint kOverlappingOutmostM2_ = 2;
    ///< index of the second outmost layer of the overlapping region

    // grid status
    const DefUint kFlagExist_ = 1 << 2;
    const DefUint kFlagCoarse2Fine0_ = 1 << 4;
    const DefUint kFlagCoarse2FineM1_ = 1 << 6;
    const DefUint kFlagFine2Coarse0_ = 1 << 7;
    const DefUint kFlagFine2CoarseM1_ = 1 << 8;

    // grid type
    /*the method is used to create grid instance for all refinement levels.
      Give different methods for each grid is possible by using other functions
      instead of SetGridInfoInterfaceAllLevels(EGridNodeType node_type) during
      initilization in function SetGridParameters(). */
    std::vector<std::shared_ptr<GridInfoInterface>> vec_ptr_grid_info_;
    void CreateSameGridInstanceForAllLevel(
        GridInfoCreatorInterface* ptr_grid_creator);

    void DefaultInitialization(const DefSizet max_level);
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
    void SetNumberOfExtendLayerForGrid(const DefSizet i_level,
        const GeometryInfoInterface& geo_info,
        std::vector<DefLUint>* const num_extend_inner_layer_neg,
        std::vector<DefLUint>* const num_extend_inner_layer_pos,
        std::vector<DefLUint>* const num_extend_outer_layer_neg,
        std::vector<DefLUint>* const num_extend_outer_layer_pos);
    
    virtual void PrintGridInfo(void) const = 0;
    virtual void SetGridParameters(void) = 0;

    virtual SFBitsetAuxInterface* GetSFBitsetAuxPtr() = 0;
    virtual std::vector<DefReal> GetDomainDxArrAsVec() const = 0;

    virtual bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefUint>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) = 0;

    // node related functions
    virtual void FindAllNeighboursSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbours) const = 0;
    virtual void FindAllNodesInACellAtLowerLevel(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const = 0;
    virtual DefSFBitset NodeAtNLowerLevel(
        const DefSizet n_level, const DefSFBitset& biset_in) const = 0;
    virtual DefLUint CalNumOfBackgroundNode() const = 0;
    virtual bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const = 0;
    virtual void TraverseBackgroundForPartition(
        const std::vector<DefLUint>& rank_load,
        const DefMap<DefUint>& background_occupied,
        std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const = 0;

    // generate grid
    void GenerateInitialMeshBasedOnGeoSerial(
        const std::vector<std::shared_ptr<GeometryInfoInterface>>&
        vec_geo_info);
    void GenerateGridFromHighToLowLevelSerial(const std::vector<std::shared_ptr
        <GeometryInfoInterface>>&vec_geo_info,
        std::vector<DefMap<DefUint>>* const ptr_sfbitset_one_lower_level);
    void InstantiateGridNodeAllLevelSerial(
        const std::vector<DefMap<DefUint>>& sfbitset_one_lower_level);

    template<InterfaceInfoHasType InterfaceInfo>
    //template<typename Type, typename InterfaceInfo>
    DefSizet CheckExistanceOfTypeAtGivenLevel(
        const DefSizet i_level, const std::string& type,
        const std::vector<std::shared_ptr<InterfaceInfo>>&
        vec_ptr_interface) const;

    virtual ~GridManagerInterface() {}

 protected:
    //// Generate grid
    const DefUint kFlag0_ = 0;  // flag to initialize status of a node
    const DefUint kFlagTrackingCell_ = 1;  // flag to initialize status of a node
    // node on the outmost interface
    const DefUint kFlagGridInterfaceOutermost_ = 1 << 2;

    // functions to generate grid (Serial)
    int FloodFillForInAndOut(const DefSFBitset& sfbitset_start,
        const DefMap<DefUint>& map_nodes_exist,
        DefMap<DefUint>* const ptr_map_nodes_outside,
        DefMap<DefUint>* const ptr_map_nodes_inside) const;
    virtual void OverlapLayerFromHighToLow(
        const DefMap<DefUint>& layer_high_level,
        DefMap<DefUint>* const ptr_layer_low_level) = 0;
    virtual void GenerateBackgroundGrid(
        const DefMap<DefUint>& map_occupied) = 0;

    virtual bool CheckCoincideBackground(const DefSizet i_level,
        const DefSFBitset& bitset_higher,
        DefSFBitset* const ptr_bitset) const = 0;

 private:
    const DefUint imax_flood_fill = 90000;
    ///<  maximum iteration for flood fill

    void ExtendGivenNumbOfLayer(DefSizet i_level,
        const std::vector<DefLUint> num_extend_neg,
        const std::vector<DefLUint> num_extend_pos,
        const DefMap<DefUint>& map_start_layer,
        DefMap<DefUint>* const ptr_map_exist,
        DefMap<DefUint>* const ptr_map_outmost) const;
    void FindOverlappingLayersBasedOnOutermostCoarse(
        const DefMap<DefUint>& layer_coarse_0,
        const DefMap<DefUint>& sfbitset_exist,
        DefMap<DefUint>* const ptr_layer_coarse_m1,
        DefMap<DefUint>* const ptr_layer_fine_0,
        DefMap<DefUint>* const ptr_layer_fine_m1,
        DefMap<DefUint>* const ptr_layer_fine_m2);
    
    virtual void ResetExtendLayerBasedOnDomainSize(
        const DefSizet i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefLUint>* const ptr_vec_extend_neg,
        std::vector<DefLUint>* const ptr_vec_extend_pos) const = 0;
    virtual void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefSizet i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const = 0;
    virtual void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets) const = 0;
    virtual void IdentifyInterfaceForACell(const DefUint flag_interface,
        const DefSFBitset bitset_in, const DefMap<DefUint>& node_exist,
        DefMap<DefUint>* const ptr_inner_layer,
        DefMap<DefUint>* const ptr_mid_layer,
        DefMap<DefUint>* const ptr_outer_layer) = 0;

    // functions to generate grid (Serial)
    virtual void GenerateGridNodeNearTrackingNode(const DefSizet i_level,
        const std::pair<ECriterionType, DefSizet>& tracking_grid_key,
        DefMap<DefUint>* const ptr_map_node_temp) const = 0;
    virtual void IdentifyTypeOfLayerByFloodFill(
        const DefSizet i_level, const DefSizet i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefUint>& map_nodes_exist,
        DefMap<DefUint>* const ptr_map_nodes_outside,
        DefMap<DefUint>* const ptr_map_nodes_inside) const = 0;
    virtual void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk) const = 0;
    virtual void ExtendOneLayerGrid(
        const DefMap<DefUint>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefUint>* const ptr_map_output_layer,
        DefMap<DefUint>* const ptr_map_exist,
        DefMap<DefUint>* const ptr_outmost_layer,
        std::vector<DefMap<DefUint>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefUint>>* const ptr_vector_boundary_max) const = 0;
    virtual void FindInterfaceBetweenGrid(
        const DefSizet i_level, const DefMap<DefUint>& map_outmost_layer,
        DefMap<DefUint>* const map_exist,
        DefMap<DefUint>* const ptr_interface_outmost,
        DefMap<DefUint>* const ptr_layer_lower_level,
        DefMap<DefUint>* const ptr_layer_lower_level_outer) = 0;

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
    std::array<DefLUint, 2> k0IntOffset_ = {1, 1};
    ///< number of offset background node
    std::array<DefReal, 2> k0DomainSize_{};  ///< domian size
    std::array<DefReal, 2> k0DomainDx_{};  ///< grid space
    /* offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (kXintOffset * kDomianDx),
    and the default value is 1. */
    std::array<DefReal, 2> k0RealOffset_{};  ///< kXintOffset * kDomianDx
    std::array<DefLUint, 2> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0IntOffset_ is included in k0MaxIndexOfBackgroundNode_

    void PrintGridInfo(void) const override;
    void SetGridParameters(void) override;
    SFBitsetAuxInterface* GetSFBitsetAuxPtr() override {
        return dynamic_cast<SFBitsetAuxInterface*>(this);
    }
    std::vector<DefReal> GetDomainDxArrAsVec() const override final {
        return { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex] };
    };

    // node related function
    void FindAllNeighboursSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbours) const override;
    void FindAllNodesInACellAtLowerLevel(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const override;
    DefSFBitset NodeAtNLowerLevel(
        const DefSizet n_level, const DefSFBitset& biset_in) const override;
    DefLUint CalNumOfBackgroundNode() const override {
        return (k0MaxIndexOfBackgroundNode_[0] - k0IntOffset_[0] + 1) *
            (k0MaxIndexOfBackgroundNode_[1] - k0IntOffset_[1] + 1);
    }
    bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const override;
    void TraverseBackgroundForPartition(
        const std::vector<DefLUint>& rank_load,
        const DefMap<DefUint>& background_occupied,
        std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const override;

    GridManager2D() {
        k0GridDims_ = 2;
        k0NumNeighbours_ = 9;
    }

protected:
   void GenerateBackgroundGrid(
        const DefMap<DefUint>& map_occupied) override;
   void OverlapLayerFromHighToLow(
       const DefMap<DefUint>& layer_high_level,
       DefMap<DefUint>* const ptr_layer_low_level) override;

    bool CheckCoincideBackground(const DefSizet i_level,
        const DefSFBitset& bitset_higher,
        DefSFBitset* const ptr_bitset) const override;
    void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefSizet i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const override;
    void ResetExtendLayerBasedOnDomainSize(
        const DefSizet i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefLUint>* const ptr_vec_extend_neg,
        std::vector<DefLUint>* const ptr_vec_extend_pos) const override;
    void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
       std::vector<DefSFBitset>* const ptr_bitsets)  const override;
    void IdentifyInterfaceForACell(const DefUint flag_interface,
        const DefSFBitset bitset_in, const DefMap<DefUint>& node_exist,
        DefMap<DefUint>* const ptr_inner_layer,
        DefMap<DefUint>* const ptr_mid_layer,
        DefMap<DefUint>* const ptr_outer_layer) override;
    bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefUint>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets)override;
    // functions to generate grid (Serial)
    void GenerateGridNodeNearTrackingNode(const DefSizet i_level,
        const std::pair<ECriterionType, DefSizet>& tracking_grid_key,
        DefMap<DefUint>* const ptr_map_node_temp) const override;
    void IdentifyTypeOfLayerByFloodFill(
        const DefSizet i_level, const DefSizet i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefUint>& map_nodes_exist,
        DefMap<DefUint>* const ptr_map_nodes_outside,
        DefMap<DefUint>* const ptr_map_nodes_inside) const override;
    void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk) const override;
    void ExtendOneLayerGrid(
        const DefMap<DefUint>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefUint>* const ptr_map_output_layer,
        DefMap<DefUint>* const ptr_map_exist,
        DefMap<DefUint>* const ptr_outmost_layer,
        std::vector<DefMap<DefUint>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefUint>>* const ptr_vector_boundary_max) 
        const override;
    void FindInterfaceBetweenGrid(
        const DefSizet i_level, const DefMap<DefUint>& map_outmost_layer,
        DefMap<DefUint>* const map_exist,
        DefMap<DefUint>* const ptr_interface_outmost,
        DefMap<DefUint>* const ptr_layer_lower_level,
        DefMap<DefUint>* const ptr_layer_lower_level_outer) override;
private:
    DefUint kFlagCurrentNodeXNeg_ = 1, kFlagCurrentNodeXPos_ = 1 << 1,
        kFlagCurrentNodeYNeg_ = 1 << 2, kFlagCurrentNodeYPos_ = 1 << 3;
    DefUint FindAllNeighboursWithSpecifiedDirection(const DefSFBitset bitset_in,
        std::array<bool, 2>& bool_neg, std::array<bool, 2>& bool_pos,
        std::vector <DefSFBitset>* const ptr_vec_neigbours) const;
    void IdentifyInterfaceNodeOnEdge(const DefUint flag_interface,
        const std::array<DefSFBitset, 2>& arr_bitset_lower,
        const DefSFBitset bitset_mid_higher,
        const DefMap<DefUint>& node_exist,
        const std::array<DefMap<DefUint>* const, 3>& arr_ptr_layer);
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
class GridManager3D :public  GridManagerInterface, public SFBitsetAux3D
 {
 public:
    std::array<DefLUint, 3> k0IntOffset_ = { 1, 1, 1};
    ///< number of offset background node
    std::array<DefReal, 3> k0DomainSize_{};  ///< domian size
    std::array<DefReal, 3> k0DomainDx_{};  ///< grid space
    /* offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (kXintOffset * kDomianDx),
    and the default value is 1. */
    std::array<DefReal, 3>k0RealOffset_{};  ///< kXintOffset * kDomianDx
    std::array<DefLUint, 3> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0IntOffset_ is included in k0MaxIndexOfBackgroundNode_

    void PrintGridInfo(void) const override;
    void SetGridParameters(void) override;

    SFBitsetAuxInterface* GetSFBitsetAuxPtr() override {
        return dynamic_cast<SFBitsetAuxInterface*>(this);
    }
    std::vector<DefReal> GetDomainDxArrAsVec() const override final {
        return { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex],
        k0DomainDx_[kZIndex] };
    };

    // node related function
    void FindAllNeighboursSFBitset(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbours) const override;
    void FindAllNodesInACellAtLowerLevel(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const override;
    DefSFBitset NodeAtNLowerLevel(
        const DefSizet n_level, const DefSFBitset& biset_in) const override;
    DefLUint CalNumOfBackgroundNode() const override {
        return (k0MaxIndexOfBackgroundNode_[0] - k0IntOffset_[0] + 1) *
            (k0MaxIndexOfBackgroundNode_[1] - k0IntOffset_[1] + 1)*
            (k0MaxIndexOfBackgroundNode_[2] - k0IntOffset_[2] + 1);
    }
    bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const override;
    void TraverseBackgroundForPartition(
        const std::vector<DefLUint>& rank_load,
        const DefMap<DefUint>& background_occupied,
        std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max) const override;

    GridManager3D() {
        k0GridDims_ = 3;
        k0NumNeighbours_ = 27;
    }


 protected:
    void GenerateBackgroundGrid(
         const DefMap<DefUint>& map_occupied) override;
    void OverlapLayerFromHighToLow(
        const DefMap<DefUint>& layer_high_level,
        DefMap<DefUint>* const ptr_layer_low_level) override;

    bool CheckCoincideBackground(const DefSizet i_level,
         const DefSFBitset& bitset_higher,
         DefSFBitset* const ptr_bitset) const override;
    void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefSizet i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const override;
    void ResetExtendLayerBasedOnDomainSize(
        const DefSizet i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefLUint>* const ptr_vec_extend_neg,
        std::vector<DefLUint>* const ptr_vec_extend_pos) const override;
    void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets)  const override;
    void IdentifyInterfaceForACell(const DefUint flag_interface,
        const DefSFBitset bitset_in, const DefMap<DefUint>& node_exist,
        DefMap<DefUint>* const ptr_inner_layer,
        DefMap<DefUint>* const ptr_mid_layer,
        DefMap<DefUint>* const ptr_outer_layer) override;
    bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefUint>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) override;
    // functions to generate grid (Serial)
    void GenerateGridNodeNearTrackingNode(const DefSizet i_level,
        const std::pair<ECriterionType, DefSizet>& tracking_grid_key,
        DefMap<DefUint>* const ptr_map_node_temp) const override;
    void IdentifyTypeOfLayerByFloodFill(
        const DefSizet i_level, const DefSizet i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefUint>& map_nodes_exist,
        DefMap<DefUint>* const ptr_map_nodes_outside,
        DefMap<DefUint>* const ptr_map_nodes_inside) const override;
    void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk)const  override;
    void ExtendOneLayerGrid(
        const DefMap<DefUint>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefUint>* const ptr_map_output_layer,
        DefMap<DefUint>* const ptr_map_exist,
        DefMap<DefUint>* const ptr_outmost_layer,
        std::vector<DefMap<DefUint>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefUint>>* const ptr_vector_boundary_max)
        const override;
    void FindInterfaceBetweenGrid(
        const DefSizet i_level, const DefMap<DefUint>& map_outmost_layer,
        DefMap<DefUint>* const map_exist,
        DefMap<DefUint>* const ptr_interface_outmost,
        DefMap<DefUint>* const ptr_layer_lower_level,
        DefMap<DefUint>* const ptr_layer_lower_level_outer) override;
 private:
     DefUint kFlagCurrentNodeXNeg_ = 1, kFlagCurrentNodeXPos_ = 1 << 1,
         kFlagCurrentNodeYNeg_ = 1 << 2, kFlagCurrentNodeYPos_ = 1 << 3,
         kFlagCurrentNodeZNeg_ = 1 << 4, kFlagCurrentNodeZPos_ = 1 << 5;
     DefUint FindAllNeighboursWithSpecifiedDirection(const DefSFBitset bitset_in,
         std::array<bool, 3>& bool_neg, std::array<bool, 3>& bool_pos,
         std::vector <DefSFBitset>* const ptr_vec_neigbours) const;
     void IdentifyInterfaceNodeOnEdge(const DefUint flag_interface,
         const std::array<DefSFBitset, 2>& arr_bitset_lower,
         const DefSFBitset bitset_mid_higher,
         const DefMap<DefUint>& node_exist,
         const std::array<DefMap<DefUint>* const, 3>& arr_ptr_layer);
     void IdentifyInterfaceNodeDiagonal(const DefUint flag_interface,
         const std::array<DefSFBitset, 4>& arr_bitset_lower,
         const DefSFBitset bitset_center_higher,
         const DefMap<DefUint>& node_exist,
         const std::array<DefMap<DefUint>* const, 3>& arr_ptr_layer);
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
*             with a sepecified node type.
* @return  if return 0, then the node type does not exist in the vector 
*          elements; else return the (indice + 1) of the element with
*          the given type.
*/
template<InterfaceInfoHasType InterfaceInfo>
//template<typename Type, typename InterfaceInfo>
DefSizet GridManagerInterface::CheckExistanceOfTypeAtGivenLevel(
    const DefSizet i_level, const std::string& type,
    const std::vector<std::shared_ptr<InterfaceInfo>>&
    vec_ptr_interface) const {
    DefSizet iter_count = 0;
    if (vec_ptr_interface.empty()) {
        return 0;
    }
    else {
        for (const auto& iter : vec_ptr_interface) {
            ++iter_count;
            if ((iter->i_level_ == i_level)
                && (iter->node_type_ == type)) {
                return iter_count;
            }          
        }
    }
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
