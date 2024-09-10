//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief class to manage grid at all refinement levels.
*/

#ifndef SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
#define SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
#include <concepts>
#include <array>
#include <set>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include <variant>
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
enum class EDomainBoundaryType{
    kCubic = 0,
};
template<typename InterfaceInfo>
concept InterfaceInfoHasType = requires(
    InterfaceInfo type_interface) {
    type_interface.node_type_;
};
/**
* @class TimeSteppingSchemeInterface
* @brief class used to manage time stepping scheme.
*/
class TimeSteppingSchemeInterface {
 public:
    std::vector<DefInt> k0TimeSteppingOrder_;
    /**< a vector to store the order of grid info will be called in one background time step*/
    virtual DefReal GetCurrentTimeStep(const DefInt i_level,
        const DefAmrLUint time_step_background, const DefInt time_step_level) const = 0;
    TimeSteppingSchemeInterface() {}
    explicit TimeSteppingSchemeInterface(const DefInt max_level) {}
    virtual ~TimeSteppingSchemeInterface() {}
};
/**
* @class TimeSteppingSchemeInterface
* @brief class used to manage multi time stepping scheme from coarse to fine level.
*/
class MultiTimeSteppingC2F : public TimeSteppingSchemeInterface {
 public:
    DefReal GetCurrentTimeStep(const DefInt i_level,
        const DefAmrLUint time_step_background, const DefInt time_step_level) const override;
    explicit MultiTimeSteppingC2F(const DefInt max_level);
};
/**
* @struct NodeBitStatus
* @brief struct declaration for node status (constants).
*/
struct NodeBitStatus {
    inline static DefInt kNodeStatus0_ = 0;
    inline static DefInt kNodeStatusCoarse2Fine0_ = 1 << 1;  ///< flag indicating node on the outmost fine layer
    inline static DefInt kNodeStatusCoarse2FineGhost_ = 1 << 2;
    inline static DefInt kNodeStatusFine2Coarse0_ = 1 << 3;  ///< flag indicating node on outmost coarse layer
    inline static DefInt kNodeStatusFine2CoarseGhost_ = 1 << 4;
    inline static DefInt kNodeStatusMpiPartitionOuter_ = 1 << 5;
    inline static DefInt kNodeStatusMpiPartitionInner_ = 1 << 6;
};
/**
* @class GridManagerInterface
* @brief class used to manage adaptive mesh refinement.
* @note parameters with prefix 'k0'will be considered
*       as constants after project initialization.
*/
class GridManagerInterface{
 public:
    /**<Constant: index corresponding to x, y, and z in vectors */
    // 2 or 3
    DefInt  k0GridDims_ = 0;  ///< dimension
    // 9 for 2D, 27 for 3D
    DefInt  k0NumNeighbors_ = 0;  ///< number of neighboring nodes
    // 0 refers to uniform mesh
    DefInt  k0MaxLevel_ = 0;  ///< maximum levels of refinement

    /* this setup is used to avoid issues caused by rapid change of
    refinement level in a small region. Number of layers is counted
    at the lower refinement level.*/
    DefInt k0IntExtendMin_ = 3;  ///< a minimum number of extending layers
    DefInt k0IntExtendRemoved_ = 1;
    ///< a number of extending layers when inside nodes will be removed
    /* number of nodes extended from position of given criterion,
    e.g. solid boundary and free surface*/

    // grid status
    const DefInt kFlagSize0_ = 0;  // flag initialize size as 0

    //  o     o     o     o  // coarse grid
    //
    //  o  x  o  x  o  x  o  // Coarse2Fine         Fine2CoarseOutmost
    //  x  x  x  x  x  x  x  //                     Fine2CoarseGhost
    //  o  x  o  x  o  x  o  // Coarse2FineOutmost  Fine2Coarse
    //  x  x  x  x  x  x  x  // find grid

    // computational domain related information
    DefSFBitset k0SFBitsetDomainMin_ = 0;
    DefSFBitset k0SFBitsetDomainMax_ = ~0;

    // creators to instantiate classes
    std::vector<std::unique_ptr<TrackingGridInfoCreatorInterface>> vec_ptr_tracking_info_creator_;

    std::vector<std::shared_ptr<SolverInterface>> vec_ptr_solver_;
    std::vector<std::shared_ptr<GridInfoInterface>> vec_ptr_grid_info_;

    // manager related functions
    void DefaultInitialization(const DefInt max_level);
    virtual void PrintGridInfo(void) const = 0;
    virtual void SetGridParameters(void) = 0;
    template<InterfaceInfoHasType InterfaceInfo>
    DefSizet CheckExistenceOfTypeAtGivenLevel(
        const DefInt i_level, const std::string& type,
        const std::vector<std::shared_ptr<InterfaceInfo>>& vec_ptr_interface) const;

    // get mesh parameters for 2D and 3D
    virtual SFBitsetAuxInterface* GetSFBitsetAuxPtr() = 0;
    virtual std::vector<DefReal> GetDomainDxArrAsVec() const = 0;
    virtual std::vector<DefAmrLUint> GetMinIndexOfBackgroundNodeArrAsVec() const = 0;
    virtual std::vector<DefAmrLUint> GetMaxIndexOfBackgroundNodeArrAsVec() const = 0;
    virtual void CalDomainBoundsAtGivenLevel(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_domain_min, std::vector<DefSFBitset>* const ptr_domain_max) const = 0;

    // setup dimension related members
    virtual void SetDomainSize(const std::vector<DefReal>& domain_size) = 0;
    virtual void SetDomainGridSize(const std::vector<DefReal>& domain_grid_size) = 0;

    // function to check nodes should be on domain boundaries
    void (GridManagerInterface::*ptr_func_insert_domain_boundary_)(const int flag_node,
        const DefSFBitset& bitset_in, GridInfoInterface* const ptr_grid_info) const = nullptr;
    EDomainBoundaryType domain_boundary_type_ = EDomainBoundaryType::kCubic;
    void InsertCubicDomainBoundary(const int flag_node, const DefSFBitset& bitset_in,
        GridInfoInterface* const ptr_grid_info) const;

    // node related functions
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
    bool InstantiateGridNode(const DefSFBitset& bitset_in,
        GridInfoInterface* const ptr_grid_info);
    virtual int CheckNodeOnDomainBoundary(const DefInt i_level, const DefSFBitset& sfbitset_in) const = 0;
    virtual bool CheckNodeNotOutsideDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::vector<DefSFCodeToUint>& sfbitset_min,
        const std::vector<DefSFCodeToUint>& sfbitset_max) const = 0;
    virtual bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) const = 0;
    virtual bool NodesBelongToOneSurfAtHigherLevel(const DefSFBitset bitset_in,
        const DefInt dir_norm,
        const DefMap<DefInt>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) const = 0;
    virtual void GridFindAllNeighborsVir(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const = 0;
    virtual void FindAllNodesInACellAtOneLevelLower(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const = 0;
    virtual DefSFBitset NodeAtNLowerLevel(
        const DefInt n_level, const DefSFBitset& sfbitset_in) const = 0;
    virtual DefAmrLUint CalNumOfBackgroundNode() const = 0;
    virtual bool CheckBackgroundOffset(const DefSFBitset& sfbitset_in) const = 0;

    // generate grid (initialization)
    /*the method is used to create grid instance for all refinement levels.
      Give different methods for each grid is possible by using other functions*/
    void CreateSameGridInstanceForAllLevel(const std::shared_ptr<SolverInterface>& ptr_solver,
        const GridInfoCreatorInterface& grid_creator);
    void CreateTrackingGridInstanceForAGeo(const DefInt i_geo, const GeometryInfoInterface& geo_info);
    void SetNumberOfExtendLayerForGrid(const DefInt i_level,
        const GeometryInfoInterface& geo_info,
        std::vector<DefAmrLUint>* const num_extend_inner_layer_neg,
        std::vector<DefAmrLUint>* const num_extend_inner_layer_pos,
        std::vector<DefAmrLUint>* const num_extend_outer_layer_neg,
        std::vector<DefAmrLUint>* const num_extend_outer_layer_pos);
    void GenerateGridFromHighToLowLevelSerial(const std::vector<std::shared_ptr
        <GeometryInfoInterface>>&vec_geo_info,
        std::vector<DefMap<DefInt>>* const ptr_sfbitset_one_lower_level);
    void InstantiateGridNodeAllLevel(const DefSFBitset sfbitset_min, const DefSFBitset sfbitset_max,
        const std::vector<DefMap<DefInt>>& sfbitset_one_lower_level);

    virtual ~GridManagerInterface() {}
    GridManagerInterface() {
        static_assert(sizeof(DefSFBitset) == sizeof(DefSFCodeToUint),
            "size of DefSFBitset and DefSFCodeToU must be the same");
        this->ptr_func_insert_domain_boundary_ = &GridManagerInterface::InsertCubicDomainBoundary;
    }

 protected:
    //// Generate grid
    const DefInt kFlag0_ = 0;  // flag to initialize status of a node

    // functions to generate grid
    void InstantiateOverlapLayerOfRefinementInterface(
        const std::vector<DefMap<DefInt>>& sfbitset_one_lower_level);
    void InstantiatePeriodicNodes(const bool bool_min, const DefInt i_dim, const DefInt i_level,
        const DefInt num_ghost_layer, const int boundary_min, const DefSFBitset sfbitset_in,
        const DefSFCodeToUint code_min_background_level, const DefSFCodeToUint code_max_background_level,
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_outer_layer,
        DefMap<DefInt>* const ptr_boundary_corner_y_min,
        DefMap<DefInt>* const ptr_boundary_corner_y_max,
        DefMap<DefInt>* const ptr_boundary_corner_z_min,
        DefMap<DefInt>* const ptr_boundary_corner_z_max);
    int MarkRefinementInterface(const DefInt i_level, const DefMap<DefInt>& sfbitset_one_lower_level,
        const DefMap<DefInt>& sfbitset_extra);
    void IdentifyInterfaceNodeOnEdgeInnermost(
        const std::array<DefSFBitset, 2>& arr_bitset_lower,
        const DefSFBitset bitset_mid_higher,
        const SFBitsetAuxInterface& sfbitset_aux,
        const DefMap<DefInt>& node_coarse_interface_innermost,
        const DefMap<std::unique_ptr<GridNode>>& node_current,
        const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer,
        DefMap<DefInt>* const ptr_node_coarse_interface_outer);
    void IdentifyInterfaceNodeOnEdge(
        const std::array<DefSFBitset, 2>& arr_bitset_lower,
        const DefSFBitset bitset_mid_higher,
        const SFBitsetAuxInterface& sfbitset_aux,
        const DefMap<DefInt>& node_outmost,
        const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer);
    void IdentifyInterfaceNodeOnEdge(
        const std::array<DefSFBitset, 2>& arr_bitset_lower,
        const DefSFBitset bitset_mid_higher,
        const SFBitsetAuxInterface& sfbitset_aux,
        const DefMap<DefInt>& node_coarse_interface_previous,
        const DefMap<DefInt>& node_coarse_interface_inter,
        const DefMap<std::unique_ptr<GridNode>>& node_current,
        const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer,
        DefMap<DefInt>* const ptr_coarse_outer);
    virtual void OverlapLayerFromHighToLow(
        const DefMap<DefInt>& layer_high_level,
        DefMap<DefInt>* const ptr_layer_low_level) = 0;
    virtual void InstantiateBackgroundGrid(const DefSFCodeToUint code_min,
    const DefSFCodeToUint code_max, const DefMap<DefInt>& map_occupied) = 0;
    virtual bool CheckCoincideBackground(const DefInt i_level,
        const DefSFBitset& bitset_higher, DefSFBitset* const ptr_bitset) const = 0;
    // function FloodFillForInAndOut only runs on rank 0
    int FloodFillForInAndOut(const DefSFBitset& sfbitset_start,
        const DefMap<DefInt>& map_nodes_exist,
        DefMap<DefInt>* const ptr_map_nodes_outside,
        DefMap<DefInt>* const ptr_map_nodes_inside) const;

 private:
    const DefInt K0IMaxFloodFill_ = static_cast<DefInt>(90000);
    ///<  maximum iteration for flood fill

    void ExtendGivenNumbOfLayer(DefInt i_level,
        const std::vector<DefAmrLUint> num_extend_neg,
        const std::vector<DefAmrLUint> num_extend_pos,
        const DefMap<DefInt>& map_start_layer,
        DefMap<DefInt>* const ptr_map_exist,
        DefMap<DefInt>* const ptr_map_outmost) const;
    void FindOverlappingLayersBasedOnOutermostCoarse(
        const DefMap<DefInt>& layer_coarse_0,
        const DefMap<DefInt>& sfbitset_exist,
        DefMap<DefInt>* const ptr_layer_coarse_m1,
        DefMap<DefInt>* const ptr_layer_fine_0,
        DefMap<DefInt>* const ptr_layer_fine_m1,
        DefMap<DefInt>* const ptr_layer_fine_m2);

    virtual void ResetExtendLayerBasedOnDomainSize(
        const DefInt i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefAmrLUint>* const ptr_vec_extend_neg,
        std::vector<DefAmrLUint>* const ptr_vec_extend_pos) const = 0;
    virtual void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const = 0;
    virtual void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets) const = 0;
    virtual void IdentifyInnermostInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
        const DefMap<DefInt>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_node_coarse_interface_outer) = 0;
    virtual void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface,
        const DefMap<DefInt>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer) = 0;
    virtual void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface_previous,
        const DefMap<DefInt>& node_coarse_interface_inner,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_outer_coarse_layer) = 0;

    // functions to generate grid (Serial)
    virtual void GenerateGridNodeNearTrackingNode(const DefInt i_level,
        const std::pair<ECriterionType, DefInt>& tracking_grid_key,
        DefMap<DefInt>* const ptr_map_node_tmp) const = 0;
    virtual void IdentifyTypeOfLayerByFloodFill(
        const DefInt i_level, const DefInt i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefInt>& map_nodes_exist,
        DefMap<DefInt>* const ptr_map_nodes_outside,
        DefMap<DefInt>* const ptr_map_nodes_inside) const = 0;
    virtual void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk) const = 0;
    virtual void ExtendOneLayerGrid(
        const DefMap<DefInt>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefInt>* const ptr_map_output_layer,
        DefMap<DefInt>* const ptr_map_exist,
        DefMap<DefInt>* const ptr_outmost_layer,
        std::vector<DefMap<DefInt>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefInt>>* const ptr_vector_boundary_max) const = 0;
    virtual void FindOutmostLayerForFineGrid(
        const DefInt i_level, const DefMap<DefInt>& map_outmost_layer,
        DefMap<DefInt>* const map_exist,
        DefMap<DefInt>* const ptr_interface_outmost,
        DefMap<DefInt>* const ptr_layer_lower_level,
        DefMap<DefInt>* const ptr_layer_lower_level_outer) = 0;

    // used for mpi only
 public:
    void InstantiateGridNodeAllLevelMpi(const int i_rank,
        const DefInt num_partition_inner_layer, const DefInt num_partition_outer_layer,
        const std::vector<DefSFCodeToUint>& vec_code_min, const std::vector<DefSFCodeToUint>& vec_code_max,
        const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefMap<DefInt>>& sfbitset_one_lower_level,
        const std::vector<DefMap<DefInt>>& sfbitset_ghost_one_lower_level,
        const DefMap<DefInt>& sfbitset_partition_interface,
        std::vector<std::map<int, DefMap<DefInt>>>* const ptr_mpi_inner_layer,
        std::vector<DefMap<DefInt>>* const ptr_mpi_outer_layer);
    void SetUpRanksSentForMpiInnerLayer(int i_rank,
        const DefSFCodeToUint interface_code_background_level, const std::vector<DefSFBitset>& vec_nodes_in_region,
        const std::vector<DefSFCodeToUint>& ull_max, const DefMap<DefInt>& mpi_inner_layer_tmp,
        std::map<int, DefMap<DefInt>>* const ptr_mpi_inner_layer);
    void SetUpMpiInnerLayerForHigherLevel(int i_rank,
        const DefInt i_level, const DefAmrLUint num_partition_outer_layer,
        const DefSFCodeToUint code_min_background_level, const DefSFCodeToUint code_max_background_level,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const std::vector<DefSFBitset>& domain_min_n_level, const std::vector<DefSFBitset>& domain_max_n_level,
        const std::vector<DefSFCodeToUint>& ull_max, const DefMap<DefInt>& partition_interface_level,
        const DefMap<DefInt>& mpi_outer_layer,
        const DefMap<std::unique_ptr<GridNode>>& map_grid_node_higher_level,
        const SFBitsetAuxInterface& sfbitset_aux, std::map<int, DefMap<DefInt>>* const ptr_mpi_inner_layer,
        DefMap<std::unique_ptr<GridNode>> *const ptr_map_grid_node);

 protected:
    void InstantiateOverlapLayerOfRefinementInterfaceMpi(
        const DefInt i_level, const DefInt num_partition_outer_layer,
        const DefInt flag_outmost_refinement_n_outer_mpi,
        const DefSFCodeToUint code_min, const DefSFCodeToUint code_max,
        const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefInt>& map_sfbitset_one_lower_level,
        const DefMap<DefInt>& sfbitset_partition_interface_background,
        DefMap<DefInt>* const ptr_outer_layer_current_level,
        DefMap<DefInt>* const ptr_outer_layer_lower_level);
    void InstantiateDomainBoundaryForMpi(const DefAmrLUint num_partition_outer_layer,
        const DefSFCodeToUint code_min_background_level, const DefSFCodeToUint code_max_background_level,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const DefMap<DefInt>& sfbitset_one_lower_level,
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_outer_layer);

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
    void PrintGridInfo(void) const override;
    void SetGridParameters(void) override;
    SFBitsetAuxInterface* GetSFBitsetAuxPtr() override {
        return dynamic_cast<SFBitsetAuxInterface*>(this);
    }
    std::vector<DefReal> GetDomainDxArrAsVec() const final {
        return { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex] };
    };
    std::vector<DefAmrLUint> GetMinIndexOfBackgroundNodeArrAsVec() const final {
        return { k0MinIndexOfBackgroundNode_[kXIndex], k0MinIndexOfBackgroundNode_[kYIndex] };
    };
    std::vector<DefAmrLUint> GetMaxIndexOfBackgroundNodeArrAsVec() const final {
        return { k0MaxIndexOfBackgroundNode_[kXIndex], k0MaxIndexOfBackgroundNode_[kYIndex] };
    };
    void CalDomainBoundsAtGivenLevel(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_domain_min, std::vector<DefSFBitset>* const ptr_domain_max) const override;

    // setup dimension related members
    void SetDomainSize(const std::vector<DefReal>& domain_size) override;
    void SetDomainGridSize(const std::vector<DefReal>& domain_grid_size) override;

    // node related function
    int CheckNodeOnDomainBoundary(const DefInt i_level, const DefSFBitset& sfbitset_in) const override;
    bool CheckNodeNotOutsideDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::vector<DefSFCodeToUint>& sfbitset_min,
        const std::vector<DefSFCodeToUint>& sfbitset_max) const override;
    bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) const override;
    bool NodesBelongToOneSurfAtHigherLevel(const DefSFBitset bitset_in,
        const DefInt dir_norm,
        const DefMap<DefInt>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) const override;
    void GridFindAllNeighborsVir(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const override;
    void FindAllNodesInACellAtOneLevelLower(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const override;
    DefSFBitset NodeAtNLowerLevel(
        const DefInt n_level, const DefSFBitset& biset_in) const override;
    DefAmrLUint CalNumOfBackgroundNode() const override {
        return (k0MaxIndexOfBackgroundNode_[0] - k0MinIndexOfBackgroundNode_[0] + 1) *
            (k0MaxIndexOfBackgroundNode_[1] - k0MinIndexOfBackgroundNode_[1] + 1);
    }
    bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const override;

    GridManager2D() {
        k0GridDims_ = 2;
        k0NumNeighbors_ = 9;
    }

 protected:
    void InstantiateBackgroundGrid(const DefSFCodeToUint code_min,
    const DefSFCodeToUint code_max, const DefMap<DefInt>& map_occupied) override;
    void OverlapLayerFromHighToLow(
        const DefMap<DefInt>& layer_high_level,
        DefMap<DefInt>* const ptr_layer_low_level) override;

    bool CheckCoincideBackground(const DefInt i_level,
        const DefSFBitset& bitset_higher,
        DefSFBitset* const ptr_bitset) const override;
    void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const override;
    void ResetExtendLayerBasedOnDomainSize(
        const DefInt i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefAmrLUint>* const ptr_vec_extend_neg,
        std::vector<DefAmrLUint>* const ptr_vec_extend_pos) const override;
    void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets)  const override;
    void IdentifyInnermostInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
        const DefMap<DefInt>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_node_coarse_interface_outer) override;
    void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface,
        const DefMap<DefInt>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer) override;
    void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface_previous,
        const DefMap<DefInt>& node_coarse_interface_inner,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_outer_coarse_layer) override;
    // functions to generate grid (Serial)
    void GenerateGridNodeNearTrackingNode(const DefInt i_level,
        const std::pair<ECriterionType, DefInt>& tracking_grid_key,
        DefMap<DefInt>* const ptr_map_node_tmp) const override;
    void IdentifyTypeOfLayerByFloodFill(
        const DefInt i_level, const DefInt i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefInt>& map_nodes_exist,
        DefMap<DefInt>* const ptr_map_nodes_outside,
        DefMap<DefInt>* const ptr_map_nodes_inside) const override;
    void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk) const override;
    void ExtendOneLayerGrid(
        const DefMap<DefInt>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefInt>* const ptr_map_output_layer,
        DefMap<DefInt>* const ptr_map_exist,
        DefMap<DefInt>* const ptr_outmost_layer,
        std::vector<DefMap<DefInt>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefInt>>* const ptr_vector_boundary_max)
        const override;
    void FindOutmostLayerForFineGrid(
        const DefInt i_level, const DefMap<DefInt>& map_outmost_layer,
        DefMap<DefInt>* const map_exist,
        DefMap<DefInt>* const ptr_interface_outmost,
        DefMap<DefInt>* const ptr_layer_lower_level,
        DefMap<DefInt>* const ptr_layer_lower_level_outer) override;

 private:
    DefInt kFlagCurrentNodeXNeg_ = 1, kFlagCurrentNodeXPos_ = 1 << 1,
        kFlagCurrentNodeYNeg_ = 1 << 2, kFlagCurrentNodeYPos_ = 1 << 3;
    DefInt FindAllNeighborsWithSpecifiedDirection(const DefSFBitset bitset_in,
        const std::array<bool, 2>& bool_neg, const std::array<bool, 2>& bool_pos,
        std::vector <DefSFBitset>* const ptr_vec_neighbors) const;

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
    void PrintGridInfo(void) const override;
    void SetGridParameters(void) override;

    SFBitsetAuxInterface* GetSFBitsetAuxPtr() override {
        return dynamic_cast<SFBitsetAuxInterface*>(this);
    }
    std::vector<DefReal> GetDomainDxArrAsVec() const final {
        return { k0DomainDx_[kXIndex], k0DomainDx_[kYIndex], k0DomainDx_[kZIndex] };
    };
    std::vector<DefAmrLUint> GetMinIndexOfBackgroundNodeArrAsVec() const final {
        return { k0MinIndexOfBackgroundNode_[kXIndex], k0MinIndexOfBackgroundNode_[kYIndex],
         k0MinIndexOfBackgroundNode_[kZIndex] };
    };
    std::vector<DefAmrLUint> GetMaxIndexOfBackgroundNodeArrAsVec() const final {
        return { k0MaxIndexOfBackgroundNode_[kXIndex], k0MaxIndexOfBackgroundNode_[kYIndex],
         k0MaxIndexOfBackgroundNode_[kZIndex] };
    };
    void CalDomainBoundsAtGivenLevel(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_domain_min, std::vector<DefSFBitset>* const ptr_domain_max) const override;

    // setup dimension related members
    void SetDomainSize(const std::vector<DefReal>& domain_size) override;
    void SetDomainGridSize(const std::vector<DefReal>& domain_grid_size) override;

    // node related function
    int CheckNodeOnDomainBoundary(const DefInt i_level, const DefSFBitset& sfbitset_in) const override;
    bool CheckNodeNotOutsideDomainBoundary(const DefSFBitset& sfbitset_in,
        const std::vector<DefSFCodeToUint>& sfbitset_min,
        const std::vector<DefSFCodeToUint>& sfbitset_max) const override;
    bool NodesBelongToOneCell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) const override;
    bool NodesBelongToOneSurfAtHigherLevel(const DefSFBitset bitset_in,
        const DefInt dir_norm,
        const DefMap<DefInt>& node_exist,
        std::vector<DefSFBitset>* const ptr_bitsets) const override;
    void GridFindAllNeighborsVir(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_neighbors) const override;
    void FindAllNodesInACellAtOneLevelLower(
        const std::vector<DefSFBitset> bitset_cell,
        std::vector<DefSFBitset>* const ptr_bitset_all) const override;
    DefSFBitset NodeAtNLowerLevel(
        const DefInt n_level, const DefSFBitset& sfbitset_in) const override;
    DefAmrLUint CalNumOfBackgroundNode() const override {
        return (k0MaxIndexOfBackgroundNode_[0] - k0MinIndexOfBackgroundNode_[0] + 1)
         * (k0MaxIndexOfBackgroundNode_[1] - k0MinIndexOfBackgroundNode_[1] + 1)
         * (k0MaxIndexOfBackgroundNode_[2] - k0MinIndexOfBackgroundNode_[2] + 1);
    }
    bool CheckBackgroundOffset(const DefSFBitset& bitset_in) const override;

    GridManager3D() {
        k0GridDims_ = 3;
        k0NumNeighbors_ = 27;
    }

 protected:
    void InstantiateBackgroundGrid(const DefSFCodeToUint code_min,
        const DefSFCodeToUint code_max, const DefMap<DefInt>& map_occupied) override;
    void OverlapLayerFromHighToLow(
        const DefMap<DefInt>& layer_high_level,
        DefMap<DefInt>* const ptr_layer_low_level) override;
    bool CheckCoincideBackground(const DefInt i_level,
         const DefSFBitset& bitset_higher,
         DefSFBitset* const ptr_bitset) const override;
    void ComputeSFBitsetOnBoundaryAtGivenLevel(const DefInt i_level,
        std::vector<DefSFBitset>* const ptr_vec_bitset_min,
        std::vector<DefSFBitset>* const ptr_vec_bitset_max) const override;
    void ResetExtendLayerBasedOnDomainSize(
        const DefInt i_level, const DefSFBitset& sfbitset_in,
        std::vector<DefAmrLUint>* const ptr_vec_extend_neg,
        std::vector<DefAmrLUint>* const ptr_vec_extend_pos) const override;
    void FindCornersForNeighbourCells(const DefSFBitset bitset_in,
        std::vector<DefSFBitset>* const ptr_bitsets) const override;
    void IdentifyInnermostInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
        const DefMap<DefInt>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer,
        DefMap<DefInt>* const ptr_node_coarse_interface_outer) override;
    void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface,
        const DefMap<DefInt>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer) override;
    void IdentifyInterfaceForACell(const DefSFBitset bitset_in,
        const DefMap<DefInt>& node_coarse_interface_previous,
        const DefMap<DefInt>& node_coarse_interface_inner,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_current,
        const DefMap<std::unique_ptr<GridNode>>& node_exist_lower,
        DefMap<DefInt>* const ptr_inner_layer, DefMap<DefInt>* const ptr_mid_layer,
        DefMap<DefInt>* const ptr_outer_layer, DefMap<DefInt>* const ptr_outer_coarse_layer) override;
    // functions to generate grid (Serial)
    void GenerateGridNodeNearTrackingNode(const DefInt i_level,
        const std::pair<ECriterionType, DefInt>& tracking_grid_key,
        DefMap<DefInt>* const ptr_map_node_tmp) const override;
    void IdentifyTypeOfLayerByFloodFill(
        const DefInt i_level, const DefInt i_geo,
        const std::vector<DefReal> flood_fill_start_point,
        const DefMap<DefInt>& map_nodes_exist,
        DefMap<DefInt>* const ptr_map_nodes_outside,
        DefMap<DefInt>* const ptr_map_nodes_inside) const override;
    void PushBackSFBitsetInFloodFill(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk)const  override;
    void ExtendOneLayerGrid(
        const DefMap<DefInt>& map_start_layer,
        const std::vector<DefSFBitset>& vec_bitset_min,
        const std::vector<DefSFBitset>& vec_bitset_max,
        const std::vector<bool>& bool_extend_neg,
        const std::vector<bool>& bool_extend_pos,
        DefMap<DefInt>* const ptr_map_output_layer,
        DefMap<DefInt>* const ptr_map_exist,
        DefMap<DefInt>* const ptr_outmost_layer,
        std::vector<DefMap<DefInt>>* const ptr_vector_boundary_min,
        std::vector<DefMap<DefInt>>* const ptr_vector_boundary_max)
        const override;
    void FindOutmostLayerForFineGrid(
        const DefInt i_level, const DefMap<DefInt>& map_outmost_layer,
        DefMap<DefInt>* const map_exist,
        DefMap<DefInt>* const ptr_interface_outmost,
        DefMap<DefInt>* const ptr_layer_lower_level,
        DefMap<DefInt>* const ptr_layer_lower_level_outer) override;

 private:
    DefInt kFlagCurrentNodeXNeg_ = 1, kFlagCurrentNodeXPos_ = 1 << 1,
        kFlagCurrentNodeYNeg_ = 1 << 2, kFlagCurrentNodeYPos_ = 1 << 3,
        kFlagCurrentNodeZNeg_ = 1 << 4, kFlagCurrentNodeZPos_ = 1 << 5;
    DefInt FindAllNeighborsWithSpecifiedDirection(const DefSFBitset bitset_in,
        const std::array<bool, 3>& bool_neg, const std::array<bool, 3>& bool_pos,
        std::vector <DefSFBitset>* const ptr_vec_neigbours) const;
    void IdentifyInterfaceNodeDiagonal(
        const std::array<DefSFBitset, 4>& arr_bitset_lower,
        const DefSFBitset bitset_center_higher,
        const DefMap<DefInt>& node_coarse_interface,
        const std::array<DefMap<DefInt>* const, 3>& arr_ptr_layer);

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
    const DefInt i_level, const std::string& type,
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
#endif  // SOURCE_AMR_PROJECT_GRID_GRID_MANAGER_H_
