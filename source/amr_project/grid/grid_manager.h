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
class GridManagerInterface {
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
    DefUint k0IntExtendMin_ = 2;  ///< a minimum number of extending layers
    DefUint k0IntExtendRemoved_ = 1;
    ///< a number of extending layers when inside nodes will be removed
    /* number of nodes extended from position of given criterion,
    e.g. solid boundary and free surface*/
    std::vector<DefLUint>  k0IntExtend_,  ///< number of extened layers
        k0XIntExtendPositive_, k0XIntExtendNegative_,
        k0YIntExtendPositive_, k0YIntExtendNegative_;
    /* number of layer extended inside the geometry
     at (i_level - 1) refienemnt level*/
    std::vector<DefLUint>  k0IntInnerExtend_;
    ///< number of extened layers inside the geometry

    // grid type
    /*the method used be used for all grid levels to create grid instance.
      Give different methods for each grid is possible by using other functions
      instead of SetGridInfoInterfaceAllLevels(EGridNodeType node_type) during
      initilization in function SetGridParameters(). */
    ///< method chosen to be solved on the grid
    std::vector<std::shared_ptr<GridInfoInterface>> vec_ptr_grid_info_;
    void CreateGridInstance(const DefUint i_level,
        std::shared_ptr<SolverInterface> ptr_sovler);

    void DefaultInitialization(const DefSizet max_level);
    int CheckIfPointOutsideDomain(
        const std::array<DefReal, 2>& coordinate_min,
        const std::array<DefReal, 2>& coordinate_max,
        const std::array<DefReal, 2>& domain_min,
        const std::array<DefReal, 2>& domain_max);
    int CheckIfPointOutsideDomain(
        const std::array<DefReal, 3>& coordinate_min,
        const std::array<DefReal, 3>& coordinate_max,
        const std::array<DefReal, 3>& domain_min,
        const std::array<DefReal, 3>& domain_max);


    virtual void DefaultInitializationDims(const DefSizet max_level) = 0;
    virtual void SetGridParameters(void) = 0;

    // generate grid
    void GenerateGridSerial(const std::vector<std::shared_ptr
         <GeometryInfoInterface>>& vec_geo_info);

    template<InterfaceInfoHasType InterfaceInfo>
    //template<typename Type, typename InterfaceInfo>
    DefSizet CheckExistanceOfTypeAtGivenLevel(
        const DefSizet i_level, const std::string& type,
        const std::vector<std::shared_ptr<InterfaceInfo>>& vec_ptr_interface);

 private:
    //// Generate grid
    const DefUint kFlag0_ = 0;  // flag to initilize status of a node
    // node with criterion at lower level
    const DefUint kFlagBitLowerLevel_ = 1;
    // node inside enclosed shape will be removed
    const DefUint kFlagBitEnclosedRemove_ = 1 << 2;
    // node inside enclosed shape will be extended at given number of layers
    const DefUint kFlagBitEnclosedExtend_ = 1 << 3;

    const DefUint imax_flood_fill = 90000;  
    ///<  maximum iteration for flood fill

    void AddNodesInstanceBasedOnLowerLevel(
        const DefSFBitset& sfbitset_lower_level,
        const DefMap<DefUint>& map_sfbitset_at_lower_level,
        std::shared_ptr<GridInfoInterface> ptr_grid_info);

    //// functions to generate grid (Serial)
    //void SearchingForNodesBasedOnGeometry(
    //    const DefSizet i_level, const DefSizet i_geo,
    //    const std::shared_ptr<criterion::GeometryInfoInterface> ptr_geo_info,
    //    DefMap<DefUint>* const ptr_map_sfbitset_at_lower_level,
    //    std::vector<DefMap<DefUint>>* const ptr_vec_map_tracking_nodes);
    //void FindTrackingNodesBaseOnGeoVertics(const DefSizet i_geo,
    //    const DefSizet i_level, const DefSizet i_geo_level,
    //    const std::vector<criterion::GeometryCoordinate>& vec_geo_coordi,
    //    DefMap<DefUint>* const ptr_map_sfbitset_at_lower_level,
    //    DefMap<DefUint>* const ptr_map_tracking_node_temp);
    //void IdentifyTypeOfLayerByFloodFill(
    //    const DefSizet i_level, const DefSizet i_geo,
    //    const std::shared_ptr<criterion::GeometryInfoInterface> ptr_geo_info,
    //    const DefMap<DefUint>& map_sfbitset,
    //    DefMap<DefUint>* const ptr_map_nodes_uncolored,
    //    DefMap<DefUint>* const ptr_map_nodes_colored);
    //void ClassifyLayersNearGeometry(const DefSizet i_level,
    //    const DefSizet i_geo, const DefUint flag_color,
    //    const std::shared_ptr<criterion::GeometryInfoInterface> ptr_geo_info,
    //    const std::vector<DefReal>& vec_origin,
    //    DefMap<DefUint>* const map_nodes_exist,
    //    DefMap<DefUint>* const ptr_map_nodes_uncolored,
    //    DefMap<DefUint>* const ptr_map_nodes_colored);
    //void FloodFillOneLayer(const DefSFBitset& sfbitset_start,
    //    const DefMap<DefUint>& map_nodes_exist,
    //    DefMap<DefUint>* const ptr_map_nodes_uncolored,
    //    DefMap<DefUint>* const ptr_map_nodes_colored);
    //void GenerateNodeLayerByLayer(const DefSizet i_level,
    //    const std::vector<DefUint>& vec_number_of_extend_layer_neg,
    //    const std::vector<DefUint>& vec_number_of_extend_layer_pos,
    //    const DefMap<DefUint>& map_nodes_starting_layer,
    //    DefMap<DefUint>* const ptr_map_nodes_exist,
    //    DefMap<DefUint>* const ptr_map_two_outmost_layer_at_lower_level);

    //void ExtendSameNumberOfLayer(const DefSizet i_level,
    //    const DefMap<DefUint>& map_tracking_node_temp,
    //    DefMap<DefUint>* const ptr_map_sfbitset_at_lower_level);


    /* begin: functions depend on k0GridDims_, using ptr to call */
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
    void FindTrackingNodesContainingAGeometryCoordinate2D(
            const std::vector<DefReal>& grid_space,
            const std::vector<DefReal>& geo_coordi,
            std::vector<DefSFBitset>* ptr_vec_sfbitset);
    void ResetExtendLayerBasedOnDomainSize2D(
        const DefSizet i_level, const DefSFBitset& bitset_in,
        std::vector<DefLUint>* ptr_vec_extend_neg,
        std::vector<DefLUint>* ptr_vec_extend_pos);
    DefSFBitset FindStartingPointForFloodFill2D(
        const DefSizet i_level, const DefSizet i_geo,
        const std::vector<DefReal>& vec_origin,
        const DefMap<DefUint>& map_nodes_exist);
    void PushBackSFBitsetInFloodFill2D(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk);
    bool AddNodesAroundAGivenNode2D(
        const DefSFBitset& sfbitset_in,
        const std::vector<bool>& vec_bool_need_add_neg,
        const std::vector<bool>& vec_bool_need_add_pos,
        DefMap<DefLUint>* const ptr_map_nodes_exist,
        DefMap<DefLUint>* const ptr_map_nodes_aded);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef DEBUG_DISABLE_3D_FUNCTIONS
    void FindTrackingNodesContainingAGeometryCoordinate3D(
        const std::vector<DefReal>& grid_space,
        const std::vector<DefReal>& geo_coordi,
        std::vector<DefSFBitset>* ptr_vec_sfbitset);
    void ResetExtendLayerBasedOnDomainSize3D(
        const DefSizet i_level, const DefSFBitset& bitset_in,
        std::vector<DefLUint>* ptr_vec_extend_neg,
        std::vector<DefLUint>* ptr_vec_extend_pos);
    DefSFBitset FindStartingPointForFloodFill3D(
        const DefSizet i_level, const DefSizet i_geo,
        const std::vector<DefReal>& vec_origin,
        const DefMap<DefUint>& map_nodes_exist);
    void PushBackSFBitsetInFloodFill3D(const DefSFBitset& bitset_in,
        std::vector<DefSFBitset>* const ptr_vec_code_stk);
    bool AddNodesAroundAGivenNode3D(
        const DefSFBitset& sfbitset_in,
        const std::vector<bool>& vec_bool_need_add_neg,
        const std::vector<bool>& vec_bool_need_add_pos,
        DefMap<DefLUint>* const ptr_map_nodes_exist,
        DefMap<DefLUint>* const ptr_map_nodes_aded);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    /* end:  functions depend on k0GridDims_ */

#ifdef DEBUG_UNIT_TEST
 private:
// gtest to access private member functions
    //FRIEND_TEST(GridManagerGeneration2D, SearchingForNodesBasedOnGeo);
    //FRIEND_TEST(GridManagerGeneration3D, SearchingForNodesBasedOnGeo);
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
    std::array<DefReal, 2> k0IntOffest_{};  ///< number of offset background node
    std::array<DefReal, 2> k0DomainSize_{};  ///< domian size
    std::array<DefReal, 2> k0DomainDx_{};  ///< grid space
    /* offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (kXintOffest * kDomianDx),
    and the default value is 1. */
    std::array<DefReal, 2> k0RealOffest_{};  ///< kXintOffest * kDomianDx
    std::array<DefReal, 2> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0IntOffest_ is included in k0MaxIndexOfBackgroundNode_

    void DefaultInitializationDims(const DefSizet max_level) override;
    void SetGridParameters(void) override;
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GridManager3D :public  GridManagerInterface, public SFBitsetAux3D {
public:
    std::array<DefReal, 3> k0IntOffest_{};  ///< number of offset background node
    std::array<DefReal, 3> k0DomainSize_{};  ///< domian size
    std::array<DefReal, 3> k0DomainDx_{};  ///< grid space
    std::vector<DefLUint> k0ZIntExtendPositive_, k0ZIntExtendNegative_;
    /* offset to avoid exceeding the boundary limits when searching nodes.
    The offset distance is (kXintOffest * kDomianDx),
    and the default value is 1. */
    std::array<DefReal, 3>k0RealOffest_{};  ///< kXintOffest * kDomianDx
    std::array<DefReal, 3> k0MaxIndexOfBackgroundNode_{};
    ///< the maximum index of background nodes in each direction*/
    // k0IntOffest_ is included in k0MaxIndexOfBackgroundNode_

    void DefaultInitializationDims(const DefSizet max_level) override;
    void SetGridParameters(void) override;
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
/**
* @brief   function to check if node type exists at a given level.
* @param[in]  i_level level of refinement.
* @param[in]  type the given type.
* @param[in]  vec_ptr_interface pointer to vector of information
*             with a sepecified node type.
* @return  if return 0, then the node type does not exsit in the vector 
*          elements; else return the (indice + 1) of the element with
*          the given type.
*/
template<InterfaceInfoHasType InterfaceInfo>
//template<typename Type, typename InterfaceInfo>
DefSizet GridManagerInterface::CheckExistanceOfTypeAtGivenLevel(
    const DefSizet i_level, const std::string& type,
    const std::vector<std::shared_ptr<InterfaceInfo>>& vec_ptr_interface) {
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
