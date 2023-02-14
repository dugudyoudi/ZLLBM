//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/grid_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define parent class for nodes
*/

#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_INFO_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_GRID_INFO_H_
#include <vector>
#include <array>
#include <map>
#include <set>
#include <memory>
#include <string>
#include "../defs_libs.h"
#include "criterion/criterion_numerates.h"
#include "./solver_info_interface.h"
namespace rootproject {
namespace amrproject{
/**
* @struct TrackingNode
* @brief structure to nodes information for traking movement
*/
struct TrackingNode {
public:
    std::set<std::pair<DefSizet, DefSizet>> set_point_index;
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
};
/**
* @struct GhostNode
* @brief structure to store ghost node information
*/
struct GhostNode {
public:
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
};
/**
* @struct GridNode
* @brief struct used to store information of a node
* @date  2022-6-4
*/
struct GridNode {
public:
    DefUint flag_status_;
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
};
/**
* @class InterfaceLayerInfo
* @brief class used to store information of
* @date  2023-1-5
*/
class InterfaceLayerInfo {
public:
    // number of extended based on the interface grid
    std::vector<DefLUint> k0ExtendOuterNeg, k0ExtendOuterPos;
    ///< number of extened layers outside the geometry
    std::vector<DefLUint> k0ExtendInnerNeg, k0ExtendInnerPos;
    ///< number of extened layers inside the geometry

    std::vector<DefMap<DefUint>> vec_outer_coarse2fine_,
        vec_outer_fine2coarse_;
    std::vector<DefMap<DefUint>> vec_inner_coarse2fine_,
        vec_inner_fine2coarse_;
};
/**
* @class TrackingGridInfoInterface
* @brief class used to store tracking grid information at i_level_ of refinement
* @date  2022-6-4
*/
class TrackingGridInfoInterface {
public:
    DefReal computational_cost_ = 1.;
    std::string node_type_;
    EGridExtendType grid_extend_type_ = EGridExtendType::kSameInAllDirections;

    // number of extended based on the trakcing grid
    std::vector<DefLUint> k0ExtendOuterNeg, k0ExtendOuterPos;
    ///< number of extened layers outside the geometry
    std::vector<DefLUint> k0ExtendInnerNeg, k0ExtendInnerPos;
    ///< number of extened layers inside the geometry

    // information of TrackingNode
    DefMap<TrackingNode> map_tracking_node_{};
    DefSizet k0NumIntForEachNode_ = 1;
    DefSizet k0NumRealForEachNode_ = 0;
    TrackingNode k0TrackNodeInstance_;
};
/**
* @class TrackingGridInfoCreatorInterface
* @brief abstract class used to generate TrackingGridInfo instance.
*/
class TrackingGridInfoCreatorInterface {
public:
    virtual std::shared_ptr<TrackingGridInfoInterface>
        CreateTrackingGridInfo() = 0;
    virtual ~TrackingGridInfoCreatorInterface() {}
};
/**
* @class GhostGridInfoInterface
* @brief class used to store ghost grid information at i_level_ of refinement
* @note
* @date  2022-6-4
*/
class GhostGridInfoInterface {
public:
    // information of grid at each level of refinement
    DefSizet i_level_ = 0;
    const DefUint kCountIndex_ = 1;
    DefReal computational_cost_ = 1.;
    std::string node_type_;

    // information of GhostNode
    DefMap<GhostNode> map_ghost_node_{};
    DefSizet k0NumIntForEachNode_ = 1;
    DefSizet k0NumRealForEachNode_ = 0;
    GridNode k0GhostNodeInstance_;
    ///< instance for a ghost node with preset vector sizes
    virtual ~GhostGridInfoInterface() {}
protected:
    virtual void InitialGhostNode(const DefSFBitset& bitset_in) = 0;
};
/**
* @class GhostGridInfoCreatorInterface
* @brief abstract class used to generate GhostGridInfo instance.
*/
class GhostGridInfoCreatorInterface {
public:
    virtual std::shared_ptr<GhostGridInfoInterface>
        CreateGhostGridInfo() = 0;
    virtual ~GhostGridInfoCreatorInterface() {}
};
/**
* @class GridInfoInterface
* @brief class used to store grid information at i_level_ of refinement
* @note
* @date  2022-6-4
*/
class GridInfoInterface {
public:
    // information of grid at each level of refinement
    DefSizet i_level_ = 0;

    DefReal computational_cost_ = 1.;
    std::string node_type_;
    std::vector<DefReal> grid_space_;
    std::shared_ptr<SolverInterface> ptr_solver_ = nullptr;

    std::map<std::pair<ECriterionType, DefSizet>,
        std::shared_ptr<TrackingGridInfoInterface>>
        map_ptr_tracking_grid_info_;
    std::shared_ptr<GhostGridInfoInterface>
        ptr_ghost_grid_info_;

    // interface between grid of differet refinement levels
    DefSizet k0NumFine2CoarseLayer_ = 3;
    DefSizet k0NumCoarse2FineLayer_ = 2;
    std::map<std::pair<ECriterionType, DefSizet>,
        std::shared_ptr<InterfaceLayerInfo>>
        map_ptr_interface_layer_info_;

    // information of GridNode
    DefMap<GridNode> map_grid_node_{};
    DefMap<DefUint> map_grid_count_exist_{};
    DefSizet k0NumIntForEachNode_ = 0;
    DefSizet k0NumRealForEachNode_ = 0;
    virtual void set_number_of_vec_elements() = 0;

    GridNode k0GridNodeInstance_;
    ///< instance for insert node with preset vector sizes
    virtual void InitialGridNode(const DefSFBitset& bitset_in) = 0;

    virtual ~GridInfoInterface() {}
};
/**
* @class GridInfoCreatorInterface
* @brief abstract class used to generate GhostGridInfo instance.
*/
class GridInfoCreatorInterface {
public:
    virtual std::shared_ptr<GridInfoInterface>
        CreateGridInfo() = 0;
    virtual ~GridInfoCreatorInterface() {}
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_SOURCE_GRID_GRID_INFO_H_
