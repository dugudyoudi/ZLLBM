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
* @note
*/
struct TrackingNode {
private:
    DefUint count_exist = 1;
    DefUint status = 0;
    std::vector<std::array<DefSizet, 2>> geo_ref;
public:
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
};
/**
* @struct GhostNode
* @brief structure to store ghost node information
* @note
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
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
};
/**
* @class TrackingGridInfoInterface
* @brief class used to store tracking grid information at i_level_ of refinement
* @note
* @date  2022-6-4
*/
class TrackingGridInfoInterface {
public:
    DefReal computational_cost_ = 1.;
    std::string node_type_;

    // information of TrackingNode
    DefMap<TrackingNode> map_tracking_node_{};
    DefSizet k0NumIntForEachNode_ = 1;
    DefSizet k0NumRealForEachNode_ = 0;
protected:
    virtual void FindTrackingNodeNearGeo(const DefSFBitset& bitset_in) = 0;
};
/**
* @class TrackingGridInfoCreatorInterface
* @brief abstract class used to generate TrackingGridInfo instance.
*/
class TrackingGridInfoCreatorInterface {
public:
    virtual std::shared_ptr<TrackingGridInfoInterface>
        CreateTrackingGridInfo() = 0;
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
protected:
    GhostNode ghost_node_instance_;
    ///< an instance for a ghost node with preset vector sizes
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

    std::map<std::pair<ECriteriolType, DefSizet>,
        std::shared_ptr<TrackingGridInfoInterface>>
        map_ptr_tracking_grid_info_;
    std::shared_ptr<GhostGridInfoInterface>
        ptr_ghost_grid_info_;

    // information of GridNode
    DefMap<GridNode> map_grid_node_{};
    DefSizet k0NumIntForEachNode_ = 0;
    DefSizet k0NumRealForEachNode_ = 0;
    virtual void set_number_of_vec_elements() = 0;

    // overlapping region between differen refinement level
    const DefUint kOverlappingOutmostM1_ = 1;
    ///< index of the outmost layer of the overlapping region 
    const DefUint kOverlappingOutmostM2_ = 2;
    ///< index of the second outmost layer of the overlapping region 
    std::vector<DefMap<DefUint>> vec_map_coarse2fine_;
    std::vector<DefMap<DefUint>> vec_map_fine2coarse_;

    GridNode grid_node_instance_;
    ///< an instance for insert node with preset vector sizes
    virtual void InitialGridNode(const DefSFBitset& bitset_in) = 0;
};
/**
* @class GridInfoCreatorInterface
* @brief abstract class used to generate GhostGridInfo instance.
*/
class GridInfoCreatorInterface {
public:
    virtual std::shared_ptr<GridInfoInterface>
        CreateGridInfo() = 0;
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_SOURCE_GRID_GRID_INFO_H_
