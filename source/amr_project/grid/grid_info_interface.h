//  Copyright (c) 2021 - 2023, Zhengliang Liu
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
#include <utility>
#include <string>
#include "../defs_libs.h"
#include "criterion/criterion_numerates.h"
#include "./solver_info_interface.h"
#include "grid/sfbitset_aux.h"
namespace rootproject {
namespace amrproject {
class OutputDataFormat;
class Base64Utility;
/**
* @struct TrackingNode
* @brief structure to nodes information for tracking movement
*/
struct TrackingNode {
 public:
    std::set<std::pair<DefAmrIndexUint, DefSizet>> set_point_index;
    std::vector<DefInt> vec_int{};
    std::vector<DefReal> vec_real{};
};
/**
* @struct GhostNode
* @brief structure to store ghost node information
*/
struct GhostNode{
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
    DefAmrUint flag_status_ = 0;
    virtual ~GridNode() {}
};
/**
* @class InterfaceLayerInfo
* @brief class used to store information of refinement interface layer
*/
class InterfaceLayerInfo {
 public:
    // number of extended based on the interface grid
    std::vector<DefAmrIndexLUint> k0ExtendOuterNeg_, k0ExtendOuterPos_;
    ///< number of extened layers outside the geometry
    std::vector<DefAmrIndexLUint> k0ExtendInnerNeg_, k0ExtendInnerPos_;
    ///< number of extened layers inside the geometry

    std::vector<DefMap<DefAmrUint>> vec_outer_coarse2fine_, vec_outer_fine2coarse_;
    std::vector<DefMap<DefAmrUint>> vec_inner_coarse2fine_, vec_inner_fine2coarse_;
};
/**
* @class TrackingGridInfoInterface
* @brief class used to store tracking grid information at i_level_ of refinement
*/
class TrackingGridInfoInterface {
 public:
    DefAmrUint computational_cost_ = 1;
    std::string node_type_;
    EGridExtendType grid_extend_type_ = EGridExtendType::kSameInAllDirections;

    // number of extended based on the tracking grid
    std::vector<DefAmrIndexLUint> k0ExtendOuterNeg_, k0ExtendOuterPos_;
    ///< number of extened layers outside the geometry
    std::vector<DefAmrIndexLUint> k0ExtendInnerNeg_, k0ExtendInnerPos_;
    ///< number of extened layers inside the geometry

    // index of creators in corresponding vector
    DefAmrIndexUint k0IndexCreator = 0;

    // information of TrackingNode
    DefMap<TrackingNode> map_tracking_node_{};
    DefAmrIndexUint k0NumIntForEachNode_ = 1;
    DefAmrIndexUint k0NumRealForEachNode_ = 0;
    TrackingNode k0TrackNodeInstance_;
};
/**
* @class TrackingGridInfoCreatorInterface
* @brief abstract class used to generate TrackingGridInfo instance.
*/
class TrackingGridInfoCreatorInterface {
 public:
    virtual std::shared_ptr<TrackingGridInfoInterface> CreateTrackingGridInfo() const = 0;
    virtual ~TrackingGridInfoCreatorInterface() {}
};
/**
* @class GhostGridInfoInterface
* @brief class used to store ghost grid information at i_level_ of refinement
*/
class GhostGridInfoInterface {
 public:
    // information of grid at each level of refinement
    DefAmrIndexUint i_level_ = 0;
    const DefAmrUint kCountIndex_ = 1;
    DefAmrUint computational_cost_ = 1;
    std::string node_type_;

    // information of GhostNode
    DefMap<GhostNode> map_ghost_node_{};
    DefAmrIndexUint k0NumIntForEachNode_ = 1;
    DefAmrIndexUint k0NumRealForEachNode_ = 0;
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
*/
class GridInfoInterface {
 public:
    // level of refinement
    DefAmrIndexUint i_level_ = 0;

    std::map<std::string, void*> kMemberNames_;
    void SetMemberVariable(const std::string& member_name, int value);

    DefAmrUint computational_cost_ = 1;
    std::string node_type_;
    std::vector<DefReal> grid_space_;
    std::shared_ptr<SolverInterface> ptr_solver_ = nullptr;

    std::map<std::pair<ECriterionType, DefAmrIndexUint>,
        std::shared_ptr<TrackingGridInfoInterface>> map_ptr_tracking_grid_info_;
    std::shared_ptr<GhostGridInfoInterface> ptr_ghost_grid_info_;

    // interface between grid of different refinement levels
    DefAmrIndexUint k0NumFine2CoarseLayer_ = 3;
    DefAmrIndexUint k0NumCoarse2FineLayer_ = 2;
    std::map<std::pair<ECriterionType, DefAmrIndexUint>,
        std::shared_ptr<InterfaceLayerInfo>> map_ptr_interface_layer_info_;

    // information of GridNode
    DefMap<std::unique_ptr<GridNode>> map_grid_node_{};
    DefMap<DefAmrUint> map_grid_count_exist_{};
    DefAmrIndexUint k0NumIntForEachNode_ = 0;
    DefAmrIndexUint k0NumRealForEachNode_ = 0;

    // domain boundary
    std::vector<DefSFBitset> k0VecBitsetDomainMin_, k0VecBitsetDomainMax_;
    ///< space filling codes of bounds for computational domain
    std::vector<DefMap<DefAmrIndexUint>> domain_boundary_min_, domain_boundary_max_;
    ///< map storing spacing filling codes of bounds (min and max in each coordinate) for computational domain

    // output related
    virtual void SetupOutputVariables() {}
    virtual void WriteOutputScalarAndVectors(FILE* const fp, const bool bool_binary,
        const Base64Utility& base64_instance,
        const OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index) const {}

    // node type
    virtual std::unique_ptr<GridNode> GridNodeCreator() {
        return std::make_unique<GridNode>();
    }
    virtual void SetPointerToCurrentNodeType() {}

    virtual void InitialGridInfo() = 0;
    virtual ~GridInfoInterface() {}
};
/**
* @class GridInfoCreatorInterface
* @brief abstract class used to generate GhostGridInfo instance.
*/
class GridInfoCreatorInterface {
 public:
    virtual std::shared_ptr<GridInfoInterface>
        CreateGridInfo() const = 0;
    virtual ~GridInfoCreatorInterface() {}
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_SOURCE_GRID_GRID_INFO_H_
