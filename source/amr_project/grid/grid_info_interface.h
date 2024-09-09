//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/grid_info_interface.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define classes to store grid information
*/

#ifndef SOURCE_AMR_PROJECT_GRID_GRID_INFO_INTERFACE_H_
#define SOURCE_AMR_PROJECT_GRID_GRID_INFO_INTERFACE_H_
#include <vector>
#include <array>
#include <map>
#include <set>
#include <memory>
#include <utility>
#include <string>
#include <functional>
#include "../defs_libs.h"
#include "criterion/criterion_numerates.h"
#include "./solver_info_interface.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_enumerates.h"
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
    ///< count of this tracking node relies on how many criterion points
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

    // functions need to be override for interpolation
    virtual void InterpolationAdditionAssignCoefficient(const GridNode& node_in, const DefReal coefficient) {
        // *this += node_in * coefficient;
    }
};
/**
* @struct OutputNodeVariableInfoInterface
* @brief class to store variables need to output for each node
*/
struct OutputNodeVariableInfoInterface {
    DefAmrIndexUint variable_dims_;
    std::string output_name_;
    virtual void WriteNodeVariable(const GridNode grid_node){}
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
    virtual std::shared_ptr<TrackingGridInfoInterface> CreateTrackingGridInfo() const {
        return std::make_shared<TrackingGridInfoInterface>();
    }
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
    SFBitsetAuxInterface* ptr_sfbitset_aux_ = nullptr;

    std::map<std::pair<ECriterionType, DefAmrIndexUint>,
        std::shared_ptr<TrackingGridInfoInterface>> map_ptr_tracking_grid_info_;
    std::shared_ptr<GhostGridInfoInterface> ptr_ghost_grid_info_;

    // interface between grid of different refinement levels
    DefAmrIndexUint k0NumFine2CoarseLayer_ = 3;  // (k0NumCoarse2FineLayer_)*2 - 1
    DefAmrIndexUint k0NumCoarse2FineLayer_ = 2;
    DefAmrIndexUint k0NumFine2CoarseGhostLayer_ = k0NumFine2CoarseLayer_/2 + 1;
    DefAmrIndexUint k0NumCoarse2FineGhostLayer_ = k0NumCoarse2FineLayer_/2;
    std::map<std::pair<ECriterionType, DefAmrIndexUint>,
        std::shared_ptr<InterfaceLayerInfo>> map_ptr_interface_layer_info_;

    // information of GridNode
    DefMap<std::unique_ptr<GridNode>> map_grid_node_{};
    DefMap<DefAmrUint> map_grid_count_exist_{};
    DefAmrIndexUint k0NumIntForEachNode_ = 0;
    DefAmrIndexUint k0NumRealForEachNode_ = 0;
    virtual void SetNodeVariablesAsZeros(GridNode* const ptr_node) {}  // will be called in interpolation

    // domain boundary related
    // noting that kFlagInsideDomain indicates nodes are not on the domain boundary
    static constexpr int kFlagInsideDomain_ = 0,  kFlagOutsideDomain_ = -1,
        kFlagXMinBoundary_ = 1, kFlagXMaxBoundary_ = 2, kFlagYMinBoundary_ = 4,
        kFlagYMaxBoundary_ = 8, kFlagZMinBoundary_ = 16, kFlagZMaxBoundary_ = 32;
    std::vector<DefSFBitset> k0VecBitsetDomainMin_, k0VecBitsetDomainMax_;
    ///< space filling codes of bounds for computational domain
    std::vector<DefMap<DefAmrIndexUint>> domain_boundary_min_, domain_boundary_max_;
    ///< map storing spacing filling codes of bounds (min and max in each coordinate) for computational domain
    int CheckIfNodeOutsideCubicDomain(const DefAmrIndexUint dims,
        const DefSFBitset& bitset_in, const SFBitsetAuxInterface& sfbitset_aux) const;
    void CheckNodesOnCubicPeriodicBoundary(const DefAmrIndexUint dims, const DefSFBitset& bitset_in,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const SFBitsetAuxInterface& sfbitset_aux, std::vector<DefSFBitset>* const ptr_nodes_periodic) const;

    // output related
    std::vector<std::unique_ptr<OutputNodeVariableInfoInterface>> output_variables_;
    virtual void SetupOutputVariables() {}
    virtual void WriteOutputScalarAndVectors(FILE* const fp, const bool bool_binary,
        const Base64Utility& base64_instance,
        const OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index) const {}

    // node type
    virtual void InitialNotComputeNodeFlag() {}
    virtual std::unique_ptr<GridNode> GridNodeCreator() {
        return std::make_unique<GridNode>();
    }
    virtual void SetPointerToCurrentNodeType() {}

    // time marching related
    virtual void SetUpGridAtBeginningOfTimeStep(
        const DefAmrIndexUint time_step, GridManagerInterface* const ptr_grid_manager) {}

    // communication between grid of different refinement levels
    virtual int TransferInfoFromCoarseGrid(const SFBitsetAuxInterface& sfbitset_aux,
        const DefAmrUint node_flag, const GridInfoInterface& grid_info_coarse) {return -1;}
    virtual int TransferInfoToCoarseGrid(const SFBitsetAuxInterface& sfbitset_aux,
        const DefAmrUint node_flag, GridInfoInterface* const ptr_grid_info_coarse) {return -1;}

    // debug
    virtual void DebugWrite() {}
    virtual void DebugWriteNode(const GridNode& node) {}

    // interpolation
 public:
    EInterpolationMethod interp_method_ = EInterpolationMethod::kLinear;
    DefAmrIndexLUint max_interp_length_ = 2;
    ///< the maximum half length of a cubic region used for interpolation
    std::function<int(const DefAmrIndexLUint, const DefAmrIndexLUint, const DefAmrUint, const DefSFBitset&,
        const amrproject::SFBitsetAuxInterface&, const std::vector<DefSFBitset>&,
        const DefMap<std::unique_ptr<GridNode>>& nodes_fine, const amrproject::GridInfoInterface& coarse_grid_info,
        const DefMap<std::unique_ptr<GridNode>>& nodes_coarse, GridNode* const ptr_node)> func_node_interp_;
    void ChooseInterpolationMethod();
    virtual void NodeInfoCoarse2fine(const GridNode& coarse_node, GridNode* const ptr_fine_node) const {}
    virtual void NodeInfoFine2Coarse(const GridNode& fine_node, GridNode* const ptr_coarse_node) const {}
    // linear interpolation
    int InterpolationLinear2D(const DefAmrIndexLUint region_length,
        const DefAmrUint flag_not_for_interp_coarse, const DefSFBitset& sfbitset_in,
        const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefSFBitset>& sfbitset_region,
        const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<GridNode>>& nodes_coarse,
        GridNode* const ptr_node);
    int InterpolationLinear3D(const DefAmrIndexLUint region_length,
        const DefAmrUint flag_not_for_interp_coarse, const DefSFBitset& sfbitset_in,
        const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefSFBitset>& sfbitset_region,
        const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<GridNode>>& nodes_coarse,
        GridNode* const ptr_node);
    // lagrangian interpolation
 private:
    struct LagrangianCoeff {
        std::vector<DefReal> coeff0, coeff1;
    };
    std::map<DefAmrIndexLUint, LagrangianCoeff> lagrangian_coefficients_;

 public:
    const LagrangianCoeff& CalculateLagrangianInterpCoeff(const DefAmrIndexLUint interp_half_length);
    int InterpolationLagrangian2D(const DefAmrIndexLUint interpolation_length,
        const DefAmrIndexLUint region_length, const DefAmrUint flag_not_for_interp_coarse,
        const DefSFBitset& sfbitset_in, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefSFBitset>& sfbitset_region, const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info,
        const DefMap<std::unique_ptr<GridNode>>& nodes_coarse, GridNode* const ptr_node);
    int InterpolationLagrangian3D(const DefAmrIndexLUint interpolation_length,
        const DefAmrIndexLUint region_length, const DefAmrUint flag_not_for_interp_coarse,
        const DefSFBitset& sfbitset_in, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefSFBitset>& sfbitset_region, const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info,
        const DefMap<std::unique_ptr<GridNode>>& nodes_coarse, GridNode* const ptr_node);

    // mpi related
 public:
    DefMap<std::unique_ptr<GridNode>> interp_nodes_outer_layer_;
    int AddGhostNodesForInterpolation(const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefAmrUint>& refinement_interface,
        const DefMap<std::unique_ptr<GridNode>>& map_nodes_lower);
    virtual bool CheckIfPeriodicDomainRequired(const DefAmrIndexUint dims,
        std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const {
        ptr_periodic_min->assign(dims, false);
        ptr_periodic_max->assign(dims, false);
        return false;
    }
    virtual int SizeOfGridNodeInfoForMpiCommunication() const {return 0;}
    virtual int CopyNodeInfoToBuffer(const DefMap<DefAmrIndexUint>& map_nodes, char* const ptr_buffer) {return 0;}
    virtual int CopyInterpolationNodeInfoToBuffer(const GridInfoInterface& grid_info_lower,
        const DefMap<DefAmrIndexUint>& map_nodes, char* const ptr_buffer) {return 0;}
    virtual int ReadInterpolationNodeInfoFromBuffer(
        const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) {return 0;}
    virtual void ComputeInfoInMpiLayers(
        const std::map<int, DefMap<DefAmrIndexUint>>& map_inner_nodes,
        const DefMap<DefAmrIndexUint>& map_outer_nodes) {}
    virtual void ReadNodeInfoFromBuffer(const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) {}
    void SetPeriodicBoundaryAsPartitionInterface(DefAmrIndexUint dims, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<bool>& bool_periodic_min, const std::vector<bool>& bool_periodic_max,
        DefMap<DefAmrIndexUint>* const ptr_partition_interface);

    // general purpose functions
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
#endif  // SOURCE_AMR_PROJECT_GRID_GRID_INFO_INTERFACE_H_
