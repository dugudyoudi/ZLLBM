//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
class GridManagerInterface;
class MpiManager;
class CriterionManager;
/**
* @struct TrackingNode
* @brief structure to nodes information for tracking movement
*/
struct TrackingNode {
 public:
    // the first element (DefInt) refers to the level of  criterion points
    std::set<std::pair<DefInt, DefSizet>> set_point_index{};
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
*/
struct GridNode {
 public:
    DefInt flag_status_ = 0;
    virtual ~GridNode() {}

    // functions need to be override for interpolation
    virtual void InterpolationAdditionAssignCoefficient(const GridNode& /*node_in*/, const DefReal /*coefficient*/) {
        // *this += node_in * coefficient;
    }
    virtual void CopyVariablesToBuffer(const GridNode& /*node_ref*/, char* const /*ptr_node_buffer*/) const {}
    virtual void ReadVariablesFromBuffer(const char* /*ptr_node_buffer*/, const GridNode* const /*ptr_node*/) {}
    virtual int CopyANodeToBufferForMpi(char* const ptr_node_buffer) const {
        std::memcpy(ptr_node_buffer, &flag_status_, sizeof(DefInt));
        int offset = sizeof(DefInt);
        return offset;
    }
    virtual int ReadANodeFromBufferForMpi(const char* ptr_node_buffer) {
        std::memcpy(&flag_status_, ptr_node_buffer, sizeof(DefInt));
        int offset = sizeof(DefInt);
        return offset;
    }
    virtual int CopyAInterpNodeToBufferForMpi(char* const ptr_node_buffer) const {
        return CopyANodeToBufferForMpi(ptr_node_buffer);
    }
    virtual int ReadAInterpNodeFromBufferForMpi(const char* ptr_node_buffer) {
        return ReadANodeFromBufferForMpi(ptr_node_buffer);
    }
    virtual int CopyANodeToBufferForCheckpoint(char* const ptr_node_buffer) const {
        return CopyANodeToBufferForMpi(ptr_node_buffer);
    }
    virtual int ReadANodeFromBufferForCheckpoint(const char* ptr_node_buffer) {
        return ReadANodeFromBufferForMpi(ptr_node_buffer);
    }
};
/**
* @struct OutputNodeVariableInfoInterface
* @brief class to store variables need to output for each node
*/
struct OutputNodeVariableInfoInterface {
    DefInt variable_dims_;
    std::string output_name_;
    virtual void WriteNodeVariable(const GridNode& grid_node){}
    virtual ~OutputNodeVariableInfoInterface() {}
};
/**
* @struct DomainInfo
* @brief class to store information of computational domain at current level
*/
struct DomainInfo {
 public:
    std::vector<DefReal> grid_space_;  ///< grid space at current grid level
    bool bool_periodic_domain_ = false;
    std::vector<bool> periodic_min_, periodic_max_;
    std::vector<DefSFBitset> domain_min_n_level_, domain_max_n_level_;
};
/**
* @class InterfaceLayerInfo
* @brief class used to store information of refinement interface layer
*/
class InterfaceLayerInfo {
 public:
    // number of extended based on the interface grid
    std::vector<DefAmrLUint> k0ExtendOuterNeg_, k0ExtendOuterPos_;
    ///< number of extened layers outside the geometry
    std::vector<DefAmrLUint> k0ExtendInnerNeg_, k0ExtendInnerPos_;
    ///< number of extened layers inside the geometry

    std::vector<DefMap<DefInt>> vec_outer_coarse2fine_, vec_outer_fine2coarse_;
    std::vector<DefMap<DefInt>> vec_inner_coarse2fine_, vec_inner_fine2coarse_;
};
/**
* @class TrackingGridInfoInterface
* @brief class used to store tracking grid information at i_level_ of refinement
*/
class TrackingGridInfoInterface {
 public:
    DefInt computational_cost_ = 1;
    std::string node_type_;
    EGridExtendType grid_extend_type_ = EGridExtendType::kSameInAllDirections;

    // number of extended based on the tracking grid
    std::vector<DefAmrLUint> k0ExtendOuterNeg_, k0ExtendOuterPos_;
    ///< number of extened layers outside the geometry
    std::vector<DefAmrLUint> k0ExtendInnerNeg_, k0ExtendInnerPos_;
    ///< number of extened layers inside the geometry

    // index of creators in corresponding vector
    DefInt k0IndexCreator = 0;

    // information of TrackingNode
    DefMap<TrackingNode> map_tracking_node_{};
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
    DefInt i_level_ = 0;
    const DefInt kCountIndex_ = 1;
    DefInt computational_cost_ = 1;
    std::string node_type_;

    // information of GhostNode
    DefMap<GhostNode> map_ghost_node_{};
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
 protected:
    DefInt i_level_ = 0;  ///< level of refinement
    std::map<std::string, void*> kMemberNames_;
    DefInt computational_cost_ = 1;
    std::string node_type_;
    std::vector<DefReal> grid_space_;
    std::weak_ptr<SolverInterface> ptr_solver_;
    SFBitsetAuxInterface* ptr_sfbitset_aux_ = nullptr;

    std::shared_ptr<GhostGridInfoInterface> ptr_ghost_grid_info_;

    // interface between grid of different refinement levels
    DefInt k0NumFine2CoarseLayer_ = 3;
    DefInt k0NumCoarse2FineLayer_ = 2;
    DefInt k0NumFine2CoarseGhostLayer_ = k0NumFine2CoarseLayer_/2 + 1;
    DefInt k0NumCoarse2FineGhostLayer_ = k0NumCoarse2FineLayer_/2;

 public:
    // set and get protected members
    DefInt GetGridLevel() const {return i_level_;}
    DefInt GetComputationalCost() const {return computational_cost_;}
    const std::vector<DefReal>& GetGridSpace() const {return grid_space_;}
    DefInt GetNumFine2CoarseLayer() const {return k0NumFine2CoarseLayer_;}
    DefInt GetNumCoarse2FineLayer() const {return k0NumCoarse2FineLayer_;}
    DefInt GetNumFine2CoarseGhostLayer() const {return k0NumFine2CoarseGhostLayer_;}
    DefInt GetNumCoarse2FineGhostLayer() const {return k0NumCoarse2FineGhostLayer_;}
    SFBitsetAuxInterface* GetPtrSFBitsetAux() const {
        if (ptr_sfbitset_aux_ == nullptr) {
            amrproject::LogManager::LogError("pointer to SFBitsetAux is nullptr");
        }
        return ptr_sfbitset_aux_;
    }
    const std::string& GetNodeType() const {return node_type_;}
    std::weak_ptr<SolverInterface> GetPtrToSolver() const {
        if (ptr_solver_.expired()) {
            amrproject::LogManager::LogError("pointer to solver is empty");
        }
        return ptr_solver_;
    }

    void SetGridLevel(const DefInt i_level);
    void SetComputationalCost(const DefInt computational_cost);
    void SetGridSpace(const std::vector<DefReal>& grid_space);
    void SetPtrSFBitsetAux(SFBitsetAuxInterface* const ptr_sfbitset_aux);
    void SetPtrSolver(const std::weak_ptr<SolverInterface>& ptr_solver);
    void SetMemberVariable(const std::string& member_name, int value);
    void SetNodeType(const std::string& node_type);
    void SetNumFine2CoarseLayer(const DefInt num_fine2coarse_layer);
    void SetNumCoarse2FineLayer(const DefInt num_coarse2fine_layer);
    DomainInfo GetDomainInfo() const;

    // simulation related
    virtual void AdvancingAtCurrentTime(const ETimeSteppingScheme time_scheme,
        const DefInt time_step_level, const DefReal time_step_current,
        MpiManager* const ptr_mpi_manager, CriterionManager* const ptr_criterion_manager) {}

    // layers for refinement
    std::map<std::pair<ECriterionType, DefInt>,
        std::shared_ptr<TrackingGridInfoInterface>> map_ptr_tracking_grid_info_;
    std::map<std::pair<ECriterionType, DefInt>,
        std::shared_ptr<InterfaceLayerInfo>> map_ptr_interface_layer_info_;

    // information of GridNode
    DefMap<std::unique_ptr<GridNode>> map_grid_node_{};
    DefMap<DefInt> map_grid_count_exist_{};
    virtual void SetNodeVariablesAsZeros(GridNode* const ptr_node) {}  // will be called in interpolation

    // domain boundary related
    // noting that kFlagInsideDomain indicates nodes are not on the domain boundary
    static constexpr DefInt kFlagInsideDomain_ = 0,  kFlagOutsideDomain_ = -1,
        kFlagXMinBoundary_ = 1, kFlagXMaxBoundary_ = 2, kFlagYMinBoundary_ = 4,
        kFlagYMaxBoundary_ = 8, kFlagZMinBoundary_ = 16, kFlagZMaxBoundary_ = 32;
    std::vector<DefSFBitset> k0VecBitsetDomainMin_, k0VecBitsetDomainMax_;
    ///< space filling codes of bounds for computational domain
    std::vector<DefMap<DefInt>> domain_boundary_min_, domain_boundary_max_;
    ///< map storing spacing filling codes of bounds (min and max in each coordinate) for computational domain
    int CheckIfNodeOutsideCubicDomain(const DefInt dims,
        const DefSFBitset& bitset_in, const SFBitsetAuxInterface& sfbitset_aux) const;
    void CheckNodesOnCubicPeriodicBoundary(const DefInt dims, const DefSFBitset& bitset_in,
        const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const SFBitsetAuxInterface& sfbitset_aux, std::vector<DefSFBitset>* const ptr_nodes_periodic) const;

    // output related
    std::vector<std::unique_ptr<OutputNodeVariableInfoInterface>> output_variables_;
    virtual void SetupOutputVariables() {}
    virtual void WriteOutputScalarAndVectors(FILE* const /*fp*/, const bool /*bool_binary*/,
        const Base64Utility& /*base64_instance*/,
        const OutputDataFormat& /*output_data_format*/,
        const DefMap<DefSizet>& /*map_node_index*/) const {}

    // node type
    virtual void InitialNotComputeNodeFlag() {}
    virtual std::unique_ptr<GridNode> GridNodeCreator() const {
        return std::make_unique<GridNode>();
    }
    virtual void SetPointerToCurrentNodeType() {}

    // time marching related
    virtual void SetUpGridAtBeginningOfTimeStep(const DefInt /*time_step*/) {}

    // communication between grid of different refinement levels
    virtual int TransferInfoFromCoarseGrid(const SFBitsetAuxInterface& /*sfbitset_aux*/,
        const DefInt /*node_flag*/, const GridInfoInterface& /*grid_info_coarse*/) {return -1;}
    virtual int TransferInfoToCoarseGrid(const SFBitsetAuxInterface& /*sfbitset_aux*/,
        const DefInt /*node_flag*/, GridInfoInterface* const /*ptr_grid_info_coarse*/) {return -1;}

    // parent grid manager
 protected:
    GridManagerInterface* ptr_parent_grid_manager_ = nullptr;

 public:
    void SetPtrToParentGridManager(GridManagerInterface* ptr_grid_manager) {
        ptr_parent_grid_manager_ = ptr_grid_manager;}
    GridManagerInterface* GetPtrToParentGridManager() const {
        if (ptr_parent_grid_manager_ == nullptr) {
            amrproject::LogManager::LogError("pointer to parent grid manager is nullptr");
        }
        return ptr_parent_grid_manager_;
    }

    // interpolation
 public:
    EInterpolationMethod interp_method_ = EInterpolationMethod::kLagrangian;
    DefAmrLUint max_interp_length_ = 2;
    ///< the maximum half length of a cubic region used for interpolation
    std::function<int(const DefAmrLUint, const DefAmrLUint, const DefInt, const DefSFBitset&,
        const amrproject::SFBitsetAuxInterface&, const std::vector<DefSFBitset>&,
        const DefMap<std::unique_ptr<GridNode>>& nodes_fine, const amrproject::GridInfoInterface& coarse_grid_info,
        const DefMap<std::unique_ptr<GridNode>>& nodes_coarse, GridNode* const ptr_node)> func_node_interp_;
    void ChooseInterpolationMethod(const DefInt dims);
    virtual void NodeInfoCoarse2fine(const GridNode& coarse_node, GridNode* const ptr_fine_node) const {}
    virtual void NodeInfoFine2Coarse(const GridNode& fine_node, GridNode* const ptr_coarse_node) const {}
    // linear interpolation
    int InterpolationLinear2D(const DefAmrLUint region_length,
        const DefInt flag_not_for_interp_coarse, const DefSFBitset& sfbitset_in,
        const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefSFBitset>& sfbitset_region,
        const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<GridNode>>& nodes_coarse,
        GridNode* const ptr_node);
    int InterpolationLinear3D(const DefAmrLUint region_length,
        const DefInt flag_not_for_interp_coarse, const DefSFBitset& sfbitset_in,
        const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefSFBitset>& sfbitset_region,
        const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<GridNode>>& nodes_coarse,
        GridNode* const ptr_node);
    // lagrangian interpolation
 private:
    struct LagrangianCoeff {
        std::vector<DefReal> coeff0, coeff1;
    };
    std::map<DefAmrLUint, LagrangianCoeff> lagrangian_coefficients_;

 public:
    const LagrangianCoeff& CalculateLagrangianInterpCoeff(const DefAmrLUint interp_half_length);
    int InterpolationLagrangian2D(const DefAmrLUint interpolation_length,
        const DefAmrLUint region_length, const DefInt flag_not_for_interp_coarse,
        const DefSFBitset& sfbitset_in, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefSFBitset>& sfbitset_region, const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info,
        const DefMap<std::unique_ptr<GridNode>>& nodes_coarse, GridNode* const ptr_node);
    int InterpolationLagrangian3D(const DefAmrLUint interpolation_length,
        const DefAmrLUint region_length, const DefInt flag_not_for_interp_coarse,
        const DefSFBitset& sfbitset_in, const SFBitsetAuxInterface& sfbitset_aux,
        const std::vector<DefSFBitset>& sfbitset_region, const DefMap<std::unique_ptr<GridNode>>& nodes_fine,
        const GridInfoInterface& coarse_grid_info,
        const DefMap<std::unique_ptr<GridNode>>& nodes_coarse, GridNode* const ptr_node);

    // mpi related
 public:
    DefMap<std::unique_ptr<GridNode>> interp_nodes_outer_layer_;
    std::vector<int> vec_num_interp_nodes_receive_;
    std::map<int, DefMap<DefInt>> interp_nodes_inner_layer_;
    int AddGhostNodesForInterpolation(const std::vector<bool>& periodic_min, const std::vector<bool>& periodic_max,
        const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefInt>& refinement_interface,
        const DefMap<std::unique_ptr<GridNode>>& map_nodes_lower);
    virtual bool CheckIfPeriodicDomainRequired(const DefInt dims,
        std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const {
        ptr_periodic_min->assign(dims, false);
        ptr_periodic_max->assign(dims, false);
        return false;
    }
    virtual int GetSizeOfGridNodeInfoForMpiCommunication() const {return sizeof(DefInt);}
    virtual int GetSizeOfGridNodeInfoForCheckPoint() const {return GetSizeOfGridNodeInfoForMpiCommunication();}
    virtual void CopyNodeInfoToBuffer(
        const std::function<void(const GridNode& node_ref, char* const)>& func_copy_buffer,
        const DefMap<DefInt>& map_nodes, char* const ptr_buffer) const;
    virtual void ReadNodeInfoFromBuffer(
        const std::function<void(const char*,  GridNode* const ptr_node)>& func_read_buffer,
        const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer);
    virtual void CopyInterpolationNodeInfoToBuffer(
        const std::function<void(const GridNode& node_ref, char* const)>& func_copy_buffer,
        const GridInfoInterface& grid_info_lower,
        const DefMap<DefInt>& map_nodes, char* const ptr_buffer);
    virtual void ReadInterpolationNodeInfoFromBuffer(
        const std::function<void(const char*,  GridNode* const ptr_node)>& func_read_buffer,
        const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer);
    virtual void ComputeInfoInMpiLayers(
        const std::map<int, DefMap<DefInt>>& map_inner_nodes, const DefMap<DefInt>& map_outer_nodes) {}
    virtual void ComputeInfoInInterpMpiLayers(
        const std::map<int, DefMap<DefInt>>& map_intper_nodes) {}
    void RemoveUnnecessaryF2CNodesOnMpiOuterLayer(const DefSFCodeToUint code_min,
        const DefSFCodeToUint code_max, const DefInt num_outer_layer, DefMap<DefInt>* const ptr_mpi_outer_layer);

    // general purpose functions
    virtual void InitialGridInfo(const DefInt dims) = 0;
    virtual ~GridInfoInterface() {}
};
/**
* @class GridInfoCreatorInterface
* @brief abstract class used to generate GhostGridInfo instance.
*/
class GridInfoCreatorInterface {
 public:
    virtual std::shared_ptr<GridInfoInterface> CreateGridInfo() const = 0;
    virtual ~GridInfoCreatorInterface() {}
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_GRID_GRID_INFO_INTERFACE_H_
