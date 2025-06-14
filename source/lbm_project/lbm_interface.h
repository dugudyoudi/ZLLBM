//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_interface.h
* @author Zhengliang Liu
* @brief define classes to manage LBM models.
* @date  2023-9-30
*/
#ifndef SOURCE_LBM_PROJECT_LBM_INTERFACE_H_
#define SOURCE_LBM_PROJECT_LBM_INTERFACE_H_
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include "../defs_libs.h"
#include "./auxiliary_inline_func.h"
#include "./solver_info_interface.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
#include "boundary_conditions/lbm_boundary_conditions.h"
#include "collision_model/lbm_collision_models.h"
#include "les_models/lbm_les_models.h"
namespace rootproject {
namespace lbmproject {
class SolverLbmInterface;
/**
* @brief structure to store node information for LBM simulation
*/
struct GridNodeLbm : public amrproject::GridNode {
    DefReal rho_ = 1.;
    const int kSizeReal = sizeof(DefReal);
    std::vector<DefReal> velocity_{};
    std::vector<DefReal> f_{}, f_collide_{};
    std::vector<DefReal> force_{};
    GridNodeLbm() {}
    GridNodeLbm(const DefReal rho0, const std::vector<DefReal>& velocity0,
        const std::vector<DefReal>& f, const std::vector<DefReal>& f_collide)
        : rho_(rho0), velocity_(velocity0), f_(f), f_collide_(f_collide) {
            velocity_.shrink_to_fit();
            f_.shrink_to_fit();
            f_collide_.shrink_to_fit();
    }
    GridNodeLbm(const DefReal rho0, const std::vector<DefReal>& velocity0,
        const std::vector<DefReal>& force0,
        const std::vector<DefReal>& f, const std::vector<DefReal>& f_collide)
        : rho_(rho0), velocity_(velocity0), f_(f), f_collide_(f_collide), force_(force0) {
            velocity_.shrink_to_fit();
            force_.shrink_to_fit();
            f_.shrink_to_fit();
            f_collide_.shrink_to_fit();
    }
    int CopyANodeToBufferForCheckpoint(char* const ptr_node_buffer) const override {
        std::memcpy(ptr_node_buffer, &flag_status_, sizeof(DefInt));
        int offset = sizeof(DefInt);
        std::memcpy(ptr_node_buffer + offset, f_.data(), f_.size() * kSizeReal);
        offset += static_cast<int>(f_.size()) * kSizeReal;
        if (!force_.empty()) {
            std::memcpy(ptr_node_buffer + offset, force_.data(), force_.size() * kSizeReal);
            offset += static_cast<int>(force_.size()) * kSizeReal;
        }
        return offset;
    }
    int ReadANodeFromBufferForCheckpoint(const char* ptr_node_buffer) override {
        std::memcpy(&flag_status_, ptr_node_buffer, sizeof(DefInt));
        int offset = sizeof(DefInt);
        std::memcpy(f_.data(), ptr_node_buffer + offset, f_.size() * kSizeReal);
        offset += static_cast<int>(f_.size()) * kSizeReal;
        if (!force_.empty()) {
            std::memcpy(force_.data(), ptr_node_buffer + offset, force_.size() * kSizeReal);
            offset += static_cast<int>(force_.size()) * kSizeReal;
        }
        return offset;
    }

    // compulsory function from the GridNode base class
    int CopyANodeToBufferForMpi(char* const ptr_node_buffer) const override {
        std::memcpy(ptr_node_buffer, &rho_, kSizeReal);
        int offset = kSizeReal;
        std::memcpy(ptr_node_buffer + offset, velocity_.data(), velocity_.size() * kSizeReal);
        offset += static_cast<int>(velocity_.size()) * kSizeReal;
        std::memcpy(ptr_node_buffer + offset, f_.data(), f_.size() * kSizeReal);
        offset += static_cast<int>(f_.size()) * kSizeReal;
        if (!force_.empty()) {
            std::memcpy(ptr_node_buffer + offset, force_.data(), force_.size() * kSizeReal);
            offset += static_cast<int>(force_.size()) * kSizeReal;
        }
        return offset;
    }
    int ReadANodeFromBufferForMpi(const char* ptr_node_buffer) override {
        std::memcpy(&rho_, ptr_node_buffer, kSizeReal);
        int offset = kSizeReal;
        std::memcpy(velocity_.data(), ptr_node_buffer + offset, velocity_.size() * kSizeReal);
        offset += static_cast<int>(velocity_.size()) * kSizeReal;
        std::memcpy(f_.data(), ptr_node_buffer + offset, f_.size() * kSizeReal);
        offset += static_cast<int>(f_.size()) * kSizeReal;
        if (!force_.empty()) {
            std::memcpy(force_.data(), ptr_node_buffer + offset, force_.size() * kSizeReal);
            offset += static_cast<int>(force_.size()) * kSizeReal;
        }
        return offset;
    }

    // GridNodeLbm& operator=(const GridNodeLbm& node_r) {
    //     this->rho_ = node_r.rho_;
    //     this->velocity_.resize(node_r.velocity_.size());
    //     std::copy(node_r.velocity_.begin(), node_r.velocity_.end(), this->velocity_.begin());
    //     this->f_collide_.resize(node_r.f_collide_.size());
    //     std::copy(node_r.f_collide_.begin(), node_r.f_collide_.end(), this->f_collide_.begin());
    //     return *this;
    // }
    // GridNodeLbm& operator+=(const GridNodeLbm& node_r) {
    //     this->rho_ += node_r.rho_;
    //     this->velocity_.resize(node_r.velocity_.size());
    //     std::transform(node_r.velocity_.begin(), node_r.velocity_.end(),
    //         this->velocity_.begin(), this->velocity_.begin(), std::plus<DefReal>());
    //     this->f_collide_.resize(node_r.f_collide_.size());
    //     std::transform(node_r.f_collide_.begin(), node_r.f_collide_.end(),
    //         this->f_collide_.begin(), this->f_collide_.begin(), std::plus<DefReal>());
    //     return *this;
    // }
    // GridNodeLbm operator+(const GridNodeLbm& node_r) const {
    //     GridNodeLbm node_result;
    //     node_result.rho_ = this->rho_ + node_r.rho_;
    //     node_result.velocity_.resize(node_r.velocity_.size());
    //     std::transform(node_r.velocity_.begin(), node_r.velocity_.end(),
    //         this->velocity_.begin(), node_result.velocity_.begin(), std::plus<DefReal>());
    //     node_result.f_collide_.resize(node_r.f_collide_.size());
    //     std::transform(node_r.f_collide_.begin(), node_r.f_collide_.end(),
    //         this->f_collide_.begin(), node_result.f_collide_.begin(), std::plus<DefReal>());
    //     return node_result;
    // }
    // GridNodeLbm operator-(const GridNodeLbm& node_r) const {
    //     GridNodeLbm node_result;
    //     node_result.rho_ = this->rho_ - node_r.rho_;
    //     node_result.velocity_.resize(node_r.velocity_.size());
    //     std::transform(node_r.velocity_.begin(), node_r.velocity_.end(),
    //         this->velocity_.begin(), node_result.velocity_.begin(), std::minus<DefReal>());
    //     node_result.f_collide_.resize(node_r.f_collide_.size());
    //     std::transform(node_r.f_collide_.begin(), node_r.f_collide_.end(),
    //         this->f_collide_.begin(), node_result.f_collide_.begin(), std::minus<DefReal>());
    //     return node_result;
    // }
    // GridNodeLbm operator*(const GridNodeLbm& node_r) const {
    //     GridNodeLbm node_result;
    //     node_result.rho_ = this->rho_ * node_r.rho_;
    //     node_result.velocity_.resize(node_r.velocity_.size());
    //     std::transform(node_r.velocity_.begin(), node_r.velocity_.end(),
    //         this->velocity_.begin(), node_result.velocity_.begin(), std::multiplies<DefReal>());
    //     node_result.f_collide_.resize(node_r.f_collide_.size());
    //     std::transform(node_r.f_collide_.begin(), node_r.f_collide_.end(),
    //         this->f_collide_.begin(), node_result.f_collide_.begin(), std::multiplies<DefReal>());
    //     return node_result;
    // }
    // GridNodeLbm operator/(const GridNodeLbm& node_r) const {
    //     GridNodeLbm node_result;
    //     node_result.rho_ = this->rho_ / node_r.rho_;
    //     node_result.velocity_.resize(node_r.velocity_.size());
    //     std::transform(node_r.velocity_.begin(), node_r.velocity_.end(),
    //         this->velocity_.begin(), node_result.velocity_.begin(), std::divides<DefReal>());
    //     node_result.f_collide_.resize(node_r.f_collide_.size());
    //     std::transform(node_r.f_collide_.begin(), node_r.f_collide_.end(),
    //         this->f_collide_.begin(), node_result.f_collide_.begin(), std::divides<DefReal>());
    //     return node_result;
    // }
    // GridNodeLbm fabs() const {
    //     GridNodeLbm node_result;
    //     node_result.rho_ = std::fabs(this->rho_);
    //     node_result.velocity_.resize(this->velocity_.size());
    //     std::transform(this->velocity_.begin(), this->velocity_.end(),
    //         node_result.velocity_.begin(), [](DefReal num) { return std::abs(num); });
    //     node_result.f_collide_.resize(this->f_collide_.size());
    //     std::transform(this->f_collide_.begin(), this->f_collide_.end(),
    //         node_result.f_collide_.begin(), [](DefReal num) { return std::abs(num); });
    //     return node_result;
    // }
    // GridNodeLbm operator*(const DefReal real_r) const {
    //     GridNodeLbm node_result;
    //     node_result.rho_ = this->rho_ * real_r;
    //     node_result.velocity_.resize(this->velocity_.size());
    //     std::transform(this->velocity_.begin(), this->velocity_.end(),
    //         node_result.velocity_.begin(), [real_r](DefReal num) { return real_r * num; });
    //     node_result.f_collide_.resize(this->f_collide_.size());
    //     std::transform(this->f_collide_.begin(), this->f_collide_.end(),
    //         node_result.f_collide_.begin(), [real_r](DefReal num) { return real_r * num; });
    //     return node_result;
    // }

    void InterpolationAdditionAssignCoefficient(const GridNode& node_in, const DefReal coefficient) override;
};
struct OutputLBMNodeVariableInfo : public amrproject::OutputNodeVariableInfoInterface {
    std::function<std::vector<DefReal>(GridNodeLbm* const)> func_get_var_;
};
/**
* @brief  class to store LBM grid information at each refinement level
*/
class GridInfoLbmInteface : public amrproject::GridInfoInterface {
 protected:
    std::vector<DefReal> f_ini_, f_collide_ini_;
    // convert pointer to current type
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_grid_nodes_ = nullptr;

 public:
    void InitialGridInfo(const DefInt dims) override;

    void SetPointerToCurrentNodeType() override;
    DefMap<std::unique_ptr<GridNodeLbm>>* GetPtrToLbmGrid();

    void AdvancingAtCurrentTime(const amrproject::ETimeSteppingScheme time_scheme,
        const DefInt time_step_level, const DefReal time_step_current,
        amrproject::MpiManager* const ptr_mpi_manager,
        amrproject::CriterionManager* const ptr_criterion_manager) override;

    // pointer to boundary conditions adopted for domain boundary
    std::map<amrproject::EDomainBoundaryDirection, std::unique_ptr<BoundaryConditionLbmInterface>>
        domain_boundary_condition_;

    // set and get functions
    std::unique_ptr<amrproject::GridNode> GridNodeCreator() const override;
    void SetNodeVariablesAsZeros(amrproject::GridNode* const ptr_node) override;
    void InitialNotComputeNodeFlag() override;
    DefInt GetNodeFlagNotStream() const { return NodeFlagNotStream_;}
    DefInt GetNodeFlagNotCollision() const { return NodeFlagNotCollision_;}
    int GetSizeOfGridNodeInfoForMpiCommunication() const override;
    int GetSizeOfGridNodeInfoForCheckPoint() const override;

 protected:
    // domain boundaries
    bool CheckIfPeriodicDomainRequired(const DefInt dims,
        std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const override;
    void ComputeDomainBoundaryCondition();
    void ComputeDomainBoundaryConditionForANode(int flag_node_boundary, const DefSFBitset bitset_in);

    // time marching related
    void SetUpGridAtBeginningOfTimeStep(const DefInt time_step) override;
    void UpdateCriterion(const DefReal time_step_current,
        const amrproject::MpiManager& ptr_mpi_manager, amrproject::CriterionManager*const ptr_criterion_manager);

    // communication between grid of different refinement levels
    int TransferInfoFromCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const DefInt node_flag_not_interp, const amrproject::GridInfoInterface& grid_info_coarse) override;
    int TransferInfoToCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const DefInt node_flag, amrproject::GridInfoInterface* const ptr_grid_info_coarse) override;

    // transfer node information between fine and coarse grids
    DefInt NodeFlagNotStream_ = 0, NodeFlagNotCollision_ = 0;
    void NodeInfoCoarse2fine(const amrproject::GridNode& coarse_node,
        amrproject::GridNode* const ptr_fine_node) const override;
    void NodeInfoFine2Coarse(const amrproject::GridNode& fine_node,
        amrproject::GridNode* const ptr_coarse_node) const override;

    // output related
    void SetupOutputVariables() override;
    void WriteOutputScalarAndVectors(FILE* const fp, const bool bool_binary,
        const amrproject::Base64Utility& base64_instance,
        const amrproject::OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index) const override;
    std::set<std::string> output_variable_name_ = {"rho", "velocity"};
    int OutputOneVariable(FILE* const fp, const bool bool_binary,
        const amrproject::Base64Utility& base64_instance,
        const OutputLBMNodeVariableInfo& output_info,
        const amrproject::OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index) const;

    // mpi related
    void ComputeInfoInMpiLayers(const std::map<int, DefMap<DefInt>>& map_inner_nodes,
        const DefMap<DefInt>& map_outer_nodes) override;
    void ComputeInfoInInterpMpiLayers(
        const std::map<int, DefMap<DefInt>>& map_intper_nodes) override;
    virtual void ComputeNodeInfoBeforeMpiCommunication(const DefSFBitset sfbitset_in,
        const SolverLbmInterface& lbm_solver);
    virtual void ComputeNodeInfoAfterMpiCommunication(const DefSFBitset sfbitset_in,
        const SolverLbmInterface& lbm_solver);

};
class GridInfoLbmIntefaceCreator :
    public amrproject::GridInfoCreatorInterface {
 public:
    std::shared_ptr<amrproject::GridInfoInterface>
        CreateGridInfo() const override {
        return std::make_shared<GridInfoLbmInteface>();
    };
};
class SolverLbmInterface :public amrproject::SolverInterface {
 protected:
    DefInt k0NumForces_ = 0;
    DefReal k0LbmViscosity_ = 1.;  // the lbm viscosity at background level, where dt_lbm is 1
    DefReal k0Rho_ = 1.;  // default density
    std::vector<DefReal> k0Velocity_, k0Force_;  // velocity and forces for initialization
    std::vector<DefReal> k0ConstForce_;

 public:
    // get and set functions
    void SetDefaultViscosity(const DefReal viscosity_in) { k0LbmViscosity_ = viscosity_in;}
    void SetDefaultDensity(const DefReal rho_in) { k0Rho_ = rho_in;}
    void SetDefaultVelocity(const std::vector<DefReal>& velocity_in) { k0Velocity_ = velocity_in;}
    void SetDefaultForce(const std::vector<DefReal>& force_in);
    void SetConstantForce(const std::vector<DefReal>& force_in);
    DefReal GetDefaultViscosity() const { return k0LbmViscosity_;}
    DefReal GetDefaultDensity() const { return k0Rho_;}
    const std::vector<DefReal>& GetDefaultVelocity() const { return k0Velocity_;}
    const std::vector<DefReal>& GetDefaultForce() const { return k0Force_;}
    const std::vector<DefReal>& GetConstantForce() const { return k0ConstForce_;}
    DefInt GetNumForces() const {
        if (bool_forces_) {
            return k0NumForces_;
        } else {
            return 0;
        }
    }
    void SetNumForces(const DefInt num_forces) {
        bool_forces_ = true;
        k0NumForces_ = num_forces;
        k0Force_.resize(num_forces);
    }

    bool bool_forces_ = false;
    const DefInt k0NumQ_;
    const DefInt k0NumQInOneDirection_ = 0;   ///< number of Q indices in each direction
    const std::vector<DefReal> k0Cx_, k0Cy_, k0Cz_;  ///< directions of particle velocity
    const std::vector<DefReal> k0Weights_;
    // noting that the nth element in k0QIndicesPos_ must be the inverse index of the nth element in k0QIndicesNeg_
    const std::vector<std::vector<DefInt>> k0QIndicesNeg_, k0QIndicesPos_;
    static constexpr DefReal kCs_Reciprocal_ = SqrtConstexpr(3.), kCs_Sq_Reciprocal_ = 3.,
        kCs_ = 1. / SqrtConstexpr(3.), kCs_Sq_ = 1./ 3.;           ///< Lattice sound speed related
    std::string GetSolverMethod() final {
        return  "Solver LBM D" + std::to_string(k0SolverDims_)
            + "Q" + std::to_string(k0NumQ_) + " Model is active.";
    }

    std::function<void(const DefReal,
        LbmCollisionOptInterface* const, GridNodeLbm* const)> func_collision_node_ = nullptr;
    void SolverInitial() override;
    void RunSolverOnGivenGrid(const amrproject::ETimeSteppingScheme time_scheme,
        const DefInt time_step_level, const DefReal time_step_current,
        const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) override;
    void InformationFromGridOfDifferentLevel(
        const DefInt time_step_level, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) override;
    void ReadAndSetupSolverParameters(amrproject::InputParser* const ptr_input_parser) override;

    virtual void InitialModelDependencies() = 0;
    void SetFunctionForNodeCollision();
    template <typename Node>
    void CollisionForGivenNodes(const DefInt i_level, const DefInt flag_not_collide,
        const DefMap<Node>& node_for_computation, GridInfoLbmInteface* const ptr_lbm_grid_nodes_info);
    template <typename Node>
    void StreamInForGivenNodes(const DefInt flag_not_stream, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const DefMap<Node>& node_for_computation, GridInfoLbmInteface* const ptr_lbm_grid_nodes_info);
    virtual void StreamInForAGivenNode(const DefSFBitset sfbitset_in,
        const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const = 0;
    virtual void StreamOutForAGivenNode(const DefSFBitset sfbitset_in,
        const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const = 0;
    virtual void Stream(const DefInt flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const;
    void Collision(const DefInt flag_not_compute, GridInfoLbmInteface* const ptr_grid_info) const;

    // function for calculating initial distribution functions based on macroscopic variables
    void SetInitialDisFuncBasedOnReferenceMacros(
        std::vector<DefReal>* const ptr_f, std::vector<DefReal>* const ptr_f_collide) const;
    void CalInitialDisFuncAsFeq(const DefReal rho, const std::vector<DefReal>& velocity,
        std::vector<DefReal>* const ptr_f, std::vector<DefReal>* const ptr_f_collide) const;

    void SetDefault2DFunctions();
    void SetDefault3DFunctions();
    void ResizeModelRelatedVectors();

    // function to calculate the equilibrium distribution functions
    bool k0BoolCompressible_ = true;  ///< use compressible assumption for feq and velocity computation
    std::function<void(const DefReal rho,
        const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq)> func_cal_feq_;
    void CalFeq2DCompressible(const DefReal rho,
        const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const;
    void CalFeq3DCompressible(const DefReal rho,
        const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const;
    void CalFeq2DIncompressible(const DefReal rho,
        const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const;
    void CalFeq3DIncompressible(const DefReal rho,
        const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) const;

    // function to calculate the LBM forcing term
    virtual std::vector<DefReal> GetAllForcesForANode(const GridNodeLbm& lbm_node) const;
    void (SolverLbmInterface::*ptr_func_cal_force_iq_)(const GridNodeLbm& node,
        const std::vector<DefReal>& force, std::vector<DefReal>* const ptr_forcing_term) const = nullptr;
    void CalForcingTerming2D(const GridNodeLbm& node, const std::vector<DefReal>& force,
        std::vector<DefReal>* const ptr_forcing_term) const;
    void CalForcingTerming3D(const GridNodeLbm& node, const std::vector<DefReal>& force,
        std::vector<DefReal>* const ptr_forcing_term) const;

    // boundary conditions
    void SetDomainBoundaryCondition(const amrproject::EDomainBoundaryDirection which_boundary,
        const ELbmBoundaryConditionScheme which_boundary_condition,
        std::map<amrproject::EDomainBoundaryDirection, std::unique_ptr<BoundaryConditionLbmInterface>>* const
        ptr_boundary_condition) const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryBounceBackCreator() const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryPeriodicCreator() const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryNonEqExtraCreator() const;

    // pointer to function calculate macroscopic variables
    // declare with std::function rather than pointer to function is used for implementations in derived class
    std::function<void(const std::vector<DefReal>&,
        DefReal* const, std::vector<DefReal>* const)> func_macro_without_force_;
    std::function<void(const DefReal dt_lbm, const std::vector<DefReal>&, const std::vector<DefReal>&,
        DefReal* const, std::vector<DefReal>* const)> func_macro_with_force_;
    std::function<void(const DefReal dt_lbm, const GridNodeLbm& lbm_node,
        DefReal* const, std::vector<DefReal>* const)> func_macro_;

    virtual ~SolverLbmInterface() {}

 protected:
    void CalMacro2DCompressible(const std::vector<DefReal>& f,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacro2DIncompressible(const std::vector<DefReal>& f,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce2DCompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
        const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce2DIncompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
        const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacro3DCompressible(const std::vector<DefReal>& f,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacro3DIncompressible(const std::vector<DefReal>& f,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce3DCompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
        const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce3DIncompressible(const DefReal dt_lbm, const std::vector<DefReal>& f,
        const std::vector<DefReal>& force, DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;

    SolverLbmInterface();
    SolverLbmInterface(const DefInt num_q, const DefInt num_q_in_one_direction,
        const std::vector<DefReal>& cx, const std::vector<DefReal>& cy, const std::vector<DefReal>& cz,
        const std::vector<DefReal>& weights,
        const std::vector<std::vector<DefInt>>& q_indices_neg, const std::vector<std::vector<DefInt>>& q_indices_pos)
        : k0NumQ_(num_q), k0NumQInOneDirection_(num_q_in_one_direction),
        k0Cx_(cx), k0Cy_(cy), k0Cz_(cz), k0Weights_(weights),
        k0QIndicesNeg_(q_indices_neg), k0QIndicesPos_(q_indices_pos) {
            ptr_grid_info_creator_ = std::make_unique<GridInfoLbmIntefaceCreator>();
    }

// LES models
 protected:
    bool bool_les_model_ = false;
    std::unique_ptr<LesModelInterface> ptr_les_model_ = std::make_unique<LesModelSmagorinsky>();

 public:
    virtual void SetLesModel(const std::string& model_type);
    virtual void SetLesModel(const ELbmLesModelType model_type);

//  collision models
 protected:
    std::map<DefInt, std::unique_ptr<LbmCollisionOptInterface>> collision_operators_;
    std::vector<ELbmCollisionOperatorType> collision_type_;

 public:
    void ReadCollisionOperator(const std::string&  operator_type);
    void ReadCollisionOperator(const std::vector<std::string>& operator_types);
    void InstantiateCollisionOperator(const DefInt i_level);
    void InstantiateCollisionOperator(const DefInt i_level, const ELbmCollisionOperatorType collision_operator_type);
    void SetCollisionOperator(const DefInt i_level, const ELbmCollisionOperatorType collision_operator_type);
    LbmCollisionOptInterface& GetCollisionOperator(DefInt i_level) const;
    // for MRT
    virtual std::vector<std::vector<DefReal>> GetMrtMMatrix() const = 0;
    virtual std::vector<std::vector<DefReal>> GetMrtImMatrix() const = 0;
    virtual std::vector<std::vector<DefReal>> InitialMrtSMatrix(const DefReal tau) const = 0;
    virtual void UpdateMrtSMatrix(const DefReal tau, std::vector<std::vector<DefReal>>* const ptr_diag_s) const = 0;
    virtual std::vector<std::vector<DefReal>> InitialMrtDMatrix(const DefReal tau) const = 0;
    virtual void UpdateMrtDMatrix(const DefReal tau, std::vector<std::vector<DefReal>>* const ptr_diag_d) const = 0;
};
/**
 * @brief function to conduct collision process for give nodes,
 *          initially designed for separating lbm simulation and mpi communication if possible.
 * @param[in] flag_not_collide flag indicating the node do not need collide.
 * @param[in] node_for_computation nodes to be computed.
 * @param[in] ptr_lbm_grid_nodes_info  pointer to class storing lbm grid nodes information.
 */
template <typename Node>
void SolverLbmInterface::CollisionForGivenNodes(const DefInt i_level, const DefInt flag_not_collide,
    const DefMap<Node>& node_for_computation, GridInfoLbmInteface* const ptr_lbm_grid_nodes_info) {
    LbmCollisionOptInterface& collision_opt = GetCollisionOperator(i_level);
    const DefReal dt_lbm = collision_opt.GetDtLbm();
    DefMap<std::unique_ptr<GridNodeLbm>>& lbm_grid_nodes = *ptr_lbm_grid_nodes_info->GetPtrToLbmGrid();

    if (ptr_lbm_grid_nodes_info->GetPtrToLbmGrid() != nullptr) {
        for (auto& iter_node : node_for_computation) {
            if (lbm_grid_nodes.find(iter_node.first) != lbm_grid_nodes.end()) {
                GridNodeLbm& lbm_node = *lbm_grid_nodes.at(iter_node.first).get();
                if (!(lbm_node.flag_status_ & flag_not_collide)) {
                    this->func_collision_node_(dt_lbm, &collision_opt, &lbm_node);
                }
            }
        }
    }
}
/**
 * @brief function to conduct stream process for give nodes,
 *          initially designed for separating lbm simulation and mpi communication if possible.
 * @param[in] flag_not_stream flag indicating the node do not need stream.
 * @param[in] node_for_computation nodes to be computed.
 * @param[in] ptr_lbm_grid_nodes_info  pointer to class storing lbm grid nodes information.
 */
template <typename Node>
void SolverLbmInterface::StreamInForGivenNodes(const DefInt flag_not_stream,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    const DefMap<Node>& node_for_computation, GridInfoLbmInteface* const ptr_lbm_grid_nodes_info) {
    DefMap<std::unique_ptr<GridNodeLbm>>& lbm_grid_nodes = *ptr_lbm_grid_nodes_info->GetPtrToLbmGrid();
    if (ptr_lbm_grid_nodes_info->GetPtrToLbmGrid() != nullptr) {
        for (auto& iter_node : node_for_computation) {
            if (lbm_grid_nodes.find(iter_node.first) != lbm_grid_nodes.end()) {
                GridNodeLbm& lbm_node = *lbm_grid_nodes.at(iter_node.first).get();
                if (!(lbm_node.flag_status_ & flag_not_stream)) {
                    StreamInForAGivenNode(iter_node.first, sfbitset_aux, &lbm_grid_nodes);
                }
            }
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_LBM_INTERFACE_H_
