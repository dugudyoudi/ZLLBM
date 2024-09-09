//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_interface.h
* @author Zhengliang Liu
* @brief define classes to manage LBM models.
* @date  2023-9-30
*/
#ifndef SOURCE_LBM_PROJECT_LBM_INTERFACE_H_
#define ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include "../defs_libs.h"
#include "auxiliary_inline_func.h"
#include "solver_info_interface.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
#include "immersed_boundary/immersed_boundary.h"
#include "boundary_conditions/lbm_boundary_conditions.h"
namespace rootproject {
namespace lbmproject {
class SolverLbmInterface;
/**
* @brief enumerate LBM collision operator type
*/
enum class ELbmCollisionOperatorType {
    kUndefined = 0,
    kLbmSrt = 1,
    kLbmMrt = 2,
};
/**
* @brief structure to store node information for LBM simulation
*/
struct GridNodeLbm : public amrproject::GridNode {
    std::vector<DefReal> f_collide_{}, f_{};
    DefReal rho_ = 1.;
    std::vector<DefReal> velocity_{};
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
        : rho_(rho0), velocity_(velocity0), force_(force0), f_(f), f_collide_(f_collide) {
            velocity_.shrink_to_fit();
            force_.shrink_to_fit();
            f_.shrink_to_fit();
            f_collide_.shrink_to_fit();
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
/**
* @brief interface class to manage LBM collision model 
*/
class LbmCollisionOptInterface {
 public:
    DefReal dt_lbm_, viscosity_lbm_, tau_collision_c2f_, tau_collision_f2c_, tau_stream_c2f_, tau_stream_f2c_;
    bool cal_tau_each_node_ = false;
    SolverLbmInterface* ptr_lbm_solver_ = nullptr;
    virtual void CalRelaxationTime() = 0;
    virtual void CalRelaxationTimeNode(const GridNodeLbm& node) = 0;
    virtual void CollisionOperator(const SolverLbmInterface& lbm_solver, GridNodeLbm* const ptr_node) const = 0;
    virtual void CalRelaxationTimeRatio();
    virtual void PostCollisionCoarse2Fine(const DefReal dt_lbm, const std::vector<DefReal>& feq,
        const std::vector<DefReal>& f_collide_coarse, std::vector<DefReal>* const ptr_f_collide_fine);
    virtual void PostCollisionFine2Coarse(const DefReal dt_lbm, const std::vector<DefReal>& feq,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse);
    virtual void PostStreamCoarse2Fine(const DefReal dt_lbm, const std::vector<DefReal>& feq,
        const std::vector<DefReal>& f_coarse, std::vector<DefReal>* const ptr_f_fine);
    virtual void PostStreamFine2Coarse(const DefReal dt_lbm, const std::vector<DefReal>& feq,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse);
    virtual void PostCollisionCoarse2FineForce(const DefReal dt_lbm, const std::vector<DefReal>& feq,
        const SolverLbmInterface& lbm_solver,
        DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
        const GridNodeLbm& node_coarse, std::vector<DefReal>* const ptr_f_collide_fine);
    virtual void PostCollisionFine2CoarseForce(const DefReal dt_lbm, const std::vector<DefReal>& feq,
        const SolverLbmInterface& lbm_solver,
        DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse);

    explicit LbmCollisionOptInterface(const SolverLbmInterface& lbm_solver);
    virtual ~LbmCollisionOptInterface() {}
};
/**
* @brief interface class to manage the SRT collision model
*/
class LbmStrCollisionOpt : public LbmCollisionOptInterface {
 public:
    DefReal tau_srt_;
    void CalRelaxationTime() override;
    void CalRelaxationTimeNode(const GridNodeLbm& node) override;
    void CollisionOperator(const SolverLbmInterface& _solver, GridNodeLbm* const ptr_node) const override;
    explicit LbmStrCollisionOpt(const SolverLbmInterface& lbm_solver)
        : LbmCollisionOptInterface(lbm_solver) {}
};
/**
* @brief interface class to manage the SRT collision model with forcing term
*/
class LbmStrForceCollisionOpt : public LbmStrCollisionOpt {
 public:
    void CollisionOperator(const SolverLbmInterface& lbm_solver, GridNodeLbm* const ptr_node) const override;
    explicit LbmStrForceCollisionOpt(const SolverLbmInterface& lbm_solver)
        : LbmStrCollisionOpt(lbm_solver) {}
};
struct OutputLBMNodeVariableInfo : public amrproject::OutputNodeVariableInfoInterface {
    std::function<std::vector<DefReal>(GridNodeLbm* const)> func_get_var_;
};
/**
* @brief  class to store LBM grid information at each refinement level
*/
class GridInfoLbmInteface : public amrproject::GridInfoInterface {
 public:
    bool bool_forces_ = false;
    std::vector<DefReal> f_ini_, f_collide_ini_;
    void InitialGridInfo() override;

    // convert pointer to current type
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_grid_nodes_ = nullptr;
    void SetPointerToCurrentNodeType() override;
    DefMap<std::unique_ptr<GridNodeLbm>>* GetPointerToLbmGrid();

    // node related
    void InitialNotComputeNodeFlag() override;
    std::unique_ptr<amrproject::GridNode> GridNodeCreator() override;
    void SetNodeVariablesAsZeros(amrproject::GridNode* const ptr_node) override;

    //  collision models
    ELbmCollisionOperatorType k0CollisionOperatorType_ = ELbmCollisionOperatorType::kLbmSrt;
    std::unique_ptr<LbmCollisionOptInterface> ptr_collision_operator_;
    void SetCollisionOperator();
    virtual std::unique_ptr<LbmCollisionOptInterface> LbmSrtCollisionOptCreator(
        const SolverLbmInterface& lbm_solver) const {
        return std::make_unique<LbmStrCollisionOpt>(lbm_solver);
    }
    virtual std::unique_ptr<LbmCollisionOptInterface> LbmSrtForceCollisionOptCreator(
        const SolverLbmInterface& lbm_solver) const {
        return std::make_unique<LbmStrForceCollisionOpt>(lbm_solver);
    }

    void DebugWrite() override;
    void DebugWriteNode(const amrproject::GridNode& node) override;

    // domain boundaries
    /**< pointer to boundary conditions adopted for domain boundary */
    std::map<ELbmBoundaryType, std::unique_ptr<BoundaryConditionLbmInterface>> domain_boundary_condition_;
    bool CheckIfPeriodicDomainRequired(const DefAmrIndexUint dims,
        std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const override;
    void ComputeDomainBoundaryCondition();
    void ComputeDomainBoundaryConditionForANode(int flag_node_boundary, const DefSFBitset bitset_in);

    // time marching related
    void SetUpGridAtBeginningOfTimeStep(const DefAmrIndexUint time_step,
        amrproject::GridManagerInterface* const ptr_grid_manager) override;

    // communication between grid of different refinement levels
    int TransferInfoFromCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const DefAmrUint node_flag_not_interp, const amrproject::GridInfoInterface& grid_info_coarse) override;
    int TransferInfoToCoarseGrid(const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        const DefAmrUint node_flag, amrproject::GridInfoInterface* const ptr_grid_info_coarse) override;

    // transfer node information between fine and coarse grids
    DefAmrUint NodeFlagNotStream_ = 0, NodeFlagNotCollision_ = 0;
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
    int k0SizeOfAllDistributionFunctions_ = 0;
    int SizeOfGridNodeInfoForMpiCommunication() const override;
    int CopyNodeInfoToBuffer(const DefMap<DefAmrIndexUint>& map_nodes, char* const ptr_buffer) override;
    int CopyInterpolationNodeInfoToBuffer(const GridInfoInterface& grid_info_lower,
        const DefMap<DefAmrIndexUint>& map_nodes, char* const ptr_buffer) override;
    int ReadInterpolationNodeInfoFromBuffer(
        const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) override;
    void ComputeInfoInMpiLayers(
        const std::map<int, DefMap<DefAmrIndexUint>>& map_inner_nodes,
        const DefMap<DefAmrIndexUint>& map_outer_nodes) override;
    void ReadNodeInfoFromBuffer(const DefSizet buffer_size, const std::unique_ptr<char[]>& buffer) override;
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
 public:
    DefAmrIndexUint k0NumQ_ = 0;
    std::vector<DefReal> k0Cx_, k0Cy_, k0Cz_;  ///< directions of particle velocity
    std::vector<DefReal> k0Weights_;
    DefAmrIndexUint k0NumQInOneDirection_ = 0;   ///< number of Q indices in one direction
    // noting that the nth element in k0QIndicesPos_ must be the inverse index of the nth element in k0QIndicesNeg_
    std::vector<std::vector<DefAmrIndexUint>> k0QIndicesNeg_, k0QIndicesPos_;
    DefReal k0LbmViscosity_ = 0.;  // the lbm viscosity at background level, where dt_lbm is 1
    DefReal k0Rho_ = 1.;  // default density
    std::vector<DefReal> k0Velocity_, k0Force_;  // default velocity and forces
    static constexpr DefReal kCs_Reciprocal_ = SqrtConstexpr(3.), kCs_Sq_Reciprocal_ = 3.,
        kCs_ = 1. / SqrtConstexpr(3.), kCs_Sq_ = 1./ 3.;           ///< Lattice sound speed related
    std::string GetSolverMethod() final {
        return  "Solver LBM D" + std::to_string(k0SolverDims_)
            + "Q" + std::to_string(k0NumQ_) + " Model is active.";
    }

    SolverLbmInterface() {
        ptr_grid_info_creator_ = std::make_unique<GridInfoLbmIntefaceCreator>();
    }
    void SolverInitial() override;
    void RunSolverOnGrid(const amrproject::ETimeSteppingScheme time_scheme,
        const DefAmrIndexUint time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) override;
    void CallDomainBoundaryCondition(const amrproject::ETimeSteppingScheme time_scheme,
        const DefAmrIndexUint time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) override;
    int InformationFromGridOfDifferentLevel(
        const amrproject::ETimingInOneStep timing, const amrproject::ETimeSteppingScheme time_scheme,
        const DefAmrIndexUint time_step_current, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) override;

    virtual void InitialModelDependencies() = 0;
    virtual void Stream(const DefAmrUint flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const;
    virtual void StreamInForAGivenNode(const DefSFBitset sfbitset_in,
        const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const = 0;
    virtual void StreamOutForAGivenNode(const DefSFBitset sfbitset_in,
        const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const = 0;

    void Collision(const DefAmrUint flag_not_compute, GridInfoLbmInteface* const ptr_grid_info) const;

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

    // function to calculate the LBM body force term
    DefReal (SolverLbmInterface::*ptr_func_cal_force_iq_)(
        const int iq, const GridNodeLbm& node) const = nullptr;
    DefReal CalForceIq2D(const int iq, const GridNodeLbm& node) const;
    DefReal CalForceIq3D(const int iq, const GridNodeLbm& node) const;

    // boundary conditions
    void SetDomainBoundaryCondition(const ELbmBoundaryType which_boundary,
        const ELbmBoundaryConditionScheme which_boundary_condition,
        std::map<ELbmBoundaryType, std::unique_ptr<BoundaryConditionLbmInterface>>* const ptr_boundary_condition) const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryBounceBackCreator() const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryPeriodicCreator() const;

    // pointer to function calculate macroscopic variables
    // declare with std::function rather than pointer to function is used for implementations in derived class
    std::function<void(const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity)> func_macro_without_force_;
    std::function<void(const DefReal dt_lbm, const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity)> func_macro_with_force_;

 protected:
    void CalMacro2DCompressible(const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacro2DIncompressible(const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce2DCompressible(const DefReal dt_lbm, const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce2DIncompressible(const DefReal dt_lbm, const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacro3DCompressible(const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacro3DIncompressible(const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce3DCompressible(const DefReal dt_lbm, const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
    void CalMacroForce3DIncompressible(const DefReal dt_lbm, const GridNodeLbm& node,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
