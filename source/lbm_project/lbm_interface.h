//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_interface.h
* @author Zhengliang Liu
* @brief define classes to manage LBM models.
* @date  2023-9-30
*/
#ifndef ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
#define ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include "../defs_libs.h"
#include "auxiliary_inline_func.h"
#include "solver_info_interface.h"
#include "grid/sfbitset_aux.h"
#include "grid/grid_info_interface.h"
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
    std::vector<DefReal> f_collide_, f_;
    DefReal rho_;
    std::vector<DefReal> velocity_;
    std::vector<DefReal> force_;
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
};
/**
* @brief interface class to manage LBM collision model 
*/
class LbmCollisionOptInterface {
 public:
    DefReal dt_lbm_, viscosity_lbm_;
    bool cal_tau_each_node_ = false;
    virtual void CalRelaxationTime() = 0;
    virtual void CalRelaxationTimeNode(const GridNodeLbm& node) = 0;
    virtual void CollisionOperator(const SolverLbmInterface& lbm_solver, GridNodeLbm* const ptr_node) const = 0;
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
struct OutputLBMNodeVariableInfo {
    std::function<std::vector<DefReal>(const GridNodeLbm&)> func_get_var_;
    DefAmrIndexUint variable_dims_;
    std::string output_name_;
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
    DefMap<std::unique_ptr<GridNodeLbm>>* ptr_lbm_grid_ = nullptr;
    void SetPointerToCurrentNodeType() override;
    DefMap<std::unique_ptr<GridNodeLbm>>* GetPointerToLbmGrid();

    // node related
    std::unique_ptr<amrproject::GridNode> GridNodeCreator() override;

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

    // domain boundaries
    /**< pointer to boundary conditions adopted for domain boundary */
    std::map<ELbmBoundaryType, std::unique_ptr<BoundaryConditionLbmInterface>> domain_boundary_condition_;
    void ComputeDomainBoundaryCondition();

    // output related
    void SetupOutputVariables() override;
    void WriteOutputScalarAndVectors(FILE* const fp, const bool bool_binary,
        const amrproject::Base64Utility& base64_instance,
        const amrproject::OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index) const override;
    std::set<std::string> output_variable_name_ = {"rho", "velocity"};
    std::vector<OutputLBMNodeVariableInfo> output_variables_;
    int OutputOneVariable(FILE* const fp, const bool bool_binary,
        const amrproject::Base64Utility& base64_instance,
        const OutputLBMNodeVariableInfo& output_info,
        const amrproject::OutputDataFormat& output_data_format,
        const DefMap<DefSizet>& map_node_index) const;
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
    static constexpr DefReal kCsReciprocal = SqrtConstexpr(3.), kCsSqReciprocal = 3.,
        kCs = 1. / SqrtConstexpr(3.), kCsSq = 1./ 3.;           ///< Lattice sound speed related
    std::string GetSolverMethod() final {
        return  "Solver LBM D" + std::to_string(k0SolverDims_)
            + "Q" + std::to_string(k0NumQ_) + " Model is active.";
    }

    SolverLbmInterface() {
        ptr_grid_info_creator_ = std::make_unique<GridInfoLbmIntefaceCreator>();
    }
    void SolverInitial() override;
    void SetNodeFlagForSolver() override;
    void RunSolver(const DefReal time_step_current,
        const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        amrproject::GridInfoInterface* const ptr_grid_info) override;

    virtual void InitialModelDependencies() = 0;
    virtual void Stream(const DefAmrUint flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
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

    // function to conduct conversion for distribution functions between fine and coarse grids
    void (SolverLbmInterface::*ptr_func_coarse2fine_)(
        const DefReal dt_lbm, const std::vector<std::vector<DefReal>>& matrix_c2f,
        const GridNodeLbm& node_coarse, GridNodeLbm* const ptr_node_fine) = nullptr;
    void Coarse2FineSrt(const DefReal dt_lbm, const std::vector<std::vector<DefReal>>& matrix_c2f,
        const GridNodeLbm& node_coarse, GridNodeLbm* const ptr_node_fine);
    void (SolverLbmInterface::*ptr_func_fine2coarse_)(
        const DefReal dt_lbm, const std::vector<std::vector<DefReal>>& matrix_f2c,
        const GridNodeLbm& node_coarse, GridNodeLbm* const ptr_node_fine) = nullptr;
    void Fine2CoarseSrt(const DefReal dt_lbm, const std::vector<std::vector<DefReal>>& matrix_f2c,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse);

    // boundary conditions
    void SetDomainBoundaryCondition(const ELbmBoundaryType which_boundary,
        const ELbmBoundaryConditionScheme which_boundary_condition,
        std::map<ELbmBoundaryType, std::unique_ptr<BoundaryConditionLbmInterface>>* const ptr_boundary_condition) const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryBounceBackCreator() const;
    virtual std::unique_ptr<BoundaryConditionLbmInterface> BoundaryPeriodicCreator() const;

 protected:
    DefAmrUint NodeFlagNotCompute_;
    inline DefReal Square(DefReal x) const {return x * x;}

    void CalMacro2DCompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacro2DIncompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroForce2DCompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroForce2DIncompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacro3DCompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacro3DIncompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroForce3DCompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroForce3DIncompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;

    std::function<void(const DefReal, GridNodeLbm* const)> func_macro_, func_macro_force_;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_LBM_INTERFACE_H_
