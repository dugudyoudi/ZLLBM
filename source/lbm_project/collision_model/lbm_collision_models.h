//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_collision_models.h
* @author Zhengliang Liu
* @brief define classes to manage LBM collision models.
* @date  2023-10-30
*/
#ifndef SOURCE_LBM_PROJECT_COLLISION_MODEL_LBM_COLLISION_MODELS_H_
#define SOURCE_LBM_PROJECT_COLLISION_MODEL_LBM_COLLISION_MODELS_H_
#include <vector>
#include "../defs_libs.h"
namespace rootproject {
namespace lbmproject {
struct GridNodeLbm;
class SolverLbmInterface;
/**
* @brief enumerate LBM collision operator type
*/
enum class ELbmCollisionOperatorType {
    kLbmSrt,
    kLbmMrt,
    kUndefined  // need to be the last which is used to setup integrated test using gtest
};
/**
* @brief interface class to manage LBM collision model 
*/
class LbmCollisionOptInterface {
 protected:
    DefReal tau0_, tau_eff_;
    DefReal dt_lbm_, viscosity_lbm_, tau_collision_c2f_, tau_collision_f2c_, tau_stream_c2f_, tau_stream_f2c_;

 public:
    DefReal GetDtLbm() const { return dt_lbm_; }
    DefReal GetViscosityLbm() const { return viscosity_lbm_; }
    DefReal GetRelaxationTime() const { return tau0_; }
    DefReal GetEffectiveRelaxationTime() const { return tau0_; }

    bool cal_tau_each_node_ = false;
    virtual void CalRelaxationTime(const DefReal vis_lbm, const SolverLbmInterface& lbm_solver);
    virtual void SetEffectiveRelaxationTime(const DefReal tau_eff, const SolverLbmInterface& lbm_solver);
    virtual void SetEffectiveRelaxationTimeForForcing(const DefReal tau_eff, const SolverLbmInterface& lbm_solver);
    virtual void CalRelaxationTimeRatio(const DefReal vis_lbm, const SolverLbmInterface& lbm_solver);
    virtual void PostCollisionCoarse2Fine(const std::vector<DefReal>& feq,
        const std::vector<DefReal>& f_collide_coarse, std::vector<DefReal>* const ptr_f_collide_fine) const;
    virtual void PostCollisionFine2Coarse(const std::vector<DefReal>& feq,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) const;
    virtual void PostStreamCoarse2Fine(const std::vector<DefReal>& feq,
        const std::vector<DefReal>& f_coarse, std::vector<DefReal>* const ptr_f_fine) const;
    virtual void PostStreamFine2Coarse(const std::vector<DefReal>& feq,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) const;
    virtual void PostCollisionCoarse2FineForce(const std::vector<DefReal>& feq,
        const SolverLbmInterface& lbm_solver,
        DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
        const GridNodeLbm& node_coarse, std::vector<DefReal>* const ptr_f_collide_fine) const;
    virtual void PostCollisionFine2CoarseForce(const std::vector<DefReal>& feq,
        const SolverLbmInterface& lbm_solver,
        DefReal (SolverLbmInterface::*ptr_func_cal_force_iq)(const int, const GridNodeLbm&) const,
        const GridNodeLbm& node_fine, GridNodeLbm* const ptr_node_coarse) const;
    virtual void CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
        const std::vector<DefReal>& force, GridNodeLbm* const ptr_node) const = 0;
    virtual void CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
        GridNodeLbm* const ptr_node) const = 0;

    LbmCollisionOptInterface(const DefInt i_level, const SolverLbmInterface& lbm_solver);
    virtual ~LbmCollisionOptInterface() {}
};
/**
* @brief interface class to manage the SRT collision model
*/
class LbmStrCollisionOpt : public LbmCollisionOptInterface {
 public:
    void CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
        const std::vector<DefReal>& force, GridNodeLbm* const ptr_node) const override;
    void CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
        GridNodeLbm* const ptr_node) const override;
    LbmStrCollisionOpt(const DefInt i_level, const SolverLbmInterface& lbm_solver)
        : LbmCollisionOptInterface(i_level, lbm_solver) {
        CalRelaxationTime(viscosity_lbm_, lbm_solver);
        CalRelaxationTimeRatio(viscosity_lbm_, lbm_solver);
    }
};
/**
* @brief interface class to manage the SRT collision model
*/
class LbmMrtCollisionOpt : public LbmCollisionOptInterface {
 public:
    void CalRelaxationTime(const DefReal vis_lbm, const SolverLbmInterface& lbm_solver) override;
    void SetEffectiveRelaxationTime(const DefReal tau_eff, const SolverLbmInterface& lbm_solver) override;
    void SetEffectiveRelaxationTimeForForcing(const DefReal tau_eff, const SolverLbmInterface& lbm_solver) override;
    void CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
        const std::vector<DefReal>& force, GridNodeLbm* const ptr_node) const override;
    void CollisionOperator(const DefInt num_q, const std::vector<DefReal>& feq,
        GridNodeLbm* const ptr_node) const override;
    LbmMrtCollisionOpt(const DefInt i_level, const SolverLbmInterface& lbm_solver);

 protected:
    void SetImSMMatrix(const DefReal relax_tau, const SolverLbmInterface& lbm_solver);
    void SetImDMMatrix(const DefReal relax_tau, const SolverLbmInterface& lbm_solver);
    std::vector<std::vector<DefReal>> matrix_m_, matrix_im_, diag_s_, matrix_im_s_m_;
    std::vector<std::vector<DefReal>> diag_d_, matrix_im_d_m_;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_COLLISION_MODEL_LBM_COLLISION_MODELS_H_
