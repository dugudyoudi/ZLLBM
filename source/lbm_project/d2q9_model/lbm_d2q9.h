//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d2q9.cpp
* @author Zhengliang Liu
* @brief define class used for LBM D2Q9 model.
* @date  2023-9-30
*/
#ifndef SOURCE_LBM_PROJECT_D2Q9_MODEL_LBM_D2Q9_H_
#define SOURCE_LBM_PROJECT_D2Q9_MODEL_LBM_D2Q9_H_
#include <array>
#include <memory>
#include <vector>
#include "d2q9_model/boundary_d2q9.h"
#include "./lbm_interface.h"
namespace rootproject {
namespace lbmproject {
class SolverLbmD2Q9 :public SolverLbmInterface {
 public:
    // f
    static constexpr DefInt kFX0Y0Z0_ = 0, kFXnY0Z0_ = 3, kFXpY0Z0_ = 1, kFX0YnZ0_ = 4, kFX0YpZ0_ = 2,
        kFXnYnZ0_ = 7, kFXnYpZ0_ = 6, kFXpYnZ0_ = 8, kFXpYpZ0_ = 5;
        /**< indices of distribution functions*/
    static constexpr std::array<std::array<DefReal, 9>, 9> kMatrixMMrt_ = {{
        {  1.,  1.,  1.,  1.,  1., 1.,  1.,  1.,  1. },
        { -4., -1., -1., -1., -1., 2.,  2.,  2.,  2. },
        {  4., -2., -2., -2., -2., 1.,  1.,  1.,  1. },
        {  0.,  1.,  0., -1.,  0., 1., -1., -1.,  1. },
        {  0., -2.,  0.,  2.,  0., 1., -1., -1.,  1. },
        {  0.,  0.,  1.,  0., -1., 1.,  1., -1., -1. },
        {  0.,  0., -2.,  0.,  2., 1.,  1., -1., -1. },
        {  0.,  1., -1.,  1., -1., 0.,  0.,  0.,  0. },
        {  0.,  0.,  0.,  0.,  0., 1., -1.,  1., -1. }
        }};
    static constexpr std::array<std::array<DefReal, 9>, 9> kMatrixImMrt_ = {{
        { 1. / 9., -1. / 9. ,  1. / 9. ,  0.     ,  0.      ,  0.     ,  0.      ,  0.     ,  0.      },
        { 1. / 9., -1. / 36., -1. / 18.,  1. / 6., -1. / 6. ,  0.     ,  0.      ,  1. / 4.,  0.      },
        { 1. / 9., -1. / 36., -1. / 18.,  0.     ,  0.      ,  1. / 6., -1. / 6. , -1. / 4.,  0.      },
        { 1. / 9., -1. / 36., -1. / 18., -1. / 6.,  1. / 6. ,  0.     ,  0.      ,  1. / 4.,  0.      },
        { 1. / 9., -1. / 36., -1. / 18.,  0.     ,  0.      , -1. / 6.,  1. / 6. , -1. / 4.,  0.      },
        { 1. / 9.,  1. / 18.,  1. / 36.,  1. / 6.,  1. / 12.,  1. / 6.,  1. / 12.,  0.     ,  1. / 4. },
        { 1. / 9.,  1. / 18.,  1. / 36., -1. / 6., -1. / 12.,  1. / 6.,  1. / 12.,  0.     , -1. / 4. },
        { 1. / 9.,  1. / 18.,  1. / 36., -1. / 6., -1. / 12., -1. / 6., -1. / 12.,  0.     ,  1. / 4. },
        { 1. / 9.,  1. / 18.,  1. / 36.,  1. / 6.,  1. / 12., -1. / 6., -1. / 12.,  0.     , -1. / 4. }
        }};
    void InitialModelDependencies() final;
    void Stream(const DefInt flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
         DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;
    void StreamOutForAGivenNode(const DefSFBitset sfbitset_in, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;
    void StreamInForAGivenNode(const DefSFBitset sfbitset_in, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;

    // std::unique_ptr<BoundaryConditionLbmInterface> BoundaryBounceBackCreator() const override {
    //     return std::make_unique<BoundaryBounceBackD2Q9>();
    // }

    SolverLbmD2Q9() : SolverLbmInterface(9, 3, InitCx(), InitCy(), {}, InitWeights(), IniIndexNeg(), IniIndexPos()) {
        k0SolverDims_ = 2;
    }

 private:
    static std::vector<DefReal> InitCx() {
        return { 0., 1., 0., -1., 0., 1., -1., -1., 1. };
    }
    static std::vector<DefReal> InitCy() {
        return { 0., 0., 1., 0., -1., 1., 1., -1., -1. };
    }
    static std::vector<DefReal> InitWeights() {
        return { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
    }
    static std::vector<std::vector<DefInt>> IniIndexNeg() {
        return {{ kFXnY0Z0_, kFXnYnZ0_, kFXnYpZ0_ }, {kFX0YnZ0_, kFXnYnZ0_, kFXpYnZ0_}};
    }
    static std::vector<std::vector<DefInt>> IniIndexPos() {
        return {{ kFXpY0Z0_, kFXpYnZ0_, kFXpYpZ0_ }, {kFX0YpZ0_, kFXpYpZ0_, kFXnYpZ0_}};
    }

 protected:
    void CalMacroD2Q9Incompressible(const std::vector<DefReal>& f,
        DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) const;
};
class SolverCreatorLbmD2Q9 final :public amrproject::SolverCreatorInterface {
 public:
    std::shared_ptr<amrproject::SolverInterface> CreateSolver() const override;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_D2Q9_MODEL_LBM_D2Q9_H_

