//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d2q9.cpp
* @author Zhengliang Liu
* @brief define class used for LBM D2Q9 model.
* @date  2023-9-30
*/
#ifndef ROOTPROJECT_SOURCE_LBM_LBM_D2Q9_H_
#define ROOTPROJECT_SOURCE_LBM_LBM_D2Q9_H_
#include <array>
#include <memory>
#include "d2q9_model/boundary_d2q9.h"
#include "lbm_interface.h"
namespace rootproject {
namespace lbmproject {
class SolverLbmD2Q9 :public SolverLbmInterface {
 public:
    // f
    static constexpr DefAmrIndexUint kFX0Y0Z0 = 0, kFXnY0Z0 = 3, kFXpY0Z0 = 1, kFX0YnZ0 = 4, kFX0YpZ0 = 2,
        kFXnYnZ0 = 7, kFXnYpZ0 = 6, kFXpYnZ0 = 8, kFXpYpZ0 = 5;
        /**< indices of distribution functions*/
    static constexpr std::array<std::array<DefReal, 9>, 9> kMatrixMMrt = {{
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
    static constexpr std::array<std::array<DefReal, 9>, 9> kMatrixImMrt = {{
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
    void Stream(const DefAmrUint flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
         DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;
    void StreamForAGivenNode(const DefSFBitset sfbitset_in, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;

    // std::unique_ptr<BoundaryConditionLbmInterface> BoundaryBounceBackCreator() const override {
    //     return std::make_unique<BoundaryBounceBackD2Q9>();
    // }

 protected:
    void CalMacroD2Q9Compressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroD2Q9Incompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroForceD2Q9Compressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
    void CalMacroForceD2Q9Incompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const;
};
class SolverCreatorLbmD2Q9 final :public amrproject::SolverCreatorInterface {
 public:
    std::shared_ptr<amrproject::SolverInterface> CreateSolver() const override;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_LBM_D2Q9_H_
