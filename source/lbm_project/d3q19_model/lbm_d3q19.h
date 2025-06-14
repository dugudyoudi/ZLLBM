//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d3q19.cpp
* @author Zhengliang Liu
* @brief define class used for LBM D3Q19 model.
* @date  2023-11-26
*/
#ifndef SOURCE_LBM_PROJECT_D3Q19_MODEL_LBM_D3Q19_H_
#define SOURCE_LBM_PROJECT_D3Q19_MODEL_LBM_D3Q19_H_
#include <array>
#include <memory>
#include <vector>
#include "./lbm_interface.h"
namespace rootproject {
namespace lbmproject {
// lattice 0:  (0  0  0)
// lattice 1 : (1  0  0)
// lattice 2 : (-1  0  0)
// lattice 3 : (0  1  0)
// lattice 4 : (0 - 1  0)
// lattice 5 : (0  0  1)
// lattice 6 : (0  0 - 1)

// lattice 7 : (+1 + 1  0)
// lattice 8 : (-1 + 1  0)
// lattice 9 : (+1 - 1  0)
// lattice 10: (-1 - 1  0)

// lattice 11: (+1  0 + 1)
// lattice 12: (-1  0 + 1)
// lattice 13: (+1  0 - 1)
// lattice 14: (-1  0 - 1)

// lattice 15: (0 + 1 + 1)
// lattice 16: (0 - 1 + 1)
// lattice 17: (0 + 1 - 1)
// lattice 18: (0 - 1 - 1)
class SolverLbmD3Q19 :public SolverLbmInterface {
 public:
    // f
    static constexpr DefInt kFX0Y0Z0_ = 0,
        kFXnY0Z0_ = 2, kFXpY0Z0_ = 1,  kFX0YnZ0_ = 4, kFX0YpZ0_ = 3, kFX0Y0Zn_ = 6, kFX0Y0Zp_ = 5,
        kFXnYnZ0_ = 10, kFXnYpZ0_ = 8, kFXpYnZ0_ = 9, kFXpYpZ0_ = 7, kFXnY0Zn_ = 14, kFXnY0Zp_ = 12,
        kFXpY0Zn_ = 13, kFXpY0Zp_ = 11, kFX0YnZn_ = 18, kFX0YnZp_ = 16, kFX0YpZn_ = 17, kFX0YpZp_ = 15;
        /**< indices of distribution functions*/

    void InitialModelDependencies() final;
    void Stream(const DefInt flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
         DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;
    void StreamOutForAGivenNode(const DefSFBitset sfbitset_in, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;
    void StreamInForAGivenNode(const DefSFBitset sfbitset_in, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
        DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const final;

    SolverLbmD3Q19() : SolverLbmInterface(19, 5, InitCx(), InitCy(), InitCz(),
        InitWeights(), IniIndexNeg(), IniIndexPos()) {
        k0SolverDims_ = 3;
        name_ = "LbmD3Q19";
        solver_type_ = "LbmD3Q19";
    }

    // for MRT
    std::vector<std::vector<DefReal>> GetMrtMMatrix() const final;
    std::vector<std::vector<DefReal>> GetMrtImMatrix() const final;
    std::vector<std::vector<DefReal>> InitialMrtSMatrix(const DefReal tau) const final;
    void UpdateMrtSMatrix(const DefReal tau, std::vector<std::vector<DefReal>>* const ptr_diag_s) const final;
    std::vector<std::vector<DefReal>> InitialMrtDMatrix(const DefReal tau) const final;
    void UpdateMrtDMatrix(const DefReal tau, std::vector<std::vector<DefReal>>* const ptr_diag_d) const final;

 private:
    static std::vector<DefReal> InitCx() {
        return { 0., 1., -1., 0.,  0., 0.,  0., 1., -1.,  1., -1., 1., -1.,  1., -1., 0., 0.,   0.,  0.};
    }
    static std::vector<DefReal> InitCy() {
        return { 0., 0.,  0., 1., -1., 0.,  0., 1.,  1., -1., -1., 0.,  0.,  0.,  0., 1., -1.,  1., -1.};
    }
    static std::vector<DefReal> InitCz() {
        return { 0., 0.,  0., 0.,  0., 1., -1., 0.,  0.,  0.,  0., 1.,  1., -1., -1., 1.,  1., -1., -1.};
    }
    static std::vector<DefReal> InitWeights() {
        return { 1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
            1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36. };
    }
    static std::vector<std::vector<DefInt>> IniIndexNeg() {
        return {{kFXnY0Z0_, kFXnYnZ0_, kFXnYpZ0_, kFXnY0Zn_, kFXnY0Zp_},
        {kFX0YnZ0_, kFXnYnZ0_, kFXpYnZ0_, kFX0YnZn_, kFX0YnZp_},
        {kFX0Y0Zn_, kFXnY0Zn_, kFXpY0Zn_, kFX0YnZn_, kFX0YpZn_}};
    }
    static std::vector<std::vector<DefInt>> IniIndexPos() {
        return {{kFXpY0Z0_, kFXpYpZ0_, kFXpYnZ0_, kFXpY0Zp_, kFXpY0Zn_},
            {kFX0YpZ0_, kFXpYpZ0_, kFXnYpZ0_, kFX0YpZp_, kFX0YpZn_},
            {kFX0Y0Zp_, kFXpY0Zp_, kFXnY0Zp_, kFX0YpZp_, kFX0YnZp_}};
    }

    // for MRT
    static constexpr std::array<std::array<DefReal, 19>, 19> kMatrixMMrt_ = {{
        { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., },
        { -30., -11., -11., -11., -11., -11., -11., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., },
        { 12., -4., -4., -4., -4., -4., -4., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., },
        { 0., 1., -1., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0., },
        { 0., -4., 4., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0., },
        { 0., 0., 0., 1., -1., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 1., -1., 1., -1., },
        { 0., 0., 0., -4., 4., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 1., -1., 1., -1., },
        { 0., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., 1., 1., -1., -1., 1., 1., -1., -1., },
        { 0., 0., 0., 0., 0., -4., 4., 0., 0., 0., 0., 1., 1., -1., -1., 1., 1., -1., -1., },
        { 0., 2., 2., -1., -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2., },
        { 0., -4., -4., 2., 2., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2., },
        { 0., 0., 0., 1., 1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0., },
        { 0., 0., 0., -2., -2., 2., 2., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0., },
        { 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., 0., 0., 0., 0., },
        { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., },
        { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., },
        { 0., 0., 0., 0., 0., 0., 0., 1., -1., 1., -1., -1., 1., -1., 1., 0., 0., 0., 0., },
        { 0., 0., 0., 0., 0., 0., 0., -1., -1., 1., 1., 0., 0., 0., 0., 1., -1., 1., -1., },
        { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., -1., -1., 1., 1., }
        }};
  static constexpr std::array<std::array<DefReal, 19>, 19> kMatrixImMrt_ = {{
        {1. / 19, - 5.  / 399,    1. /  21,  0. ,   0.   ,   0. ,   0.   ,   0.,    0.   ,   0.     ,
            0.     ,   0.     ,   0.     ,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19, - 11. / 2394, - 1. /  63,  0.1, - 0.1  ,   0. ,   0.   ,   0.,    0.   ,   1. / 18,
        - 1. / 18,   0.     ,   0.     ,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19, - 11. / 2394, - 1. /  63, -0.1,   0.1  ,   0. ,   0.   ,   0.,    0.   ,   1. / 18,
        - 1. / 18,   0.     ,   0.     ,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19, - 11. / 2394, - 1. /  63,  0. ,   0.   ,   0.1, - 0.1  ,   0.,    0.   , - 1. / 36,
            1. / 36,   1. / 12, - 1. / 12,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19, - 11. / 2394, - 1. /  63,  0. ,   0.   , - 0.1,   0.1  ,   0.,    0.   , - 1. / 36,
            1. / 36,   1. / 12, - 1. / 12,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19, - 11. / 2394, - 1. /  63,  0. ,   0.   ,   0. ,   0.   ,   0.1, - 0.1  , - 1. / 36,
            1. / 36, - 1. / 12,   1. / 12,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19, - 11. / 2394, - 1. /  63,  0. ,   0.   ,   0. ,   0.   , - 0.1,   0.1  , - 1. / 36,
            1. / 36, - 1. / 12,   1. / 12,   0.  ,   0.  ,   0.  ,   0.   ,   0.   ,   0.},
        {1. / 19,    4. / 1197,   1. / 252,  0.1,   0.025,   0.1,   0.025,   0.,    0.   ,   1. / 36,
            1. / 72,   1. / 12,   1. / 24,   0.25,   0.  ,   0.  ,   0.125, - 0.125,   0.},
        {1. / 19,    4. / 1197,   1. / 252, -0.1, - 0.025,   0.1,   0.025,   0.,    0.   ,   1. / 36,
            1. / 72,   1. / 12,   1. / 24, - 0.25,   0.  ,   0.  , - 0.125, - 0.125,   0.},
        {1. / 19,    4. / 1197,   1. / 252,  0.1,   0.025, - 0.1, - 0.025,   0.,    0.   ,   1. / 36,
            1. / 72,   1. / 12,   1. / 24, - 0.25,   0.  ,   0.  ,   0.125,   0.125,   0.},
        {1. / 19,    4. / 1197,   1. / 252, -0.1, - 0.025, - 0.1, - 0.025,   0.,    0.   ,   1. / 36,
            1. / 72,   1. / 12,   1. / 24,   0.25,   0.  ,   0.  , - 0.125,   0.125,   0.},
        {1. / 19,    4. / 1197,   1. / 252,  0.1,   0.025,   0. ,   0.   ,   0.1,   0.025,   1. / 36,
            1. / 72, - 1. / 12, - 1. / 24,   0.  ,   0.  ,   0.25, - 0.125,   0.   ,   0.125},
        {1. / 19,    4. / 1197,   1. / 252, -0.1, - 0.025,   0. ,   0.   ,   0.1,   0.025,   1. / 36,
            1. / 72, - 1. / 12, - 1. / 24,   0.  ,   0.  , - 0.25,   0.125,   0.   ,   0.125},
        {1. / 19,    4. / 1197,   1. / 252,  0.1,   0.025,   0. ,   0.   , - 0.1, - 0.025,   1. / 36,
            1. / 72, - 1. / 12, - 1. / 24,   0.  ,   0.  , - 0.25, - 0.125,   0.   , - 0.125},
        {1. / 19,    4. / 1197,   1. / 252, -0.1, - 0.025,   0. ,   0.   , - 0.1, - 0.025,   1. / 36,
            1. / 72, - 1. / 12, - 1. / 24,   0.  ,   0.  ,   0.25,   0.125,   0.   , - 0.125},
        {1. / 19,    4. / 1197,   1. / 252,  0. ,   0.   ,   0.1,   0.025,   0.1,   0.025, - 1. / 18,
        - 1. / 36,   0.     ,   0.     ,   0.  ,   0.25,   0.  ,   0.   ,   0.125, - 0.125},
        {1. / 19,    4. / 1197,   1. / 252,  0. ,   0.   , - 0.1, - 0.025,   0.1,   0.025, - 1. / 18,
        - 1. / 36,   0.     ,   0.     ,   0.  , - 0.25,   0.  ,   0.   , - 0.125, - 0.125},
        {1. / 19,    4. / 1197,   1. / 252,  0. ,   0.   ,   0.1,   0.025, - 0.1, - 0.025, - 1. / 18,
        - 1. / 36,   0.     ,   0.     ,   0.  , - 0.25,   0.  ,   0.   ,   0.125,   0.125},
        {1. / 19,    4. / 1197,   1. / 252,  0. ,   0.   , - 0.1, - 0.025, - 0.1, - 0.025, - 1. / 18,
        - 1. / 36,   0.     ,   0.     ,   0.  ,   0.25,   0.  ,   0.   , - 0.125,   0.125}
        }};
        std::array<DefReal, 19> k0VectorSMrt_ =  {
            0., 1.5, 1.4, 0., 1.2, 0., 1.2, 0., 1.2, 0., 1.4, 0., 1.4, 0., 0., 0., 1.98, 1.98, 1.98};
        // sk[0] = 0.0; sk[3] = 0.0; sk[5] = 0.0; sk[7] = 0.0;
        // sk[1] = 1.5;  // omega_e
        // sk[2] = 1.4;  // omega_epsilon
        // sk[4] = 1.2; sk[6] = 1.2; sk[8] = 1.2;  // omega_q
        // sk[9] = tau; sk[11] = tau; sk[13] = tau; sk[14] = tau; sk[15] = tau;  // omega_nu
        // sk[10] = 1.4;  sk[12] = 1.4;  // omega_pi
        // sk[16] = 1.98; sk[17] = 1.98; sk[18] = 1.98;  // omega_m
};
class SolverCreatorLbmD3Q19 final :public amrproject::SolverCreatorInterface {
 public:
    std::shared_ptr<amrproject::SolverInterface> CreateSolver() const override;
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_D3Q19_MODEL_LBM_D3Q19_H_
