//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_boundary_conditions.h
* @author Zhengliang Liu
* @brief define classes to manage LBM models.
* @date  2023-9-30
*/
#ifndef SOURCE_LBM_PROJECT_BOUNDARY_CONDITIONS_LBM_BOUNDARY_CONDITIONS_H_
#define SOURCE_LBM_PROJECT_BOUNDARY_CONDITIONS_LBM_BOUNDARY_CONDITIONS_H_
#include <vector>
#include <array>
#include <map>
#include "grid/grid_enumerates.h"
namespace rootproject {
namespace lbmproject {
class SolverLbmInterface;
/**
* @brief enumerate boundary conditions
*/
enum class ELbmBoundaryConditionScheme {
    kUndefined = 0,
    kBounceBack = 1,
    kPeriodic = 2,
    kNonEqExtrapolation = 3
};
class GridInfoLbmInteface;
/**
* @brief interface class to manage LBM boundary conditions 
*/
class BoundaryConditionLbmInterface {
 public:
    ELbmBoundaryConditionScheme boundary_scheme_ = ELbmBoundaryConditionScheme::kUndefined;
    void GetBoundaryNInverseIndices(const amrproject::EDomainBoundaryDirection boundary_dir,
        const SolverLbmInterface& lbm_solver,
        std::vector<DefInt>* const ptr_indices,
        std::vector<DefInt>* const ptr_inverse_indices) const;
    virtual void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const pr_grid_info) const = 0;
    virtual void SetValues(const std::vector<DefReal> values) = 0;
    virtual ~BoundaryConditionLbmInterface() {}
};
/**
* @brief  class to manage 2D bounce back boundary conditions 
*/
class BoundaryBounceBack2D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override;

 protected:
    std::array<DefReal, 2> boundary_velocity_ = {0., 0.};
};
/**
* @brief  class to manage 3D bounce back boundary conditions 
*/
class BoundaryBounceBack3D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override;

 protected:
    std::array<DefReal, 3> boundary_velocity_ = {0., 0., 0.};
};
/**
* @brief  class to manage 2D periodic boundary conditions 
*/
class BoundaryPeriodic2D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override {}
};
/**
* @brief  class to manage 3D periodic boundary conditions 
*/
class BoundaryPeriodic3D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override {}
};
/**
* @brief  class to manage 2D non-equilibrium extrapolation boundary conditions 
*/
class BoundaryNonEqExtrapolation2D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override;

 protected:
    std::array<DefReal, 2> boundary_velocity_ = {0., 0.};
};
/**
* @brief  class to manage 3D non-equilibrium extrapolation boundary conditions 
*/
class BoundaryNonEqExtrapolation3D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const amrproject::EDomainBoundaryDirection boundary_dir,
        const DefMap<DefInt>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override;

protected:
    std::array<DefReal, 3> boundary_velocity_ = {0., 0., 0.};
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_BOUNDARY_CONDITIONS_LBM_BOUNDARY_CONDITIONS_H_
