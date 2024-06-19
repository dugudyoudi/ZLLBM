//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_boundary_conditions.h
* @author Zhengliang Liu
* @brief define classes to manage LBM models.
* @date  2023-9-30
*/
#ifndef ROOTPROJECT_SOURCE_LBM_BOUNDARY_CONDITIONS_H_
#define ROOTPROJECT_SOURCE_LBM_BOUNDARY_CONDITIONS_H_
#include <vector>
#include <array>
#include <map>
#include "../../defs_libs.h"
namespace rootproject {
namespace lbmproject {
class SolverLbmInterface;
/**
* @brief enumerate boundary types
*/
enum class ELbmBoundaryType {
    kUndefined = 0,
    kBoundaryXMin = 1,
    kBoundaryXMax = 2,
    kBoundaryYMin = 3,
    kBoundaryYMax = 4,
    kBoundaryZMin = 5,
    kBoundaryZMax = 6
};
/**
* @brief enumerate boundary conditions
*/
enum class ELbmBoundaryConditionScheme {
    kUndefined = 0,
    kBounceBack = 1,
    kPeriodic = 2
};
class GridInfoLbmInteface;
/**
* @brief interface class to manage LBM boundary conditions 
*/
class BoundaryConditionLbmInterface {
 public:
    ELbmBoundaryConditionScheme boundary_scheme_ = ELbmBoundaryConditionScheme::kUndefined;
    void GetBoundaryNInverseIndices(const ELbmBoundaryType boundary_type,
        const SolverLbmInterface& lbm_solver,
        std::vector<DefAmrIndexUint>* const ptr_indices,
        std::vector<DefAmrIndexUint>* const ptr_inverse_indices) const;
    virtual void CalBoundaryCondition(const ELbmBoundaryType boundary_type,
        const DefMap<DefAmrIndexUint>& boundary_nodes,
        GridInfoLbmInteface* const pr_grid_info) const = 0;
    virtual void SetValues(const std::vector<DefReal> values) = 0;
    virtual ~BoundaryConditionLbmInterface() {}
};
/**
* @brief  class to manage 2D bounce back boundary conditions 
*/
class BoundaryBounceBack2D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const ELbmBoundaryType boundary_type,
        const DefMap<DefAmrIndexUint>& boundary_nodes,
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
    void CalBoundaryCondition(const ELbmBoundaryType boundary_type,
        const DefMap<DefAmrIndexUint>& boundary_nodes,
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
    void CalBoundaryCondition(const ELbmBoundaryType boundary_type,
        const DefMap<DefAmrIndexUint>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override {}
};
/**
* @brief  class to manage 3D periodic boundary conditions 
*/
class BoundaryPeriodic3D : public BoundaryConditionLbmInterface {
 public:
    void CalBoundaryCondition(const ELbmBoundaryType boundary_type,
        const DefMap<DefAmrIndexUint>& boundary_nodes,
        GridInfoLbmInteface* const ptr_grid_info) const override;
    void SetValues(const std::vector<DefReal> values) override {}
};
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_LBM_BOUNDARY_CONDITIONS_H_
