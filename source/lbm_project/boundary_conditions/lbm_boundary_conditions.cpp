//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_bounce_back.cpp
* @author Zhengliang Liu
* @brief functions used to implement bounce back boundary condition.
* @date  2023-11-6
*/
#include <memory>
#include <string>
#include "boundary_conditions/lbm_boundary_conditions.h"
#include "io/log_write.h"
#include "lbm_interface.h"
namespace rootproject {
namespace lbmproject {
 /**
 * @brief function to get indices of distribution functions need to compute at the boundary and their inverses.
 * @param[in] which_boundary reuse enum boundary type to indicate condition implemented on which boundary.
 * @param[in] which_boundary_condition implemented boundary condition.
 * @param[out] ptr_boundary_condition pointer to boundary conditions.
 */   
void SolverLbmInterface::SetDomainBoundaryCondition(
    const ELbmBoundaryType which_boundary, const ELbmBoundaryConditionScheme which_boundary_condition,
    std::map<ELbmBoundaryType, std::unique_ptr<BoundaryConditionLbmInterface>>* const ptr_boundary_condition) const {
    switch (which_boundary_condition) {
    case ELbmBoundaryConditionScheme::kBounceBack: {
            if (ptr_boundary_condition->find(which_boundary) == ptr_boundary_condition->end()) {
                ptr_boundary_condition->insert({which_boundary, BoundaryBounceBackCreator()});
                ptr_boundary_condition->at(which_boundary)->boundary_scheme_ = ELbmBoundaryConditionScheme::kBounceBack;
            }
        }
        break;
    case ELbmBoundaryConditionScheme::kPeriodic: {
            if (ptr_boundary_condition->find(which_boundary) == ptr_boundary_condition->end()) {
                ptr_boundary_condition->insert({which_boundary, BoundaryPeriodicCreator()});
                ptr_boundary_condition->at(which_boundary)->boundary_scheme_ = ELbmBoundaryConditionScheme::kPeriodic;
            }
        }
        break;
    default:
        break;
    }
}
/**
 * @brief function to create instance for bounce back boundary condition.
 */
std::unique_ptr<BoundaryConditionLbmInterface> SolverLbmInterface::BoundaryBounceBackCreator() const {
    if (k0SolverDims_ == 2) {
        return std::make_unique<BoundaryBounceBack2D>();
    } else if (k0SolverDims_ == 3) {
        return std::make_unique<BoundaryBounceBack3D>();
    } else {
        amrproject::LogManager::LogError("Dimension of the solver should be 2 or 3"
            " for instantiating bounce back boundary condition in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return nullptr;
    }
}
/**
 * @brief function to create instance for periodic boundary condition.
 */
std::unique_ptr<BoundaryConditionLbmInterface> SolverLbmInterface::BoundaryPeriodicCreator() const {
    if (k0SolverDims_ == 2) {
        return std::make_unique<BoundaryPeriodic2D>();
    } else if (k0SolverDims_ == 3) {
        return std::make_unique<BoundaryPeriodic3D>();
    } else {
        amrproject::LogManager::LogError("Dimension of the solver should be 2 or 3"
            " for instantiating periodic boundary condition in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return nullptr;
    }
}
/**
 * @brief function to get indices of distribution functions need to compute at the boundary and their inverses.
 * @param[in] boundary_type type of the boundary.
 * @param[in] lbm_solver reference to LBM solver.
 * @param[out] ptr_indices indices of distribution functions need to compute.
 * @param[out] ptr_inverse_indices inverse indices.
 */
void BoundaryConditionLbmInterface::GetBoundaryNInverseIndices(
    const ELbmBoundaryType boundary_type,
    const SolverLbmInterface& lbm_solver,
    std::vector<DefAmrIndexUint>* const ptr_indices,
    std::vector<DefAmrIndexUint>* const ptr_inverse_indices) const {
    switch (boundary_type) {
        case ELbmBoundaryType::kBoundaryXMin:
            // nodes at negative x boundary need to compute distribution functions at positive x
            *ptr_indices = lbm_solver.k0QIndicesPos_.at(kXIndex);
            *ptr_inverse_indices = lbm_solver.k0QIndicesNeg_.at(kXIndex);
            ptr_indices->shrink_to_fit();
            ptr_inverse_indices->shrink_to_fit();
            break;
        case ELbmBoundaryType::kBoundaryXMax:
            *ptr_indices = lbm_solver.k0QIndicesNeg_.at(kXIndex);
            *ptr_inverse_indices = lbm_solver.k0QIndicesPos_.at(kXIndex);
            ptr_indices->shrink_to_fit();
            ptr_inverse_indices->shrink_to_fit();
            break;
        case ELbmBoundaryType::kBoundaryYMin:
            *ptr_indices = lbm_solver.k0QIndicesPos_.at(kYIndex);
            *ptr_inverse_indices = lbm_solver.k0QIndicesNeg_.at(kYIndex);
            ptr_indices->shrink_to_fit();
            ptr_inverse_indices->shrink_to_fit();
            break;
        case ELbmBoundaryType::kBoundaryYMax:
            *ptr_indices = lbm_solver.k0QIndicesNeg_.at(kYIndex);
            *ptr_inverse_indices = lbm_solver.k0QIndicesPos_.at(kYIndex);
            ptr_indices->shrink_to_fit();
            ptr_inverse_indices->shrink_to_fit();
            break;
        case ELbmBoundaryType::kBoundaryZMin:
            if (lbm_solver.k0SolverDims_ == 2) {
                amrproject::LogManager::LogError("Dimension for LBM solver is 2, does not support "
                    " boundary in z direction in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            *ptr_indices = lbm_solver.k0QIndicesPos_.at(kZIndex);
            *ptr_inverse_indices = lbm_solver.k0QIndicesNeg_.at(kZIndex);
            break;
        case ELbmBoundaryType::kBoundaryZMax:
            if (lbm_solver.k0SolverDims_ == 2) {
                amrproject::LogManager::LogError("Dimension for LBM solver is 2, does not support "
                    " boundary in z direction in "
                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            }
            *ptr_indices = lbm_solver.k0QIndicesNeg_.at(kZIndex);
            *ptr_inverse_indices = lbm_solver.k0QIndicesPos_.at(kZIndex);
            break;
        default:
            amrproject::LogManager::LogError("Boundary type is not defined in "
                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            break;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject