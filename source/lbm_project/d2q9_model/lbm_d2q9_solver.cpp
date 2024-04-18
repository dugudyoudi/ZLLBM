//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d2q9_solver.cpp
* @author Zhengliang Liu
* @brief alternative functions for LBM simulation using D2Q9 model.
* @date  2023-9-30
*/
#include <string>
#include "lbm_d2q9.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to perform propagation from other nodes to the given node in the LBM simulation.
 * @param[in] flag_not_compute flag indicating whether to compute or not.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD2Q9::Stream(const DefAmrUint flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    if (ptr_map_grid_nodes != nullptr) {
        DefSFBitset sfbitset_tmp, sfbitset_tmp1;
        const amrproject::SFBitsetAux2D sfbitset_aux2d = dynamic_cast<const amrproject::SFBitsetAux2D&>(sfbitset_aux);
        for (auto& iter_node : *ptr_map_grid_nodes) {
            if (iter_node.second->flag_status_ & flag_not_compute) {
            } else {
                // f(0, 0)
                iter_node.second->f_.at(kFX0Y0Z0) = iter_node.second->f_collide_.at(kFX0Y0Z0);
                // f(-x, 0)
                sfbitset_tmp = sfbitset_aux2d.FindXNeg(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXnY0Z0) = iter_node.second->f_collide_.at(kFXnY0Z0);
                }
                // f(-x, -y)
                sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYnZ0) = iter_node.second->f_collide_.at(kFXnYnZ0);
                }
                // f(-x, +y)
                sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYpZ0) = iter_node.second->f_collide_.at(kFXnYpZ0);
                }
                // f(+x, 0)
                sfbitset_tmp = sfbitset_aux2d.FindXPos(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXpY0Z0) = iter_node.second->f_collide_.at(kFXpY0Z0);
                }
                // f(+x, -y)
                sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYnZ0) = iter_node.second->f_collide_.at(kFXpYnZ0);
                }
                // f(+x, +y)
                sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYpZ0) = iter_node.second->f_collide_.at(kFXpYpZ0);
                }
                // f(0, -y)
                sfbitset_tmp = sfbitset_aux2d.FindYNeg(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YnZ0) = iter_node.second->f_collide_.at(kFX0YnZ0);
                }
                // f(0, +y)
                sfbitset_tmp = sfbitset_aux2d.FindYPos(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YpZ0) = iter_node.second->f_collide_.at(kFX0YpZ0);
                }
            }
        }
    }
}
/**
 * @brief function to calculate macroscopic variables based on D2Q9 distribution functions.
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param ptr_node pointer to a grid node storing LBM related information.
 */
void SolverLbmD2Q9::CalMacroD2Q9Compressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const {
    ptr_node->rho_ = ptr_node->f_[kFX0Y0Z0] + ptr_node->f_[kFXnY0Z0] + ptr_node->f_[kFXpY0Z0] + ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFX0YpZ0] + ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXnYpZ0]
        + ptr_node->f_[kFXpYnZ0] + ptr_node->f_[kFXpYpZ0];
    ptr_node->velocity_[kXIndex] = (ptr_node->f_[kFXpY0Z0] - ptr_node->f_[kFXnY0Z0]
        + ptr_node->f_[kFXpYpZ0] - ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kYIndex] = (ptr_node->f_[kFX0YpZ0] - ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFXpYpZ0] + ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] - ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kXIndex]/=ptr_node->rho_;
    ptr_node->velocity_[kYIndex]/=ptr_node->rho_;
}
/**
 * @brief function to calculate macroscopic variables based on D2Q9 distribution functions (incompressible).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param ptr_node pointer to a grid node storing LBM related information.
 */
void SolverLbmD2Q9::CalMacroD2Q9Incompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const {
    ptr_node->rho_ = ptr_node->f_[kFX0Y0Z0] + ptr_node->f_[kFXnY0Z0] + ptr_node->f_[kFXpY0Z0] + ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFX0YpZ0] + ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXnYpZ0]
        + ptr_node->f_[kFXpYnZ0] + ptr_node->f_[kFXpYpZ0];
    ptr_node->velocity_[kXIndex] = (ptr_node->f_[kFXpY0Z0] - ptr_node->f_[kFXnY0Z0]
        + ptr_node->f_[kFXpYpZ0] - ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kYIndex] = (ptr_node->f_[kFX0YpZ0] - ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFXpYpZ0] + ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] - ptr_node->f_[kFXpYnZ0]);
}
/**
 * @brief function to calculate macroscopic variables based on D2Q9 distribution functions (with forcing term).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param ptr_node pointer to a grid node storing LBM related information.
 */
void SolverLbmD2Q9::CalMacroForceD2Q9Compressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const {
    ptr_node->rho_ = ptr_node->f_[kFX0Y0Z0] + ptr_node->f_[kFXnY0Z0] + ptr_node->f_[kFXpY0Z0] + ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFX0YpZ0] + ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXnYpZ0]
        + ptr_node->f_[kFXpYnZ0] + ptr_node->f_[kFXpYpZ0];
    ptr_node->velocity_[kXIndex] = (ptr_node->f_[kFXpY0Z0] - ptr_node->f_[kFXnY0Z0]
        + ptr_node->f_[kFXpYpZ0] - ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kYIndex] = (ptr_node->f_[kFX0YpZ0] - ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFXpYpZ0] + ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] - ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kXIndex] += 0.5 * ptr_node->force_[kXIndex] * dt_lbm;
    ptr_node->velocity_[kYIndex] += 0.5 * ptr_node->force_[kYIndex] * dt_lbm;
    ptr_node->velocity_[kXIndex]/=ptr_node->rho_;
    ptr_node->velocity_[kYIndex]/=ptr_node->rho_;
}
/**
 * @brief function to calculate macroscopic variables based on D2Q9 distribution functions
 *        (with forcing term and incompressible).
 * @param[in] dt_lbm time spacing of LBM at current refinement level.
 * @param ptr_node pointer to a grid node storing LBM related information.
 */
void SolverLbmD2Q9::CalMacroForceD2Q9Incompressible(const DefReal dt_lbm, GridNodeLbm* const ptr_node) const {
    ptr_node->rho_ = ptr_node->f_[kFX0Y0Z0] + ptr_node->f_[kFXnY0Z0] + ptr_node->f_[kFXpY0Z0] + ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFX0YpZ0] + ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXnYpZ0]
        + ptr_node->f_[kFXpYnZ0] + ptr_node->f_[kFXpYpZ0];
    ptr_node->velocity_[kXIndex] = (ptr_node->f_[kFXpY0Z0] - ptr_node->f_[kFXnY0Z0]
        + ptr_node->f_[kFXpYpZ0] - ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] + ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kYIndex] = (ptr_node->f_[kFX0YpZ0] - ptr_node->f_[kFX0YnZ0]
        + ptr_node->f_[kFXpYpZ0] + ptr_node->f_[kFXnYpZ0] - ptr_node->f_[kFXnYnZ0] - ptr_node->f_[kFXpYnZ0]);
    ptr_node->velocity_[kXIndex] += 0.5 * ptr_node->force_[kXIndex] * dt_lbm;
    ptr_node->velocity_[kYIndex] += 0.5 * ptr_node->force_[kYIndex] * dt_lbm;
}
}  // end namespace lbmproject
}  // end namespace rootproject