//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d3q19_solver.cpp
* @author Zhengliang Liu
* @brief alternative functions for LBM simulation using D3Q19 model.
* @date  2023-9-30
*/
#include "d3q19_model/lbm_d3q19.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to perform streaming step in the LBM simulation.
 * @param[in] flag_not_compute flag indicating whether to compute or not.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD3Q19::Stream(const DefInt flag_not_compute, const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    if (ptr_map_grid_nodes != nullptr) {
        DefSFBitset sfbitset_tmp, sfbitset_tmp1;
        const amrproject::SFBitsetAux3D sfbitset_aux3d = dynamic_cast<const amrproject::SFBitsetAux3D&>(sfbitset_aux);
        for (auto& iter_node : *ptr_map_grid_nodes) {
            if (iter_node.second->flag_status_ & flag_not_compute) {
            } else {
                iter_node.second->f_.at(kFX0Y0Z0) = iter_node.second->f_collide_.at(kFX0Y0Z0);
                // f(-x, 0, 0)
                sfbitset_tmp = sfbitset_aux3d.FindXPos(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXnY0Z0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXnY0Z0);
                }
                // f(-x, -y, 0)
                sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXnYnZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYnZ0);
                }
                // f(-x, +y, 0)
                sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXnYpZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYpZ0);
                }
                // f(-x, 0, -z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXnY0Zn) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnY0Zn);
                }
                // f(-x, 0, +z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXnY0Zp) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnY0Zp);
                }
                // f(+x, 0, 0)
                sfbitset_tmp = sfbitset_aux3d.FindXNeg(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXpY0Z0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXpY0Z0);
                }
                // f(+x, -y, 0)
                sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXpYnZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYnZ0);
                }
                // f(+x, +y, 0)
                sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXpYpZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYpZ0);
                }
                // f(+x, 0, -z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXpY0Zn) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpY0Zn);
                }
                // f(+x, 0, +z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFXpY0Zp) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpY0Zp);
                }
                // f(0, -y, 0)
                sfbitset_tmp = sfbitset_aux3d.FindYPos(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0YnZ0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YnZ0);
                }
                // f(0, -y, -z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0YnZn) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YnZn);
                }
                // f((0, -y, +z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0YnZp) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YnZp);
                }
                // f(0, +y, 0)
                sfbitset_tmp = sfbitset_aux3d.FindYNeg(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0YpZ0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YpZ0);
                }
                // f(0, +y, -z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0YpZn) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YpZn);
                }
                // f((0, +y, +z)
                sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
                if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0YpZp) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YpZp);
                }
                // f(0, 0, -z)
                sfbitset_tmp = sfbitset_aux3d.FindZPos(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0Y0Zn) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0Y0Zn);
                }
                // f(0, 0, +z)
                sfbitset_tmp = sfbitset_aux3d.FindZNeg(iter_node.first);
                if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                    iter_node.second->f_.at(kFX0Y0Zp) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0Y0Zp);
                }
            }
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
