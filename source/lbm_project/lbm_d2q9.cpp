//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_d2q9.cpp
* @author Zhengliang Liu
* @brief functions used for LBM D2Q9 model.
* @date  2023-9-30
*/
#include <string>
#include "lbm_d2q9.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to create instance and set constants and pointer to functions for the D2Q9 model.
 * @return shared pointer to the instance of solver using D2Q9 model.
 */
std::shared_ptr<amrproject::SolverInterface> SolverCreatorLbmD2Q9::CreateSolver() {
    std::shared_ptr<SolverLbmD2Q9> ptr_temp = std::make_shared<SolverLbmD2Q9>();
    ptr_temp->k0SolverDims_ = 2;
    ptr_temp->k0NumQ_ = 9;
    ptr_temp->node_type_ = "LbmD2Q9";
    ptr_temp->k0Cx_ = { 0., 1., 0., -1., 0., 1., -1., -1., 1. };
    ptr_temp->k0Cy_ = { 0., 0., 1., 0., -1., 1., 1., -1., -1. };
    ptr_temp->k0Weights_ = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
    ptr_temp->ptr_func_cal_feq_ = &SolverLbmD2Q9::CalFeq2D;
    return ptr_temp;
}
void SolverLbmD2Q9::InitialSetIndices() {
    if (k0SolverDims_ != 2) {
        amrproject::LogManager::LogError("k0SolverDims_ for LBM D2Q9 Model should be"
            " 2 rather than: " + std::to_string(k0SolverDims_)
             + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
 * @brief function to perform streaming step in the LBM simulation.
 * @param flag_not_compute flag indicating whether to compute or not.
 * @param sfbitset_aux2d class to manage functions of spacing filling code related manipulations.
 * @param ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD2Q9::Stream(const DefAmrUint flag_not_compute, const amrproject::SFBitsetAux2D& sfbitset_aux2d,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    DefSFBitset sfbitset_tmp, sfbitset_tmp1;
    for (auto& iter_node : *ptr_map_grid_nodes) {
        if (iter_node.second->flag_status_ & flag_not_compute) {
        } else {
            // f(-x, 0)
            sfbitset_tmp = sfbitset_aux2d.FindXPos(iter_node.first);
            if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFXnY0Z0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXnY0Z0);
            }
            // f(-x, -y)
            sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
            if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFXnYnZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYnZ0);
            }
            // f(-x, +y)
            sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
            if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFXnYpZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYpZ0);
            }
            // f(+x, 0)
            sfbitset_tmp = sfbitset_aux2d.FindXNeg(iter_node.first);
            if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFXpY0Z0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXpY0Z0);
            }
            // f(+x, -y)
            sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
            if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFXpYnZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYnZ0);
            }
            // f(+x, +y)
            sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
            if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFXpYpZ0) = ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYpZ0);
            }
            // f(0, -y)
            sfbitset_tmp = sfbitset_aux2d.FindYPos(iter_node.first);
            if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFX0YnZ0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YnZ0);
            }
            // f(0, +y)
            sfbitset_tmp = sfbitset_aux2d.FindYNeg(iter_node.first);
            if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
                iter_node.second->f_.at(kFX0YpZ0) = ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YpZ0);
            }
        }
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject