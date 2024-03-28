//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_solver_d2q9.cpp
* @author Zhengliang Liu
* @brief mandatory functions for LBM simulation using D2Q9 model.
* @date  2023-9-30
*/
#include <string>
#include "d2q9_model/lbm_d2q9.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to create instance and set constants and pointer to functions for the D2Q9 model.
 * @return shared pointer to the instance of solver using D2Q9 model.
 */
std::shared_ptr<amrproject::SolverInterface> SolverCreatorLbmD2Q9::CreateSolver() const {
    std::shared_ptr<SolverLbmD2Q9> ptr_temp = std::make_shared<SolverLbmD2Q9>();
    ptr_temp->k0SolverDims_ = 2;
    ptr_temp->k0NumQ_ = 9;
    ptr_temp->solver_type_ = "LbmD2Q9";
    ptr_temp->k0Cx_ = { 0., 1., 0., -1., 0., 1., -1., -1., 1. };
    ptr_temp->k0Cy_ = { 0., 0., 1., 0., -1., 1., 1., -1., -1. };
    ptr_temp->k0Weights_ = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };
    ptr_temp->k0NumQInOneDirection_ = 3;
    ptr_temp->k0QIndicesNeg_ = {{ptr_temp->kFXnY0Z0, ptr_temp->kFXnYnZ0, ptr_temp->kFXnYpZ0},
        {ptr_temp->kFX0YnZ0, ptr_temp->kFXnYnZ0, ptr_temp->kFXpYnZ0}};
    ptr_temp->k0QIndicesPos_ = {{ptr_temp->kFXpY0Z0, ptr_temp->kFXpYpZ0, ptr_temp->kFXpYnZ0},
        {ptr_temp->kFX0YpZ0, ptr_temp->kFXpYpZ0, ptr_temp->kFXnYpZ0}};
    ptr_temp->ResizeModelRelatedVectors();
    return ptr_temp;
}
void SolverLbmD2Q9::InitialModelDependencies() {
    if (k0SolverDims_ != 2) {
        amrproject::LogManager::LogError("k0SolverDims_ for LBM D2Q9 Model should be"
            " 2 rather than: " + std::to_string(k0SolverDims_)
             + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (k0BoolCompressible_) {
        this->func_macro_ = [this](const DefReal dt_lbm, GridNodeLbm* const ptr_node) {
            this->CalMacroD2Q9Compressible(dt_lbm, ptr_node);
        };
        this->func_macro_force_ = [this](const DefReal dt_lbm, GridNodeLbm* const ptr_node) {
            this->CalMacroForceD2Q9Compressible(dt_lbm, ptr_node);
        };
    } else {
        this->func_macro_ = [this](const DefReal dt_lbm, GridNodeLbm* const ptr_node) {
            this->CalMacroD2Q9Incompressible(dt_lbm, ptr_node);
        };
        this->func_macro_force_ = [this](const DefReal dt_lbm, GridNodeLbm* const ptr_node) {
            this->CalMacroForceD2Q9Incompressible(dt_lbm, ptr_node);
        };
    }
}
/**
 * @brief function to perform streaming step for a node in the LBM simulation.
 * @param[in] sfbitset_in space filling code of the current node.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD2Q9::StreamForAGivenNode(const DefSFBitset sfbitset_in,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    const amrproject::SFBitsetAux2D sfbitset_aux2d = dynamic_cast<const amrproject::SFBitsetAux2D&>(sfbitset_aux);
    DefSFBitset sfbitset_tmp, sfbitset_tmp1;
    // f(0, 0)
    ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0) = ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0);
    // f(-x, 0)
    sfbitset_tmp = sfbitset_aux2d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXnY0Z0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Z0);
    }
    // f(-x, -y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYnZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYnZ0);
    }
    // f(-x, +y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYpZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYpZ0);
    }
    // f(+x, 0)
    sfbitset_tmp = sfbitset_aux2d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXpY0Z0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Z0);
    }
    // f(+x, -y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYnZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYnZ0);
    }
    // f(+x, +y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYpZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYpZ0);
    }
    // f(0, -y)
    sfbitset_tmp = sfbitset_aux2d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YnZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZ0);
    }
    // f(0, +y)
    sfbitset_tmp = sfbitset_aux2d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YpZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZ0);
    }
}
/**
 * @brief function to perform streaming step in the LBM simulation.
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
}  // end namespace lbmproject
}  // end namespace rootproject
