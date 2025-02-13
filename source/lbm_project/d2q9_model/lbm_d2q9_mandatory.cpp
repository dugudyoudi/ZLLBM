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
    std::shared_ptr<SolverLbmD2Q9> ptr_tmp = std::make_shared<SolverLbmD2Q9>();
    ptr_tmp->SetSolverDims(2);
    ptr_tmp->SetSolverType("LbmD2Q9");
    ptr_tmp->ResizeModelRelatedVectors();
    return ptr_tmp;
}
void SolverLbmD2Q9::InitialModelDependencies() {
    if (k0SolverDims_ != 2) {
        amrproject::LogManager::LogError("k0SolverDims_ for LBM D2Q9 Model should be"
            " 2 rather than: " + std::to_string(k0SolverDims_)
             + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
    if (!k0BoolCompressible_) {
        this->func_macro_without_force_ = [this](const std::vector<DefReal>& f,
            DefReal* const ptr_rho, std::vector<DefReal>* const ptr_velocity) {
            this->CalMacroD2Q9Incompressible(f, ptr_rho, ptr_velocity);
        };
    }
}
/**
 * @brief function to perform propagation of a node to others in the LBM simulation.
 * @param[in] sfbitset_in space filling code of the current node.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD2Q9::StreamOutForAGivenNode(const DefSFBitset sfbitset_in,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    const amrproject::SFBitsetAux2D sfbitset_aux2d = dynamic_cast<const amrproject::SFBitsetAux2D&>(sfbitset_aux);
    DefSFBitset sfbitset_tmp, sfbitset_tmp1;
    // f(0, 0)
    ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0_)
        = ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0_);
    // f(-x, 0)
    sfbitset_tmp = sfbitset_aux2d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXnY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Z0_);
    }
    // f(-x, -y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYnZ0_);
    }
    // f(-x, +y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYpZ0_);
    }
    // f(+x, 0)
    sfbitset_tmp = sfbitset_aux2d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXpY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Z0_);
    }
    // f(+x, -y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYnZ0_);
    }
    // f(+x, +y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYpZ0_);
    }
    // f(0, -y)
    sfbitset_tmp = sfbitset_aux2d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZ0_);
    }
    // f(0, +y)
    sfbitset_tmp = sfbitset_aux2d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZ0_);
    }
}
/**
 * @brief function to perform propagation from other nodes to the given node in the LBM simulation.
 * @param[in] sfbitset_in space filling code of the current node.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD2Q9::StreamInForAGivenNode(const DefSFBitset sfbitset_in,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    const amrproject::SFBitsetAux2D sfbitset_aux2d = dynamic_cast<const amrproject::SFBitsetAux2D&>(sfbitset_aux);
    DefSFBitset sfbitset_tmp, sfbitset_tmp1;
    ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0_) =
        ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0_);
    // f(-x, 0)
    sfbitset_tmp = sfbitset_aux2d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXnY0Z0_);
    }
    // f(-x, -y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYnZ0_);
    }
    // f(-x, +y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYpZ0_);
    }
    // f(+x, 0)
    sfbitset_tmp = sfbitset_aux2d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXpY0Z0_);
    }
    // f(+x, -y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYnZ0_);
    }
    // f(+x, +y)
    sfbitset_tmp1 = sfbitset_aux2d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYpZ0_);
    }
    // f(0, -y)
    sfbitset_tmp = sfbitset_aux2d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YnZ0_);
    }
    // f(0, +y)
    sfbitset_tmp = sfbitset_aux2d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YpZ0_);
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
