//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_solver_d3q19.cpp
* @author Zhengliang Liu
* @brief mandatory functions for LBM simulation using D2Q9 model.
* @date  2023-9-30
*/
#include <string>
#include <memory>
#include "d3q19_model/lbm_d3q19.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to create instance and set constants and pointer to functions for the D3Q19 model.
 * @return shared pointer to the instance of solver using D3Q19 model.
 */
std::shared_ptr<amrproject::SolverInterface> SolverCreatorLbmD3Q19::CreateSolver() const {
    std::shared_ptr<SolverLbmD3Q19> ptr_tmp = std::make_shared<SolverLbmD3Q19>();
    ptr_tmp->SetSolverType("LbmD3Q19");
    ptr_tmp->ResizeModelRelatedVectors();
    return ptr_tmp;
}
void SolverLbmD3Q19::InitialModelDependencies() {
    if (k0SolverDims_ != 3) {
        amrproject::LogManager::LogError("k0SolverDims_ for LBM D3Q19 Model should be"
            " 3 rather than: " + std::to_string(k0SolverDims_)
             + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }
}
/**
 * @brief function to perform propagation of a node to others in the LBM simulation.
 * @param[in] sfbitset_in space filling code of the current node.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD3Q19::StreamOutForAGivenNode(const DefSFBitset sfbitset_in,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    const amrproject::SFBitsetAux3D sfbitset_aux3d = dynamic_cast<const amrproject::SFBitsetAux3D&>(sfbitset_aux);
    DefSFBitset sfbitset_tmp, sfbitset_tmp1;
    // f(-x, 0, 0)
    ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0_) =
        ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0_);
    // f(-x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXnY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Z0_);
    }
    // f(-x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYnZ0_);
    }
    // f(-x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYpZ0_);
    }
    // f(-x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnY0Zn_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Zn_);
    }
    // f(-x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnY0Zp_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Zp_);
    }
    // f(+x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXpY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Z0_);
    }
    // f(+x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYnZ0_);
    }
    // f(+x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYpZ0_);
    }
    // f(+x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpY0Zn_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Zn_);
    }
    // f(+x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpY0Zp_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Zp_);
    }
    // f(0, -y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZ0_);
    }
    // f(0, -y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YnZn_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZn_);
    }
    // f((0, -y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YnZp_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZp_);
    }
    // f(0, +y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZ0_);
    }
    // f(0, +y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YpZn_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZn_);
    }
    // f((0, +y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YpZp_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZp_);
    }
    // f(0, 0, -z)
    sfbitset_tmp = sfbitset_aux3d.FindZNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0Y0Zn_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Zn_);
    }
    // f(0, 0, +z)
    sfbitset_tmp = sfbitset_aux3d.FindZPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0Y0Zp_) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Zp_);
    }
}
/**
 * @brief function to perform propagation from other nodes to the given node in the LBM simulation.
 * @param[in] sfbitset_in space filling code of the current node.
 * @param[in]  sfbitset_aux class to manage functions of spacing filling code related manipulations.
 * @param[out] ptr_map_grid_nodes pointer to grid nodes for LBM simulation.
 */
void SolverLbmD3Q19::StreamInForAGivenNode(const DefSFBitset sfbitset_in,
    const amrproject::SFBitsetAuxInterface& sfbitset_aux,
    DefMap<std::unique_ptr<GridNodeLbm>>* const ptr_map_grid_nodes) const {
    const amrproject::SFBitsetAux3D sfbitset_aux3d = dynamic_cast<const amrproject::SFBitsetAux3D&>(sfbitset_aux);
    DefSFBitset sfbitset_tmp, sfbitset_tmp1;
    ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0_) =
        ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0_);
    // f(-x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXnY0Z0_);
    }
    // f(-x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYnZ0_);
    }
    // f(-x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYpZ0_);
    }
    // f(-x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Zn_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnY0Zn_);
    }
    // f(-x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Zp_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnY0Zp_);
    }
    // f(+x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Z0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXpY0Z0_);
    }
    // f(+x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpYnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYnZ0_);
    }
    // f(+x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpYpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYpZ0_);
    }
    // f(+x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Zn_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpY0Zn_);
    }
    // f(+x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Zp_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpY0Zp_);
    }
    // f(0, -y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YnZ0_);
    }
    // f(0, -y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZn_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YnZn_);
    }
    // f((0, -y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZp_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YnZp_);
    }
    // f(0, +y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZ0_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YpZ0_);
    }
    // f(0, +y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZn_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YpZn_);
    }
    // f((0, +y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZp_) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YpZp_);
    }
    // f(0, 0, -z)
    sfbitset_tmp = sfbitset_aux3d.FindZPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Zn_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0Y0Zn_);
    }
    // f(0, 0, +z)
    sfbitset_tmp = sfbitset_aux3d.FindZNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Zp_) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0Y0Zp_);
    }
}

}  // end namespace lbmproject
}  // end namespace rootproject
