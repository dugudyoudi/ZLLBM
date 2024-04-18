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
    std::shared_ptr<SolverLbmD3Q19> ptr_temp = std::make_shared<SolverLbmD3Q19>();
    ptr_temp->k0SolverDims_ = 3;
    ptr_temp->k0NumQ_ = 19;
    ptr_temp->solver_type_ = "LbmD3Q19";

    ptr_temp->k0Cx_ = { 0., 1., -1., 0.,  0., 0.,  0., 1., -1.,  1., -1., 1., -1.,  1., -1., 0., 0.,   0.,  0.};
    ptr_temp->k0Cy_ = { 0., 0.,  0., 1., -1., 0.,  0., 1.,  1., -1., -1., 0.,  0.,  0.,  0., 1., -1.,  1., -1.};
    ptr_temp->k0Cz_ = { 0., 0.,  0., 0.,  0., 1., -1., 0.,  0.,  0.,  0., 1.,  1., -1., -1., 1.,  1., -1., -1.};
    ptr_temp->k0Weights_ = { 1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
        1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36. };
    ptr_temp->k0NumQInOneDirection_ = 5;
    ptr_temp->k0QIndicesNeg_ = {{ptr_temp->kFXnY0Z0, ptr_temp->kFXnYnZ0, ptr_temp->kFXnYpZ0,
        ptr_temp->kFXnY0Zn, ptr_temp->kFXnY0Zp},
        {ptr_temp->kFX0YnZ0, ptr_temp->kFXnYnZ0, ptr_temp->kFXpYnZ0,
        ptr_temp->kFX0YnZn, ptr_temp->kFX0YnZp},
        {ptr_temp->kFX0Y0Zn, ptr_temp->kFXnY0Zn, ptr_temp->kFXpY0Zn,
        ptr_temp->kFX0YnZn, ptr_temp->kFX0YpZn}};
    ptr_temp->k0QIndicesPos_ = {{ptr_temp->kFXpY0Z0, ptr_temp->kFXpYpZ0, ptr_temp->kFXpYnZ0,
        ptr_temp->kFXpY0Zp, ptr_temp->kFXpY0Zn},
        {ptr_temp->kFX0YpZ0, ptr_temp->kFXpYpZ0, ptr_temp->kFXnYpZ0,
        ptr_temp->kFX0YpZp, ptr_temp->kFX0YpZn},
        {ptr_temp->kFX0Y0Zp, ptr_temp->kFXpY0Zp, ptr_temp->kFXnY0Zp,
        ptr_temp->kFX0YpZp, ptr_temp->kFX0YnZp}};
    ptr_temp->ResizeModelRelatedVectors();
    return ptr_temp;
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
    ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0) = ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0);
    // f(-x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXnY0Z0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Z0);
    }
    // f(-x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYnZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYnZ0);
    }
    // f(-x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnYpZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnYpZ0);
    }
    // f(-x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnY0Zn) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Zn);
    }
    // f(-x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXnY0Zp) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXnY0Zp);
    }
    // f(+x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFXpY0Z0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Z0);
    }
    // f(+x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYnZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYnZ0);
    }
    // f(+x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpYpZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpYpZ0);
    }
    // f(+x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpY0Zn) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Zn);
    }
    // f(+x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFXpY0Zp) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFXpY0Zp);
    }
    // f(0, -y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YnZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZ0);
    }
    // f(0, -y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YnZn) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZn);
    }
    // f((0, -y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YnZp) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YnZp);
    }
    // f(0, +y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0YpZ0) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZ0);
    }
    // f(0, +y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YpZn) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZn);
    }
    // f((0, +y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp1)->f_.at(kFX0YpZp) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0YpZp);
    }
    // f(0, 0, -z)
    sfbitset_tmp = sfbitset_aux3d.FindZNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0Y0Zn) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Zn);
    }
    // f(0, 0, +z)
    sfbitset_tmp = sfbitset_aux3d.FindZPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_tmp)->f_.at(kFX0Y0Zp) =
            ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Zp);
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
    ptr_map_grid_nodes->at(sfbitset_in)->f_collide_.at(kFX0Y0Z0) =
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Z0);
    // f(-x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Z0) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXnY0Z0);
    }
    // f(-x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnYnZ0) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYnZ0);
    }
    // f(-x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnYpZ0) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnYpZ0);
    }
    // f(-x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Zn) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnY0Zn);
    }
    // f(-x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXnY0Zp) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXnY0Zp);
    }
    // f(+x, 0, 0)
    sfbitset_tmp = sfbitset_aux3d.FindXNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Z0) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFXpY0Z0);
    }
    // f(+x, -y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpYnZ0) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYnZ0);
    }
    // f(+x, +y, 0)
    sfbitset_tmp1 = sfbitset_aux3d.FindYNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpYpZ0) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpYpZ0);
    }
    // f(+x, 0, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Zn) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpY0Zn);
    }
    // f(+x, 0, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFXpY0Zp) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFXpY0Zp);
    }
    // f(0, -y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZ0) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YnZ0);
    }
    // f(0, -y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZn) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YnZn);
    }
    // f((0, -y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YnZp) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YnZp);
    }
    // f(0, +y, 0)
    sfbitset_tmp = sfbitset_aux3d.FindYNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZ0) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0YpZ0);
    }
    // f(0, +y, -z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZPos(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZn) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YpZn);
    }
    // f((0, +y, +z)
    sfbitset_tmp1 = sfbitset_aux3d.FindZNeg(sfbitset_tmp);
    if (ptr_map_grid_nodes->find(sfbitset_tmp1) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0YpZp) =
            ptr_map_grid_nodes->at(sfbitset_tmp1)->f_collide_.at(kFX0YpZp);
    }
    // f(0, 0, -z)
    sfbitset_tmp = sfbitset_aux3d.FindZPos(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Zn) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0Y0Zn);
    }
    // f(0, 0, +z)
    sfbitset_tmp = sfbitset_aux3d.FindZNeg(sfbitset_in);
    if (ptr_map_grid_nodes->find(sfbitset_tmp) != ptr_map_grid_nodes->end()) {
        ptr_map_grid_nodes->at(sfbitset_in)->f_.at(kFX0Y0Zp) =
            ptr_map_grid_nodes->at(sfbitset_tmp)->f_collide_.at(kFX0Y0Zp);
    }
}

}  // end namespace lbmproject
}  // end namespace rootproject
