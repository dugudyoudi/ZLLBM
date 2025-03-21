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
    // f(0, 0, 0)
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
/**
 * @brief function to get constant matrix for MRT collision.
 * @return vectors storing the constant matrix defined in current class.
 */
std::vector<std::vector<DefReal>> SolverLbmD3Q19::GetMrtMMatrix() const {
    std::vector<std::vector<DefReal>> matrix;
    for (const auto& row : kMatrixMMrt_) {
        matrix.emplace_back(row.begin(), row.end());
    }
    return matrix;
}
/**
 * @brief function to get constant inverse matrix for MRT collision.
 * @return vectors storing the constant inverse matrix defined in current class.
 */
std::vector<std::vector<DefReal>> SolverLbmD3Q19::GetMrtImMatrix() const {
    std::vector<std::vector<DefReal>> matrix;
    for (const auto& row : kMatrixImMrt_) {
        matrix.emplace_back(row.begin(), row.end());
    }
    return matrix;
}
/**
 * @brief function to initialize diagonal matrix for non-equilibrium distribution functions in MRT collision.
 * @param[in] relax_tau reciprocal of relaxation time.
 * @return     diagonal matrix for non-equilibrium distribution functions.
 */
std::vector<std::vector<DefReal>> SolverLbmD3Q19::InitialMrtSMatrix(const DefReal relax_tau) const {
    std::vector<std::vector<DefReal>> diag_s(k0NumQ_);
    std::array<DefReal, 19> sk(k0VectorSMrt_);
    sk[9] = relax_tau; sk[11] = relax_tau; sk[13] = relax_tau; sk[14] = relax_tau; sk[15] = relax_tau;  // omega_nu
    for (DefInt j = 0; j < k0NumQ_; ++j) {
        diag_s.at(j).resize(k0NumQ_);
        for (DefInt i = 0; i < k0NumQ_; ++i) {
            if (i == j) {
                diag_s.at(j)[i] = sk[i];
            } else {
                diag_s.at(j)[i] = 0.;
            }
        }
    }
    return diag_s;
}
/**
 * @brief function to initialize diagonal matrix for forcing term in MRT collision.
 * @param[in] relax_tau reciprocal of relaxation time.
 * @return     diagonal matrix for forcing term.
 */
std::vector<std::vector<DefReal>> SolverLbmD3Q19::InitialMrtDMatrix(const DefReal relax_tau) const {
    std::vector<std::vector<DefReal>> diag_d(k0NumQ_);
    std::array<DefReal, 19> sk(k0VectorSMrt_);
    sk[9] = relax_tau; sk[11] = relax_tau; sk[13] = relax_tau; sk[14] = relax_tau; sk[15] = relax_tau;  // omega_nu
    for (DefInt j = 0; j < k0NumQ_; ++j) {
        diag_d.at(j).resize(k0NumQ_);
        for (DefInt i = 0; i < k0NumQ_; ++i) {
            if (i == j) {
                diag_d.at(j)[i] = 1. - 0.5 * sk[i];
            } else {
                diag_d.at(j)[i] = 0.;
            }
        }
    }
    return diag_d;
}
/**
 * @brief function to update elements related to relaxation time in diagonal matrix for non-equilibrium distribution functions in MRT collision.
 * @param[in]    relax_tau reciprocal of relaxation time.
 * @param[out]   ptr_diag_s pointer to matrix for non-equilibrium distribution functions.
 */
void SolverLbmD3Q19::UpdateMrtSMatrix(
    const DefReal relax_tau, std::vector<std::vector<DefReal>>* const ptr_diag_s) const {
    const std::array<DefInt, 5> tau_index = {9, 11, 13, 14, 15};
    for (DefInt i = 0; i < 5; ++i) {
        ptr_diag_s->at(tau_index[i]).at(tau_index[i]) = relax_tau;
    }
}
/**
 * @brief function to update elements related to relaxation time in diagonal matrix for forcing term in MRT collision.
 * @param[in]    relax_tau reciprocal of relaxation time.
 * @param[out]   ptr_diag_d pointer to matrix for forcing term.
 */
void SolverLbmD3Q19::UpdateMrtDMatrix(
    const DefReal relax_tau, std::vector<std::vector<DefReal>>* const ptr_diag_d) const {
    const std::array<DefInt, 5> tau_index = {9, 11, 13, 14, 15};
    for (DefInt i = 0; i < 5; ++i) {
        ptr_diag_d->at(tau_index[i]).at(tau_index[i]) = 1. - 0.5 * relax_tau;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
