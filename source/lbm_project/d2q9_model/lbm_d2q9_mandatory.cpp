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
/**
 * @brief function to get constant matrix for MRT collision.
 * @return vectors storing the constant matrix defined in current class.
 */
std::vector<std::vector<DefReal>> SolverLbmD2Q9::GetMrtMMatrix() const {
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
std::vector<std::vector<DefReal>> SolverLbmD2Q9::GetMrtImMatrix() const {
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
std::vector<std::vector<DefReal>> SolverLbmD2Q9::InitialMrtSMatrix(const DefReal relax_tau) const {
    std::vector<std::vector<DefReal>> diag_s(k0NumQ_);
    std::array<DefReal, 9> sk(k0VectorSMrt_);
    sk[7] = relax_tau;
    sk[8] = relax_tau;
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
std::vector<std::vector<DefReal>> SolverLbmD2Q9::InitialMrtDMatrix(const DefReal relax_tau) const {
    std::vector<std::vector<DefReal>> diag_d(k0NumQ_);
    std::array<DefReal, 9> sk(k0VectorSMrt_);
    sk[7] = relax_tau;
    sk[8] = relax_tau;
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
void SolverLbmD2Q9::UpdateMrtSMatrix(
    const DefReal relax_tau, std::vector<std::vector<DefReal>>* const ptr_diag_s) const {
    ptr_diag_s->at(7).at(7) = relax_tau;
    ptr_diag_s->at(8).at(8) = relax_tau;
}
/**
 * @brief function to update elements related to relaxation time in diagonal matrix for forcing term in MRT collision.
 * @param[in]    relax_tau reciprocal of relaxation time.
 * @param[out]   ptr_diag_d pointer to matrix for forcing term.
 */
void SolverLbmD2Q9::UpdateMrtDMatrix(
    const DefReal relax_tau, std::vector<std::vector<DefReal>>* const ptr_diag_d) const {
    ptr_diag_d->at(7).at(7) = 1. - 0.5 * relax_tau;
    ptr_diag_d->at(8).at(8) = 1. - 0.5 * relax_tau;
}
}  // end namespace lbmproject
}  // end namespace rootproject
