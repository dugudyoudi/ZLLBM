//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_les_smagorinsky.cpp
* @author Zhengliang Liu
* @brief functions used for Smagorinsky LES model.
* @date  2025-2-20
*/
#include "./lbm_interface.h"
#include "les_models/lbm_les_models.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to calculate relaxation time based on Smagorinsky LES model.
 * @param[in] dt_lbm  time spacing of current level.
 * @param[in] tau_0 relaction time computed from universal parameters, e.g. physical viscosity.
 * @param[in] feq equilibrium distribution function. 
 * @param[in] node reference to a LBM node.
 * @param[in] lbm_solver reference to the LBM solver.
 * @return calculated relaxation time based on LES model.
 */
DefReal LesModelSmagorinsky ::CalSgsRelaxationTimeWithoutForce(const DefReal dt_lbm,
    const DefReal tau_0, const std::vector<DefReal>& feq,
    const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) const {
    const DefInt num_q = lbm_solver.k0NumQ_;
    DefReal pi_xx = 0., pi_xy = 0., pi_yy = 0., pi_zz = 0., pi_xz = 0., pi_yz = 0.;
    if (lbm_solver.GetSolverDim() == 2) {
        for (DefInt iq = 1; iq < num_q; ++iq) {
            pi_xx += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cx_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_xy += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_yy += lbm_solver.k0Cy_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
        }
    } else {
        for (DefInt iq = 1; iq < num_q; ++iq) {
            pi_xx += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cx_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_xy += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_yy += lbm_solver.k0Cy_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_zz += lbm_solver.k0Cz_.at(iq) * lbm_solver.k0Cz_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_xz += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cz_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_yz += lbm_solver.k0Cy_.at(iq) * lbm_solver.k0Cz_.at(iq) * (node.f_.at(iq) - feq.at(iq));
        }
    }
    DefReal pi_neq2 = pi_xx * pi_xx + pi_yy * pi_yy +pi_zz * pi_zz
        + 2. * (pi_xy * pi_xy + pi_xz * pi_xz + pi_yz * pi_yz);
    return 0.5 * (sqrt(tau_0 * tau_0 + k0Smagorinsky * sqrt(2.*pi_neq2)) - tau_0);
}
/**
 * @brief function to calculate relaxation time based on Smagorinsky LES model.
 * @param[in] dt_lbm  time spacing of current level.
 * @param[in] tau_0 relaction time computed from universal parameters, e.g. physical viscosity.
 * @param[in] feq equilibrium distribution function.
 * @param[in] force body force acting on the node.
 * @param[in] node reference to a LBM node.
 * @param[in] lbm_solver reference to the LBM solver.
 * @return calculated relaxation time based on LES model.
 */
DefReal LesModelSmagorinsky ::CalSgsRelaxationTimeWithForce(const DefReal dt_lbm,
    const DefReal tau_0, const std::vector<DefReal>& feq, const std::vector<DefReal>& force,
    const GridNodeLbm& node, const SolverLbmInterface& lbm_solver) const {
    const DefInt num_q = lbm_solver.k0NumQ_;
    DefReal pi_xx = 0., pi_xy = 0., pi_yy = 0., pi_zz = 0., pi_xz = 0., pi_yz = 0.;
    if (lbm_solver.GetSolverDim() == 2) {
        for (DefInt iq = 1; iq < num_q; ++iq) {
            pi_xx += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cx_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_xy += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_yy += lbm_solver.k0Cy_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
        }
        pi_xx += dt_lbm * force.at(kXIndex) * node.velocity_[kXIndex];
        pi_yy += dt_lbm * force.at(kYIndex) * node.velocity_[kYIndex];
        pi_xy += 0.5 * dt_lbm * (force.at(kXIndex) * node.velocity_[kXIndex]
            + force.at(kYIndex) * node.velocity_[kYIndex]);
    } else {
        for (DefInt iq = 1; iq < num_q; ++iq) {
            pi_xx += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cx_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_xy += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_yy += lbm_solver.k0Cy_.at(iq) * lbm_solver.k0Cy_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_zz += lbm_solver.k0Cz_.at(iq) * lbm_solver.k0Cz_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_xz += lbm_solver.k0Cx_.at(iq) * lbm_solver.k0Cz_.at(iq) * (node.f_.at(iq) - feq.at(iq));
            pi_yz += lbm_solver.k0Cy_.at(iq) * lbm_solver.k0Cz_.at(iq) * (node.f_.at(iq) - feq.at(iq));
        }
        pi_xx += dt_lbm * force.at(kXIndex) * node.velocity_[kXIndex];
        pi_yy += dt_lbm * force.at(kYIndex) * node.velocity_[kYIndex];
        pi_xy += 0.5 * dt_lbm * (force.at(kXIndex) * node.velocity_[kXIndex]
            + force.at(kYIndex) * node.velocity_[kYIndex]);
        pi_zz += dt_lbm * force.at(kZIndex) * node.velocity_[kZIndex];
        pi_xz += 0.5 * dt_lbm * (force.at(kXIndex) * node.velocity_[kXIndex]
        + force.at(kZIndex) * node.velocity_[kZIndex]);
        pi_yz += 0.5 * dt_lbm * (force.at(kYIndex) * node.velocity_[kYIndex]
        + force.at(kZIndex) * node.velocity_[kZIndex]);
    }
    DefReal pi_neq2 = pi_xx * pi_xx + pi_yy * pi_yy +pi_zz * pi_zz
        + 2. * (pi_xy * pi_xy + pi_xz * pi_xz + pi_yz * pi_yz);
    return 0.5 * (sqrt(tau_0 * tau_0 + k0Smagorinsky * sqrt(2.*pi_neq2)) - tau_0);
}
}  // end namespace lbmproject
}  // end namespace rootproject
