//    Copyright (c) 2022, Zhengliang Liu
//    All rights reserved

/**
* @file test_grid.cpp
* @author Zhengliang Liu
* @date  2022-5-12
* @brief  test grid generation and related MPI functions
*/
#include "amr_manager.h"
#include "d2q9_model/lbm_d2q9.h"
using namespace rootproject;
using namespace rootproject::amrproject;
amrproject::AmrManager* ptr_amr_instance_;
DefReal u_max_ = 0.1, dx_ = 0.2, max_domain_height_ = 1., max_domain_length_ = 1.;
DefAmrIndexLUint max_t_ = 1000;
std::vector<DefReal> u_analytical_;
DefReal u_analytical_sum_;
void SetTestDependentParameters(DefAmrIndexUint max_level,
    const amrproject::SolverCreatorInterface& solver_creator) {
    ptr_amr_instance_->DefaultInitialization(2, max_level);
    ptr_amr_instance_->ptr_io_manager_->bool_binary_ = false;
    ptr_amr_instance_->ptr_io_manager_->vtk_instance_.vtk_ghost_cell_option_ =
        amrproject::EVtkWriterGhostCellOption::kPartitionMultiBlock;

    // grid related parameters //
    ptr_amr_instance_->ptr_grid_manager_->SetDomainSize({ max_domain_length_, max_domain_height_ });
    ptr_amr_instance_->ptr_grid_manager_->SetDomainGridSize({ dx_ });
    // end grid related parameters //

    ptr_amr_instance_->SetupParameters();
    // set grid node type and solver at all levels the same
    ptr_amr_instance_->AddSolverToGridManager(solver_creator);
    ptr_amr_instance_->SetDependentInfoForAllLevelsTheSame(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_solver_.at(0));

    for (auto& iter_grid : ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_) {
        lbmproject::GridInfoLbmInteface& grid_ref
            = *dynamic_cast<lbmproject::GridInfoLbmInteface*>(iter_grid.get());
        grid_ref.bool_forces_ = false;
        lbmproject::SolverLbmInterface& solver_ref
            = *dynamic_cast<lbmproject::SolverLbmInterface*>(grid_ref.ptr_solver_.get());
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryXNeg,
            lbmproject::ELbmBoundaryConditionScheme::kPeriodic, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryXPos,
            lbmproject::ELbmBoundaryConditionScheme::kPeriodic, &grid_ref.domain_boundary_condition_);
    }


    lbmproject::SolverLbmInterface& solver_ref = *dynamic_cast<lbmproject::SolverLbmInterface*>(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_solver_.at(0).get());
    solver_ref.k0BoolCompressible_ = true;
    DefReal tau_0 = sqrt(3. / 16) + 0.5;  // BGK gives the exact solution
    solver_ref.k0LbmViscosity_ = sqrt(3. / 16) * lbmproject::SolverLbmInterface::kCs_Sq_;
    DefReal lbm_height = max_domain_height_ / dx_ + 1.;  // bounce back wall at 0.5 distance to the node
    DefSizet num_probes = static_cast<DefSizet>(max_domain_height_ / dx_ + 1. + kEps);
    u_analytical_.resize(num_probes);
    u_analytical_sum_ = 0.;
    for (auto i_probe = 0; i_probe < num_probes; ++i_probe) {
        u_analytical_.at(i_probe) = u_max_ * ((i_probe + 0.5) / lbm_height);
        u_analytical_sum_ += Square(u_analytical_.at(i_probe));
    }

    ptr_amr_instance_->SetupSolverForGrids();

    ptr_amr_instance_->InitializeMesh();

}
void SetDomainBoundaryOtherThanPeriodic(
    lbmproject::ELbmBoundaryConditionScheme boundary_condition_type) {
    for (auto& iter_grid : ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_) {
        lbmproject::GridInfoLbmInteface& grid_ref
            = *dynamic_cast<lbmproject::GridInfoLbmInteface*>(iter_grid.get());
        lbmproject::SolverLbmInterface& solver_ref
            = *dynamic_cast<lbmproject::SolverLbmInterface*>(grid_ref.ptr_solver_.get());
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryYNeg,
            boundary_condition_type, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryYPos,
            boundary_condition_type, &grid_ref.domain_boundary_condition_);
        // set y boundary velocity to u_max
        grid_ref.domain_boundary_condition_.at(lbmproject::ELbmBoundaryType::kBoundaryYPos)
            ->SetValues({ u_max_, 0. });
    }
}
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    ptr_amr_instance_ = amrproject::AmrManager::GetInstance();
    ptr_amr_instance_->program_name_ = "test_tmp";

    DefAmrIndexUint dims = 2;   // dimension
    DefAmrIndexUint max_refinement_level = 0;  // maximum refinement level

    lbmproject::SolverCreatorLbmD2Q9 solver_creator = lbmproject::SolverCreatorLbmD2Q9();
    SetTestDependentParameters(0, solver_creator);
    SetDomainBoundaryOtherThanPeriodic(lbmproject::ELbmBoundaryConditionScheme::kBounceBack);

    for (DefAmrIndexLUint it = 0; it < max_t_; ++it) {
        ptr_amr_instance_->TimeMarching(it);
    }

    DefReal error_sum = 0., u_calculated;
    amrproject::SFBitsetAux2D& sfbitset_aux = *dynamic_cast<amrproject::SFBitsetAux2D*>(
        ptr_amr_instance_->ptr_grid_manager_.get());
    lbmproject::GridInfoLbmInteface& grid_ref = *dynamic_cast<lbmproject::GridInfoLbmInteface*>(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_.at(0).get());
    DefSFBitset code_tmp = sfbitset_aux.SFBitsetEncoding({ DefAmrIndexLUint(max_domain_length_ / dx_ + kEps) + 1, 1 });
    DefAmrIndexUint i_y = 0;
    while (grid_ref.ptr_lbm_grid_nodes_->find(code_tmp) != grid_ref.ptr_lbm_grid_nodes_->end()) {
        u_calculated = grid_ref.ptr_lbm_grid_nodes_->at(code_tmp)->velocity_[kXIndex];
        error_sum += Square(u_calculated - u_analytical_.at(i_y));
        code_tmp = sfbitset_aux.FindYPos(code_tmp);
        ++i_y;
    }
    error_sum = sqrt(error_sum) / sqrt(u_analytical_sum_);

    MPI_Finalize();
}
