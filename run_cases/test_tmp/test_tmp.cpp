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
const DefReal u_max_ = 0.1, dx_ = 0.2, max_domain_height_ = 1., max_domain_length_ = 1.6;
std::vector<DefReal> u_analytical_;
DefReal u_analytical_sum_;
void SetTestDependentParameters(
    const amrproject::SolverCreatorInterface& solver_creator) {
    ptr_amr_instance_->ptr_io_manager_->bool_binary_ = false;
    ptr_amr_instance_->ptr_io_manager_->vtk_instance_.vtk_ghost_cell_option_ =
        amrproject::EVtkWriterGhostCellOption::kPartitionMultiBlock;

    // grid related parameters //
    ptr_amr_instance_->ptr_grid_manager_->SetDomainSize({max_domain_length_, max_domain_height_});
    ptr_amr_instance_->ptr_grid_manager_->SetDomainGridSize({dx_});
    // end grid related parameters //

    ptr_amr_instance_->SetupParameters();
    // set grid node type and solver at all levels the same
    ptr_amr_instance_->AddSolverToGridManager(solver_creator);
    ptr_amr_instance_->SetDependentInfoForAllLevelsTheSame(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_solver_.at(0));

    for (auto& iter_grid : ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_) {
        lbmproject::GridInfoLbmInteface& grid_ref
            = *dynamic_cast<lbmproject::GridInfoLbmInteface*>(iter_grid.get());
        lbmproject::SolverLbmInterface& solver_ref
            = *dynamic_cast<lbmproject::SolverLbmInterface*>(grid_ref.GetPtrSolver());
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryXMin,
            lbmproject::ELbmBoundaryConditionScheme::kPeriodic, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryXMax,
            lbmproject::ELbmBoundaryConditionScheme::kPeriodic, &grid_ref.domain_boundary_condition_);
    }

    lbmproject::SolverLbmInterface& solver_ref = *dynamic_cast<lbmproject::SolverLbmInterface*>(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_solver_.at(0).get());
    solver_ref.k0BoolCompressible_ = true;
    solver_ref.SetDefaultViscosity(sqrt(3./16) * lbmproject::SolverLbmInterface::kCs_Sq_);
    DefReal lbm_height = max_domain_height_ /dx_ + 1.;  // bounce back wall at 0.5 distance to the node
    DefSizet num_probes = static_cast<DefSizet>(max_domain_height_/dx_ + 1. + kEps);
    u_analytical_.resize(num_probes);
    u_analytical_sum_ = 0.;
    for (auto i_probe = 0; i_probe < num_probes; ++i_probe) {
        u_analytical_.at(i_probe)  = u_max_* ((i_probe + 0.5)/lbm_height);
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
            = *dynamic_cast<lbmproject::SolverLbmInterface*>(grid_ref.GetPtrSolver());
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryYMin,
            boundary_condition_type, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(lbmproject::ELbmBoundaryType::kBoundaryYMax,
            boundary_condition_type, &grid_ref.domain_boundary_condition_);
        // set y boundary velocity to u_max
        grid_ref.domain_boundary_condition_.at(lbmproject::ELbmBoundaryType::kBoundaryYMax)
            ->SetValues({u_max_, 0.});
    }
}
void CalFeqLinear(const DefReal rho,
    const std::vector<DefReal>& velocity, std::vector<DefReal>* const ptr_feq) {
    lbmproject::SolverLbmInterface& solver_ref = *dynamic_cast<lbmproject::SolverLbmInterface*>(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_solver_.at(0).get());
    ptr_feq->resize(solver_ref.k0NumQ_);
    DefReal c_uv = 0.;
    for (int iq = 0; iq < solver_ref.k0NumQ_; ++iq) {
        c_uv = velocity.at(kXIndex) * solver_ref.k0Cx_.at(iq)
            + velocity.at(kYIndex) * solver_ref.k0Cy_.at(iq);
        ptr_feq->at(iq) = solver_ref.k0Weights_.at(iq) * (rho + 3. * c_uv);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    ptr_amr_instance_ = amrproject::AmrManager::GetInstance();
    DefInt dims = 2;   // dimension
    DefInt max_refinement_level = 1;  // maximum refinement level
    ptr_amr_instance_->DefaultInitialization(dims, max_refinement_level);

    // geometry related parameters //
    ptr_amr_instance_->ptr_grid_manager_->vec_ptr_tracking_info_creator_.push_back(
        std::make_unique<amrproject::TrackingGridInfoCreatorInterface>());
    amrproject::GeometryInfoOriginCreator geo_creator;
    ptr_amr_instance_->ptr_criterion_manager_->vec_ptr_geometries_.push_back(
        geo_creator.CreateGeometryInfo(dims));
    amrproject::GeometryInfoOrigin* ptr_geo_tmp =
        dynamic_cast<amrproject::GeometryInfoOrigin*>(ptr_amr_instance_->
        ptr_criterion_manager_->vec_ptr_geometries_.at(0).get());
    ptr_geo_tmp->ptr_geo_shape_ = std::make_unique<amrproject::GeoShapeDefaultLine2D>();
    amrproject::GeoShapeDefaultLine2D* ptr_line =
        dynamic_cast<amrproject::GeoShapeDefaultLine2D*>(ptr_geo_tmp->ptr_geo_shape_.get());
    ptr_line->start_point_ = { 0., 0};
    ptr_line->end_point_ = { max_domain_length_, 0};
    ptr_geo_tmp->ptr_tracking_grid_info_creator_ =
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_tracking_info_creator_.at(0).get();
    ptr_geo_tmp->geometry_center_ = { max_domain_length_/2, max_domain_height_/2 };
    ptr_geo_tmp->k0XIntExtendPositive_ = { 4, 4 };
    ptr_geo_tmp->k0XIntExtendNegative_ = { 4, 4 };
    ptr_geo_tmp->k0YIntExtendPositive_ = { 4, 4 };
    ptr_geo_tmp->k0YIntExtendNegative_ = { 4, 4 };
    ptr_geo_tmp->k0IntInnerExtend_ = { 2, 2 };
    ///< number of extened layers

    ptr_geo_tmp->grid_extend_type_ = amrproject::EGridExtendType::kInAndOut;
    ptr_geo_tmp->i_level_ = max_refinement_level;
    /* used for generating predefined geometries, number of input parameters
    is based on the type of geometry_shape_*/
    DefReal dx = dx_ / DefReal(std::pow(2, max_refinement_level));
    ptr_geo_tmp->SetOffset({2*dx, 2*dx, 2*dx});

    lbmproject::SolverCreatorLbmD2Q9 solver_creator = lbmproject::SolverCreatorLbmD2Q9();
    SetTestDependentParameters(solver_creator);
    SetDomainBoundaryOtherThanPeriodic(lbmproject::ELbmBoundaryConditionScheme::kBounceBack);

    lbmproject::SolverLbmInterface& solver_ref = *dynamic_cast<lbmproject::SolverLbmInterface*>(
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_solver_.at(0).get());
        solver_ref.func_cal_feq_ =  [](const DefReal rho, const std::vector<DefReal>& velocity,
        std::vector<DefReal>* const ptr_feq) {
        CalFeqLinear(rho, velocity, ptr_feq);
    };
    for (auto& iter :
        ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_.at(0)->map_grid_node_) {
        lbmproject::GridNodeLbm & lbm_node = *dynamic_cast<lbmproject::GridNodeLbm*>(iter.second.get());
        CalFeqLinear(1, {u_max_, 0}, &lbm_node.f_collide_);
        CalFeqLinear(1, {u_max_, 0}, &lbm_node.f_);
        lbm_node.velocity_ = {u_max_, 0};
    }


    // // check if nodes for sending and receiving are correct, only used for debug purposes
    // for (int i_level = 0; i_level< 2; ++i_level) {
    //     ptr_amr_instance_->ptr_mpi_manager_->CheckMpiNodesCorrespondence(
    //         *ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_.at(i_level).get());
    // }


    // if (ptr_amr_instance_->ptr_mpi_manager_->rank_id_ == 1) {
    //     for (auto& iter :
    //     ptr_amr_instance_->ptr_grid_manager_->vec_ptr_grid_info_.at(1)->map_grid_node_) {
    //     lbmproject::GridNodeLbm & lbm_node = *dynamic_cast<lbmproject::GridNodeLbm*>(iter.second.get());
    //     CalFeqLinear(1, {u_max_, 0}, &lbm_node.f_);
    //     CalFeqLinear(1, {u_max_, 0}, &lbm_node.f_collide_);
    // }
    // }


    for (DefAmrLUint it = 0; it < 1; ++it) {
        ptr_amr_instance_->TimeMarching(it);
    }
    MPI_Finalize();
}
