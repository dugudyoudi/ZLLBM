//    Copyright (c) 2021 - 2025, Zhengliang Liu
//    All rights reserved

/**
* @file cylinder2d.cpp
* @author Zhengliang Liu
* @brief  flow past a 2D cylinder using IB-LBM.
*/
#include "./amr_manager.h"
#include "d2q9_model/lbm_d2q9.h"
#include "immersed_boundary/immersed_boundary.h"
#include "immersed_boundary/geometry_ib_shape.h"
using namespace rootproject;
int main(int argc, char* argv[]) {
    // case related settings
    lbmproject::SolverCreatorLbmD2Q9 solver_creator = lbmproject::SolverCreatorLbmD2Q9();
    const lbmproject::ELbmCollisionOperatorType collision_type = lbmproject::ELbmCollisionOperatorType::kLbmSrt;
    DefAmrLUint max_t = 20000;
    DefInt max_refinement_level = 2;
    std::vector<DefReal> domain_size = {30, 20};
    DefReal grid_size = 0.1;
    std::array<DefReal, 2> center = {10, 10};
    DefReal diameter = 1.;
    std::vector<DefReal> velocity = {0.1, 0};;
    DefReal Re = 100;
    lbmproject::ELbmBoundaryConditionScheme boundary_condition =
        lbmproject::ELbmBoundaryConditionScheme::kNonEqExtrapolation;

    // startup project manager
    amrproject::AmrManager* ptr_amr_instance;
    DefInt dims = 2;
    ptr_amr_instance = amrproject::AmrManager::GetInstance();
    ptr_amr_instance->StartupAmrManager(argc, argv);
    ptr_amr_instance->StartupInitialization(dims);
    ptr_amr_instance->ptr_grid_manager_->SetMaxLevel(max_refinement_level);
    ptr_amr_instance->program_name_ = "cylinder2d";
    ptr_amr_instance->ptr_io_manager_->bool_binary_ = false;

    // set geometry parameters
    lbmproject::GeometryInfoImmersedBoundaryCreator geo_creator;
    ptr_amr_instance->ptr_criterion_manager_->vec_ptr_geometries_.push_back(
        geo_creator.CreateGeometryInfo(dims));
    std::shared_ptr<lbmproject::GeometryInfoImmersedBoundary> ptr_geo_tmp =
        std::dynamic_pointer_cast<lbmproject::GeometryInfoImmersedBoundary>(ptr_amr_instance->
        ptr_criterion_manager_->vec_ptr_geometries_.at(0));
    ptr_geo_tmp->SetWriteIBForce(true);
    ptr_geo_tmp->ptr_geo_shape_ = std::make_unique<lbmproject::GeoShapeIBCircle2D>(ptr_geo_tmp);
    lbmproject::GeoShapeIBCircle2D* ptr_cylinder =
        dynamic_cast<lbmproject::GeoShapeIBCircle2D*>(ptr_geo_tmp->ptr_geo_shape_.get());
    ptr_cylinder->center_ = center;
    ptr_cylinder->radius_ = diameter/2.;
    ptr_geo_tmp->SetGridExtendType(amrproject::EGridExtendType::kInAndOut);
    ptr_geo_tmp->SetLevel(max_refinement_level);
    ptr_geo_tmp->SetXExtendPositive({ 16, 16, 16 });
    ptr_geo_tmp->SetXExtendNegative({ 8, 8 , 8 });
    ptr_geo_tmp->SetYExtendPositive({ 8, 8, 8 });
    ptr_geo_tmp->SetYExtendNegative({ 8, 8, 8 });
    ptr_geo_tmp->SetInnerExtend({ 8, 8 });
    ptr_geo_tmp->SetStencilDis(1);

    // set grid parameters
    ptr_amr_instance->ptr_grid_manager_->vec_ptr_tracking_info_creator_.push_back(
        std::make_unique<amrproject::TrackingGridInfoCreatorInterface>());
    ptr_geo_tmp->SetPtrTrackingGridInfoCreator(
        ptr_amr_instance->ptr_grid_manager_->vec_ptr_tracking_info_creator_.at(0).get());
    ptr_amr_instance->ptr_grid_manager_->SetDomainSize(domain_size);
    ptr_amr_instance->ptr_grid_manager_->SetDomainGridSize({grid_size});
    ptr_amr_instance->SetupInputDependentParameters();
    ptr_amr_instance->AddSolverToGridManager(solver_creator);
    ptr_amr_instance->SetSameSolverDependentInfoForAllGrids(
        ptr_amr_instance->ptr_grid_manager_->vec_ptr_solver_.at(0));

    // set solver parameters
    lbmproject::SolverLbmInterface& solver_ref = *dynamic_cast<lbmproject::SolverLbmInterface*>(
            ptr_amr_instance->ptr_grid_manager_->vec_ptr_solver_.at(0).get());
    for (auto& iter_grid : ptr_amr_instance->ptr_grid_manager_->vec_ptr_grid_info_) {
        lbmproject::GridInfoLbmInteface& grid_ref
            = *dynamic_cast<lbmproject::GridInfoLbmInteface*>(iter_grid.get());
        solver_ref.SetCollisionOperator(grid_ref.GetGridLevel(), collision_type);
        solver_ref.SetNumForces(dims);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryXMin,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryXMax,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryYMin,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryYMax,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        for (auto& iter_boundary : grid_ref.domain_boundary_condition_) {
            iter_boundary.second->SetValues(velocity);
        }
    }
    const DefReal nu = velocity.at(0)* diameter/grid_size/Re;
    solver_ref.SetDefaultViscosity(nu);
    solver_ref.SetDefaultVelocity(velocity);

    // initialize grids
    ptr_amr_instance->InitializeAllSolvers();
    ptr_amr_instance->SetupMesh(0);

    // check mesh
    ptr_amr_instance->CheckMeshAfterInitialization();

    for (DefAmrLUint it = 0; it < max_t; ++it) {
        ptr_amr_instance->TimeMarching(it);
    }

    ptr_amr_instance->FinalizeSimulation();

    return 0;
}

