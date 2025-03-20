//    Copyright (c) 2021 - 2024, Zhengliang Liu
//    All rights reserved

/**
* @file sphere3d.cpp
* @author Zhengliang Liu
* @date  2025-1-29
* @brief  flow past a 3D shpere using IB-LBM.
*/
#include "./amr_manager.h"
#include "d3q19_model/lbm_d3q19.h"
#include "immersed_boundary/immersed_boundary.h"
#include "immersed_boundary/geometry_ib_shape.h"
using namespace rootproject;
int main(int argc, char* argv[]) {
    // case related settings
    const DefInt dims = 3;
    std::string input_name = "sphere3d.in";

    lbmproject::SolverCreatorLbmD3Q19 solver_creator = lbmproject::SolverCreatorLbmD3Q19();
    lbmproject::ELbmBoundaryConditionScheme boundary_condition =
        lbmproject::ELbmBoundaryConditionScheme::kNonEqExtrapolation;

    // startup project manager
    amrproject::AmrManager* ptr_amr_instance;
    ptr_amr_instance = amrproject::AmrManager::GetInstance();
    ptr_amr_instance->program_name_ = "sphere3d";

    ptr_amr_instance->StartupAmrManager(argc, argv);
    ptr_amr_instance->StartupInitialization(dims);
    ptr_amr_instance->AddSolverToGridManager(solver_creator);

    // read input file
    DefInt max_t;
    lbmproject::GeoIBTypeReader geo_ib_reader;
    amrproject::InputParser input_parser(input_name);
    ptr_amr_instance->BroadCastInputParse(&input_parser);
    input_parser.GetValue<DefInt>("max_time_step", &max_t);
    ptr_amr_instance->ptr_grid_manager_->ReadAndSetupGridParameters(&input_parser);
    ptr_amr_instance->ptr_criterion_manager_->ReadAndSetGeoParametersBasedOnShape(
        dims, ptr_amr_instance->ptr_grid_manager_->GetMaxLevel(), &input_parser, geo_ib_reader);
    ptr_amr_instance->ptr_grid_manager_->vec_ptr_solver_.at(0)->ReadAndSetupSolverParameters(&input_parser);

    if (ptr_amr_instance->ptr_mpi_manager_ != nullptr) {
        if (ptr_amr_instance->ptr_mpi_manager_->GetRankId() == 0) {
            input_parser.PrintUnusedParameters();
        }
    } else {
        input_parser.PrintUnusedParameters();
    }

    // use default tracking node type for all geometries
    ptr_amr_instance->ptr_grid_manager_->vec_ptr_tracking_info_creator_.push_back(
        std::make_unique<amrproject::TrackingGridInfoCreatorInterface>());
    for (auto& iter_geo : ptr_amr_instance->ptr_criterion_manager_->vec_ptr_geometries_) {
        iter_geo->SetPtrTrackingGridInfoCreator(
            ptr_amr_instance->ptr_grid_manager_->vec_ptr_tracking_info_creator_.at(0).get());
    }

    ptr_amr_instance->SetupDependentParameters();

    lbmproject::SolverLbmInterface& solver_ref = *dynamic_cast<lbmproject::SolverLbmInterface*>(
        ptr_amr_instance->ptr_grid_manager_->vec_ptr_solver_.at(0).get());
    ptr_amr_instance->SetSameSolverDependentInfoForAllGrids(
        ptr_amr_instance->ptr_grid_manager_->vec_ptr_solver_.at(0));


    // set boundary conditions
    for (auto& iter_grid : ptr_amr_instance->ptr_grid_manager_->vec_ptr_grid_info_) {
        lbmproject::GridInfoLbmInteface& grid_ref
            = *dynamic_cast<lbmproject::GridInfoLbmInteface*>(iter_grid.get());
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryXMin,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryXMax,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryYMin,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryYMax,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryZMin,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        solver_ref.SetDomainBoundaryCondition(amrproject::EDomainBoundaryDirection::kBoundaryZMax,
            boundary_condition, &grid_ref.domain_boundary_condition_);
        for (auto& iter_boundary : grid_ref.domain_boundary_condition_) {
            iter_boundary.second->SetValues(solver_ref.GetDefaultVelocity());
        }
    }

    // initialize grids
    ptr_amr_instance->InitializeAllSolvers();
    ptr_amr_instance->InitializeMesh();

    // check mesh
    ptr_amr_instance->CheckMeshAfterInitialization();

    for (DefInt it = 0; it < max_t; ++it) {
        ptr_amr_instance->TimeMarching(it);
    }

    ptr_amr_instance->FinalizeSimulation();

    return 0;
}

