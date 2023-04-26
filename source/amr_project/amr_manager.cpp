//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include <filesystem>
#include "io/log_write.h"
#include "auxiliary_inline_func.h"
#include "amr_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
void AmrManager::LoadModules(DefUint dims) {
#ifdef ENABLE_MPI
    ptr_mpi_manager_ = std::make_unique<MpiManager>();
#endif
    if (dims == 2) {
        ptr_grid_manager_ = std::make_unique<GridManager2D>();
    } else {
        ptr_grid_manager_ = std::make_unique<GridManager3D>();
    }
    ptr_io_manager_ = std::make_unique<IoManager>();
    ptr_criterion_manager_ = std::make_unique<CriterionManager>();
}
/**
* @brief      function to set default parameters for all modules.
* @param[in]  argc    number of inputs from command line.
* @param[in]  argv    inputs from command line.
*/
void AmrManager::DefaultInitialization(DefUint dims, DefSizet max_level,
    int argc, char* argv[]) {
    //name of the current executable
    std::filesystem::path file_name =
        std::filesystem::path(argv[0]).filename();
    file_name.replace_extension();
    program_name_ = file_name.generic_string();

    LoadModules(dims);

    // mpi settings
#ifdef ENABLE_MPI
    ptr_mpi_manager_->StartupMpi(argc, argv);
#endif  // ENABLE_MPI

    ptr_io_manager_->DefaultInitialization();
    ptr_grid_manager_->DefaultInitialization(max_level);
}
/**
* @brief   function to set parameters according to inputs for all modules
*/
void AmrManager::SetupParameters() {
    // setup grid parametres
    ptr_grid_manager_->SetGridParameters();

#ifdef ENABLE_MPI
    // setup mpi parametres
    ptr_mpi_manager_->SetMpiParameters();
#endif

    ptr_io_manager_->SetIoParameters();
}
/**
* @brief   function to initialize simulation
*/
void AmrManager::InitializeSimulation() {
    int rank_id = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI

    if (rank_id == 0) {
        /*ptr_criterion_manager_->InitialGeometrySerial(
            ptr_grid_manager_->k0RealOffset_);*/
        ptr_grid_manager_->GenerateInitialMeshBasedOnGeoSerial(
            ptr_criterion_manager_->vec_ptr_geometries_);
    }
}
/**
* @brief   function to create the same type of grid
*          instance for all levels of refinement.
* @param[in]  node_type  type of grid node.
* @note  just for convience. Grid instance can be specified for each level.
*/
void AmrManager::SetTheSameLevelDependentInfoForAllLevels(
    SolverCreatorInterface* const ptr_solver_creator) {
    std::shared_ptr<SolverInterface> ptr_solver =
        ptr_solver_creator->CreateSolver();
    ptr_grid_manager_->CreateSameGridInstanceForAllLevel(
        ptr_solver->ptr_grid_info_creator_.get());
}
/**
* @brief   function to finalize simulation
*/
void AmrManager::FinalizeSimulation() {
    //ptr_io_manager_->OutputFlowfield(
    //    program_name_, ptr_grid_manager_.get(),
    //    ptr_criterion_manager_.get());
}
}  // end amrproject
}  // end namespace rootproject
