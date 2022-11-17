//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include "io/log_write.h"
#include "auxiliary_inline_func.h"
#include "amr_manager.h"
namespace rootproject {
namespace amrproject {
void AmrManager::LoadModules() {
#ifdef ENABLE_MPI
    ptr_mpi_manager = std::make_shared<mpi::MpiManager> ptr_mpi_manager();
#endif
    //ptr_grid_manager_ = std::make_shared<grid::GridManagerInterface>();
    ptr_io_manager_ = std::make_shared<io::IoManager>();
    ptr_criterion_manager_ = std::make_shared<criterion::CriterionManager>();
}
/**
* @brief      function to set default parameters for all modules.
* @param[in]  argc    number of inputs from command line.
* @param[in]  argv    inputs from command line.
* @note
*/
void AmrManager::DefaultInitialization(DefUint dims, DefSizet max_level,
    int argc, char* argv[]) {

    LoadModules();

    // mpi settings
#ifdef ENABLE_MPI
    ptr_mpi_manager_->StartupMpi(argc, argv);
#endif  // ENABLE_MPI

    ptr_io_manager_->DefaultInitialization();
    ptr_grid_manager_->DefaultInitialization(max_level);
    ptr_criterion_manager_->k0GeoDims_ = dims;
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
void AmrManager::InitialSimulation() {
    SetupParameters();
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManagerrank_id_;
#endif  // ENABLE_MPI

    if (rank_id == 0) {
        /*ptr_criterion_manager_->InitialGeometrySerial(
            ptr_grid_manager_->k0RealOffest_);*/
        //ptr_grid_manager_->GenerateGridSerial(
        //    ptr_criterion_manager_->vec_ptr_geometries_);
    }
}
/**
* @brief   function to create the same type of grid
*          instance for all levels of refinement.
* @param[in]  node_type  type of grid node.
* @note  just for convience. Grid instance can be specified for each level.
*/
void AmrManager::SetTheSameLevelDependentInfoForAllLevels(
    std::shared_ptr<SolverCreatorInterface> ptr_solver_creator) {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::GetInstance()->rank_id_;
#endif  // ENABLE_MPI
    std::shared_ptr<SolverInterface> ptr_sovler =
        ptr_solver_creator->CreateSolver();
    for (DefSizet i_level = 0; i_level < ptr_grid_manager_->k0MaxLevel_ + 1;
        ++i_level) {
        ptr_grid_manager_->vec_ptr_grid_info_.emplace_back(
            ptr_sovler->ptr_grid_info_creator_->CreateGridInfo());
        grid::GridInfoInterface& grid_ref =
            *(ptr_grid_manager_->vec_ptr_grid_info_).back();
        grid_ref.i_level_ = i_level;
        // link sovler to grid at i_level of refinement 
        grid_ref.ptr_solver_ = ptr_sovler;
        // set grid space in all directions
        grid_ref.grid_space_ = std::vector<DefReal>(
            ptr_grid_manager_->k0GridDims_, 0.);
        //for (DefUint idim = 0; idim < ptr_grid_manager_->k0GridDims_; ++idim) {
        //    grid_ref.grid_space_.at(idim) =
        //        ptr_grid_manager_->k0DomainDx_.at(idim) /
        //        static_cast<DefReal>(TwoPowerN(i_level));
        //}
        // set computational cost for each node 2^i_level
        grid_ref.computational_cost_ =
            static_cast<DefReal>(TwoPowerN(i_level));
        grid_ref.set_number_of_vec_elements();
    }
}
/**
* @brief   function to finalize simulation
*/
void AmrManager::FinalizeSimulation() {
    ptr_io_manager_->OutputFlowfield(ptr_grid_manager_,
        ptr_criterion_manager_);
}
}  // end amrproject
}  // end namespace rootproject
