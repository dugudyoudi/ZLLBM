//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file obj_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage all processes.
* @date  2022-5-16
* @note .
*/
#include <memory>
#include <filesystem>
#include <vector>
#include "io/log_write.h"
#include "auxiliary_inline_func.h"
#include "amr_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
void AmrManager::LoadModules(DefAmrIndexUint dims) {
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
* @param[in]  dims    dimension of the mesh.
* @param[in]  max_level    maxim refinement level.
* @param[in]  argc    number of inputs from command line.
* @param[in]  argv    inputs from command line.
*/
void AmrManager::DefaultInitialization(DefAmrIndexUint dims, DefAmrIndexUint max_level, int argc, char* argv[]) {
    //  name of the current executable
    std::filesystem::path file_name =
        std::filesystem::path(argv[0]).filename();
    file_name.replace_extension();
    program_name_ = file_name.generic_string();

    LoadModules(dims);

    // mpi settings
#ifdef ENABLE_MPI
    ptr_mpi_manager_->StartupMpi(argc, argv);
#endif  // ENABLE_MPI

    LogManager::LogStartTime();

    ptr_io_manager_->DefaultInitialization();
    ptr_grid_manager_->DefaultInitialization(max_level);
}
/**
* @brief   function to set parameters according to inputs for all modules
*/
void AmrManager::SetupParameters() {
    // setup grid parameters
    ptr_grid_manager_->SetGridParameters();

#ifdef ENABLE_MPI
    // setup mpi parameters
    ptr_mpi_manager_->SetMpiParameters();
#endif

    ptr_io_manager_->SetIoParameters();
}
/**
* @brief   function to initialize simulation
*/
void AmrManager::InitializeMesh() {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = ptr_mpi_manager_->rank_id_;
#endif  // ENABLE_MPI

    std::array<DefSFBitset, 2> sfbitset_bound_current;
    std::vector<DefMap<DefAmrIndexUint>> sfbitset_one_lower_level(ptr_grid_manager_->k0MaxLevel_ + 1);
    std::vector<DefReal> real_offset(ptr_grid_manager_->k0GridDims_);
    DefAmrIndexUint i_geo = 0;
    if (rank_id == 0) {
        ptr_criterion_manager_->InitialAllGeometrySerial(ptr_grid_manager_->k0GridDims_, real_offset);
    }
    for (const auto& iter_geo_info : ptr_criterion_manager_->vec_ptr_geometries_) {
        ptr_grid_manager_->CreateTrackingGridInstanceForAGeo(i_geo, *iter_geo_info);
        ++i_geo;
    }
    if (rank_id == 0) {
        ptr_grid_manager_->GenerateGridFromHighToLowLevelSerial(
         ptr_criterion_manager_->vec_ptr_geometries_, &sfbitset_one_lower_level);
        sfbitset_bound_current = {ptr_grid_manager_->k0SFBitsetDomainMin_, ptr_grid_manager_->k0SFBitsetDomainMax_};
    }

#ifdef ENABLE_MPI
    // mpi partition sending and receiving nodes
    std::vector<DefAmrIndexLUint> vec_cost;
    for (auto iter_grid : ptr_grid_manager_->vec_ptr_grid_info_) {
        vec_cost.push_back(iter_grid->computational_cost_);
    }
    std::vector<DefMap<DefAmrIndexUint>> sfbitset_one_lower_level_current_rank(ptr_grid_manager_->k0MaxLevel_ + 1),
       sfbitset_ghost_one_lower_level_current_rank(ptr_grid_manager_->k0MaxLevel_ + 1);
    ptr_mpi_manager_->sfbitset_min_current_rank_ = sfbitset_bound_current.at(0);
    ptr_mpi_manager_->sfbitset_max_current_rank_ = sfbitset_bound_current.at(1);
    DefMap<DefAmrIndexUint> partition_interface_background;
    if (ptr_grid_manager_->k0GridDims_ == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
        GridManager2D* ptr_grid_manager_2d = dynamic_cast<GridManager2D*>(ptr_grid_manager_.get());
        ptr_mpi_manager_->SendNReceiveGridInfoAtGivenLevels(ptr_grid_manager_->kFlagSize0_,
            ptr_grid_manager_->kNodeStatusCoarse2Fine0_,
            ptr_grid_manager_->k0GridDims_, ptr_grid_manager_->k0MaxLevel_,
            ptr_grid_manager_->k0SFBitsetDomainMin_, ptr_grid_manager_->k0SFBitsetDomainMax_,
            {ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_[kXIndex],
            ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_[kYIndex]},
            {ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_[kXIndex],
            ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_[kYIndex]},
            vec_cost, *ptr_grid_manager_2d, ptr_grid_manager_->vec_ptr_tracking_info_creator_,
            sfbitset_one_lower_level, &sfbitset_bound_current, &sfbitset_one_lower_level_current_rank,
            &sfbitset_ghost_one_lower_level_current_rank, &(ptr_grid_manager_->vec_ptr_grid_info_));



            // for (const auto& iter_ghost : sfbitset_ghost_one_lower_level_current_rank.at(2)) {
            //     std::vector<DefReal> indices, indices1;
            //     ptr_grid_manager_2d->SFBitsetComputeCoordinateVir(iter_ghost.first, {0.01, 0.01}, &indices);
            //     if (rank_id == 1 && std::fabs(indices[0] - 0.56) < 0.001) {
            //     ptr_grid_manager_2d->SFBitsetComputeCoordinateVir(iter_ghost.first, {0.01, 0.01}, &indices1);
            //     std::cout << indices1[0] - 0.02 << " " << indices1[1] - 0.02 << std::endl;
            //     }
            // }

        if (rank_id == 0) {
            sfbitset_one_lower_level.clear();
            sfbitset_one_lower_level.shrink_to_fit();
        }
        ptr_mpi_manager_->FindInterfaceForPartitionFromMinNMax(
            ptr_mpi_manager_->vec_sfbitset_min_all_ranks_.at(rank_id),
            ptr_mpi_manager_->vec_sfbitset_max_all_ranks_.at(rank_id),
            ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_,
            *ptr_grid_manager_2d, &partition_interface_background);

        ptr_grid_manager_->InstantiateGridNodeAllLevelMpi(rank_id, ptr_mpi_manager_->k0NumPartitionInnerLayers_,
            ptr_mpi_manager_->k0NumPartitionOuterLayers_, ptr_mpi_manager_->vec_sfbitset_min_all_ranks_,
            ptr_mpi_manager_->vec_sfbitset_max_all_ranks_, *ptr_grid_manager_2d, sfbitset_one_lower_level_current_rank,
            sfbitset_ghost_one_lower_level_current_rank, partition_interface_background,
            &ptr_mpi_manager_->mpi_communication_inner_layers_, &ptr_mpi_manager_->mpi_communication_outer_layers_);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
        GridManager3D* ptr_grid_manager_3d = dynamic_cast<GridManager3D*>(ptr_grid_manager_.get());
        ptr_mpi_manager_->SendNReceiveGridInfoAtGivenLevels(ptr_grid_manager_->kFlagSize0_,
            ptr_grid_manager_->kNodeStatusCoarse2Fine0_,
            ptr_grid_manager_->k0GridDims_, ptr_grid_manager_->k0MaxLevel_,
            ptr_grid_manager_->k0SFBitsetDomainMin_, ptr_grid_manager_->k0SFBitsetDomainMax_,
            {ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_[kXIndex],
            ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_[kYIndex],
            ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_[kZIndex]},
            {ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_[kXIndex],
            ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_[kYIndex],
            ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_[kZIndex]},
            vec_cost, *ptr_grid_manager_3d, ptr_grid_manager_->vec_ptr_tracking_info_creator_,
            sfbitset_one_lower_level, &sfbitset_bound_current, &sfbitset_one_lower_level_current_rank,
            &sfbitset_ghost_one_lower_level_current_rank, &(ptr_grid_manager_->vec_ptr_grid_info_));

        if (rank_id == 0) {
            sfbitset_one_lower_level.clear();
            sfbitset_one_lower_level.shrink_to_fit();
        }
        ptr_mpi_manager_->FindInterfaceForPartitionFromMinNMax(
            ptr_mpi_manager_->vec_sfbitset_min_all_ranks_.at(rank_id),
            ptr_mpi_manager_->vec_sfbitset_max_all_ranks_.at(rank_id),
            ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_,
            *ptr_grid_manager_3d, &partition_interface_background);

        ptr_grid_manager_->InstantiateGridNodeAllLevelMpi(rank_id, ptr_mpi_manager_->k0NumPartitionInnerLayers_,
            ptr_mpi_manager_->k0NumPartitionOuterLayers_, ptr_mpi_manager_->vec_sfbitset_min_all_ranks_,
            ptr_mpi_manager_->vec_sfbitset_max_all_ranks_, *ptr_grid_manager_3d, sfbitset_one_lower_level_current_rank,
            sfbitset_ghost_one_lower_level_current_rank, partition_interface_background,
            &ptr_mpi_manager_->mpi_communication_inner_layers_, &ptr_mpi_manager_->mpi_communication_outer_layers_);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    }

    // add nodes on both the refinement and partition interfaces
    // which are only stored in coarse to fine refinement interfaces
    const DefAmrIndexUint flag0 = ptr_grid_manager_->kFlagSize0_;
    for (DefAmrIndexUint i_level = 0; i_level < ptr_grid_manager_->k0MaxLevel_; ++i_level) {
        for (const auto & iter_interfaces :
         ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_) {
            for (const auto & iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
                for (const auto & iter_node : iter_coarse2fine) {
                    if (ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_grid_node_.find(iter_node.first)
                        == ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_grid_node_.end()) {
                        ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_grid_node_.insert(
                            {iter_node.first, ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->k0GridNodeInstance_});
                    }
                }
            }
        }
    }
#else   // mesh on rank 0 is the only one when run serially
    ptr_grid_manager_->InstantiateGridNodeAllLevel(
     sfbitset_bound_current.at(0), sfbitset_bound_current.at(1), sfbitset_one_lower_level);
#endif  // ENABLE_MPI
}
/**
* @brief   function to create the same type of grid
*          instance for all levels of refinement.
* @param[in]  node_type  type of grid node.
* @note  just for convience. Grid instance can be specified for each level.
*/
void AmrManager::SetDependentInfoForAllLevelsTheSame(
    SolverCreatorInterface* const ptr_solver_creator) {
    std::shared_ptr<SolverInterface> ptr_solver = ptr_solver_creator->CreateSolver();
    ptr_grid_manager_->CreateSameGridInstanceForAllLevel(ptr_solver->ptr_grid_info_creator_.get());
}
/**
* @brief   function to finalize simulation
*/
void AmrManager::FinalizeSimulation() {
    ptr_io_manager_->OutputFlowfield(program_name_, ptr_grid_manager_.get(), ptr_criterion_manager_.get());
#ifdef ENABLE_MPI
    ptr_mpi_manager_->FinalizeMpi();
#endif  // ENABLE_MPI
}
}  // end namespace amrproject
}  // end namespace rootproject
