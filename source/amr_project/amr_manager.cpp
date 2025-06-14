//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file amr_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage modules for amr project.
* @date  2022-5-16
* @note .
*/
#include <memory>
#include <filesystem>
#include <vector>
#include "io/log_write.h"
#include "./auxiliary_inline_func.h"
#include "./amr_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief      function to setup the amr program based on command line inputs.
* @param[in]  argc    number of inputs from command line.
* @param[in]  argv    inputs from command line.
*/
int AmrManager::StartupAmrManager(int argc, char* argv[]) {
    //  name of the current executable
    std::filesystem::path file_name = std::filesystem::path(argv[0]).filename();
    file_name.replace_extension();
    program_name_ = file_name.generic_string();

    int flag_mpi =  0;
#ifdef ENABLE_MPI
    flag_mpi = MPI_Init(&argc, &argv);
#endif
    return flag_mpi;
}
/**
* @brief     function to load modules of amr manager.
*/
void AmrManager::LoadModules(DefInt dims) {
#ifdef ENABLE_MPI
    ptr_mpi_manager_ = std::make_unique<MpiManager>();
#endif
    if (dims == 2) {
        ptr_grid_manager_ = std::make_unique<GridManager2D>();
    } else if (dims == 3) {
        ptr_grid_manager_ = std::make_unique<GridManager3D>();
    }
    ptr_io_manager_ = std::make_unique<IoManager>();
    ptr_criterion_manager_ = std::make_unique<CriterionManager>();
}
/**
* @brief      function to set default parameters for all modules.
* @param[in]  dims    dimension of the mesh.
* @param[in]  max_level    maxim refinement level.
*/
void AmrManager::StartupInitialization(DefInt dims) {
    LoadModules(dims);

    // mpi settings
#ifdef ENABLE_MPI
    ptr_mpi_manager_->SetUpMpi();
#endif  // ENABLE_MPI

    LogManager::LogStartTime();

    ptr_io_manager_->SetupOutputFormat();
}
/**
* @brief   function to set parameters dependent on input parameters for all modules
*/
void AmrManager::SetupInputDependentParameters() {
    // setup grid parameters
    ptr_grid_manager_->SetupDependentGridParameters();

    ptr_io_manager_->SetupDependentIOParameters();
}
/**
* @brief   function to initialize mesh.
*/
void AmrManager::InitializeMesh() {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = ptr_mpi_manager_->GetRankId();
#endif  // ENABLE_MPI
    std::array<DefSFBitset, 2> sfbitset_bound_current;
    const DefInt max_level = ptr_grid_manager_->GetMaxLevel();
    std::vector<DefMap<DefInt>> sfbitset_one_lower_level(max_level + 1);

    DefInt i_geo = 0;
    for (const auto& iter_geo_info : ptr_criterion_manager_->vec_ptr_geometries_) {
        ptr_grid_manager_->CreateTrackingGridInstanceForAGeo(i_geo, *iter_geo_info);
        ++i_geo;
    }
    if (rank_id == 0) {
        ptr_grid_manager_->GenerateGridFromHighToLowLevelSerial(
            ptr_criterion_manager_->vec_ptr_geometries_, &sfbitset_one_lower_level);
        sfbitset_bound_current = {ptr_grid_manager_->GetSFbitsetforDomainMin(),
            ptr_grid_manager_->GetSFbitsetforDomainMax()};
    }

#ifdef ENABLE_MPI
    if (rank_id == 0 && sfbitset_one_lower_level.size() > 1) {
        ptr_grid_manager_->ResetBackgroundGridAsMpiLayer(ptr_mpi_manager_->GetNumPartitionOuterLayers(),
            *(ptr_grid_manager_->vec_ptr_grid_info_.at(0)), sfbitset_one_lower_level.at(1),
            &sfbitset_one_lower_level.at(0));
    }
    // mpi partition sending and receiving nodes
    std::vector<DefInt> vec_cost = ptr_grid_manager_->GetNodeCostAtEachLevel();
    std::vector<DefMap<DefInt>> sfbitset_one_lower_level_current_rank(max_level + 1),
       sfbitset_ghost_one_lower_level_current_rank(max_level + 1);
    ptr_mpi_manager_->SetSFBitsetMinCurrentRank(sfbitset_bound_current.at(0));
    ptr_mpi_manager_->SetSFBitsetMaxCurrentRank(sfbitset_bound_current.at(1));
    DefMap<DefInt> partition_interface_background;
    if (ptr_grid_manager_->k0GridDims_ == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
        GridManager2D* ptr_grid_manager_2d = dynamic_cast<GridManager2D*>(ptr_grid_manager_.get());
        ptr_mpi_manager_->IniSendNReceiveGridInfoAtAllLevels(ptr_grid_manager_->kFlagSize0_,
            NodeBitStatus::kNodeStatusCoarse2Fine0_,
            ptr_grid_manager_->k0GridDims_, max_level,
            ptr_grid_manager_->GetSFbitsetforDomainMin(), ptr_grid_manager_->GetSFbitsetforDomainMax(),
            {ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_[kXIndex],
            ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_[kYIndex]},
            {ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_[kXIndex],
            ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_[kYIndex]},
            vec_cost, *ptr_grid_manager_2d, ptr_grid_manager_->vec_ptr_tracking_info_creator_,
            sfbitset_one_lower_level, &sfbitset_bound_current, &sfbitset_one_lower_level_current_rank,
            &sfbitset_ghost_one_lower_level_current_rank, &(ptr_grid_manager_->vec_ptr_grid_info_));
        if (rank_id == 0) {
            sfbitset_one_lower_level.clear();
            sfbitset_one_lower_level.shrink_to_fit();
        }
        ptr_mpi_manager_->IniFindInterfaceForPartitionFromMinNMax(
            ptr_mpi_manager_->GetSFCodeMinAllRanks().at(rank_id),
            ptr_mpi_manager_->GetSFCodeMaxAllRanks().at(rank_id),
            ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_,
            *ptr_grid_manager_2d, &partition_interface_background);

        ptr_grid_manager_->InstantiateGridNodeAllLevelMpi(rank_id, ptr_mpi_manager_->GetNumPartitionInnerLayers(),
            ptr_mpi_manager_->GetNumPartitionOuterLayers(), ptr_mpi_manager_->GetSFCodeMinAllRanks(),
            ptr_mpi_manager_->GetSFCodeMaxAllRanks(),
            *ptr_grid_manager_2d, sfbitset_one_lower_level_current_rank,
            sfbitset_ghost_one_lower_level_current_rank, partition_interface_background,
            &ptr_mpi_manager_->mpi_communication_inner_layers_, &ptr_mpi_manager_->mpi_communication_outer_layers_);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
        GridManager3D* ptr_grid_manager_3d = dynamic_cast<GridManager3D*>(ptr_grid_manager_.get());
        ptr_mpi_manager_->IniSendNReceiveGridInfoAtAllLevels(ptr_grid_manager_->kFlagSize0_,
            NodeBitStatus::kNodeStatusCoarse2Fine0_,
            ptr_grid_manager_->k0GridDims_, max_level,
            ptr_grid_manager_->GetSFbitsetforDomainMin(), ptr_grid_manager_->GetSFbitsetforDomainMax(),
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
        ptr_mpi_manager_->IniFindInterfaceForPartitionFromMinNMax(
            ptr_mpi_manager_->GetSFCodeMinAllRanks().at(rank_id),
            ptr_mpi_manager_->GetSFCodeMaxAllRanks().at(rank_id),
            ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_,
            *ptr_grid_manager_3d, &partition_interface_background);

        ptr_grid_manager_->InstantiateGridNodeAllLevelMpi(rank_id, ptr_mpi_manager_->GetNumPartitionInnerLayers(),
            ptr_mpi_manager_->GetNumPartitionOuterLayers(), ptr_mpi_manager_->GetSFCodeMinAllRanks(),
            ptr_mpi_manager_->GetSFCodeMaxAllRanks(),
            *ptr_grid_manager_3d, sfbitset_one_lower_level_current_rank,
            sfbitset_ghost_one_lower_level_current_rank, partition_interface_background,
            &ptr_mpi_manager_->mpi_communication_inner_layers_, &ptr_mpi_manager_->mpi_communication_outer_layers_);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    }
    for (DefInt i_level = 0; i_level < max_level; ++i_level) {
        // add nodes on both the refinement and partition interfaces
        // which are only stored in coarse to fine refinement interfaces
        for (const auto & iter_interfaces :
            ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_ptr_interface_layer_info_) {
            for (const auto& iter_coarse2fine : iter_interfaces.second->vec_outer_coarse2fine_) {
                for (const auto& iter_node : iter_coarse2fine) {
                    if (ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_grid_node_.find(iter_node.first)
                        == ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_grid_node_.end()) {
                        ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->map_grid_node_.insert(
                            {iter_node.first, ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->GridNodeCreator()});
                    }
                }
            }
        }
        if (max_level > 0) {
            GridInfoInterface& grid_info = *ptr_grid_manager_->vec_ptr_grid_info_.at(i_level + 1);
            ptr_mpi_manager_->SendNReceiveSFbitsetForInterpolation(i_level + 1, *grid_info.GetPtrSFBitsetAux(),
                grid_info.interp_nodes_outer_layer_, &grid_info.vec_num_interp_nodes_receive_,
                &grid_info.interp_nodes_inner_layer_);
            for (const auto& iter_layer : grid_info.interp_nodes_inner_layer_) {
                for (const auto& iter_node : iter_layer.second) {
                    if (grid_info.map_grid_node_.find(iter_node.first) != grid_info.map_grid_node_.end()) {
                        grid_info.map_grid_node_.at(iter_node.first)->flag_status_
                            |= NodeBitStatus::kNodeStatusMpiInterpInner_;
                    }
                }
            }
        }
    }

#else   // mesh on rank 0 is the only one when run serially
    ptr_mpi_manager_->SetSFBitsetMinCurrentRank(sfbitset_bound_current.at(0),
        sfbitset_bound_current.at(1));
    ptr_grid_manager_->InstantiateGridNodeAllLevel(
        sfbitset_bound_current.at(0), sfbitset_bound_current.at(1), sfbitset_one_lower_level);
#endif  // ENABLE_MPI

    for (auto& iter_geo : ptr_criterion_manager_->vec_ptr_geometries_) {
        // shape information need to be updated during initialization
        bool bool_tmp = iter_geo->GetNeedUpdateShape();
        iter_geo->SetNeedUpdateShape(true);
        iter_geo->SetupGeometryInfo(0., *ptr_mpi_manager_.get(),
            *ptr_grid_manager_->vec_ptr_grid_info_.at(iter_geo->GetLevel()));
        iter_geo->SetNeedUpdateShape(bool_tmp);
    }
}
/**
* @brief   function to set up the mesh for the next time step.
* @param[in]  time_step time step restarting from.
*/
void AmrManager::RestartFromCheckPointMesh(DefInt time_step) {
#ifdef ENABLE_MPI
    ptr_io_manager_->ReadCheckPointData(time_step, program_name_, ptr_grid_manager_.get(),
        ptr_criterion_manager_.get(), ptr_mpi_manager_.get());
#else
    logger::LogError("Restart from checkpoint is not supported in serial mode");
#endif
}
void AmrManager::SetupMesh(DefInt time_step) {
    // initial geometries
    std::array<DefReal, 3> real_offset{};
    std::vector<DefReal> domain_dx = ptr_grid_manager_->GetDomainDxArrAsVec();
    std::vector<DefAmrLUint> domain_min_index = ptr_grid_manager_->GetMinIndexOfBackgroundNodeArrAsVec();
    for (DefInt i_dims = 0; i_dims < ptr_grid_manager_->k0GridDims_; ++i_dims) {
        real_offset.at(i_dims) = domain_min_index[i_dims] * domain_dx[i_dims];
    }
    ptr_criterion_manager_->InitialAllGeometrySerial(domain_dx.at(kXIndex), real_offset);

    if (time_step == 0) {
        InitializeMesh();
    } else {
        RestartFromCheckPointMesh(time_step);
    }

    InstantiateTimeSteppingScheme();
}
/**
* @brief function to setup default grid related parameters.
* @param[in]  max_level  maximum refinement level.
*/
void AmrManager::InstantiateTimeSteppingScheme() {
    switch (k0TimeSteppingType_) {
    case ETimeSteppingScheme::kMultiSteppingC2F:
        ptr_time_stepping_scheme_ = std::make_unique<MultiTimeSteppingC2F>(
            ptr_grid_manager_->GetMaxLevel());
        break;
    default:
        LogManager::LogError("undefined time stepping type in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        break;
    }
}
/**
* @brief      function to manage grid information during one time step at the background refinement level.
* @param[in]  time_step_background   current background time step.
* @param[in]  sfbitset_aux   class manages space filling curves.
*/
void AmrManager::TimeMarching(const DefAmrLUint time_step_background) {
    current_time_step_ = time_step_background;
    // record number of time step at i_level
    std::vector<DefInt> time_step_level(ptr_grid_manager_->GetMaxLevel() + 1, 0);
    DefInt i_level;
    DefReal time_step_current;
    for (auto iter_level = ptr_time_stepping_scheme_->k0TimeSteppingOrder_.begin();
        iter_level != ptr_time_stepping_scheme_->k0TimeSteppingOrder_.end(); ++iter_level) {
        i_level = *iter_level;
        ++time_step_level[i_level];
        time_step_current = ptr_time_stepping_scheme_->GetCurrentTimeStep(
            i_level, time_step_background, time_step_level[i_level]);

        GridInfoInterface& grid_ref = *(ptr_grid_manager_->vec_ptr_grid_info_.at(i_level));
        grid_ref.SetUpGridAtBeginningOfTimeStep(time_step_level[i_level]);
        grid_ref.AdvancingAtCurrentTime(k0TimeSteppingType_, time_step_level[i_level],
            time_step_current, ptr_mpi_manager_.get(), ptr_criterion_manager_.get());
    }
#ifdef ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif  // ENABLE_MPI
}
/**
* @brief   function to initialize solvers.
*/
void AmrManager::InitializeAllSolvers() {
    for (auto& iter_solver : ptr_grid_manager_->vec_ptr_solver_) {
        iter_solver->SolverInitial();
    }
}
/**
* @brief   function to create solver instance in GridManager.
* @param[in]  ptr_solver_creator  pointer to class for instantiating a given solver.
*/
void AmrManager::AddSolverToGridManager(const SolverCreatorInterface& solver_creator) {
    ptr_grid_manager_->vec_ptr_solver_.emplace_back(solver_creator.CreateSolver());
    if (ptr_grid_manager_->k0GridDims_ != ptr_grid_manager_->vec_ptr_solver_.back()->GetSolverDim()) {
        LogManager::LogWarning("Grid dimension (" + std::to_string(ptr_grid_manager_->k0GridDims_)
            +") mismatches with solver dimension ("
            + std::to_string(ptr_grid_manager_->vec_ptr_solver_.back()->GetSolverDim())
            +") for the solver " + ptr_grid_manager_->vec_ptr_solver_.back()->GetSolverMethod());
    }
}
/**
* @brief   function to set solver dependent information as the same.
* @param[in]  ptr_solver  pointer to a solver.
*/
void AmrManager::SetSameSolverDependentInfoForAllGrids(const std::shared_ptr<SolverInterface>& ptr_solver) {
    ptr_grid_manager_->CreateSameGridInstanceForAllLevel(
        ptr_solver, *ptr_solver->ptr_grid_info_creator_.get());
}
/**
* @brief   function to finalize simulation

*/
void AmrManager::FinalizeSimulation() {
    ptr_io_manager_->OutputMeshData(program_name_, ptr_grid_manager_.get(), ptr_criterion_manager_.get());

    if (bool_write_checkpoint_) {
        if (ptr_grid_manager_->k0GridDims_ == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
            GridManager2D& grid_manager_2d = dynamic_cast<GridManager2D&>(*ptr_grid_manager_.get());
            ptr_io_manager_->WriteCheckPointData(current_time_step_, program_name_, grid_manager_2d,
                grid_manager_2d, *ptr_criterion_manager_.get(), *ptr_mpi_manager_.get());
#endif
        } else if (ptr_grid_manager_->k0GridDims_ == 3) {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
            GridManager3D& grid_manager_3d = dynamic_cast<GridManager3D&>(*ptr_grid_manager_.get());
            ptr_io_manager_->WriteCheckPointData(current_time_step_, program_name_, grid_manager_3d,
                grid_manager_3d, *ptr_criterion_manager_.get(), *ptr_mpi_manager_.get());
#endif
        }
    }

    LogManager::LogEndTime();
#ifdef ENABLE_MPI
    ptr_mpi_manager_->FinalizeMpi();
#endif  // ENABLE_MPI
}
void AmrManager::BroadCastInputParse(InputParser* const ptr_input_parser) const {
#ifdef ENABLE_MPI
    MpiManager* ptr_mpi = GetPointerToMpiManager();
    if (ptr_mpi == nullptr) {
        LogManager::LogError("Pointer to mpi manager is nullptr");
    } else {
        ptr_input_parser->BroadcastInputData(*ptr_mpi);
    }
#endif  // ENABLE_MPI
}
}  // end namespace amrproject
}  // end namespace rootproject
