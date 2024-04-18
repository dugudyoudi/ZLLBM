//  Copyright (c) 2021 - 2024, Zhengliang Liu
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
/**
* @brief      function to set program features.
* @param[in]  argc    number of inputs from command line.
* @param[in]  argv    inputs from command line.
*/
int AmrManager::SetUpProgramFeature(int argc, char* argv[]) {
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
*/
void AmrManager::DefaultInitialization(DefAmrIndexUint dims, DefAmrIndexUint max_level) {
    LoadModules(dims);

    // mpi settings
#ifdef ENABLE_MPI
    ptr_mpi_manager_->SetUpMpi();
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

    ptr_io_manager_->SetIoParameters();
}
/**
* @brief   function to initialize mesh.
*/
void AmrManager::InitializeMesh() {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = ptr_mpi_manager_->rank_id_;
#endif  // ENABLE_MPI
    std::array<DefSFBitset, 2> sfbitset_bound_current;
    std::vector<DefMap<DefAmrIndexUint>> sfbitset_one_lower_level(ptr_grid_manager_->k0MaxLevel_ + 1);
    DefAmrIndexUint i_geo = 0;
    if (rank_id == 0) {
        std::vector<DefReal> real_offset(ptr_grid_manager_->k0GridDims_);
        std::vector<DefReal> domain_dx = ptr_grid_manager_->GetDomainDxArrAsVec();
        std::vector<DefAmrIndexLUint> domain_min_index = ptr_grid_manager_->GetMinIndexOfBackgroundNodeArrAsVec();
        for (DefAmrIndexUint i_dims = 0; i_dims < ptr_grid_manager_->k0GridDims_; ++i_dims) {
            real_offset.at(i_dims) = domain_min_index[i_dims] * domain_dx[i_dims];
        }

        ptr_criterion_manager_->InitialAllGeometrySerial(ptr_grid_manager_->k0GridDims_,
            domain_dx.at(kXIndex), real_offset);
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
        ptr_mpi_manager_->IniSendNReceiveGridInfoAtAllLevels(ptr_grid_manager_->kFlagSize0_,
            NodeBitStatus::kNodeStatusCoarse2Fine0_,
            ptr_grid_manager_->k0GridDims_, ptr_grid_manager_->k0MaxLevel_,
            ptr_grid_manager_->k0SFBitsetDomainMin_, ptr_grid_manager_->k0SFBitsetDomainMax_,
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
            ptr_mpi_manager_->vec_sfcode_min_all_ranks_.at(rank_id),
            ptr_mpi_manager_->vec_sfcode_max_all_ranks_.at(rank_id),
            ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_,
            *ptr_grid_manager_2d, &partition_interface_background);

        ptr_grid_manager_->InstantiateGridNodeAllLevelMpi(rank_id, ptr_mpi_manager_->k0NumPartitionInnerLayers_,
            ptr_mpi_manager_->k0NumPartitionOuterLayers_, ptr_mpi_manager_->vec_sfcode_min_all_ranks_,
            ptr_mpi_manager_->vec_sfcode_max_all_ranks_,
            *ptr_grid_manager_2d, sfbitset_one_lower_level_current_rank,
            sfbitset_ghost_one_lower_level_current_rank, partition_interface_background,
            &ptr_mpi_manager_->mpi_communication_inner_layers_, &ptr_mpi_manager_->mpi_communication_outer_layers_);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
        GridManager3D* ptr_grid_manager_3d = dynamic_cast<GridManager3D*>(ptr_grid_manager_.get());
        ptr_mpi_manager_->IniSendNReceiveGridInfoAtAllLevels(ptr_grid_manager_->kFlagSize0_,
            NodeBitStatus::kNodeStatusCoarse2Fine0_,
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
        ptr_mpi_manager_->IniFindInterfaceForPartitionFromMinNMax(
            ptr_mpi_manager_->vec_sfcode_min_all_ranks_.at(rank_id),
            ptr_mpi_manager_->vec_sfcode_max_all_ranks_.at(rank_id),
            ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_,
            *ptr_grid_manager_3d, &partition_interface_background);

        ptr_grid_manager_->InstantiateGridNodeAllLevelMpi(rank_id, ptr_mpi_manager_->k0NumPartitionInnerLayers_,
            ptr_mpi_manager_->k0NumPartitionOuterLayers_, ptr_mpi_manager_->vec_sfcode_min_all_ranks_,
            ptr_mpi_manager_->vec_sfcode_max_all_ranks_,
            *ptr_grid_manager_3d, sfbitset_one_lower_level_current_rank,
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
                            {iter_node.first, ptr_grid_manager_->vec_ptr_grid_info_.at(i_level)->GridNodeCreator()});
                    }
                }
            }
        }
    }
#else   // mesh on rank 0 is the only one when run serially
    ptr_grid_manager_->InstantiateGridNodeAllLevel(
     sfbitset_bound_current.at(0), sfbitset_bound_current.at(1), sfbitset_one_lower_level);
#endif  // ENABLE_MPI

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
            ptr_grid_manager_->k0MaxLevel_);
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
void AmrManager::TimeMarching(const DefAmrIndexLUint time_step_background) {
    // record number of time step at i_level
    std::vector<DefAmrIndexUint> time_step_level(ptr_grid_manager_->k0MaxLevel_ + 1, 0);
    DefAmrIndexUint i_level;
    DefReal time_step_current;
    for (auto iter_level = ptr_time_stepping_scheme_->k0TimeSteppingOrder_.begin();
        iter_level != ptr_time_stepping_scheme_->k0TimeSteppingOrder_.end(); ++iter_level) {
        i_level = *iter_level;
        ++time_step_level[i_level];
        time_step_current = ptr_time_stepping_scheme_->GetCurrentTimeStep(
            i_level, time_step_background, time_step_level[i_level]);
        GridInfoInterface& grid_ref = *(ptr_grid_manager_->vec_ptr_grid_info_.at(i_level));
        grid_ref.SetUpGridAtBeginningOfTimeStep(time_step_level[i_level], ptr_grid_manager_.get());

        // use information from previous time step
        grid_ref.ptr_solver_->InformationFromGridOfDifferentLevel(
            ETimingInOneStep::kStepBegin, k0TimeSteppingType_,
            time_step_level[i_level], *ptr_grid_manager_->GetSFBitsetAuxPtr(), &grid_ref);

#ifdef ENABLE_MPI
        std::vector<MpiManager::BufferSizeInfo> send_buffer_info, receive_buffer_info;
        std::vector<std::vector<MPI_Request>> vec_vec_reqs_send, vec_vec_reqs_receive;
        std::vector<std::unique_ptr<char[]>> vec_ptr_buffer_receive, vec_ptr_buffer_send;
        ptr_mpi_manager_->SendNReceiveGridNodes(&send_buffer_info, &receive_buffer_info,
            &vec_vec_reqs_send, &vec_vec_reqs_receive, &vec_ptr_buffer_send, &vec_ptr_buffer_receive, &grid_ref);
#endif  //  ENABLE_MPI


        // update criterion information for tracking nodes
        for (auto& iter : grid_ref.map_ptr_tracking_grid_info_) {
            if (iter.first.first == ECriterionType::kGeometry) {
                ptr_criterion_manager_->vec_ptr_geometries_.at(iter.first.second)->UpdateGeometry(time_step_current);
            }
        }

        grid_ref.ptr_solver_->RunSolverOnGrid(k0TimeSteppingType_,
            time_step_level[i_level], *ptr_grid_manager_->GetSFBitsetAuxPtr(), &grid_ref);

#ifdef ENABLE_MPI
        ptr_mpi_manager_->WaitAndReadGridNodesFromBuffer(send_buffer_info,
            receive_buffer_info, vec_ptr_buffer_receive,
            &vec_vec_reqs_send, &vec_vec_reqs_receive, &grid_ref);
        MPI_Barrier(MPI_COMM_WORLD);
#endif  //  ENABLE_MPI

        // use information in current time step
        grid_ref.ptr_solver_->InformationFromGridOfDifferentLevel(
            ETimingInOneStep::kStepEnd, k0TimeSteppingType_,
            time_step_level[i_level], *ptr_grid_manager_->GetSFBitsetAuxPtr(), &grid_ref);

        grid_ref.ptr_solver_->FinalizeAtTimeStepEnd(k0TimeSteppingType_,
            time_step_level[i_level], *ptr_grid_manager_->GetSFBitsetAuxPtr(), &grid_ref);
    }
}
/**
* @brief   function to initialize solvers.
*/
void AmrManager::SetupSolverForGrids() {
    for (auto& iter_solver : ptr_grid_manager_->vec_ptr_solver_) {
        iter_solver->SolverInitial();
    }
}
/**
* @brief   function to create solver instance in GridManager.
* @param[in]  ptr_solver_creator  pointer to class for instantiating a given solver.
*/
void AmrManager::AddSolverToGridManager(const SolverCreatorInterface& solver_creator) {
    ptr_grid_manager_->vec_ptr_solver_.push_back(solver_creator.CreateSolver());
}
/**
* @brief   function to set solver dependent information as the same.
* @param[in]  ptr_solver  pointer to a solver.
*/
void AmrManager::SetDependentInfoForAllLevelsTheSame(const std::shared_ptr<SolverInterface>& ptr_solver) {
    ptr_grid_manager_->CreateSameGridInstanceForAllLevel(
        ptr_solver, *ptr_solver->ptr_grid_info_creator_.get());
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
