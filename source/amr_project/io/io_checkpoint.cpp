//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file io_checkpoint.cpp
* @author Zhengliang Liu
* @brief functions used to manage IO processes for grids.
*/
#include <string>
#include <filesystem>
#include "io/log_write.h"
#include "io/io_manager.h"
#include "criterion/criterion_manager.h"
#include "grid/grid_manager.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @brief function to write checkpoint data for the purpose of restart.
* @param[in]  time_step  time step writing checkpoint.
* @param[in]  prog_name  name of the program.
* @param[in]  sfbitset_aux  class manages space filling curves.
* @param[in]  grid_manager  class used to manage grid information.
* @param[in]  criterion_manager  class used to manage criterion information.
* @param[in]  mpi_manager  class used to manage mpi processes.
*/
void IoManager::WriteCheckPointData(const DefAmrLUint time_step, const std::string& prog_name,
    const SFBitsetAuxInterface& sfbitset_aux,
    const GridManagerInterface& grid_manager,
    const CriterionManager& criterion_manager, const MpiManager& mpi_manager) const {
    std::string folder_name = prog_name + "_checkpoint" + std::to_string(time_step);
    if (mpi_manager.GetRankId() == 0) {
        if (!std::filesystem::exists(folder_name)) {
            std::filesystem::create_directories(folder_name);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const DefInt max_level = grid_manager.GetMaxLevel();
    std::map<DefSFCodeToUint, BackgroundLoadData> cost_node_level;
    std::vector<DefAmrLUint> num_of_nodes, num_of_interface_nodes;

    DefAmrLUint cost_current_rank = mpi_manager.ComputeComputationalLoadOnEachRank(max_level,
        sfbitset_aux, grid_manager.vec_ptr_grid_info_, &cost_node_level, &num_of_nodes, &num_of_interface_nodes);
    mpi_manager.WriteCheckPointNodesAtWhichLevels(max_level, cost_current_rank,
        folder_name + "/loadings", cost_node_level);
    mpi_manager.WriteCheckPointGridNodes(folder_name + "/grid_nodes", num_of_nodes,
        cost_node_level, grid_manager.vec_ptr_grid_info_);
    mpi_manager.WriteCheckPointInterfaceNodes(folder_name + "/interface_nodes", num_of_interface_nodes,
        cost_node_level, grid_manager.vec_ptr_grid_info_);
}
/**
* @brief function to read checkpoint data for the purpose of restart.
* @param[in]  time_step  time step restarting from.
* @param[in]  prog_name  name of the program.
* @param[out]  ptr_grid_manager  pointer to class used to manage grid information.
* @param[out]  ptr_criterion_manager   pointer to class used to manage criterion information.
* @param[out]  ptr_mpi_manager   pointer to class used to manage mpi processes.
*/
void IoManager::ReadCheckPointData(const DefAmrLUint time_step, const std::string& prog_name,
    GridManagerInterface* const ptr_grid_manager,
    CriterionManager* const ptr_criterion_manager, MpiManager* const ptr_mpi_manager) const {
    std::string folder_name = prog_name + "_checkpoint" + std::to_string(time_step);

    if (!(std::filesystem::exists(folder_name) && std::filesystem::is_directory(folder_name))) {
        LogManager::LogError("Checkpoint folder does not exist: " + folder_name
            + " in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    }

    // read and calculate range of space filling code for each rank
    std::map<DefSFCodeToUint, BackgroundLoadData> cost_node_level;
    DefAmrLUint cost_all = ptr_mpi_manager->ReadCheckPointNodesAtWhichLevels(
        folder_name + "/loadings", &cost_node_level);
    std::vector<DefInt> vec_cost = ptr_grid_manager->GetNodeCostAtEachLevel();
    std::vector<DefAmrLUint> num_of_grid_nodes, num_of_interface_nodes;
    ptr_mpi_manager->ComputeMinNMaxSFbitsetForEachRank(cost_all, vec_cost, cost_node_level,
        &num_of_grid_nodes, &num_of_interface_nodes);

    const SFBitsetAuxInterface& sfbitset_aux = *ptr_grid_manager->GetPtrToSFBitsetAux();
    const DefSFCodeToUint code_min_background_level =
        sfbitset_aux.SFBitsetToSFCode(ptr_mpi_manager->GetSFBitsetMinCurrentRank()),
        code_max_background_level =
        sfbitset_aux.SFBitsetToSFCode(ptr_mpi_manager->GetSFBitsetMaxCurrentRank());
    DefMap<DefInt> partition_interface_background;
    if (ptr_grid_manager->k0GridDims_ == 2) {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
        GridManager2D* ptr_grid_manager_2d = dynamic_cast<GridManager2D*>(ptr_grid_manager);
        ptr_mpi_manager->IniFindInterfaceForPartitionFromMinNMax(
            code_min_background_level, code_max_background_level,
            ptr_grid_manager_2d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_2d->k0MaxIndexOfBackgroundNode_,
            dynamic_cast<const SFBitsetAux2D&>(sfbitset_aux), &partition_interface_background);
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
    } else {
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
        GridManager3D* ptr_grid_manager_3d = dynamic_cast<GridManager3D*>(ptr_grid_manager);
        ptr_mpi_manager->IniFindInterfaceForPartitionFromMinNMax(
            code_min_background_level, code_max_background_level,
            ptr_grid_manager_3d->k0MinIndexOfBackgroundNode_, ptr_grid_manager_3d->k0MaxIndexOfBackgroundNode_,
            dynamic_cast<const SFBitsetAux3D&>(sfbitset_aux), &partition_interface_background);
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
    }

    // read node data from checkpoint file
    std::vector<DefMap<DefInt>> mpi_interface_each_level;
    ptr_mpi_manager->ReadCheckPointGridNodes(folder_name + "/grid_nodes", num_of_grid_nodes,
        partition_interface_background, sfbitset_aux, &mpi_interface_each_level, ptr_grid_manager);
    std::vector<std::map<int, DefMap<DefInt>>> outer_mpi_nodes;
    ptr_mpi_manager->SearchForMpiOuterLayerBasedOnInterface(mpi_interface_each_level,
        sfbitset_aux, *ptr_grid_manager, &outer_mpi_nodes);
    ptr_mpi_manager->SendNReceiveCheckPointMpiLayers(outer_mpi_nodes, ptr_grid_manager);
    std::vector<InterfaceIndexForMpi> vec_interface_index;
    std::map<int, std::pair<int, DefMap<std::vector<InterfaceIndexForMpi>>>> interface_to_send;
    ptr_mpi_manager->ReadCheckPointInterfaceNodes(folder_name + "/interface_nodes", num_of_interface_nodes,
        &vec_interface_index, &ptr_grid_manager->vec_ptr_grid_info_, &interface_to_send);

    // communicate nodes in mpi communication layers
    ptr_mpi_manager->CommunicateCheckPointMpiLayers(ptr_grid_manager);
    ptr_mpi_manager->CommunicateCheckPointInterfaceLayers(
        vec_interface_index, interface_to_send, ptr_grid_manager);

    // add nodes for interpolation
    DefInt max_level = ptr_grid_manager->GetMaxLevel();
    DefInt dims = ptr_grid_manager->k0GridDims_;
    for (DefInt i_level = 1; i_level <= max_level; ++i_level) {
        GridInfoInterface& grid_info = *ptr_grid_manager->vec_ptr_grid_info_.at(i_level);
        std::vector<bool> periodic_min(dims, false), periodic_max(dims, false);;
        grid_info.CheckIfPeriodicDomainRequired(dims, &periodic_min, &periodic_max);
        const DefInt maxlayer = grid_info.GetNumFine2CoarseLayer();
        for (auto& iter_interface : grid_info.map_ptr_interface_layer_info_) {
            for (DefInt i_layer = maxlayer - grid_info.GetNumFine2CoarseGhostLayer();
                i_layer < maxlayer; ++i_layer) {
                if (iter_interface.second->vec_inner_fine2coarse_.size() > i_layer) {
                    grid_info.AddGhostNodesForInterpolation(periodic_min, periodic_max,
                        sfbitset_aux, iter_interface.second->vec_inner_fine2coarse_.at(i_layer),
                        ptr_grid_manager->vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_);
                }
                if (iter_interface.second->vec_outer_fine2coarse_.size() > i_layer) {
                    grid_info.AddGhostNodesForInterpolation(periodic_min, periodic_max,
                        sfbitset_aux, iter_interface.second->vec_outer_fine2coarse_.at(i_layer),
                        ptr_grid_manager->vec_ptr_grid_info_.at(i_level - 1)->map_grid_node_);
                }
            }
        }
    }

    for (DefInt i_level = 0; i_level < max_level; ++i_level) {
        if (max_level > 0) {
            amrproject::GridInfoInterface& grid_info =
                *ptr_grid_manager->vec_ptr_grid_info_.at(i_level+ 1);
            ptr_mpi_manager->SendNReceiveSFbitsetForInterpolation(
                i_level + 1, sfbitset_aux,
                grid_info.interp_nodes_outer_layer_, &grid_info.vec_num_interp_nodes_receive_,
                &grid_info.interp_nodes_inner_layer_);
            for (const auto& iter_layer : grid_info.interp_nodes_inner_layer_) {
                for (const auto& iter_node : iter_layer.second) {
                    if (grid_info.map_grid_node_.find(iter_node.first) != grid_info.map_grid_node_.end()) {
                        grid_info.map_grid_node_.at(iter_node.first)->flag_status_
                            |= amrproject::NodeBitStatus::kNodeStatusMpiInterpInner_;
                    }
                }
            }
        }
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
