//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_grid_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM grid interface.
*/
#include "./lbm_interface.h"
#include "io/log_write.h"
#include "io/vtk_writer.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to check domain boundary is periodic.
 * @param[in] dims dimensions of grid.
 * @param[out] ptr_periodic_min pointer to vector indicating periodic boundaries at minimum.
 * @param[out] ptr_periodic_max pointer to vector indicating periodic boundaries at maximum.
 * @return if true, at least one of the domain boundaries is periodic.
 * @note initially designed for finding periodic boundary in order to setup communication layers for them.
 */
bool GridInfoLbmInteface::CheckIfPeriodicDomainRequired(const DefInt dims,
    std::vector<bool>* const ptr_periodic_min, std::vector<bool>* const ptr_periodic_max) const {
    ptr_periodic_min->assign(dims, false);
    ptr_periodic_max->assign(dims, false);
    bool bool_has_periodic = false;
    if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryXMin)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryXMin)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_min->at(kXIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryYMin)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryYMin)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_min->at(kYIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryXMax)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryXMax)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_max->at(kXIndex) = true;
        bool_has_periodic = true;
    }
    if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryYMax)
        != domain_boundary_condition_.end()
        && domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryYMax)->boundary_scheme_
        == ELbmBoundaryConditionScheme::kPeriodic) {
        ptr_periodic_max->at(kYIndex) = true;
        bool_has_periodic = true;
    }
    if (dims == 3) {
        if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryZMin)
            != domain_boundary_condition_.end()
            && domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryZMin)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic) {
            ptr_periodic_min->at(kZIndex) = true;
            bool_has_periodic = true;
        }
        if (domain_boundary_condition_.find(amrproject::EDomainBoundaryDirection::kBoundaryZMax)
            != domain_boundary_condition_.end()
            && domain_boundary_condition_.at(amrproject::EDomainBoundaryDirection::kBoundaryZMax)->boundary_scheme_
            == ELbmBoundaryConditionScheme::kPeriodic) {
            ptr_periodic_max->at(kZIndex) = true;
            bool_has_periodic = true;
        }
    }
    return bool_has_periodic;
}
/**
 * @brief function to calculate the size of the grid node information needed for MPI communication.
 * @return size of the grid node information for MPI communication.
 */
int GridInfoLbmInteface::GetSizeOfGridNodeInfoForMpiCommunication() const {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    int size_info =  static_cast<int>(ptr_lbm_solver->k0NumQ_*sizeof(DefReal));
    size_info += static_cast<int>((1 + ptr_lbm_solver->GetSolverDim()) * sizeof(DefReal));
    if (ptr_lbm_solver->bool_forces_) {
        size_info += static_cast<int>(ptr_lbm_solver->GetNumForces()*sizeof(DefReal));
    }
    return size_info;
}
/**
 * @brief function to compute information of node in inner mpi layer before communication.
 * @param[in] sfbitset_in space filling code of the input node.
 * @param[in] lbm_solver class to manage LBM solver.
 */
void GridInfoLbmInteface::ComputeNodeInfoBeforeMpiCommunication(
    const DefSFBitset sfbitset_in, const SolverLbmInterface& lbm_solver) {
    lbm_solver.StreamInForAGivenNode(sfbitset_in, *ptr_sfbitset_aux_, ptr_lbm_grid_nodes_);
}
/**
 * @brief function to compute information of node in mpi communication layers.
 * @param map_inner_nodes container storing space filling codes of inner mpi communication layers will be sent to other ranks.
 * @param map_outer_nodes container storing space filling codes of outer mpi communication layer of the current rank.
 */
void GridInfoLbmInteface::ComputeInfoInMpiLayers(const std::map<int, DefMap<DefInt>>& map_inner_nodes,
    const DefMap<DefInt>& map_outer_nodes) {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    LbmCollisionOptInterface& collision_operator = ptr_lbm_solver->GetCollisionOperator(i_level_);
    DefReal dt_lbm = collision_operator.GetDtLbm();

    // collision for nodes in outer and inner MPI communication layers
    DefInt flag_not_collide = NodeFlagNotCollision_&(~(amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_));
    ptr_lbm_solver->CollisionForGivenNodes<DefInt>(i_level_, flag_not_collide, map_outer_nodes, this);
    DefMap<DefInt> map_one_layer_near_inner;
    std::vector<DefSFBitset> vec_neighbor;

    for (const auto& iter_layer : map_inner_nodes) {
        for (const auto& iter_node : iter_layer.second) {
            if (!(ptr_lbm_grid_nodes_->at(iter_node.first)->flag_status_&flag_not_collide)) {
                ptr_lbm_solver->func_collision_node_(
                    dt_lbm, &collision_operator, ptr_lbm_grid_nodes_->at(iter_node.first).get());
            }
            ptr_sfbitset_aux_->SFBitsetFindAllNeighborsVir(iter_node.first, &vec_neighbor);
            for (const auto& iter_neighbor : vec_neighbor) {
                if (map_one_layer_near_inner.find(iter_neighbor) == map_one_layer_near_inner.end()) {
                    map_one_layer_near_inner.insert({iter_neighbor, 0});
                }
            }
        }
    }

    // Since stream step will be performed before mpi communication, post-collision distribution
    // functions of neighboring nodes are need
    ptr_lbm_solver->CollisionForGivenNodes<DefInt>(i_level_, flag_not_collide, map_one_layer_near_inner, this);

    DefInt flag_not_stream = NodeFlagNotStream_&(~(amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_));
    for (const auto& iter_layer : map_inner_nodes) {
        ptr_lbm_solver->StreamInForGivenNodes<DefInt>(
            flag_not_stream, *ptr_sfbitset_aux_, iter_layer.second, this);
    }
}
/**
 * @brief function to compute information of node in mpi communication layers for interpolation.
 * @param map_interp_nodes container storing space filling codes of mpi communication layers for interpolation.
 * @note need to be called before ComputeInfoInMpiLayers
 */
void GridInfoLbmInteface::ComputeInfoInInterpMpiLayers(const std::map<int, DefMap<DefInt>>& map_interp_nodes) {
    SolverLbmInterface* ptr_lbm_solver = nullptr;
    if (auto ptr_tmp = ptr_solver_.lock()) {
        ptr_lbm_solver = dynamic_cast<SolverLbmInterface*>(ptr_tmp.get());
    } else {
        amrproject::LogManager::LogError("LBM solver is not created.");
    }
    LbmCollisionOptInterface& collision_operator = ptr_lbm_solver->GetCollisionOperator(i_level_);
    DefReal dt_lbm = collision_operator.GetDtLbm();

    // collision for nodes in outer and inner MPI communication layers
    DefInt flag_not_collide = NodeFlagNotCollision_&(~(amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_
        |amrproject::NodeBitStatus::kNodeStatusMpiInterpInner_));
    DefMap<DefInt> map_one_layer_near_inner;
    std::vector<DefSFBitset> vec_neighbor;

    for (const auto& iter_layer : map_interp_nodes) {
        for (const auto& iter_node : iter_layer.second) {
            if (ptr_lbm_grid_nodes_->find(iter_node.first) != ptr_lbm_grid_nodes_->end()) {
                if (!(ptr_lbm_grid_nodes_->at(iter_node.first)->flag_status_&flag_not_collide)) {
                    ptr_lbm_solver->func_collision_node_(
                        dt_lbm, &collision_operator, ptr_lbm_grid_nodes_->at(iter_node.first).get());
                }
                ptr_sfbitset_aux_->SFBitsetFindAllNeighborsVir(iter_node.first, &vec_neighbor);
                for (const auto& iter_neighbor : vec_neighbor) {
                    if (map_one_layer_near_inner.find(iter_neighbor) == map_one_layer_near_inner.end()) {
                        map_one_layer_near_inner.insert({iter_neighbor, 0});
                    }
                }
            }
        }
    }
    // Since stream step will be performed before mpi communication, post-collision distribution
    // functions of neighboring nodes are need
    ptr_lbm_solver->CollisionForGivenNodes<DefInt>(i_level_, flag_not_collide, map_one_layer_near_inner, this);

    DefInt flag_not_stream = NodeFlagNotStream_&(~(amrproject::NodeBitStatus::kNodeStatusMpiPartitionOuter_
        |amrproject::NodeBitStatus::kNodeStatusMpiPartitionInner_
        |amrproject::NodeBitStatus::kNodeStatusMpiInterpInner_));
    for (const auto& iter_layer : map_interp_nodes) {
        ptr_lbm_solver->StreamInForGivenNodes<DefInt>(
            flag_not_stream, *ptr_sfbitset_aux_, iter_layer.second, this);
    }
}
/**
 * @brief function to compute information of node in inner mpi layer after communication.
 * @param[in] sfbitset_in space filling code of the input node.
 * @param[in] lbm_solver class to manage LBM solver.
 */
void GridInfoLbmInteface::ComputeNodeInfoAfterMpiCommunication(
    const DefSFBitset sfbitset_in, const SolverLbmInterface& lbm_solver) {
}
}  // end namespace lbmproject
}  // end namespace rootproject
