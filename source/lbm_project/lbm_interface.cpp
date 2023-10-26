//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM models.
* @date  2023-9-30
*/
#include "lbm_interface.h"
#include "io/log_write.h"
namespace rootproject {
namespace lbmproject {
/**
* @brief  class to reinterpret type of grid nodes as LBM node type.
*/
void GridInfoLbmInteface::SetPointerToCurrentNodeType() {
    if (!map_grid_node_.empty()) {
        auto& first_element = map_grid_node_.begin()->second;
        if (dynamic_cast<GridNodeLbm*>(first_element.get())) {
            // The elements in map_nodes are of type GridNodeLbm,
            // assuming all nodes in map_grid_node_ are the same type.
            ptr_lbm_grid_ = reinterpret_cast<DefMap<std::unique_ptr<GridNodeLbm>>*>(&map_grid_node_);
        } else {
            std::string msg = "type of nodes stored in map_grid_node_ is not GridNodeLbm, "
            "please check if appropriate node creator is available in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__);
            amrproject::LogManager::LogError(msg);
        }
    }
}
/**
 * Get a pointer to the map store LBM nodes.
 */
DefMap<std::unique_ptr<GridNodeLbm>>* GridInfoLbmInteface::GetPointerToLbmGrid() {
    if (ptr_lbm_grid_ == nullptr) {
        SetPointerToCurrentNodeType();
    }
    return ptr_lbm_grid_;
}
void GridInfoLbmInteface::InitialGridNode(const DefSFBitset& bit_set_in) {
 
}
void SolverLbmInterface::SolverInitial() {
    InitialSetIndices();
}
/**
 * @brief function to calculate the equilibrium distribution functions for a 2D grid node.
 * @param[in] node  grid node containing LBM related information.
 * @param[out] ptr_feq pointer to the vector to store the calculated equilibrium distribution functions.
 */
void SolverLbmInterface::CalFeq2D(
    const GridNodeLbm& node, std::vector<DefReal>* const ptr_feq) const {
    ptr_feq->resize(this->k0NumQ_);
    DefReal c_uv = 0.;
    DefReal uv_sq = Square(node.velocity_.at(kXIndex)) + Square(node.velocity_.at(kYIndex));
    for (int iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = node.velocity_.at(kXIndex) * k0Cx_.at(iq) + node.velocity_.at(kYIndex) * k0Cy_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) * node.rho_ * (1. + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
/**
 * @brief function to calculate the equilibrium distribution functions for a 32D grid node.
 * @param[in] node  grid node containing LBM related information.
 * @param[out] ptr_feq pointer to the vector to store the calculated equilibrium distribution functions.
 */
void SolverLbmInterface::CalFeq3D(
    const GridNodeLbm& node, std::vector<DefReal>* const ptr_feq) const {
    ptr_feq->resize(this->k0NumQ_);
    DefReal c_uv = 0.;
    DefReal uv_sq = Square(node.velocity_.at(kXIndex)) + Square(node.velocity_.at(kYIndex))
        + Square(node.velocity_.at(kZIndex));
    for (int iq = 0; iq < this->k0NumQ_; ++iq) {
        c_uv = node.velocity_.at(kXIndex) * k0Cx_.at(iq) + node.velocity_.at(kYIndex) * k0Cy_.at(iq)
            + node.velocity_.at(kZIndex) * k0Cz_.at(iq);
        ptr_feq->at(iq) = k0Weights_.at(iq) * node.rho_ * (1. + 3. * c_uv + 4.5 * c_uv * c_uv - 1.5 * uv_sq);
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject