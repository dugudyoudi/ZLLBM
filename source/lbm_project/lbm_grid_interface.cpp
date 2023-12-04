//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_grid_interface.cpp
* @author Zhengliang Liu
* @brief functions used for manage LBM grid interface.
* @date  2023-9-30
*/
#include "lbm_interface.h"
#include "io/log_write.h"
#include "io/vtk_writer.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to create grid node.
 */
std::unique_ptr<amrproject::GridNode> GridInfoLbmInteface::GridNodeCreator() {
    const SolverLbmInterface& lbm_solver = *(std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_));
    if (!bool_forces_) {
        return std::make_unique<GridNodeLbm>(
            lbm_solver.k0Rho_, lbm_solver.k0Velocity_, f_ini_, f_collide_ini_);
    } else {
        return std::make_unique<GridNodeLbm>(
            lbm_solver.k0Rho_, lbm_solver.k0Velocity_, lbm_solver.k0Force_, f_ini_, f_collide_ini_);
    }
}
/**
 * @brief function to initializes the grid information.
 */
void GridInfoLbmInteface::InitialGridInfo() {
    SolverLbmInterface* ptr_lbm_solver = std::dynamic_pointer_cast<SolverLbmInterface>(ptr_solver_).get();
    ptr_lbm_solver->SetInitialDisFuncBasedOnReferenceMacros(&f_ini_, &f_collide_ini_);
    SetCollisionOperator();
    ptr_collision_operator_->viscosity_lbm_ = ptr_lbm_solver->k0LbmViscosity_;
    ptr_collision_operator_->dt_lbm_ = 1./ static_cast<DefReal>(TwoPowerN(i_level_));
    ptr_collision_operator_->CalRelaxationTime();
}
/**
* @brief  function to reinterpret type of grid nodes as LBM node type.
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
 * @brief function to get pointer to the map store LBM nodes.
 */
DefMap<std::unique_ptr<GridNodeLbm>>* GridInfoLbmInteface::GetPointerToLbmGrid() {
    if (ptr_lbm_grid_ == nullptr) {
        SetPointerToCurrentNodeType();
    }
    return ptr_lbm_grid_;
}
/**
 * @brief function to call assigned boundary conditions for each domain boundary.
 */
void GridInfoLbmInteface::ComputeDomainBoundaryCondition() {
    for (auto i =0; i < domain_boundary_min_.size(); ++i) {
        switch (i) {
        case kXIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXNeg)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXNeg)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryXNeg, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for x minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kYIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYNeg)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYNeg)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryYNeg, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for y minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kZIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZNeg)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZNeg)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryZNeg, domain_boundary_min_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for z minimum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        default:
            amrproject::LogManager::LogError("Type of boundary not defined for the " + std::to_string(i)
                + "th minimum domain boundary in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            break;
        }
    }
    for (auto i =0; i < domain_boundary_max_.size(); ++i) {
        switch (i) {
        case kXIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryXPos)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryXPos)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryXPos, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for x maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kYIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryYPos)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryYPos)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryYPos, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for y maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        case kZIndex: {
                if (domain_boundary_condition_.find(ELbmBoundaryType::kBoundaryZPos)
                    != domain_boundary_condition_.end()) {
                    domain_boundary_condition_.at(ELbmBoundaryType::kBoundaryZPos)->CalBoundaryCondition(
                        ELbmBoundaryType::kBoundaryZPos, domain_boundary_max_.at(i), this);
                } else {
                    amrproject::LogManager::LogError("Boundary condition for z maximum domain boundary not found "
                        + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                }
            }
            break;
        default:
            amrproject::LogManager::LogError("Type of boundary not defined for the " + std::to_string(i)
                + "th maximum domain boundary in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
            break;
        }
    }
}
/**
 * @brief function to setup output infomation on variables of LBM node.
 */
void GridInfoLbmInteface::SetupOutputVariables() {
    GetPointerToLbmGrid();  // assign pointer to ptr_lbm_grid_ if it is nullptr.
    OutputLBMNodeVariableInfo output_info_temp;
    for (const auto& iter_str : output_variable_name_) {
        if (iter_str == "rho") {
            output_info_temp.variable_dims_ = 1;
            output_info_temp.func_get_var_ = [] (const GridNodeLbm& node)->std::vector<DefReal> {
                return {node.rho_};
            };
            output_info_temp.output_name_ = "rho";
            output_variables_.push_back(output_info_temp);
        } else if (iter_str == "velocity") {
            output_info_temp.variable_dims_ = ptr_solver_->k0SolverDims_;
            output_info_temp.func_get_var_ = [] (const GridNodeLbm& node)->std::vector<DefReal> {
                return node.velocity_;
            };
            output_info_temp.output_name_ = "velocity";
            output_variables_.push_back(output_info_temp);
        }
}
}
/**
* @brief   function to write scalar and vector variables.
* @param[in]  fp   pointer to output file.
* @param[in]  base64_instance reference to class to convert data to uint8 type.
* @param[in]  bool_binary write data in binary or ascii format.
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  map_node_index indices of nodes for unstructured grid.
*/
void GridInfoLbmInteface::WriteOutputScalarAndVectors(
    FILE* const fp, const bool bool_binary,
    const amrproject::Base64Utility& base64_instance,
    const amrproject::OutputDataFormat& output_data_format,
    const DefMap<DefSizet>& map_node_index) const {
    for (const auto& iter_var : output_variables_) {
        OutputOneVariable(fp, bool_binary, base64_instance, iter_var, output_data_format, map_node_index);
    }
}
/**
* @brief function to write a scalar of each node.
* @param[in]  fp   pointer to output file.
* @param[in]  bool_binary write data in binary or ascii format.
* @param[in] output_info output infomation of a node variable.
* @param[in]  base64_instance reference to class to convert data to uint8 type.
* @param[in]  output_data_format output data (real or integer) format
* @param[in]  map_node_index indices of nodes for unstructured grid.
*/
int GridInfoLbmInteface::OutputOneVariable(
    FILE* const fp, const bool bool_binary,
    const amrproject::Base64Utility& base64_instance,
    const OutputLBMNodeVariableInfo& output_info,
    const amrproject::OutputDataFormat& output_data_format,
    const DefMap<DefSizet>& map_node_index) const {
    if (ptr_lbm_grid_ == nullptr) {
        amrproject::LogManager::LogWarning("LBM grid is not found for output in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return 1;
    }
    const DefMap<std::unique_ptr<GridNodeLbm>>& map_grid_node = *ptr_lbm_grid_;
    std::string str_format, str_temp;
    DefAmrIndexUint dims = output_info.variable_dims_;
    if (bool_binary) {
        str_format =  "binary";
    } else {
        str_format =  "ascii";
    }
    str_temp.assign("      <DataArray NumberOfComponents=\"" + std::to_string(dims)
        + "\" type=\"" + output_data_format.output_real_.format_name_
        + "\" Name=\"" + output_info.output_name_ + "\" format=\"" + str_format + "\">\n");
    fprintf_s(fp, str_temp.c_str());

    if (bool_binary) {
        std::vector<uint8_t> vec_uint8{}, vec_base64{};
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            const std::vector<DefReal>& vec_var =  output_info.func_get_var_(*map_grid_node.at(iter->first).get());
            for (DefAmrIndexUint i_dims = 0;  i_dims < dims; ++i_dims) {
                base64_instance.AddToVecChar(output_data_format.output_real_.CastType(vec_var.at(i_dims)), &vec_uint8);
            }
        }
        base64_instance.Encode(&vec_uint8, &vec_base64);
        for (const auto& iter : vec_base64) {
            fprintf_s(fp, "%c", iter);
        }
        fprintf_s(fp, "\n");
    } else {
        std::string str_format = "     "
            + output_data_format.output_real_.printf_format_;
        for (auto iter = map_node_index.begin();
            iter != map_node_index.end(); ++iter) {
            const std::vector<DefReal>& vec_var =  output_info.func_get_var_(*map_grid_node.at(iter->first).get());
            for (DefAmrIndexUint i_dims = 0; i_dims < dims; ++i_dims) {
                fprintf_s(fp, "  ");
                fprintf_s(fp, str_format.c_str(), vec_var.at(i_dims));
            }
            fprintf_s(fp, "\n");
        }
    }
    fprintf_s(fp, "      </DataArray>\n");
    return 0;
}
}  // end namespace lbmproject
}  // end namespace rootproject