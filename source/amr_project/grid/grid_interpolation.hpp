//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file grid_interpolation.hpp
* @author Zhengliang Liu
* @brief functions used for node interpolation.
* @date  2023-12-9
*/
#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_INTERPOLATION_HPP_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_INTERPOLATION_HPP_
#include <memory>
#include <vector>
#include <algorithm>
#include <string>
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function for 2D linear interpolation.
 * @param[in] region_length length of the sfbitset_coarse_region for calculating indices.
 * @param[in] flag_not_for_interp_coarse flag indicates node cannot be used for interpolation at coarse level.
 * @param[in] sfbitset_in input space filling code at fine level.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] sfbitset_coarse_region space filling codes in a region at coarse level.
 * @param[in] nodes_fine information of fine nodes.
 * @param[in] nodes_coarse information of coarse nodes.
 * @param[out] ptr_node pointer to the node to be interpolated.
 * @return 0 success; -1 failure since at least one is not found in fine or coarse, 
 * -2 node is not at the center of the region, -3 size of input region does not match the region length.
 * @note need to overload operator += Node_T, Node_T * DefReal
 */
template <typename Node_T>
requires std::is_base_of_v<GridNode, Node_T>
int GridInfoInterface::InterpolationLinear2D(const DefAmrIndexLUint region_length,
    const DefAmrUint flag_not_for_interp_coarse,
    const DefSFBitset& sfbitset_in, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefSFBitset>& sfbitset_coarse_region,
    const DefMap<std::unique_ptr<Node_T>>& nodes_fine,
    const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<Node_T>>& nodes_coarse,
    Node_T* const ptr_node) {
#ifdef DEBUG_CHECK_GRID
    // check if the input space filling code of one level lower is the same as that of the region center
    if (sfbitset_aux.SFBitsetToNLowerLevelVir(1, sfbitset_in)
        != sfbitset_coarse_region.at((region_length *2 + 1) * (region_length - 1))) {
        LogManager::LogWarning("space filling code of one level lower is not the same as that of the region center in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -2;
    } else if (sfbitset_coarse_region.size() != 4 * region_length * region_length) {
        LogManager::LogWarning("size of input space filling code does not match the region length in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -3;
    } else {
#endif
        SetNodeVariablesAsZeros(ptr_node);
        // check if x bit exist at current level
        std::array<DefReal, 2> coeffi_x, coeffi_y;
        DefSFBitset sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeXRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            coeffi_x[0] = 0.5;
            coeffi_x[1] = 0.5;
        } else {
            coeffi_x[0] = 1.;
            coeffi_x[1] = 0.;
        }
        sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeYRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            coeffi_y[0] = 0.5;
            coeffi_y[1] = 0.5;
        } else {
            coeffi_y[0] = 1.;
            coeffi_y[1] = 0.;
        }
        DefAmrIndexLUint index_x, index_y, num_elements = 2 * region_length;
        DefSFBitset sfbitset_fine;
        Node_T node_coarse2fine;
        for (auto iy = 0; iy <= 1; ++iy) {
            index_y = (iy + region_length - 1)* (2 * region_length) + region_length - 1;
            for (auto ix = 0; ix <= 1; ++ix) {
                index_x = index_y + ix;
                if (std::fabs(coeffi_x[ix] * coeffi_y[iy]) > kEps) {
                    if (nodes_coarse.find(sfbitset_coarse_region.at(index_x)) != nodes_coarse.end()&&(
                        (nodes_coarse.at(sfbitset_coarse_region.at(index_x))->flag_status_&flag_not_for_interp_coarse)
                        != flag_not_for_interp_coarse)) {
                        coarse_grid_info.NodeInfoCoarse2fine(
                            *nodes_coarse.at(sfbitset_coarse_region.at(index_x)), &node_coarse2fine);
                        *ptr_node += node_coarse2fine * (coeffi_x[ix] * coeffi_y[iy]);
                    } else {
                        sfbitset_fine =  sfbitset_aux.SFBitsetToNHigherLevelVir(1, sfbitset_coarse_region.at(index_x));
                        if (nodes_fine.find(sfbitset_fine) != nodes_fine.end()) {
                            *ptr_node += (*nodes_fine.at(sfbitset_fine).get()) * (coeffi_x[ix] * coeffi_y[iy]);
                        } else {
                            std::vector<DefReal> coordinates(2), grid_space =
                                {grid_space_.at(kXIndex), grid_space_.at(kYIndex)};
                            sfbitset_aux.SFBitsetComputeCoordinateVir(sfbitset_fine, grid_space, &coordinates);
                            LogManager::LogWarning("node ("+ std::to_string(coordinates[kXIndex]) + ", "
                                + std::to_string(coordinates[kYIndex])+ ") is not found in fine or coarse in "
                                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                            return -1;
                        }
                    }
                }
            }
        }
#ifdef DEBUG_CHECK_GRID
    }
#endif
    return 0;
}
/**
 * @brief function for 2D Lagrangian interpolation.
 * @param[in] interpolation_length length of region used for interpolation.
 * @param[in] region_length length of the sfbitset_coarse_region for calculating indices.
 * @param[in] flag_not_for_interp_coarse flag indicates node cannot be used for interpolation at coarse level.
 * @param[in] sfbitset_in input space filling code at fine level.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] sfbitset_coarse_region space filling codes in a region at coarse level.
 * @param[in] nodes_fine information of fine nodes.
 * @param[in] nodes_coarse information of coarse nodes.
 * @param[out] ptr_node pointer to the node to be interpolated.
 * @return 0 success; -1 failure since at least one is not found in fine or coarse, 
 * -2 node is not at the center of the region, -3 size of input region does not match the region length.
 * @note need to overload operator += Node_T, Node_T * DefReal
 */
template <typename Node_T>
requires std::is_base_of_v<GridNode, Node_T>
int GridInfoInterface::InterpolationLagrangian2D(const DefAmrIndexLUint interpolation_length,
    const DefAmrIndexLUint region_length, const DefAmrUint flag_not_for_interp_coarse, const DefSFBitset& sfbitset_in,
    const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefSFBitset>& sfbitset_coarse_region,
    const DefMap<std::unique_ptr<Node_T>>& nodes_fine,
    const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<Node_T>>& nodes_coarse,
    Node_T* const ptr_node) {
#ifdef DEBUG_CHECK_GRID
    // check if the input space filling code of one level lower is the same as that of the region center
    if (sfbitset_aux.SFBitsetToNLowerLevelVir(1, sfbitset_in)
        != sfbitset_coarse_region.at((region_length *2 + 1) * (region_length - 1))) {
        LogManager::LogWarning("space filling code of one level lower is not the same as that of the region center in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -2;
    } else if (sfbitset_coarse_region.size() != 4 * region_length * region_length) {
        LogManager::LogWarning("size of input space filling code does not match the region length in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -3;
    } else {
#endif
        SetNodeVariablesAsZeros(ptr_node);
        // check if x bit exist at current level
        DefAmrIndexLUint num_coeff = 2 * interpolation_length;
        std::vector<DefReal> coeffi_x(num_coeff), coeffi_y(num_coeff);
        DefSFBitset sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeXRef_.at(sfbitset_aux.kRefCurrent_);
        const LagrangianCoeff& coeff = CalculateLagrangianInterpCoeff(interpolation_length);

        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            std::copy(coeff.coeff1.begin(), coeff.coeff1.end(), coeffi_x.begin());
        } else {
            std::copy(coeff.coeff0.begin(), coeff.coeff0.end(), coeffi_x.begin());
        }
        sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeYRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            std::copy(coeff.coeff1.begin(), coeff.coeff1.end(), coeffi_y.begin());
        } else {
            std::copy(coeff.coeff0.begin(), coeff.coeff0.end(), coeffi_y.begin());
        }
        DefAmrIndexLUint index_x, index_y, num_elements = 2 * region_length;
        DefSFBitset sfbitset_fine;
        Node_T node_coarse2fine;
        for (DefAmrIndexLUint iy = 0; iy <= num_coeff; ++iy) {
            index_y = (iy + region_length - interpolation_length)* (2 * region_length)
                + region_length - interpolation_length;
            for (DefAmrIndexLUint ix = 0; ix <= num_coeff; ++ix) {
                index_x = index_y + ix;
                if (std::fabs(coeffi_x[ix] * coeffi_y[iy]) > kEps) {
                    if (nodes_coarse.find(sfbitset_coarse_region.at(index_x)) != nodes_coarse.end()&&(
                        (nodes_coarse.at(sfbitset_coarse_region.at(index_x))->flag_status_&flag_not_for_interp_coarse)
                        != flag_not_for_interp_coarse)) {
                        coarse_grid_info.NodeInfoCoarse2fine(
                            *nodes_coarse.at(sfbitset_coarse_region.at(index_x)), &node_coarse2fine);
                        *ptr_node += node_coarse2fine * (coeffi_x[ix] * coeffi_y[iy]);
                    } else {
                        sfbitset_fine =  sfbitset_aux.SFBitsetToNHigherLevelVir(1, sfbitset_coarse_region.at(index_x));
                        if (nodes_fine.find(sfbitset_fine) != nodes_fine.end()) {
                            *ptr_node += (*nodes_fine.at(sfbitset_fine).get()) * (coeffi_x[ix] * coeffi_y[iy]);
                        } else {
                            std::vector<DefReal> coordinates(2), grid_space =
                                {grid_space_.at(kXIndex), grid_space_.at(kYIndex)};
                            sfbitset_aux.SFBitsetComputeCoordinateVir(sfbitset_fine, grid_space, &coordinates);
                            LogManager::LogWarning("node ("+ std::to_string(coordinates[kXIndex]) + ", "
                                + std::to_string(coordinates[kYIndex])+ ") is not found in fine or coarse in "
                                + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                            return -1;
                        }
                    }
                }
            }
        }
#ifdef DEBUG_CHECK_GRID
    }
#endif
    return 0;
}
/**
 * @brief function for 3D linear interpolation.
 * @param[in] region_length length of the sfbitset_coarse_region for calculating indices.
 * @param[in] flag_not_for_interp_coarse flag indicates node cannot be used for interpolation at coarse level.
 * @param[in] sfbitset_in input space filling code at fine level.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] sfbitset_coarse_region space filling codes in a region at coarse level.
 * @param[in] nodes_fine information of fine nodes.
 * @param[in] nodes_coarse information of coarse nodes.
 * @param[out] ptr_node pointer to the node to be interpolated.
 * @return 0 success; -1 failure since at least one is not found in fine or coarse, 
 * -2 node is not at the center of the region, -3 size of input region does not match the region length.
 * @note need to overload operator += Node_T, Node_T * DefReal
 */
template <typename Node_T>
requires std::is_base_of_v<GridNode, Node_T>
int GridInfoInterface::InterpolationLinear3D(const DefAmrIndexLUint region_length,
    const DefAmrUint flag_not_for_interp_coarse,
    const DefSFBitset& sfbitset_in, const SFBitsetAuxInterface& sfbitset_aux,
    const std::vector<DefSFBitset>& sfbitset_coarse_region,
    const DefMap<std::unique_ptr<Node_T>>& nodes_fine,
    const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<Node_T>>& nodes_coarse,
    Node_T* const ptr_node) {
#ifdef DEBUG_CHECK_GRID
    // check if the input space filling code of one level lower is the same as that of the region center
    if (sfbitset_aux.SFBitsetToNLowerLevelVir(1, sfbitset_in)
        != sfbitset_coarse_region.at((region_length*region_length*4 +region_length*2 + 1)*(region_length - 1))) {
        LogManager::LogWarning("space filling code of one level lower is not the same as that of the region center in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -2;
    } else if (sfbitset_coarse_region.size() != 8 * region_length * region_length * region_length) {
        LogManager::LogWarning("size of input space filling code does not match the region length in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -3;
    } else {
#endif
        SetNodeVariablesAsZeros(ptr_node);
        // check if x bit exist at current level
        std::array<DefReal, 2> coeffi_x, coeffi_y, coeffi_z;
        DefSFBitset sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeXRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            coeffi_x[0] = 0.5;
            coeffi_x[1] = 0.5;
        } else {
            coeffi_x[0] = 1.;
            coeffi_x[1] = 0.;
        }
        sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeYRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            coeffi_y[0] = 0.5;
            coeffi_y[1] = 0.5;
        } else {
            coeffi_y[0] = 1.;
            coeffi_y[1] = 0.;
        }
        sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeZRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            coeffi_z[0] = 0.5;
            coeffi_z[1] = 0.5;
        } else {
            coeffi_z[0] = 1.;
            coeffi_z[1] = 0.;
        }
        DefAmrIndexLUint index_x, index_y, index_z, num_elements = 2 * region_length;
        DefSFBitset sfbitset_fine;
        Node_T node_coarse2fine;
        for (auto iz = 0; iz <= 1; ++iz) {
            index_z =(iz + region_length - 1)*(2 * region_length) + region_length - 1;
            for (auto iy = 0; iy <= 1; ++iy) {
                index_y = (index_z + iy)*(2 * region_length) + region_length - 1;
                for (auto ix = 0; ix <= 1; ++ix) {
                    index_x = index_y + ix;
                    if (std::fabs(coeffi_x[ix] * coeffi_y[iy]* coeffi_z[iz]) > kEps) {
                        if (nodes_coarse.find(sfbitset_coarse_region.at(index_x)) != nodes_coarse.end()
                            &&((nodes_coarse.at(sfbitset_coarse_region.at(index_x))->flag_status_
                            &flag_not_for_interp_coarse)!= flag_not_for_interp_coarse)) {
                            coarse_grid_info.NodeInfoCoarse2fine(
                                *nodes_coarse.at(sfbitset_coarse_region.at(index_x)), &node_coarse2fine);
                            *ptr_node += node_coarse2fine * (coeffi_x[ix] * coeffi_y[iy]* coeffi_z[iz]);
                        } else {
                            sfbitset_fine = sfbitset_aux.SFBitsetToNHigherLevelVir(
                                1, sfbitset_coarse_region.at(index_x));
                            if (nodes_fine.find(sfbitset_fine) != nodes_fine.end()) {
                                *ptr_node += (*nodes_fine.at(sfbitset_fine).get())
                                    *(coeffi_x[ix]*coeffi_y[iy]*coeffi_z[iz]);
                            } else {
                                std::vector<DefReal> coordinates(3), grid_space =
                                    {grid_space_.at(kXIndex), grid_space_.at(kYIndex), grid_space_.at(kZIndex)};
                                sfbitset_aux.SFBitsetComputeCoordinateVir(sfbitset_fine, grid_space, &coordinates);
                                LogManager::LogWarning("node ("+ std::to_string(coordinates[kXIndex]) + ", "
                                    + std::to_string(coordinates[kYIndex])+ ", "
                                    + std::to_string(coordinates[kZIndex])+ ") is not found in fine or coarse in "
                                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                                return -1;
                            }
                        }
                    }
                }
            }
        }
#ifdef DEBUG_CHECK_GRID
    }
#endif
    return 0;
}
/**
 * @brief function for 3D Lagrangian interpolation.
 * @param[in] interpolation_length length of region used for interpolation.
 * @param[in] region_length length of the sfbitset_coarse_region for calculating indices.
 * @param[in] flag_not_for_interp_coarse flag indicates node cannot be used for interpolation at coarse level.
 * @param[in] sfbitset_in input space filling code at fine level.
 * @param[in] sfbitset_aux class to manage functions for space filling code computation.
 * @param[in] sfbitset_coarse_region space filling codes in a region at coarse level.
 * @param[in] nodes_fine information of fine nodes.
 * @param[in] nodes_coarse information of coarse nodes.
 * @param[out] ptr_node pointer to the node to be interpolated.
 * @return 0 success; -1 failure since at least one is not found in fine or coarse, 
 * -2 node is not at the center of the region, -3 size of input region does not match the region length.
 * @note need to overload operator += Node_T, Node_T * DefReal
 */
template <typename Node_T>
requires std::is_base_of_v<GridNode, Node_T>
int GridInfoInterface::InterpolationLagrangian3D(const DefAmrIndexLUint interpolation_length,
    const DefAmrIndexLUint region_length,
    const DefAmrUint flag_not_for_interp_coarse, const DefSFBitset& sfbitset_in,
    const SFBitsetAuxInterface& sfbitset_aux, const std::vector<DefSFBitset>& sfbitset_coarse_region,
    const DefMap<std::unique_ptr<Node_T>>& nodes_fine,
    const GridInfoInterface& coarse_grid_info, const DefMap<std::unique_ptr<Node_T>>& nodes_coarse,
    Node_T* const ptr_node) {
#ifdef DEBUG_CHECK_GRID
    // check if the input space filling code of one level lower is the same as that of the region center
    if (sfbitset_aux.SFBitsetToNLowerLevelVir(1, sfbitset_in)
        != sfbitset_coarse_region.at((region_length*region_length*4 +region_length*2 + 1)*(region_length - 1))) {
        LogManager::LogWarning("space filling code of one level lower is not the same as that of the region center in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -2;
    } else if (sfbitset_coarse_region.size() != 8 * region_length * region_length * region_length) {
        LogManager::LogWarning("size of input space filling code does not match the region length in "
            + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
        return -3;
    } else {
#endif
        SetNodeVariablesAsZeros(ptr_node);
        // check if x bit exist at current level
        DefAmrIndexLUint num_coeff = 2 * interpolation_length;
        std::vector<DefReal> coeffi_x(num_coeff), coeffi_y(num_coeff), coeffi_z(num_coeff);
        DefSFBitset sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeXRef_.at(sfbitset_aux.kRefCurrent_);
        const LagrangianCoeff& coeff = CalculateLagrangianInterpCoeff(interpolation_length);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            std::copy(coeff.coeff1.begin(), coeff.coeff1.end(), coeffi_x.begin());
        } else {
            std::copy(coeff.coeff0.begin(), coeff.coeff0.end(), coeffi_x.begin());
        }
        sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeYRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            std::copy(coeff.coeff1.begin(), coeff.coeff1.end(), coeffi_y.begin());
        } else {
            std::copy(coeff.coeff0.begin(), coeff.coeff0.end(), coeffi_y.begin());
        }
        sfbitset_bit = sfbitset_in&sfbitset_aux.k0SFBitsetTakeZRef_.at(sfbitset_aux.kRefCurrent_);
        if ((sfbitset_bit&sfbitset_aux.k0SfBitsetCurrentLevelBits_) != 0) {
            std::copy(coeff.coeff1.begin(), coeff.coeff1.end(), coeffi_z.begin());
        } else {
            std::copy(coeff.coeff0.begin(), coeff.coeff0.end(), coeffi_z.begin());
        }
        DefAmrIndexLUint index_x, index_y, index_z, num_elements = 2 * region_length;
        DefSFBitset sfbitset_fine;
        Node_T node_coarse2fine;
        for (DefAmrIndexLUint iz = 0; iz <= num_coeff; ++iz) {
            index_z =(iz + region_length - interpolation_length)* (2 * region_length)
                + region_length - interpolation_length;
            for (DefAmrIndexLUint iy = 0; iy <= num_coeff; ++iy) {
                index_y = (index_z + iy)*(2 * region_length) + region_length - interpolation_length;
                for (DefAmrIndexLUint ix = 0; ix <= num_coeff; ++ix) {
                    index_x = index_y + ix;
                    if (std::fabs(coeffi_x[ix] * coeffi_y[iy]* coeffi_z[iz]) > kEps) {
                        if (nodes_coarse.find(sfbitset_coarse_region.at(index_x)) != nodes_coarse.end()
                            &&((nodes_coarse.at(sfbitset_coarse_region.at(index_x))->flag_status_
                            &flag_not_for_interp_coarse)!= flag_not_for_interp_coarse)) {
                            coarse_grid_info.NodeInfoCoarse2fine(
                                *nodes_coarse.at(sfbitset_coarse_region.at(index_x)), &node_coarse2fine);
                            *ptr_node += node_coarse2fine * (coeffi_x[ix] * coeffi_y[iy]* coeffi_z[iz]);
                        } else {
                            sfbitset_fine =  sfbitset_aux.SFBitsetToNHigherLevelVir(
                                1, sfbitset_coarse_region.at(index_x));
                            if (nodes_fine.find(sfbitset_fine) != nodes_fine.end()) {
                                *ptr_node += (*nodes_fine.at(sfbitset_fine).get())
                                    *(coeffi_x[ix]*coeffi_y[iy]*coeffi_z[iz]);
                            } else {
                                std::vector<DefReal> coordinates(3), grid_space =
                                    {grid_space_.at(kXIndex), grid_space_.at(kYIndex), grid_space_.at(kZIndex)};
                                sfbitset_aux.SFBitsetComputeCoordinateVir(sfbitset_fine, grid_space, &coordinates);
                                LogManager::LogWarning("node ("+ std::to_string(coordinates[kXIndex]) + ", "
                                    + std::to_string(coordinates[kYIndex])+ ", "
                                    + std::to_string(coordinates[kZIndex])+ ") is not found in fine or coarse in "
                                    + std::string(__FILE__) + " at line " + std::to_string(__LINE__));
                                return -1;
                            }
                        }
                    }
                }
            }
        }
#ifdef DEBUG_CHECK_GRID
    }
#endif
    return 0;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_INTERPOLATION_HPP_
