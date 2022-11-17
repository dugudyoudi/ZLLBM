//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file grid_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage grid related processes.
* @date  2022-6-7
* @note  functions from other managers will be called.
*/
#include <string>
#include "auxiliary_inline_func.h"
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
#include "io/log_write.h"
#include "mpi/mpi_manager.h"
namespace rootproject {
namespace amrproject {
namespace grid {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
/**
* @brief function to setup default grid related parameters.
*/
void GridManager2D::DefaultInitializationDims(const DefSizet max_level) {
    k0GridDims_ = 2;
    k0IntOffest_ = {1, 1};
    k0NumNeighbours_ = 9;
}
/**
* @brief function to setup and check grid related parameters.
*/
void GridManager2D::SetGridParameters() {
    int rank_id = 0;
#ifdef ENABLE_MPI
    rank_id = mpi::MpiManager::GetInstance()->rank_id_;
#endif  // ENABLE_MPI

    if (rank_id == 0) {
        // check if length of computational domain is given
        if (k0DomainSize_.at(kXIndex) < kEps) {
            io::LogError("Domain length in x direction (k0DomainSize_[0])"
                " should be a positive value");
        } else if (k0DomainSize_.at(kYIndex) < kEps) {
            io::LogError("Domain length in x direction (k0DomainSize_[1])"
                " should be a positive value");
        }
        // check if grid space is given
        if (k0DomainDx_.at(kXIndex) < kEps
            && k0DomainDx_.at(kYIndex) < kEps) {
            io::LogError("Grid space of x or y(k0DomianDx_)"
                " shoud be positive values");
        }
    }  // end if (rank_id == 0)

    // set grid space if not all grid spaces are given
    if (k0DomainDx_.at(kXIndex) < kEps) {
        k0DomainDx_.at(kXIndex) = k0DomainDx_.at(kYIndex);
    }
    if (k0DomainDx_.at(kYIndex) < kEps) {
        k0DomainDx_.at(kYIndex) = k0DomainDx_.at(kXIndex);
    }


    // caculate number of background nodes in each direction
    k0MaxIndexOfBackgroundNode_ = {
        static_cast<DefLUint>(k0DomainSize_[kXIndex]
        / k0DomainDx_[kXIndex]) + k0IntOffest_[kXIndex],
        static_cast<DefLUint>(k0DomainSize_[kYIndex]
        / k0DomainDx_[kYIndex]) + k0IntOffest_[kYIndex] };


    // print information of grid parameters
    if (rank_id == 0) {
        io::LogInfo("Dimension is: " + std::to_string(k0GridDims_));
        io::LogInfo("Maximum refinement level is: "
            + std::to_string(k0MaxLevel_));
        io::LogInfo("Domain size is: "
            + std::to_string(k0DomainSize_.at(kXIndex)) + " X "
            + std::to_string(k0DomainSize_.at(kYIndex)));
        io::LogInfo("Grid space dx is: "
            + std::to_string(k0DomainDx_.at(kXIndex)) + ", and dy is: "
            + std::to_string(k0DomainDx_.at(kYIndex)));

    }  // end if (rank_id == 0)
    // set offsets
    k0RealOffest_[kXIndex] = k0IntOffest_[kXIndex] * k0DomainDx_[kXIndex];
    k0RealOffest_[kYIndex] = k0IntOffest_[kYIndex] * k0DomainDx_[kYIndex];


    // check if domain size may exceed range of morton code
    /* the criterion is the maximum index for
    *  (domain size + offset distance)/(minimum grid spacing).
    *  Actually, this is a strict requirement considering high resolution
    *  might be needed near the computational domain. If high resolution
    *  region is far from the boundary, maximum index for
    *  (domain size + offset distance)/(maximum grid spacing)
    *  could be used for larger domain
    */
    // number of bits available for background mesh in one dimension
    DefSizet bit_max = kSFBitsetBit / k0GridDims_ - k0MaxLevel_;
    DefSizet index_max = TwoPowerN(bit_max);
    DefLUint scale_i_level = static_cast<DefLUint>(TwoPowerN(k0MaxLevel_));
    if (k0MaxIndexOfBackgroundNode_.at(kXIndex) > index_max) {
        io::LogError("Domain size exceeds the limist of sfbitset in"
            " x direciont, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h");
    }
    if (k0MaxIndexOfBackgroundNode_.at(kYIndex) > index_max) {
        io::LogError("Domain size exceeds the limist of sfbitset in"
            " y direciont, try to increase number of bits for "
            " storing space filling code (kSFBitsetBit) in defs_libs.h");
    }
}
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
}
}  // end namespace amrproject
}  // end namespace rootproject
