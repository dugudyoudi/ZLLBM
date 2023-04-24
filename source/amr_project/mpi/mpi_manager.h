//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file mpi.h
* @author Zhengliang Liu
* @date  2022-5-16
*/

#ifndef ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_
#define ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#include "grid/grid_manager.h"
#include "criterion/criterion_manager.h"
namespace rootproject {
namespace amrproject {
/**
* @class MpiManager
* @brief class used to manage the mpi processess.
*/
class MpiManager{
 public:
    int num_of_ranks_ = 1;  ///< total number of mpi ranks
    int rank_id_;  ///< current rank


    void StartupMpi(int argc, char* argv[]);
    void SetMpiParameters();

    void PartiteGridBySpaceFillingCurves(
        const std::vector<DefMap<DefUint>>& sfbitset_one_lower_level,
        GridManagerInterface const& grid_manager,
        std::vector<DefSFBitset>* const ptr_bitset_min,
        std::vector<DefSFBitset>* const ptr_bitset_max);

    DefSizet SerializeData(const DefMap<DefUint>& map_nodes,
        std::unique_ptr<char[]>& buffer) const;
};
} //  end namespace amrproject
}  //  end namespace rootproject
#endif  //  ENABLE_MPI
#endif  //  ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_

