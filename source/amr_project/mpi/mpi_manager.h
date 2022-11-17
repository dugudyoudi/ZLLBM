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
namespace rootproject {
namespace mpi {
/**
* @class MpiManager
* @brief class used to manage the mpi processess.
*/
class MpiManager{
 public:
    int numb_of_ranks_ = 1;  ///< total number of mpi ranks
    static int rank_id_;  ///< current rank
    DefSizet level_for_mpiblock_ = 0;  /* if 0, the i_level is
              the same as background mesh, otherwise 
              num_of_nodes_per_mpiblock_ is 2^level_for_mpiblock_*/
    DefSizet num_of_nodes_per_mpiblock_ = 0;
    /**< number of nodes in earch direction for partition blocks*/
            // partition blocks are based on Z curves, thus the number
            // of nodes is the same in all directions
    std::vector<std::array<DefSFBitset, 2>> vec_partition_morton{};
    ///< vector store sfbitsets of partition results for all mpi ranks

    void StartupMpi(int argc, char* argv[]);
    void SetMpiParameters();

 private:
    void PartiteBackground();
};
} //  end namespace mpi
}  //  end namespace rootproject
#endif  //  ENABLE_MPI
#endif  //  ROOTPROJECT_SOURCE_MPI_MPI_MANAGER_H_

