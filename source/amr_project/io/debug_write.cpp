//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file log_write.cpp
* @author Zhengliang Liu
* @date  2022-8-14
*/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include <vector>
#include "../defs_libs.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif  // ENABLE_MPI
#include "io/debug_write.h"
namespace rootproject {
namespace amrproject {
// static member
void DebugWriterManager::WriteCoordinatesInPts(const DefInt dims, const std::string& datafile_name,
    const std::vector<DefReal>& grid_offset, const std::vector<DefReal>& grid_space,
    const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefInt>& map_points) {
    FILE* fp = fopen(("pointdata_rank_" + datafile_name + ".pts").c_str(), "w");
    if (!fp) {
        int rank_id = 0;
#ifdef ENABLE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
#endif  // ENABLE_MPI
        LogManager::LogError("File on node " + std::to_string(rank_id)
            + " was not opened for writing csv data in WriteCoordinatesInPts."
            + " in "+ std::string(__FILE__) + " at line " + std::to_string(__LINE__));
    } else {
        DefSizet num_points = map_points.size();

        fprintf(fp, "%zd\n", num_points);

        std::string str_format = "%.6f  %.6f %.6f\n";
        std::vector<DefReal> coordinates;
        for (const auto& iter : map_points) {
            sfbitset_aux.SFBitsetComputeCoordinateVir(iter.first, grid_space, &coordinates);
            if (dims == 2) {
                fprintf(fp, str_format.c_str(),
                    coordinates.at(kXIndex), coordinates.at(kYIndex), 0.);
            } else {
                fprintf(fp, str_format.c_str(),
                    coordinates.at(kXIndex), coordinates.at(kYIndex), coordinates.at(kZIndex));
            }
        }
        fclose(fp);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
