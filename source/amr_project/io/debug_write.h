//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file log_write.h
* @author Zhengliang Liu
* @date  2022-8-14
*/
#ifndef SOURCE_AMR_PROJECT_IO_DEBUG_WRITE_H_
#define SOURCE_AMR_PROJECT_IO_DEBUG_WRITE_H_
#include <string>
#include <vector>
#include "grid/sfbitset_aux.h"
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
class DebugWriterManager {
 public:
    static void  WriteCoordinatesInPts(
        const DefInt dims, const std::string& datafile_name,
        const std::vector<DefReal>& grid_offset, const std::vector<DefReal>& grid_space,
        const SFBitsetAuxInterface& sfbitset_aux, const DefMap<DefInt>& map_points);
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_IO_DEBUG_WRITE_H_
