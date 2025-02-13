//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/grid_enumerates.h
* @author Zhengliang Liu
* @date  2022-5-16
* @brief define enumerate classes for grid information
*/
#ifndef SOURCE_AMR_PROJECT_GRID_GRID_ENUMERATES_H_
#define SOURCE_AMR_PROJECT_GRID_GRID_ENUMERATES_H_
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
enum class EInterpolationMethod : DefInt {
    kLinear = 1,
    kLagrangian = 2
};
/**
* @brief enumerate time stepping schemes
*/
enum class ETimeSteppingScheme {
    kMultiSteppingC2F = 1
};
/**
* @brief enumerate timing during a time step
*/
enum class ETimingInOneStep {
    kUndefined = 0,
    kStepBegin = 1,
    kStepEnd = 2
};
enum class EDomainBoundaryType{
    kCubic = 0,
};
/**
* @brief enumerate boundary directions
*/
enum class EDomainBoundaryDirection {
    kUndefined = 0,
    kBoundaryXMin = 1,
    kBoundaryXMax = 2,
    kBoundaryYMin = 3,
    kBoundaryYMax = 4,
    kBoundaryZMin = 5,
    kBoundaryZMax = 6
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_GRID_GRID_ENUMERATES_H_
