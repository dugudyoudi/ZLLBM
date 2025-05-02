//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file criterion_numerates.h
* @author Zhengliang Liu
* @date  2022-11-17
* @brief  define the class to manager numerates for criteria.
*/

#ifndef SOURCE_AMR_PROJECT_CRITERION_CRITERION_NUMERATES_H_
#define SOURCE_AMR_PROJECT_CRITERION_CRITERION_NUMERATES_H_
namespace rootproject {
namespace amrproject {
enum class ECriterionType : int {
    kUndefined = 0,
    kGeometry = 1
};
/**
* @struct EGridExtendType
* @brief enum class to store grid extending type
*/
enum class EGridExtendType : int {
    kSameInAllDirections = 1,
    kInAndOut = 2
};
///////////// geometry related //////////////
/**
* @class EGeometryCellType
* @brief enum class to store data type for geometry connection type
* @note numbering is the same as the cell types of vtk
*/
enum class EGeometryCellType : int {
    kUndefined = 0,
    kPolyLine = 4,
    kTriangle = 5
};
/**
* @class EGeometryStatus
* @brief enum class to identify geometry status whether it cause dynamic mesh adaptation
*/
enum class EGeometryStatus : int {
    kVirtual = 0,
    kStatic = 1,
    kMoving = 2
};
/**
* @class DefaultGeometryShapeType
* @brief enum class to store default geometry shape type
*/
enum class DefaultGeoShapeType : int {
    kUndefined = 0,
    kReadFromFile = 1,
    kCircle = 2,
    kCube = 3
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_CRITERION_NUMERATES_H_
