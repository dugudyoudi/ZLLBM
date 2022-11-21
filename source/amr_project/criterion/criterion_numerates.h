//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file criterion_numerates.h
* @author Zhengliang Liu
* @date  2022-11-17
* @brief  define the class to manager numerates for criteria.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_NUMERATES_H_
#define ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_NUMERATES_H_
namespace rootproject {
namespace amrproject {
enum class ECriteriolType {
    kUndefined = 0,
    kGeometry = 1
};
///////////// geometry related //////////////
/**
* @class EGeometryOutputkDataType
* @brief enum class to store data type for geometry output
* @note numbering is the same as the cell types of vtk
*/
enum class EGeometryCellType {
    kUndefined = 0,
    kPolyLine = 4,
    kTriangle = 5
};
/**
* @class EGeometryEnclosedType
* @brief enum class to store data type for shape enclose type
*/
enum class EGeometryEnclosedType {
    kOpen = 0,
    kEnclosedRemove = 1,
    kEnclosedExtend = 2
};
/**
* @class DefaultGeometryShapeType
* @brief enum class to store default geometry shape type
*/
enum class DefaultGeoShapeType {
    kUndefined = 0,
    kReadFromFile = 1,
    kCircle = 2,
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_NUMERATES_H_
