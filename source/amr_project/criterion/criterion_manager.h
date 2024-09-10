//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file criterion_mamnager.h
* @author Zhengliang Liu
* @date  2022-5-18
* @brief  define the class to manager criteria.
*/

#ifndef SOURCE_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
#define SOURCE_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
#include <vector>
#include <memory>
#include "criterion/geometry_info_connection.h"
#include "criterion/geometry_info_origin.h"
namespace rootproject {
namespace amrproject {
/**
* @class CriterionManager
* @brief class used to manage refinement criterion.
* @date  2022-5-23
*/
class CriterionManager {
 public:
    // geometry related
    DefInt numb_of_geometry_ = 0;
    std::vector<std::shared_ptr<GeometryInfoInterface>> vec_ptr_geometries_;

    void InitialAllGeometrySerial(const DefInt dims,
        const DefReal reference_dx, std::vector<DefReal> vec_real_offset);
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
