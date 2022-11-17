//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file criterion_mamnager.h
* @author Zhengliang Liu
* @date  2022-5-18
* @brief  define the class to manager criteria.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
#define ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
#include "criterion/geometry_info_connection.h"
namespace rootproject {
namespace amrproject {
namespace criterion {
/**
* @class CriterionManager
* @brief class used to manage refinement criterion.
* @date  2022-5-23
*/
class CriterionManager {
 public:
    // geometry related
    DefUint  k0GeoDims_ = 0;  ///< dimension
    DefSizet numb_of_geometry_ = 0;
    std::vector<std::shared_ptr<GeometryInfoInterface>> vec_ptr_geometries_;

    // factory design pattern to generate GeometryInfo instance
    void CreateGeometryInfo(const DefUint i_level,
         const DefUint geo_data_type);
    void InitialGeometrySerial(std::vector<DefReal> vec_real_offset);
};
}  // end namespace criterion
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
