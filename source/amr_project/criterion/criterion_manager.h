//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file criterion_mamnager.h
* @author Zhengliang Liu
* @date  2022-5-18
* @brief  define the class to manager criteria.
*/

#ifndef ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
#define ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
#include <vector>
#include <memory>
#include "criterion/geometry_info_connection.h"
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
    DefAmrIndexUint numb_of_geometry_ = 0;
    std::vector<std::shared_ptr<GeometryInfoInterface>> vec_ptr_geometries_;

    std::unique_ptr<DefaultGeoManager> ptr_default_geo_manager_;

    void InitialAllGeometrySerial(const DefAmrIndexUint dims, std::vector<DefReal> vec_real_offset);

    CriterionManager() {
        ptr_default_geo_manager_ = std::make_unique<DefaultGeoManager>();
    }
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
