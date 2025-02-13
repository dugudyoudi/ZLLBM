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
#include <string>
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


    void ReadAndSetGeoParametersBasedOnShape(const DefInt dims,  const DefInt max_level, const InputParser& input_parser,
        const GeoTypeReader& geo_reader = GeoTypeReader(), const GeoShapeReader& shape_reader = GeoShapeReader());

    void InitialAllGeometrySerial(const DefReal reference_dx, const std::array<DefReal, 3>& real_offset);
    void PushbackAGeometryBasedOnShape(const DefInt dims, const std::string& geo_type,
        const std::string& shape_id, const std::string& shape_type,
        const GeoTypeReader& geo_reader = GeoTypeReader(), const GeoShapeReader& shape_reader = GeoShapeReader());
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_CRITERION_MANAGER_H_
