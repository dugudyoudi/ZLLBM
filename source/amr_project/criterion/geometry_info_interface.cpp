//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.cpp
* @author Zhengliang Liu
* @brief functions to manage geometry information
* @date  2022-8-25
*/
#include "io/log_write.h"
#include "io/input_parser.h"
#include "criterion/geometry_info_interface.h"
#include "criterion/geometry_info_origin.h"
#include "criterion/geometry_info_connection.h"
namespace rootproject {
namespace amrproject {
/**
 * @brief function to choose grid extension type based on input string.
 * @param[in] type_string string of extension type.
 */
void GeometryInfoInterface::ChooseGridExtendType(const std::string type_string) {
    if (type_string == "in_and_out") {
        grid_extend_type_ = EGridExtendType::kInAndOut;
    } else if (type_string == "same_in_all_directions") {
        grid_extend_type_ = EGridExtendType::kSameInAllDirections;
    } else {
        LogManager::LogError("grid extend type is not supported: " + type_string);
    }
}
/**
 * @brief function to read and set geometry parameters.
 * @param[in] level default geometry level.
 * @param[in, out] ptr_geo_parameters pointer to map storing geometry parameters.
 */
void GeometryInfoInterface::ReadAndSetGeoParameters(const DefInt level,
    std::map<std::string, ParserData>* const ptr_geo_parameters) {
    DefInt geo_level = level;
    InputParser input_parser;
    input_parser.GetValue<DefInt>("geo.level", ptr_geo_parameters, &geo_level);
    SetLevel(geo_level);
    std::string extend_type;
    if (input_parser.GetValue<std::string>("geo.extend_type", ptr_geo_parameters, &extend_type)) {
        ChooseGridExtendType(extend_type);
    }
    DefInt default_extend = 6;
    std::vector<DefInt> x_extend_pos(default_extend, geo_level), x_extend_neg(default_extend, geo_level),
        y_extend_pos(default_extend, geo_level), y_extend_neg(default_extend, geo_level),
        z_extend_pos(default_extend, geo_level), z_extend_neg(default_extend, geo_level);
    input_parser.GetValue<DefInt>("geo.x_extend_pos", ptr_geo_parameters, &x_extend_pos);
    SetXExtendPositive(x_extend_pos);
    input_parser.GetValue<DefInt>("geo.x_extend_neg", ptr_geo_parameters, &x_extend_neg);
    SetXExtendNegative(x_extend_neg);
    input_parser.GetValue<DefInt>("geo.y_extend_pos", ptr_geo_parameters, &y_extend_pos);
    SetYExtendPositive(y_extend_pos);
    input_parser.GetValue<DefInt>("geo.y_extend_neg", ptr_geo_parameters, &y_extend_neg);
    SetYExtendNegative(y_extend_neg);
    if (k0GeoDim_ == 3) {
        input_parser.GetValue<DefInt>("geo.z_extend_pos", ptr_geo_parameters, &z_extend_pos);
        SetZExtendPositive(z_extend_pos);
        input_parser.GetValue<DefInt>("geo.z_extend_neg", ptr_geo_parameters, &z_extend_neg);
        SetZExtendNegative(z_extend_neg);
    }
    std::vector<DefInt> inner_extend(default_extend, geo_level);
    if (input_parser.GetValue<DefInt>("geo.inner_extend", ptr_geo_parameters, &inner_extend)) {
        SetInnerExtend(inner_extend);
    }
}
/**
 * @brief function to initialize geometry.
 * @param dx reference spatial step.
 */
void GeometryInfoInterface::InitialGeometry(const DefReal dx) {
    if (ptr_geo_shape_ != nullptr) {
        ptr_geo_shape_->InitialShape(dx);
    }
}
/**
 * @brief function to update geometry.
 * @param sum_t time step at current level.
 * @return 0 the geometry is updated, 1 the geometry is not updated
 */
void GeometryInfoInterface::UpdateGeometry(const DefReal sum_t) {
    if (ptr_geo_shape_ != nullptr && geometry_status_ == EGeometryStatus::kMoving) {
        ptr_geo_shape_->UpdateShape(sum_t);
    }
}
/**
 * @brief function to instantiate class of geometry information based on geometry type.
 * @param[in] dims dimension of the geometry.
 * @param[in] geo_type type of geometry.
 * @return share pointer of geometry information
 */
std::shared_ptr<GeometryInfoInterface> GeoTypeReader::ReadGeoType(
    const DefInt dims, const std::string& geo_type) const {
    if (geo_type == "origin") {
        return std::make_shared<GeometryInfoOrigin>(dims);
    } else if (geo_type == "connection") {
        return std::make_shared<GeometryInfoConnection>(dims);
    } else {
        return nullptr;
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
