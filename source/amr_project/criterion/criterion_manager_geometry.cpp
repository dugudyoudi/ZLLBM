//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_manager.cpp
* @author Zhengliang Liu
* @brief functions used to manage changes of geometries
* @date  2022-5-23
*/
#include "grid/grid_manager.h"
#include "io/log_write.h"
#include "io/input_parser.h"
#include "criterion/criterion_manager.h"
#include "criterion/geometry_default_shape.h"
#include "criterion/geometry_info_interface.h"
namespace rootproject {
namespace amrproject {
/**
* @brief   function to initialize geometries
* @param[in]  reference_dx reference spatial step.
* @param[in]  real_offset offsets of the mesh.
*/
void CriterionManager::InitialAllGeometrySerial(
    const DefReal reference_dx, const std::array<DefReal, 3>& real_offset) {
    numb_of_geometry_ = static_cast<DefInt>(vec_ptr_geometries_.size());
    for (auto i_geo = 0; i_geo < numb_of_geometry_; ++i_geo) {
        vec_ptr_geometries_.at(i_geo)->SetGeoIndex(i_geo);
        // assign offset distance
        vec_ptr_geometries_.at(i_geo)->SetOffset(real_offset);
        // initial geometry vertex
        vec_ptr_geometries_.at(i_geo)->InitialGeometry(reference_dx);
    }
}
/**
* @brief   function to setup geometries from input parser
* @param[in]  dims dimension of geometries.
* @param[in] default_level default refinement level for geometries.
* @param[in]  input_parser reference of input parser.
* @param[in]  geo_reader reference of geometry reader.
* @param[in]  shape_reader reference of shape reader.
*/
void CriterionManager::ReadAndSetGeoParametersBasedOnShape(const DefInt dims, const DefInt default_level,
    const InputParser& input_parser, const GeoTypeReader& geo_reader, const GeoShapeReader& shape_reader) {
    std::string key_for_this_func =  "geo.shape";
    std::string geo_type = "origin", shape_type = "default", shape_id;
    const auto input_map = input_parser.GetNestedMapInput();
    if (input_map.find(key_for_this_func) != input_map.end()) {
        for (const auto iter_name : input_map.at(key_for_this_func)) {
            if (iter_name.second.find("geo.type") != iter_name.second.end()) {
                geo_type = iter_name.second.at("geo.type");
            }
            if (iter_name.second.find("shape.type") != iter_name.second.end()) {
                shape_type = iter_name.second.at("shape.type");
            }
            if (iter_name.second.find("shape.id") != iter_name.second.end()) {
                shape_id = iter_name.second.at("shape.id");
            } else {
                LogManager::LogError("shape.id is not found for geometry: " + iter_name.first);
            }
            PushbackAGeometryBasedOnShape(dims, geo_type, shape_id, shape_type, geo_reader, shape_reader);
            std::weak_ptr<GeometryInfoInterface> ptr_geo_tmp = vec_ptr_geometries_.back();
            if (auto ptr_geo = ptr_geo_tmp.lock()) {
                ptr_geo->SetName(iter_name.first);
                if (!ptr_geo->ptr_geo_shape_->ReadAndSetGeoShapeParameters(iter_name.second)) {
                    LogManager::LogError("shape parameters are not set for geometry: " + iter_name.first);
                }
                ptr_geo->ReadAndSetGeoParameters(default_level, iter_name.second);
            }
        }
    }
}
/**
* @brief   function to instantiate a geometry based on shape.
* @param[in]  dims dimension of geometries.
* @param[in]  geo_type type of geometry.
* @param[in]  shape_id predefined id for shapes.
* @param[in]  shape_type shape type.
* @param[in]  geo_reader reference of geometry reader.
* @param[in]  shape_reader reference of shape reader.
*/
void CriterionManager::PushbackAGeometryBasedOnShape(const DefInt dims, const std::string& geo_type,
    const std::string& shape_id, const std::string& shape_type,
    const GeoTypeReader& geo_reader, const GeoShapeReader& shape_reader) {
    vec_ptr_geometries_.push_back(geo_reader.ReadGeoType(dims, geo_type));
    if (vec_ptr_geometries_.back() == nullptr) {
        LogManager::LogError("Geometry type " + geo_type + " is not supported.");
    } else {
        vec_ptr_geometries_.back()->ptr_geo_shape_ = shape_reader.ReadAndCreateGeoShape(
            dims, shape_id, vec_ptr_geometries_.back(), shape_type);
    }
}
}  // end namespace amrproject
}  // end namespace rootproject
