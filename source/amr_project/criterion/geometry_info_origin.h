//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_info_interface.h
* @author Zhengliang Liu
* @date  2023-4-25
* @brief  define classes to store geometry information with isolated points
*/
#ifndef SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_ORIGIN_H_
#define SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_ORIGIN_H_
#include <memory>
#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <utility>
#include "criterion/geometry_info_Interface.h"
#include "criterion/geometry_info_connection.h"
namespace rootproject {
namespace amrproject {
class GeometryInfoOrigin : public GeometryInfoInterface {
 public:
    void InitialGeometry(const DefReal dx) override;
    void FindTrackingNodeBasedOnGeo(
        const SFBitsetAuxInterface& sfbitset_aux, GridInfoInterface* const ptr_grid_info) override;

    DefSizet GetNumOfGeometryPoints() const final {
        return vec_vertices_.size();
    }

    explicit GeometryInfoOrigin(const DefInt dims) : GeometryInfoInterface(dims) {
        this->node_type_ = "GeometryInfoOrigin";
    }
    virtual ~GeometryInfoOrigin() {}
};
class GeometryInfoOriginCreator :public GeometryInfoCreatorInterface {
 public:
    std::shared_ptr<GeometryInfoInterface> CreateGeometryInfo(const DefInt dims) const override {
        std::shared_ptr<GeometryInfoOrigin> ptr_tmp = std::make_shared<GeometryInfoOrigin>(dims);
        return ptr_tmp;
    };
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_INFO_ORIGIN_H_
