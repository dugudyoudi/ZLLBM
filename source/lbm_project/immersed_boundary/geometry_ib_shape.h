//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_ib_shape.h
* @author Zhengliang Liu
* @date  2022-5-25
* @brief  define classes to manage default shapes for immersed boundary method.
*/

#ifndef  SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_GEOMETRY_IB_SHAPE_H_
#define  SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_GEOMETRY_IB_SHAPE_H_
#include <array>
#include <memory>
#include <string>
#include "criterion/geometry_default_shape.h"
namespace rootproject {
namespace lbmproject {
class GeoIBShapeReader : public amrproject::GeoShapeReader {
 public:
    std::unique_ptr<amrproject::GeoShapeInterface> ReadAndCreateGeoShape(
        const DefInt dims, const std::string& shape_id,
        const std::weak_ptr<amrproject::GeometryInfoInterface>& ptr_geo) const override;
    virtual ~GeoIBShapeReader() {}
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeoShapeIBCircle2D : public amrproject::GeoShapeDefaultCircle2D {
 public:
    void UpdateShape(const DefReal sum_t) override;
    explicit GeoShapeIBCircle2D(const std::weak_ptr<amrproject::GeometryInfoInterface> ptr_geo)
        : amrproject::GeoShapeDefaultCircle2D(ptr_geo) {}
};
class GeoShapeIBLine2D : public amrproject::GeoShapeDefaultLine2D {
 public:
    void UpdateShape(const DefReal sum_t) override;
    explicit GeoShapeIBLine2D(const std::weak_ptr<amrproject::GeometryInfoInterface> ptr_geo)
        : amrproject::GeoShapeDefaultLine2D(ptr_geo) {}
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeoShapeIBCube3D : public amrproject::GeoShapeDefaultCube3D {
 public:
    explicit GeoShapeIBCube3D(const std::weak_ptr<amrproject::GeometryInfoInterface> ptr_geo)
    : amrproject::GeoShapeDefaultCube3D(ptr_geo) {}
};
class GeoShapeIBQuadrilateral3D : public amrproject::GeoShapeDefaultQuadrilateral3D {
 public:
    void UpdateShape(const DefReal sum_t) override;
    explicit GeoShapeIBQuadrilateral3D(const std::weak_ptr<amrproject::GeometryInfoInterface> ptr_geo)
    : amrproject::GeoShapeDefaultQuadrilateral3D(ptr_geo) {}
};
class GeoShapeIBSphere3D : public amrproject::GeoShapeDefaultSphere3D {
 public:
    void UpdateShape(const DefReal sum_t) override;
    explicit GeoShapeIBSphere3D(const std::weak_ptr<amrproject::GeometryInfoInterface> ptr_geo)
    : amrproject::GeoShapeDefaultSphere3D(ptr_geo) {}
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_GEOMETRY_IB_SHAPE_H_
