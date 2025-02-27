//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file geometry_default_shape.h
* @author Zhengliang Liu
* @date  2022-5-25
* @brief  define classes to generate default shapes.
*/

#ifndef SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_DEFAULT_SHAPE_H_
#define SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_DEFAULT_SHAPE_H_
#include <array>
#include <limits>
#include <memory>
#include <string>
#include <map>
#include "../defs_libs.h"
#include "io/input_parser.h"
namespace rootproject {
namespace amrproject {
class  GeometryInfoInterface;
class GeoShapeInterface {
 public:
    virtual bool ReadAndSetGeoShapeParameters(std::map<std::string, ParserData>* const ptr_shape_parameters) {
        return true;
    }
    virtual void InitialShape(const DefReal dx) = 0;
    virtual void UpdateShape(const DefReal sum_t) {}
    virtual ~GeoShapeInterface() {}
    explicit GeoShapeInterface(const std::weak_ptr<GeometryInfoInterface> ptr_geo) : ptr_geo_info_(ptr_geo) {}

 protected:
    std::weak_ptr<GeometryInfoInterface> ptr_geo_info_;
};
class GeoShapeReader {
 public:
    virtual std::unique_ptr<GeoShapeInterface> ReadAndCreateGeoShape(
        const DefInt dims, const std::string& shape_id, const std::weak_ptr<GeometryInfoInterface>& ptr_geo) const;
    virtual ~GeoShapeReader() {}
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeoShapeDefaultCircle2D : public GeoShapeInterface {
 public:
    DefReal radius_ = 0.5;
    std::array<DefReal, 2> center_;
    bool ReadAndSetGeoShapeParameters(std::map<std::string, ParserData>* const ptr_shape_parameters) override;
    void InitialShape(const DefReal dx) override;
    virtual ~GeoShapeDefaultCircle2D() {}
    explicit GeoShapeDefaultCircle2D(const std::weak_ptr<GeometryInfoInterface> ptr_geo) : GeoShapeInterface(ptr_geo) {}
};
class GeoShapeDefaultLine2D : public GeoShapeInterface {
 public:
    std::array<DefReal, 2> start_point_, end_point_;
    bool ReadAndSetGeoShapeParameters(std::map<std::string, ParserData>* const ptr_shape_parameters) override;
    void InitialShape(const DefReal dx) override;
    virtual ~GeoShapeDefaultLine2D() {}
    explicit GeoShapeDefaultLine2D(const std::weak_ptr<GeometryInfoInterface> ptr_geo) : GeoShapeInterface(ptr_geo) {}
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeoShapeDefaultCube3D : public GeoShapeInterface {
 public:
    DefReal length_ = 1.;
    std::array<DefReal, 3> center_;
    bool ReadAndSetGeoShapeParameters(std::map<std::string, ParserData>* const ptr_shape_parameters) override;
    void InitialShape(const DefReal dx) override;
    virtual ~GeoShapeDefaultCube3D() {}
    explicit GeoShapeDefaultCube3D(const std::weak_ptr<GeometryInfoInterface> ptr_geo) : GeoShapeInterface(ptr_geo) {}
};
class GeoShapeDefaultQuadrilateral3D : public GeoShapeInterface {
 public:
    std::array<DefReal, 3> start_point_, neighbor_point_, diagonal_point_;
    bool ReadAndSetGeoShapeParameters(std::map<std::string, ParserData>* const ptr_shape_parameters) override;
    void InitialShape(const DefReal dx) override;
    virtual ~GeoShapeDefaultQuadrilateral3D() {}
    explicit GeoShapeDefaultQuadrilateral3D(const std::weak_ptr<GeometryInfoInterface> ptr_geo)
        : GeoShapeInterface(ptr_geo) {}

 protected:
    DefInt edge1_num_points_ = 0, edge2_num_points_ = 0;
};
class GeoShapeDefaultSphere3D : public GeoShapeInterface {
 public:
    DefReal radius_ = 0.5;
    std::array<DefReal, 3> center_;
    bool ReadAndSetGeoShapeParameters(std::map<std::string, ParserData>* const ptr_shape_parameters) override;
    void InitialShape(const DefReal dx) override;
    virtual ~GeoShapeDefaultSphere3D() {}
    explicit GeoShapeDefaultSphere3D(const std::weak_ptr<GeometryInfoInterface> ptr_geo) : GeoShapeInterface(ptr_geo) {}

 protected:
    DefInt method_generate_sphere_ = 2;
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_DEFAULT_SHAPE_H_
