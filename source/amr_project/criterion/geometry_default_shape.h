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
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
class   GeometryInfoInterface;
class GeoShapeInterface {
 public:
    GeometryInfoInterface* ptr_geo_info_ = nullptr;
    virtual void InitialShape(const DefReal dx) = 0;
    virtual void UpdateShape(const DefReal sum_t) = 0;
    virtual ~GeoShapeInterface() {}
};
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeoShapeDefaultCircle2D : public GeoShapeInterface {
 public:
    DefReal radius_ = 0.5;
    void InitialShape(const DefReal dx) override;
    void UpdateShape(const DefReal sum_t) override;
    virtual ~GeoShapeDefaultCircle2D() {}
};
class GeoShapeDefaultLine2D : public GeoShapeInterface {
 public:
    std::array<DefReal, 2> start_point_, end_point_;
    void InitialShape(const DefReal dx) override;
    void UpdateShape(const DefReal sum_t) override;
    virtual ~GeoShapeDefaultLine2D() {}
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeoShapeDefaultCubic3D : public GeoShapeInterface {
 public:
    DefReal length_ = 1.;
    void InitialShape(const DefReal dx) override;
    void UpdateShape(const DefReal sum_t) override;
    virtual ~GeoShapeDefaultCubic3D() {}
};
class GeoShapeDefaultQuadrilateral3D : public GeoShapeInterface {
 public:
    std::array<DefReal, 3> start_point_, neighbor_point_, diagonal_point_;
    void InitialShape(const DefReal dx) override;
    void UpdateShape(const DefReal sum_t) override;
    virtual ~GeoShapeDefaultQuadrilateral3D() {}

 protected:
    DefInt edge1_num_points_ = 0, edge2_num_points_ = 0;
};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_CRITERION_GEOMETRY_DEFAULT_SHAPE_H_
