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
#include "criterion/geometry_default_shape.h"
namespace rootproject {
namespace lbmproject {
#ifndef  DEBUG_DISABLE_2D_FUNCTIONS
class GeoShapeIBCircle2D : public amrproject::GeoShapeDefaultCircle2D {
 public:
};
class GeoShapeIBLine2D : public amrproject::GeoShapeDefaultLine2D {
 public:
    void UpdateShape(const DefReal sum_t) override;
};
#endif  // DEBUG_DISABLE_2D_FUNCTIONS
#ifndef  DEBUG_DISABLE_3D_FUNCTIONS
class GeoShapeIBCubic3D : public amrproject::GeoShapeDefaultCubic3D {
 public:
};
class GeoShapeIBQuadrilateral3D : public amrproject::GeoShapeDefaultQuadrilateral3D {
 public:

};
#endif  // DEBUG_DISABLE_3D_FUNCTIONS
}  // end namespace lbmproject
}  // end namespace rootproject
#endif  // SOURCE_LBM_PROJECT_IMMERSED_BOUNDARY_GEOMETRY_IB_SHAPE_H_
