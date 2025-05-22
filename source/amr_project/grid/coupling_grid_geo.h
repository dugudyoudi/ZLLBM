//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/coupling_grid_geo.h
* @author Zhengliang Liu
* @date  2022-8-24
* @brief preset class for tracking and ghost nodes
*/
#ifndef SOURCE_AMR_PROJECT_GRID_COUPLING_GRID_GEO_H_
#define SOURCE_AMR_PROJECT_GRID_COUPLING_GRID_GEO_H_
#include <memory>
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace amrproject {
class GhostGridInfoGeo : public GhostGridInfoInterface {
    void InitialGhostNode(const DefSFBitset& bitset_in) override {};
};
/**
* @class GhostGridInfoInterface
* @brief abstract class used to generate GhostGridInfo instance.
*/
class GhostGridInfoGeoCreator : public GhostGridInfoCreatorInterface {
 public:
    std::shared_ptr<GhostGridInfoInterface>
        CreateGhostGridInfo() override {
        return std::make_shared<GhostGridInfoGeo>();
    };
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_GRID_COUPLING_GRID_GEO_H_
