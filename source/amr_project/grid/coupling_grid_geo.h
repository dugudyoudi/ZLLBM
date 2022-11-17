//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file /grid/coupling_grid_geo.h
* @author Zhengliang Liu
* @date  2022-8-24
* @brief preset class for tracking and ghost nodes
*/
#ifndef ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_COUPLING_GRID_GEO_H_
#define ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_COUPLING_GRID_GEO_H_
#include "grid/grid_info_interface.h"
namespace rootproject {
namespace amrproject{
namespace grid {
class TrackingGridInfoGeo :public TrackingGridInfoInterface {
};
class TrackingGridInfoGeoCreator :public TrackingGridInfoCreatorInterface {
public:
   std::shared_ptr<TrackingGridInfoInterface>
        CreateTrackingGridInfo() override {
       return std::make_shared<TrackingGridInfoGeo>();
   };
};
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
}  // end namsapce grid
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_AMR_PROJECT_GRID_COUPLING_GRID_GEO_H_
