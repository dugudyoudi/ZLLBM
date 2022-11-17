#include "grid/grid_info_preset.h"
namespace rootproject {
namespace amrproject {
namespace grid {
std::shared_ptr<TrackingGridInfoInterface>
TrackingGridInfoGeoCreator::CreateTrackingGridInfo() {
    std::shared_ptr<TrackingGridInfoGeo> ptr_temp =
        std::make_shared<TrackingGridInfoGeo>();
    ptr_temp->node_type_ = "TrackingGeo";
    return ptr_temp;
}
}  // end namsapce grid
}  // end namespace amrproject
}  // end namespace rootproject