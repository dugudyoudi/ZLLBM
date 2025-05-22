//  Copyright (c) 2021 - 2025, Zhengliang Liu
//  All rights reserved

/**
* @file lbm_les_models.cpp
* @author Zhengliang Liu
* @brief functions used for managing LES model.
*/
#include <string>
#include "./lbm_interface.h"
namespace rootproject {
namespace lbmproject {
/**
 * @brief function to set LES model based on input key word.
 * @param[in] model_type identifier of a LES model
 */
void SolverLbmInterface::SetLesModel(const std::string& model_type) {
    if (model_type == "smagorinsky") {
        SetLesModel(ELbmLesModelType::kSmagorinsky);
    } else {
        amrproject::LogManager::LogError("LES model " + model_type + "  is not defined, "
            "please choosing from: smagorinsky");
    }
}
/**
 * @brief function to set LES model based on input enumerate.
 * @param[in] model_type identifier of a LES model
 */
void SolverLbmInterface::SetLesModel(const ELbmLesModelType model_type) {
    bool_les_model_ = true;
    switch (model_type) {
        case ELbmLesModelType::kSmagorinsky:
            ptr_les_model_ = std::make_unique<LesModelSmagorinsky>();
            break;
        default:
            amrproject::LogManager::LogError("LES model type is not defined");
            break;
    }
}
}  // end namespace lbmproject
}  // end namespace rootproject
