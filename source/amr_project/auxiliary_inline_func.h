//  Copyright (c) 2021 - 2023, Zhengliang Liu
//  All rights reserved

/**
* @file auxiliary_inline_func.h
* @author Zhengliang Liu
* @brief amp project inline functions
* @date  2022-9-3
* @note .
*/
#ifndef ROOTPROJECT_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
#define ROOTPROJECT_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_
#include <array>
#include <fstream>
#include "../defs_libs.h"
namespace rootproject {
namespace amrproject {
inline DefReal Square(DefReal x) {
    return x * x;
}
/**
* @brief function to calculate 2 power N.
* @param[in]  n        power exponent.
*/
inline DefSizet TwoPowerN(DefSizet n) {
    DefSizet two = 1;
    two <<= n;
    return two;
}
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // ROOTPROJECT_AMR_PROJECT_AUXILIARY_INLINE_FUNC_H_