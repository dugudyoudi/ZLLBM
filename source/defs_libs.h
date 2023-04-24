//  Copyright (c) 2022, Zhengliang Liu
//  All rights reserved

/**
* @file def_libs.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief contain definations and libaries will be used in all modules.
*/

#ifndef ROOTPROJECT_SOURCE_DEF_LIB_H_
#define ROOTPROJECT_SOURCE_DEF_LIB_H_
// libraries
#include <bitset>
#include <cstdint>
#include <unordered_map>
#include <iostream>
#include "./config.h"  // configuration file generated
                     // by CMAKE based on "config.h.in"
namespace rootproject {
using DefReal = double;
using DefInt = int;
using DefLInt = int;
using DefUint = unsigned int;
using DefSizet = size_t;
using DefLUint = unsigned int;
using DefTypeUint = unsigned int;
using DefSFCodeToUint = uint64_t;

// bitset to store code for space filling curves 
static constexpr DefUint kSFBitsetBit = 64;
using DefSFBitset = std::bitset<kSFBitsetBit>;
class HashFunc {
 public:
    size_t operator()(const DefSFBitset&bit_set_in) const noexcept {
        size_t hashVal =bit_set_in.to_ullong();
        return hashVal;
    }
};
template <typename value>
using DefMap = std::unordered_map<DefSFBitset, value, HashFunc>;

// numerical parameters
const DefReal kEps = 1.e-10;
/**< Constant: a small number for comparing floats*/
const DefReal kPi = 4. * atan(1.);  /**< Constant: pi*/

// globle constants
static constexpr DefUint kXIndex = 0, kYIndex = 1, kZIndex = 2;
/**< indices of vectors for x, y, z diretions*/
}  // end namespace rootproject
#endif  // ROOTPROJECT_SOURCE_DEF_LIB_H_
