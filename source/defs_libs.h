//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file defs_libs.h
* @author Zhengliang Liu
* @date  2022-5-12
* @brief contain definations and libraries will be used in all modules.
*/

#ifndef SOURCE_DEFS_LIBS_H_
#define SOURCE_DEFS_LIBS_H_
// libraries
#include <cmath>
#include <bitset>
#include <cstdint>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "./config.h"  // configuration file generated
namespace rootproject {
using DefReal = double;
using DefInt = int32_t;
using DefUint = uint32_t;
using DefSizet = size_t;
// unsigned integer for indices will be long like indices of nodes and vertices
using DefAmrLUint = uint32_t;
// unsigned integer for space filling code
using DefSFCodeToUint = uint64_t;

// bitset to store code for space filling curves
static constexpr DefInt kSFBitsetBit = 64;
using DefSFBitset = std::bitset<kSFBitsetBit>;
class HashFunc {
 public:
    DefSFCodeToUint operator()(const DefSFBitset&bit_set_in) const noexcept {
        DefSFCodeToUint hashVal = bit_set_in.to_ullong();
        return hashVal;
    }
};
template <typename value>
using DefMap = std::unordered_map<DefSFBitset, value, HashFunc>;

// numerical parameters
const DefReal kEps = DefReal(1.e-6);
/**< Constant: a small number for comparing floats*/
const DefReal kPi = DefReal(4.) * atan(DefReal(1.));  /**< Constant: pi*/

// global constants
static constexpr DefInt kXIndex = 0, kYIndex = 1, kZIndex = 2;
/**< indices of vectors for x, y, z directions*/
}  // end namespace rootproject
#endif  // SOURCE_DEFS_LIBS_H_
