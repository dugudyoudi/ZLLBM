//  Copyright (c) 2021 - 2025, Zhengliang Liu
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
#include <cstddef>
#include <vector>
#include <mpi.h>
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
/**
 * class BitField
 * @brief class to store an arbitrary number of bits using the smallest number of bytes.
 */
class BitField {
 public:
    explicit BitField(std::size_t num_bits)
        : num_bits_(num_bits), num_bytes_((num_bits + 7) / 8), data_(num_bytes_, 0) {}
    explicit BitField(DefInt num_bits)
        : BitField(static_cast<std::size_t>(num_bits)) {
        if (num_bits < 0) {
            throw std::invalid_argument("Number of bits cannot be negative");
        }
    }

    void Clear() {
        std::fill(data_.begin(), data_.end(), 0);
    }

    bool operator==(const BitField& other) const {
        return num_bits_ == other.num_bits_ && data_ == other.data_;
    }
    bool operator!=(const BitField& other) const {
        return !(*this == other);
    }

    void Set(std::size_t pos, bool value) {
        if (pos >= num_bits_) throw std::out_of_range("Bit position out of range");
        std::size_t byte_index = pos / 8;
        std::size_t bit_index = pos % 8;
        if (value) {
            data_[byte_index] |= (1 << bit_index);
        }  else {
            data_[byte_index] &= ~(1 << bit_index);
        }
    }

    bool Get(std::size_t pos) const {
        if (pos >= num_bits_) throw std::out_of_range("Bit position out of range");
        std::size_t byte_index = pos / 8;
        std::size_t bit_index = pos % 8;
        return (data_[byte_index] >> bit_index) & 1;
    }

    void* Data() { return data_.data(); }
    const void* Data() const { return data_.data(); }

    std::size_t GetNumBytes() const { return num_bytes_; }
    std::size_t GetNumBits() const { return num_bits_; }

 private:
    std::size_t num_bits_;
    std::size_t num_bytes_;
    std::vector<uint8_t> data_;
};

// numerical parameters
const DefReal kEps = DefReal(1.e-6);
/**< Constant: a small number for comparing floats*/
const DefReal kPi = DefReal(4.) * atan(DefReal(1.));  /**< Constant: pi*/

// global constants
static constexpr DefInt kXIndex = 0, kYIndex = 1, kZIndex = 2;
/**< indices of vectors for x, y, z directions*/
}  // end namespace rootproject
#endif  // SOURCE_DEFS_LIBS_H_
