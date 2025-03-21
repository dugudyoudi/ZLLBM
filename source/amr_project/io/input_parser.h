//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file input_parser.h
* @author Zhengliang Liu
* @brief define the class used to read input parameters.
* @date  2025-2-1
*/
#ifndef SOURCE_AMR_PROJECT_IO_INPUT_PARSER_H_
#define SOURCE_AMR_PROJECT_IO_INPUT_PARSER_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <string>
#include <type_traits>
namespace rootproject {
namespace amrproject {
template <typename T>
    concept ParserAllowedType = std::is_fundamental_v<T> || std::is_same_v<T, std::string>;
struct ParserData {
 public:
    bool bool_has_read_ = false;
    std::string str_value_;
};
class InputParser {
 private:
    char delimiter_ = ' ';
    std::map<std::string, ParserData> map_input_{};
    std::map<std::string, std::map<std::string, std::map<std::string, ParserData>>> nested_map_input_{};

    std::vector<std::string> Split(const std::string& str, char delimiter) const;
    std::string Trim(const std::string& str) const;
    void CheckScope(const std::string& key, const std::string& expected_scope) const;

    template <typename T>
    T CheckType(const std::string& value) const {
        if constexpr (std::is_same_v<T, bool>) {
            if (value == "true" || value == "1") {
                return true;
            } else if (value == "false" || value == "0") {
                return false;
            } else {
                throw std::runtime_error("Invalid type: expected boolean value");
            }
        } else if constexpr (std::is_integral_v<T>) {
            try {
                return static_cast<T>(std::stoll(value));
            } catch (const std::invalid_argument&) {
                throw std::runtime_error("Invalid type: expected integral number");
            }
        } else if constexpr (std::is_floating_point_v<T>) {
            try {
                return static_cast<T>(std::stod(value));
            } catch (const std::invalid_argument&) {
                throw std::runtime_error("Invalid type: expected floating-point number");
            }
        } else if constexpr (std::is_same_v<T, std::string>) {
            return value;
        } else {
            throw std::runtime_error("Unsupported type");
        }
    }

    void ParseNestedStructure(const std::string& key, const std::string& filename, std::istringstream& stream);

 public:
    bool print_values_when_read_ = false;
    std::map<std::string, ParserData>& GetMapInput()  { return map_input_; }
    std::map<std::string, std::map<std::string, std::map<std::string, ParserData>>>& GetNestedMapInput() {
        return nested_map_input_;
    }

    void PrintUnusedParameters() const;

    explicit InputParser(const std::string& filename);
    InputParser() {}

    template <ParserAllowedType T>
    std::string ValuesToOutputStr(const std::vector<T>& input) const {
        std::string result = "(";
        for (auto iter = input.begin(); iter != std::prev(input.end()); ++iter) {
            result += std::to_string(*iter) + ", ";
        }
        if (!input.empty()) {
            result += std::to_string(input.back());
        }
        result += ")";
        return result;
    }
    template <ParserAllowedType T, std::size_t N>
    std::string ValuesToOutputStr(const std::array<T, N>& input) const {
        std::string result = "(";
        for (auto iter = input.begin(); iter != std::prev(input.end()); ++iter) {
            result += std::to_string(*iter) + ", ";
        }
        if (!input.empty()) {
            result += std::to_string(input.back());
        }
        result += ")";
        return result;
    }
    template <typename T>
    bool GetNestedValue(const std::string& key, const std::string& name,
        const std::string& nested_key, T* const ptr_variable) {
        auto it = nested_map_input_.find(key);
        if (it == nested_map_input_.end()) {
            return false;
        }
        auto nested_it = it->second.find(name);
        if (nested_it == it->second.end()) {
            return false;
        }
        auto value_it = nested_it->second.find(nested_key);
        if (value_it == nested_it->second.end()) {
            return false;
        }
        *ptr_variable = CheckType<T>(value_it->second.str_value_);
        value_it->second.bool_has_read_ = true;
        return true;
    }

    // read single value
    template <typename T>
    bool GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        T* const ptr_variable) const {
        if (ptr_map_variable->find(key) == ptr_map_variable->end()) {
            return false;
        }
        *ptr_variable = CheckType<T>(ptr_map_variable->at(key).str_value_);
        ptr_map_variable->at(key).bool_has_read_ = true;
        return true;
    }
    template <typename T>
    bool GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        T* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return Read(key, ptr_map_variable, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key, T* const ptr_variable) {
        return GetValue(key, &map_input_, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key, T* const ptr_variable, const std::string& expected_scope) {
        CheckScope(key, expected_scope);
        return GetValue(key, &map_input_, ptr_variable);
    }

    // read std::vector
    template <typename T>
    bool GetValue(const std::string& key,
        std::map<std::string, ParserData>* const ptr_map_variable, std::vector<T>* const ptr_variable) const {
        if (ptr_map_variable->find(key) == ptr_map_variable->end()) {
            return false;
        }
        std::vector<std::string> tokens = Split(ptr_map_variable->at(key).str_value_, delimiter_);
        ptr_variable->resize(tokens.size());
        for (size_t i = 0; i < tokens.size(); ++i) {
            ptr_variable->at(i) = CheckType<T>(tokens.at(i));
        }
        ptr_map_variable->at(key).bool_has_read_ = true;
        return true;
    }
    template <typename T>
    bool GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        std::vector<T>* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, ptr_map_variable, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key, std::vector<T>* const ptr_variable) {
        return GetValue(key, &map_input_, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key,
        std::vector<T>* const ptr_variable, const std::string& expected_scope) {
        CheckScope(key, expected_scope);
        return GetValue(key, &map_input_, ptr_variable);
    }

    template <typename T>
    bool GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        std::vector<T>* const ptr_variable, size_t expected_elements) const {
        if (ptr_map_variable->find(key) == ptr_map_variable->end()) {
            return false;
        }
        std::vector<std::string> tokens = Split(ptr_map_variable->at(key).str_value_, delimiter_);
        if (tokens.size() != expected_elements) {
            throw std::runtime_error("Expected " + std::to_string(expected_elements) +
                " elements for " + key + ", got " + std::to_string(tokens.size()));
        }
        ptr_variable->clear();
        for (const auto& token : tokens) {
            ptr_variable->push_back(CheckType<T>(token));
        }
        ptr_map_variable->at(key).bool_has_read_ = true;
        return true;
    }
    template <typename T>
    bool GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        std::vector<T>* const ptr_variable, const std::string& expected_scope, size_t expected_elements) const {
        CheckScope(key, expected_scope);
        return GetValue(key, ptr_map_variable, ptr_variable, expected_elements);
    }
    template <typename T>
    bool GetValue(const std::string& key, std::vector<T>* const ptr_variable, size_t expected_elements) {
        return GetValue(key, &map_input_, ptr_variable, expected_elements);
    }
    template <typename T>
    bool GetValue(const std::string& key, std::vector<T>* const ptr_variable,
        const std::string& expected_scope, size_t expected_elements) {
        CheckScope(key, expected_scope);
        return GetValue(key, &map_input_, ptr_variable, expected_elements);
    }

    // read std::array
    template <typename T, std::size_t N>
    bool GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        std::array<T, N>* const ptr_variable) const {
        if (ptr_map_variable->find(key) == ptr_map_variable->end()) {
            return false;
        }

        std::vector<std::string> tokens = Split(ptr_map_variable->at(key).str_value_, delimiter_);
        if (tokens.size() != N) {
            throw std::runtime_error("Expected " + std::to_string(N) +
                " elements for " + key + ", got " + std::to_string(tokens.size()));
        }

        for (std::size_t i = 0; i < N; ++i) {
            ptr_variable->at(i) = CheckType<T>(tokens[i]);
        }
        ptr_map_variable->at(key).bool_has_read_ = true;
        return true;
    }
    template <typename T, std::size_t N>
    void GetValue(const std::string& key, std::map<std::string, ParserData>* const ptr_map_variable,
        std::array<T, N>* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, ptr_map_variable, ptr_variable);
    }
    template <typename T, std::size_t N>
    bool GetValue(const std::string& key, std::array<T, N>* const ptr_variable) {
        return GetValue(key, &map_input_, ptr_variable);
    }
    template <typename T, std::size_t N>
    void GetValue(const std::string& key,
        std::array<T, N>* const ptr_variable, const std::string& expected_scope) {
        CheckScope(key, expected_scope);
        return GetValue(key, &map_input_, ptr_variable);
    }
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_IO_INPUT_PARSER_H_
