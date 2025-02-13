//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file input_parser.h
* @author Zhengliang Liu
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
class InputParser {
 private:
    char delimiter_ = ' ';
    std::map<std::string, std::string> map_input_;
    std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> nested_map_input_;

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
    const std::map<std::string, std::string>& GetMapInput() const { return map_input_; }
    const std::map<std::string, std::map<std::string, std::map<std::string, std::string>>>& GetNestedMapInput() const {
        return nested_map_input_;
    }

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

    // read single value
    template <typename T>
    bool GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        T* const ptr_variable) const {
        if (map_variable.find(key) == map_variable.end()) {
            return false;
        }
        *ptr_variable = CheckType<T>(map_variable.at(key));
        return true;
    }
    template <typename T>
    bool GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        T* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return Read(key, map_variable, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key, T* const ptr_variable) const {
        return GetValue(key, map_input_, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key, T* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_input_, ptr_variable);
    }
    template <typename T>
    bool GetNestedValue(const std::string& key, const std::string& name,
        const std::string& nested_key, T* const ptr_variable) const {
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
        *ptr_variable = CheckType<T>(value_it->second);
        return true;
    }

    // read std::vector
    template <typename T>
    bool GetValue(const std::string& key,
        const std::map<std::string, std::string>& map_variable, std::vector<T>* const ptr_variable) const {
        if (map_variable.find(key) == map_variable.end()) {
            return false;
        }
        std::vector<std::string> tokens = Split(map_variable.at(key), delimiter_);
        ptr_variable->resize(tokens.size());
        for (auto i = 0; i < tokens.size(); ++i) {
            ptr_variable->at(i) = CheckType<T>(tokens.at(i));
        }
        return true;
    }
    template <typename T>
    bool GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        std::vector<T>* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_variable, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key, std::vector<T>* const ptr_variable) const {
        return GetValue(key, map_input_, ptr_variable);
    }
    template <typename T>
    bool GetValue(const std::string& key,
        std::vector<T>* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_input_, ptr_variable);
    }

    template <typename T>
    bool GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        std::vector<T>* const ptr_variable, size_t expected_elements) const {
        if (map_variable.find(key) == map_variable.end()) {
            return false;
        }
        std::vector<std::string> tokens = Split(map_variable[key], delimiter_);
        if (tokens.size() != expected_elements) {
            throw std::runtime_error("Expected " + std::to_string(expected_elements) +
                " elements for " + key + ", got " + std::to_string(tokens.size()));
        }
        ptr_variable->clear();
        for (const auto& token : tokens) {
            ptr_variable->push_back(CheckType<T>(token));
        }
        return true;
    }
    template <typename T>
    bool GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        std::vector<T>* const ptr_variable, const std::string& expected_scope, size_t expected_elements) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_variable, ptr_variable, expected_elements);
    }
    template <typename T>
    bool GetValue(const std::string& key, std::vector<T>* const ptr_variable, size_t expected_elements) const {
        return GetValue(key, map_input_, ptr_variable, expected_elements);
    }
    template <typename T>
    bool GetValue(const std::string& key, std::vector<T>* const ptr_variable,
        const std::string& expected_scope, size_t expected_elements) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_input_, ptr_variable, expected_elements);
    }

    // read std::array
    template <typename T, std::size_t N>
    bool GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        std::array<T, N>* const ptr_variable) const {
        if (map_variable.find(key) == map_variable.end()) {
            return false;
        }

        std::vector<std::string> tokens = Split(map_variable.at(key), delimiter_);
        if (tokens.size() != N) {
            throw std::runtime_error("Expected " + std::to_string(N) +
                " elements for " + key + ", got " + std::to_string(tokens.size()));
        }

        for (std::size_t i = 0; i < N; ++i) {
            ptr_variable->at(i) = CheckType<T>(tokens[i]);
        }
        return true;
    }
    template <typename T, std::size_t N>
    void GetValue(const std::string& key, const std::map<std::string, std::string>& map_variable,
        std::array<T, N>* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_variable, ptr_variable);
    }
    template <typename T, std::size_t N>
    bool GetValue(const std::string& key,
        std::array<T, N>* const ptr_variable) const {
        return GetValue(key, map_input_, ptr_variable);
    }
    template <typename T, std::size_t N>
    void GetValue(const std::string& key,
        std::array<T, N>* const ptr_variable, const std::string& expected_scope) const {
        CheckScope(key, expected_scope);
        return GetValue(key, map_input_, ptr_variable);
    }
};
}  // end namespace amrproject
}  // end namespace rootproject
#endif  // SOURCE_AMR_PROJECT_IO_INPUT_PARSER_H_
