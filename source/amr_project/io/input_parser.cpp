//  Copyright (c) 2021 - 2024, Zhengliang Liu
//  All rights reserved

/**
* @file input_parser.cpp
* @author Zhengliang Liu
* @brief functions used to read input parameters.
* @date  2025-2-1
*/
#include <string>
#include <vector>
#include "io/input_parser.h"
#include "io/log_write.h"
namespace rootproject {
namespace amrproject {
/**
* @brief      function to split a string into parameters
* @param[in]  str    intput string.
* @param[in] delimiter used to split the string.
*/
std::vector<std::string> InputParser::Split(const std::string& str, char delimiter) const {
    std::vector<std::string> parameters;
    std::string element;
    std::istringstream param_str(str);
    while (std::getline(param_str, element, delimiter)) {
        parameters.push_back(element);
    }
    return parameters;
}
/**
* @brief  function to trim whitespace from a string
* @param[in] str    intput string.
*/
std::string InputParser::Trim(const std::string& str) const {
    size_t comment_pos = str.find('#');
    // If '#' is found, ignore everything after it
    std::string trimmed_str = (comment_pos != std::string::npos) ? str.substr(0, comment_pos) : str;

    size_t first = trimmed_str.find_first_not_of(' ');
    if (first == std::string::npos) return "";
    size_t last = trimmed_str.find_last_not_of(' ');
    return trimmed_str.substr(first, last - first + 1);
}
/**
* @brief  function to the key prefix matches the expected scope
* @param[in] key    intput string for the key.
* @param[in] expected_scope  expected scope.
*/
void InputParser::CheckScope(const std::string& key, const std::string& expected_scope) const {
    size_t prefix_end = key.find('.');
    if (prefix_end == std::string::npos) {
        LogManager::LogError("Invalid input format: missing scope prefix");
    }
    std::string scope = key.substr(0, prefix_end);
    if (scope != expected_scope) {
        LogManager::LogError("Variable " + key + " does not belong to scope " + expected_scope);
    }
}
/**
* @brief  constructor of InputParser
* @param[in] filename  name of input file.
*/
InputParser::InputParser(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        LogManager::LogError("Could not open file: " + filename);
    }

    std::string line;
    bool inside_block = false;
    std::string block_key;
    std::string block_content;
    ParserData parser_data;
    while (std::getline(file, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;

        size_t delimiter_pos = line.find('=');

        std::string key = Trim(line.substr(0, delimiter_pos));
        std::string value = Trim(line.substr(delimiter_pos + 1));

        if (value.front() == '{') {
            inside_block = true;
            block_key = key;
            block_content.clear();
            continue;
        } else if (line.back() == '}') {
            inside_block = false;
            block_content += line.substr(0, line.size() - 1);
            std::istringstream nested_stream(block_content);
            ParseNestedStructure(block_key, filename, nested_stream);
            continue;
        }

        if (inside_block) {
            block_content += line + "\n";
        } else {
            if (map_input_.find(key) == map_input_.end()) {
                parser_data.bool_has_read_ = false;
                parser_data.str_value_ = value;
                map_input_.insert({key, parser_data});
            } else {
                LogManager::LogError("Duplicate key word: " + key + " in input file " + filename);
            }
        }
    }
}
/**
* @brief function to assign inputs in curly braces to nested map.
* @param[in] key   key outside curly braces
* @param[in] filename  name of input file.
* @param[in] stream  string stream storing inputs in curly braces.
*/
void InputParser::ParseNestedStructure(const std::string& key,
    const std::string& filename, std::istringstream& stream) {
    std::string line;
    ParserData parser_data;
    std::map<std::string, ParserData> nested_map;
    std::string name;
    while (std::getline(stream, line)) {
        line = Trim(line);
        if (line.empty()) {
            continue;
        }

        auto pos = line.find('=');
        if (pos == std::string::npos) {
            LogManager::LogError("Invalid nested line: " + line + " in input file" + filename);
        }

        std::string nested_key = Trim(line.substr(0, pos));
        parser_data.bool_has_read_ = false;
        parser_data.str_value_ = Trim(line.substr(pos + 1));

        if (nested_key == "identifier") {
            name = parser_data.str_value_;
        } else {
            nested_map[nested_key] = parser_data;
        }
    }

    if (name.empty()) {
        LogManager::LogError("The 'identifier' key is required for nested structure: " + key
            + " in input file" + filename);
    }

    // Check if the name already exists for the given key
    if (nested_map_input_[key].find(name) != nested_map_input_[key].end()) {
        LogManager::LogError("Duplicate name is not allowed for key: " + key + " with name: " + name
            + " in input file " + filename);
    }
    nested_map_input_.at(key).insert({name, nested_map});
}
/**
* @brief function to print unused parameters in input setting.
*/
void InputParser::PrintUnusedParameters() const {
    for (const auto & iter : map_input_) {
        if (!iter.second.bool_has_read_) {
            LogManager::LogWarning("Input " + iter.first + " has not been used.");
        }
    }
    for (const auto& iter : nested_map_input_) {
        for (const auto& iter_id : iter.second) {
            for (const auto& iter_key : iter_id.second) {
                if (!iter_key.second.bool_has_read_) {
                    LogManager::LogWarning("Input " + iter_key.first + " with identifier "
                        + iter_id.first + " in scope " + iter.first + " has not been used.");
                }
            }
        }
    }
}

}  // end namespace amrproject
}  // end namespace rootproject
