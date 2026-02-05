//
//  InparamYAML.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  YAML parser for input parameters
//  based on mini-yaml:
//  https://github.com/jimmiebergmann/mini-yaml

#ifndef InparamYAML_hpp
#define InparamYAML_hpp

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include "yaml/Yaml.hpp"
#pragma clang diagnostic pop

#include "bstring.hpp"

class InparamYAML {
public:
    // constructor
    InparamYAML(const std::string &name):
    mName(name), mRoot() {
        // nothing
    }
    
    // parse
    void parse(const std::string &fname);
    
    // verbose
    std::string verbose() const;
    
    // get value by keyword vector
    template <typename T>
    T get(const std::vector<std::string> &keyword) const {
        // initialize nodePtr with root
        const Yaml::Node *nodePtr = &mRoot;
        
        // recursively increase depth of nodePtr
        std::string errKey = "";
        for (int depth = 0; depth < keyword.size(); depth++) {
            // update key series for error messege
            errKey += keyword[depth];
            
            // array or object?
            if (keyword[depth].front() == '[' ||
                keyword[depth].front() == '{') {
                ////////////////////// array //////////////////////
                // check sequence type
                if (!nodePtr->IsSequence() &&
                    nodePtr->As<std::string>() != "YAML_EMPTY_SEQUENCE") {
                    throw std::runtime_error("InparamYAML::get || "
                                             "Non-array parent node. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
                
                // check and get index string
                if ((keyword[depth].front() == '[' &&
                     keyword[depth].back() != ']') ||
                    (keyword[depth].front() == '{' &&
                     keyword[depth].back() != '}')) {
                    throw std::runtime_error("InparamYAML::get || "
                                             "Unmatched array brackets. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
                std::string indexStr =
                keyword[depth].substr(1, keyword[depth].size() - 2);
                
                // get count by [?]
                if (indexStr == "?") {
                    if (depth != keyword.size() - 1) {
                        throw std::runtime_error("InparamYAML::get || "
                                                 "Array length operator [?] "
                                                 "must be the last keyword. || "
                                                 "Keyword: " + errKey + " || "
                                                 "In YAML: " + mName);
                    }
                    // first check empty
                    if (nodePtr->As<std::string>() == "YAML_EMPTY_SEQUENCE") {
                        return bstring::cast<T>("0", "InparamYAML::get");
                    }
                    return bstring::cast<T>(bstring::toString(nodePtr->Size()),
                                            "InparamYAML::get");
                }
                
                // get key instead of value
                bool returnKey = false;
                if (keyword[depth].front() == '{') {
                    if (depth != keyword.size() - 1) {
                        throw std::runtime_error("InparamYAML::get || "
                                                 "Array key operator {#} "
                                                 "must be the last keyword. || "
                                                 "Keyword: " + errKey + " || "
                                                 "In YAML: " + mName);
                    }
                    returnKey = true;
                }
                
                // get index
                int index;
                try {
                    index = bstring::cast<int>(indexStr, "InparamYAML::get");
                } catch(...) {
                    throw std::runtime_error("InparamYAML::get || "
                                             "Error parsing array index. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
                
                // first check empty
                if (nodePtr->As<std::string>() == "YAML_EMPTY_SEQUENCE") {
                    throw std::runtime_error("InparamYAML::get || "
                                             "Index out of range. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
                
                // check and get array element
                if (index < nodePtr->Size()) {
                    // access array element by iter because
                    // mini-yaml does not provide const random access
                    auto iter = nodePtr->Begin();
                    for (int imove = 0; imove < index; imove++) {
                        iter++;
                    }
                    
                    if (returnKey) {
                        // move deeper
                        nodePtr = &((*iter).second);
                        // return key
                        return bstring::cast<T>((*nodePtr->Begin()).first,
                                                "InparamYAML::get");
                    } else {
                        // update valuePtr to array element
                        nodePtr = &((*iter).second);
                    }
                } else {
                    throw std::runtime_error("InparamYAML::get || "
                                             "Index out of range. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
            } else {
                ////////////////////// object //////////////////////
                // check and get object member
                if (nodePtr->IsMap() &&
                    nodePtr->As<std::string>() != "YAML_EMPTY_SEQUENCE") {
                    // access object member by iter because
                    // mini-yaml does not provide const random access
                    auto it = nodePtr->Begin();
                    for (; it != nodePtr->End(); it++) {
                        if ((*it).first == keyword[depth]) {
                            break;
                        }
                    }
                    // key not found
                    if (it == nodePtr->End()) {
                        throw std::runtime_error("InparamYAML::get || "
                                                 "Error finding keyword. || "
                                                 "Keyword: " + errKey + " || "
                                                 "In YAML: " + mName);
                    }
                    // update valuePtr to object element
                    nodePtr = &((*it).second);
                } else {
                    throw std::runtime_error("InparamYAML::get || "
                                             "Non-object parent node. || "
                                             "Keyword: " + errKey + " || "
                                             "In YAML: " + mName);
                }
            }
            
            // update key series for error messege
            if (depth < keyword.size() - 1) {
                errKey += ":";
            }
        }
        
        // get string and use bstring to cast
        const std::string &res = nodePtr->As<std::string>();
        try {
            return bstring::cast<T>(res, "InparamYAML::get");
        } catch (...) {
            throw std::runtime_error("InparamYAML::get || "
                                     "Error casting value: " + res + " || "
                                     "Keyword: " + errKey + " || "
                                     "In YAML: " + mName);
        }
    }
    
    // get value by keyword string
    template <typename T>
    T get(const std::string &keyword) const {
        return get<T>(bstring::split(keyword, ":"));
    }
    
    // this function is disabled for safety
    // get value with default
    // template <typename T>
    // T getWithDefault(const std::string &keyword, const T &defaultVal) const {
    //     try {
    //         return get<T>(keyword);
    //     } catch (...) {
    //         return defaultVal;
    //     }
    // }
    
    // get value with string-typed options
    template <typename T>
    T getWithOptions(const std::string &keyword,
                     const std::map<std::string, T> &options) const {
        const std::string &res = get<std::string>(keyword);
        try {
            return options.at(res);
        } catch (...) {
            try {
                return bstring::cast<T>(res, "InparamYAML::getWithOptions");
            } catch (...) {
                throw std::runtime_error("InparamYAML::getWithOptions || "
                                         "Error casting value: " + res + " || "
                                         "Keyword: " + keyword + " || "
                                         "In YAML: " + mName);
            }
        }
    }
    
    // get value with string-typed limits
    template <typename T>
    T getWithLimits(const std::string &keyword,
                    const std::map<std::string, T> &limits) const {
        const std::string &res = get<std::string>(keyword);
        try {
            return limits.at(res);
        } catch (...) {
            throw std::runtime_error("InparamYAML::getWithLimits || "
                                     "Unacceptable value: " + res + " || "
                                     "Keyword: " + keyword + " || "
                                     "In YAML: " + mName);
        }
    }
    
    // get value with numerical bounds
    template <typename T>
    T getWithBounds(const std::string &keyword, T lower,
                    T upper = std::numeric_limits<T>::max()) const {
        T res = get<T>(keyword);
        if (res < lower || res > upper) {
            const std::string &range = bstring::range(lower, upper);
            throw std::runtime_error("InparamYAML::getWithBounds || "
                                     + keyword + " ∉ " + range + " || "
                                     "Keyword: " + keyword + " || "
                                     "In YAML: " + mName);
        }
        return res;
    }
    
    // get array by keyword string
    template <typename T>
    std::vector<T> getVector(const std::string &keyword) const {
        // get size
        std::vector<std::string> kwVec = bstring::split(keyword, ":");
        kwVec.push_back("[?]");
        int size = get<int>(kwVec);
        
        // read elements one by one
        std::vector<T> result;
        result.reserve(size);
        for (int index = 0; index < size; index++) {
            kwVec.back() = "[" + bstring::toString(index) + "]";
            result.push_back(get<T>(kwVec));
        }
        return result;
    }
    
    // check existence
    bool contains(const std::string &keyword) const {
        try {
            get<std::string>(keyword);
            return true;
        } catch (...) {
            return false;
        }
    }
    
private:
    const std::string mName;
    Yaml::Node mRoot;
};

#endif /* InparamYAML_hpp */
