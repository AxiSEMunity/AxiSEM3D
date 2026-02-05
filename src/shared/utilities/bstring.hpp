//
//  bstring.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  string tools

#ifndef bstring_hpp
#define bstring_hpp

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/type_index.hpp>
#pragma clang diagnostic pop

#include <map>

namespace bstring {
    // verbose format constants
    namespace fmt {
        const int gFloatPrecision = 6;
        const int gBoxWidth = 64;
        const int gBoxNSpaceEquals = 2;
        const char gBoxTitleFill = '=';
        const char gBoxSubTitleFill = '_';
        const char gWarningFill = '*';
        const char gExceptionFill = '!';
    }
    
    
    ////////////////////////// cast //////////////////////////
    // pretty type name from type
    template <typename parType>
    std::string typeName() {
        return boost::typeindex::type_id<parType>().pretty_name();
    }
    
    // pretty type name from object
    template <typename parType>
    std::string typeName(const parType &obj) {
        return boost::typeindex::type_id_runtime(obj).pretty_name();
    }
    
    // cast string to parType
    template <typename parType>
    parType cast(const std::string &in_str, const std::string &caller) {
        try {
            // specialization for bool
            std::string val = in_str;
            if (typeid(parType) == typeid(bool)) {
                boost::to_upper<std::string>(val);
                if (val == "TRUE" || val == "YES" || val == "ON") {
                    val = "1";
                }
                if (val == "FALSE" || val == "NO" || val == "OFF") {
                    val = "0";
                }
            }
            return boost::lexical_cast<parType>(val);
        } catch (...) {
            throw std::runtime_error("bstring::cast || Error casting string \""
                                     + in_str + "\" to type \""
                                     + typeName<parType>() + "\" || "
                                     "Called by " + caller);
        }
    }
    
    // print parType to string
    template <typename parType>
    std::string toString(const parType &par,
                         int precision = fmt::gFloatPrecision) {
        // use stream
        std::stringstream ss;
        ss.precision(precision);
        ss << par;
        // bool
        if (typeid(parType) == typeid(bool)) {
            return (ss.str() == "1") ? "true" : "false";
        }
        return ss.str();
    }
    
    // cast range to string
    template <typename T>
    std::string range(T min, T max, char fbrac = '[', char bbrac=']',
                      int precision = fmt::gFloatPrecision) {
        std::stringstream ss;
        ss.precision(precision);
        ss << fbrac << min << ", " << max << bbrac;
        return ss.str();
    }
    
    
    ////////////////////////// string tools //////////////////////////
    // replace
    std::string replace(const std::string &in,
                        const std::string &from, const std::string &to);
    
    // split
    std::vector<std::string> split(const std::string &in,
                                   const std::string &sep);
    
    // read lines from a file
    std::vector<std::string> readLines(const std::string &fname,
                                       const std::string &caller);
    
    // read all from a file
    std::string readAll(const std::string &fname,
                        const std::string &caller);
    
    // get current datetime, format YYYY-MM-DDTHH:mm:ss
    std::string currentDateTime();
    
    
    ////////////////////////// verbose boxes //////////////////////////
    // filled string
    std::string filled(int length, char fill = ' ');
    
    // ****** title ******
    std::string boxTitle(const std::string &title,
                         char fill = fmt::gBoxTitleFill,
                         int length = fmt::gBoxWidth);
    // subtitle___________
    std::string boxSubTitle(int indent, const std::string &subTitle,
                            char fill = fmt::gBoxSubTitleFill,
                            int length = fmt::gBoxWidth);
    // *******************
    std::string boxBaseline(char fill = fmt::gBoxTitleFill,
                            int length = fmt::gBoxWidth);
    
    // key = value1
    //       value2
    //       ...
    template <typename T>
    std::string boxEquals(int indent, int keyWidth,
                          const std::string &key, const std::vector<T> &values,
                          const std::string &equalSign = "=",
                          bool oneLine = false,
                          int nspaces = fmt::gBoxNSpaceEquals,
                          int precision = fmt::gFloatPrecision) {
        std::stringstream ss;
        ss.precision(precision);
        ss << filled(indent);
        ss.width(keyWidth);
        ss << std::left << key;
        ss << filled(nspaces) << equalSign << filled(nspaces);
        if (oneLine) {
            if (values.size() == 0) {
                ss << "[]\n";
            } else {
                std::stringstream ssv;
                for (T val: values) {
                    ssv << val << ", ";
                }
                ss << "[" << ssv.str().substr(0, ssv.str().size() - 2) << "]\n";
            }
        } else {
            int totalIndent = 0;
            for (T val: values) {
                ss << filled(totalIndent) << toString(val) << "\n";
                totalIndent = indent + keyWidth + nspaces * 2 + 1;
            }
            if (totalIndent == 0) {
                ss << "\n";
            }
        }
        return ss.str();
    }
    
    // key = value
    template <typename T>
    std::string boxEquals(int indent, int keyWidth,
                          const std::string &key, const T &value,
                          const std::string &equalSign = "=",
                          int nspaces = fmt::gBoxNSpaceEquals,
                          int precision = fmt::gFloatPrecision) {
        return boxEquals(indent, keyWidth, key, std::vector<T>({value}),
                         equalSign, false, nspaces, precision);
    }
    
    // const char *
    std::string boxEquals(int indent, int keyWidth,
                          const std::string &key, const char *value,
                          const std::string &equalSign = "=",
                          int nspaces = fmt::gBoxNSpaceEquals,
                          int precision = fmt::gFloatPrecision);
    
    // std::map
    template <typename T>
    std::string boxEquals(int indent, int keyWidth,
                          const std::map<std::string, T> &keyValue,
                          const std::string &equalSign = "=",
                          int nspaces = fmt::gBoxNSpaceEquals,
                          int precision = fmt::gFloatPrecision) {
        std::stringstream ss;
        for (auto it = keyValue.begin(); it != keyValue.end(); ++it) {
            ss << boxEquals(indent, keyWidth, it->first, it->second,
                            equalSign, nspaces, precision);
        }
        return ss.str();
    }
    
    ////////////////////////// error //////////////////////////
    // format warning
    std::string warning(const std::string &what,
                        const std::string &head = "WARNING FROM AxiSEM3D",
                        char fill = fmt::gWarningFill);
    
    // format exception
    std::string exception(const std::exception &e);
}

#endif /* bstring_hpp */
