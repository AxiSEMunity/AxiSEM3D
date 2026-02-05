//
//  InparamYAML.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  YAML parser for input parameters
//  based on mini-yaml:
//  https://github.com/jimmiebergmann/mini-yaml

#include "InparamYAML.hpp"
#include "mpi.hpp"

// parse
void InparamYAML::parse(const std::string &fname) {
    // read input file
    std::vector<std::string> lines;
    if (mpi::root()) {
        lines = bstring::readLines(fname, "InparamYAML::parse");
    }
    mpi::bcast(lines);
    
    // remove comments
    for (std::string &line: lines) {
        line = line.substr(0, line.find('#'));
    }
    
    // change array [] into "-" because mini-yaml does not support []
    for (std::string &line: lines) {
        int front = (int)line.find('[');
        if (front == std::string::npos) {
            continue;
        }
        // find ']'
        int back = (int)line.find(']', front);
        if (back == std::string::npos) {
            throw std::runtime_error("InparamYAML::parse || "
                                     "Array brackets [] does not match. || "
                                     "Input yaml file: || " + fname);
        }
        // split
        const std::string &str = line.substr(front + 1, back - front - 1);
        const std::vector<std::string> &vec = bstring::split(str, ",");
        if (vec.front() == "") {
            // empty sequence
            line = line.substr(0, front) + "YAML_EMPTY_SEQUENCE";
        } else {
            // write header in this line
            std::string converted = line.substr(0, front) + "\n";
            // write elements with -
            for (const std::string &elem: vec) {
                if (elem.find(":") != std::string::npos) {
                    throw std::runtime_error("InparamYAML::parse || "
                                             "Array given in [] must contain "
                                             "non-object elements. || "
                                             "Input yaml file: || " + fname);
                }
                // use location of '[' as indent
                converted += bstring::filled(front) + "- " + elem + "\n";
            }
            line = converted;
        }
    }
    
    // merge lines
    std::string strArr;
    for (const std::string &line: lines) {
        strArr += line + "\n";
    }
    
    // parse input
    try {
        Yaml::Parse(mRoot, strArr);
    } catch (...) {
        throw std::runtime_error("InparamYAML::parse || "
                                 "Error parsing YAML by mini-yaml || "
                                 "Input yaml file: || " + fname);
    }
}

// verbose
std::string InparamYAML::verbose() const {
    // get contents
    std::string contents;
    Yaml::Serialize(mRoot, contents);
    
    // add indent
    contents = "  " + bstring::replace(contents, "\n", "\n  ");
    contents = contents.substr(0, contents.size() - 2);
    
    // add subtitle
    std::stringstream ss;
    ss << bstring::boxSubTitle(0, mName);
    ss << contents;
    return ss.str();
}
