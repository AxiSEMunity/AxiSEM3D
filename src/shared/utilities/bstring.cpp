//
//  bstring.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  string tools

#include "bstring.hpp"
#include <fstream>

extern "C" {
#include <time.h>
};

namespace bstring {
    ////////////////////////// string tools //////////////////////////
    // replace
    std::string replace(const std::string &in,
                        const std::string &from, const std::string &to) {
        std::string out(in);
        boost::replace_all(out, from, to);
        return out;
    }
    
    // split
    std::vector<std::string> split(const std::string &in,
                                   const std::string &sep) {
        const std::string &orig = boost::trim_copy(in);
        std::vector<std::string> substrs;
        boost::split(substrs, orig, boost::is_any_of(sep),
                     boost::token_compress_on);
        for (std::string &substr: substrs) {
            substr = boost::trim_copy(substr);
        }
        return substrs;
    }
    
    // read lines from a file
    std::vector<std::string> readLines(const std::string &fname,
                                       const std::string &caller) {
        std::ifstream fs(fname);
        if (!fs) {
            throw std::runtime_error("bstring::readLines || "
                                     "Error opening input file: || "
                                     + fname + "|| Called by " + caller);
        }
        std::vector<std::string> lines;
        std::string line = "";
        while (getline(fs, line)) {
            lines.push_back(line);
        }
        fs.close();
        return lines;
    }
    
    // read all from a file
    std::string readAll(const std::string &fname,
                        const std::string &caller) {
        const std::vector<std::string> &lines = readLines(fname, caller);
        std::stringstream ss;
        for (const std::string &line: lines) {
            ss << line << "\n";
        }
        return ss.str();
    }
    
    // get current datetime, format YYYY-MM-DDTHH:mm:ss
    std::string currentDateTime() {
        time_t now = time(0);
        struct tm tstruct = *localtime(&now);
        char buf[32] = "";
        strftime(buf, sizeof(buf), "%Y-%m-%dT%X", &tstruct);
        return std::string(buf);
    }
    
    
    ////////////////////////// verbose boxes //////////////////////////
    // filled string
    std::string filled(int length, char fill) {
        std::stringstream ss;
        ss.width(length);
        ss.fill(fill);
        ss << "";
        return ss.str();
    }
    
    // ****** title ******
    std::string boxTitle(const std::string &title, char fill, int length) {
        int tsize = (int)title.size() + 2;
        length = std::max(length, tsize + 2); // at least one *
        int nright = (length - tsize) / 2;
        int nleft = length - tsize - nright;
        std::stringstream ss;
        ss << filled(nleft, fill) << " " + title + " ";
        ss << filled(nright, fill) << "\n";
        return ss.str();
    }
    
    // subtitle___________
    std::string boxSubTitle(int indent, const std::string &subTitle,
                            char fill, int length) {
        std::stringstream ss;
        int nfill = length - indent - (int)subTitle.size();
        ss << filled(indent) << subTitle << filled(nfill, fill) << "\n";
        return ss.str();
    }
    
    // *******************
    std::string boxBaseline(char fill, int length) {
        return filled(length, fill) + "\n";
    }
    
    // const char *
    std::string boxEquals(int indent, int keyWidth,
                          const std::string &key, const char *value,
                          const std::string &equalSign,
                          int nspaces, int precision) {
        return boxEquals(indent, keyWidth, key, std::string(value),
                         equalSign, nspaces, precision);
    }
    
    
    ////////////////////////// error //////////////////////////
    // format warning
    std::string warning(const std::string &what,
                        const std::string &head, char fill) {
        // split
        std::vector<std::string> strs = split(what, "|");
        // unknown source
        if (strs.size() <= 1) {
            strs.insert(strs.begin(), "Unknown");
        }
        // source
        std::string src = "FROM: " + strs[0] + "\n";
        // message
        std::string msg = "";
        for (int istr = 1; istr < strs.size(); istr++) {
            msg += (istr == 1) ? "WHAT: " : "      ";
            msg += strs[istr] + "\n";
        }
        // time
        std::string time = "TIME: " + currentDateTime() + "\n";
        // box
        std::stringstream ss;
        ss << boxTitle(head, fill);
        ss << src << msg << time;
        ss << boxBaseline(fill) << "\n\n";
        return ss.str();
    }
    
    // format exception
    std::string exception(const std::exception &e) {
        return warning(e.what(), "AxiSEM3D ABORTED UPON RUNTIME EXCEPTION",
                       fmt::gExceptionFill);
    }
}
