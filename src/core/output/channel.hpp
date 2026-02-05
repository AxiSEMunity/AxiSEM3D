//
//  channel.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  channels for wavefield output

#ifndef channel_hpp
#define channel_hpp

#include <vector>
#include <map>
#include <string>

namespace channel {
    // coordinate system for wavefield output
    // spz: s, phi, z (solver intrinsic)
    // RTZ: radial, transverse, vertical
    // ENZ: east, north, vertical
    // xyz: x, y, z in source-centered
    enum class WavefieldCS {spz, RTZ, ENZ, xyz};
    const std::map<WavefieldCS, std::string> WavefieldCS_Str = {
        {WavefieldCS::spz, "spz"},
        {WavefieldCS::RTZ, "RTZ"},
        {WavefieldCS::ENZ, "ENZ"},
        {WavefieldCS::xyz, "xyz"}};
    
    // solid
    namespace solid {
        // allowed field types
        enum class FieldType {Displ, Nabla, Strain, Curl, Stress};
        
        // allowed channels
        const std::map<int,
        std::tuple<std::string, std::vector<std::string>, FieldType, int>>
        gChannelMap = {
            // displ
            {1, {"U1", {"U1", "U"}, FieldType::Displ, 0}},
            {2, {"U2", {"U2", "U"}, FieldType::Displ, 1}},
            {3, {"U3", {"U3", "U"}, FieldType::Displ, 2}},
            {4, {"U_NORM", {"U_NORM", "|U|"}, FieldType::Displ, -1}},
            // nabla
            {5, {"G11", {"G11", "G"}, FieldType::Nabla, 0}},
            {6, {"G12", {"G12", "G"}, FieldType::Nabla, 1}},
            {7, {"G13", {"G13", "G"}, FieldType::Nabla, 2}},
            {8, {"G21", {"G21", "G"}, FieldType::Nabla, 3}},
            {9, {"G22", {"G22", "G"}, FieldType::Nabla, 4}},
            {10, {"G23", {"G23", "G"}, FieldType::Nabla, 5}},
            {11, {"G31", {"G31", "G"}, FieldType::Nabla, 6}},
            {12, {"G32", {"G32", "G"}, FieldType::Nabla, 7}},
            {13, {"G33", {"G33", "G"}, FieldType::Nabla, 8}},
            {14, {"G_I1", {"G_I1", "Gii"}, FieldType::Nabla, -1}},
            // strain
            {15, {"E11", {"E11", "E"}, FieldType::Strain, 0}},
            {16, {"E22", {"E22", "E"}, FieldType::Strain, 1}},
            {17, {"E33", {"E33", "E"}, FieldType::Strain, 2}},
            {18, {"E23", {"E23", "E", "E32"}, FieldType::Strain, 3}},
            {19, {"E31", {"E31", "E", "E13"}, FieldType::Strain, 4}},
            {20, {"E12", {"E12", "E", "E21"}, FieldType::Strain, 5}},
            {21, {"E_I1", {"E_I1", "Eii"}, FieldType::Strain, -1}},
            {22, {"E_J2", {"E_J2"}, FieldType::Strain, -2}},
            // curl
            {23, {"R1", {"R1", "R"}, FieldType::Curl, 0}},
            {24, {"R2", {"R2", "R"}, FieldType::Curl, 1}},
            {25, {"R3", {"R3", "R"}, FieldType::Curl, 2}},
            {26, {"R_NORM", {"R_NORM", "|R|"}, FieldType::Curl, -1}},
            // stress
            {27, {"S11", {"S11", "S"}, FieldType::Stress, 0}},
            {28, {"S22", {"S22", "S"}, FieldType::Stress, 1}},
            {29, {"S33", {"S33", "S"}, FieldType::Stress, 2}},
            {30, {"S23", {"S23", "S", "S32"}, FieldType::Stress, 3}},
            {31, {"S31", {"S31", "S", "S13"}, FieldType::Stress, 4}},
            {32, {"S12", {"S12", "S", "S21"}, FieldType::Stress, 5}},
            {33, {"S_I1", {"S_I1", "Sii"}, FieldType::Stress, -1}},
            {34, {"S_J2", {"S_J2"}, FieldType::Stress, -2}},
        };
        
        // channel options
        struct ChannelOptions {
            // constructor
            ChannelOptions(WavefieldCS wcs,
                           const std::vector<std::string> &userChannels);
            
            // wavefield cs
            const WavefieldCS mWCS;
            
            // verified channels
            std::vector<int> mStdChannels;
            
            // buffer
            bool mNeedBufferU = false;
            bool mNeedBufferG = false;
            bool mNeedBufferE = false;
            bool mNeedBufferR = false;
            bool mNeedBufferS = false;
        };
    }
    
    // fluid
    namespace fluid {
        // allowed field types
        enum class FieldType {Chi, Displ, Pressure, Delta};
        
        // allowed channels
        const std::map<int,
        std::tuple<std::string, std::vector<std::string>, FieldType, int>>
        gChannelMap = {
            // chi
            {1, {"X", {"X"}, FieldType::Chi, 0}},
            // displ
            {2, {"U1", {"U1", "U"}, FieldType::Displ, 0}},
            {3, {"U2", {"U2", "U"}, FieldType::Displ, 1}},
            {4, {"U3", {"U3", "U"}, FieldType::Displ, 2}},
            {5, {"U_NORM", {"U_NORM", "|U|"}, FieldType::Displ, -1}},
            // pressure
            {6, {"P", {"P"}, FieldType::Pressure, 0}},
            // delta
            // Delta output is disabled below because the stiffness term
            // on a fluid point does NOT equal to Delta but Delta * Jacobian
            // {7, {"D", {"D"}, FieldType::Delta, 0}},
        };
        
        // channel options
        struct ChannelOptions {
            // constructor
            ChannelOptions(WavefieldCS wcs,
                           const std::vector<std::string> &userChannels);
            
            // wavefield cs
            const WavefieldCS mWCS;
            
            // verified channels
            std::vector<int> mStdChannels;
            
            // buffer
            bool mNeedBufferX = false;
            bool mNeedBufferU = false;
            bool mNeedBufferP = false;
            bool mNeedBufferD = false;
        };
    }
}

#endif /* channel_hpp */
