//
//  vicinity.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/19/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element-GLL vicinity

#define next4(n) ((n == 3) ? 0 : (n + 1))
#define prev4(n) ((n == 0) ? 3 : (n - 1))

#include "vicinity.hpp"
#include "metis.hpp"
#include <map>

namespace vicinity {
    /////////////////////// internal functions ///////////////////////
    void sharedDOF(const eigen::IMatX4_RM &connectivity, int elemA, int elemB,
                   int &nodeOrEdgeA, int &nodeOrEdgeB, bool &sharingEdge) {
        // get node lists
        const auto &nodeListA = connectivity.row(elemA);
        const auto &nodeListB = connectivity.row(elemB);
        
        // find first common node
        bool found = false;
        int nA = -1, nB = -1;
        for (nA = 0; nA < 4; nA++) {
            for (nB = 0; nB < 4; nB++) {
                if (nodeListA(nA) == nodeListB(nB)) {
                    found = true;
                    break;
                }
            }
            if (found) {
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("vicinity:: sharedDOF || Impossible.");
        }
        
        // check node or edge
        if (nodeListA[next4(nA)] == nodeListB[prev4(nB)]) {
            sharingEdge = true;
            nodeOrEdgeA = nA;
            nodeOrEdgeB = prev4(nB);
        } else if (nodeListA[prev4(nA)] == nodeListB[next4(nB)]) {
            sharingEdge = true;
            nodeOrEdgeA = prev4(nA);
            nodeOrEdgeB = nB;
        } else {
            sharingEdge = false;
            nodeOrEdgeA = nA;
            nodeOrEdgeB = nB;
        }
    }
    
    
    /////////////////////// external functions ///////////////////////
    // connectivity to element-GLL vicinity
    int connectivityToElementGLL(const eigen::IMatX4_RM &connectivity,
                                 eigen::IMatXN_RM &elementGLL,
                                 std::vector<eigen::IColX> &neighbours) {
        using namespace constants;
        
        // get neighbours from metis
        metis::formNeighbourhood(connectivity, 1, neighbours);
        
        // init element-GLL vicinity
        int nelem = (int)connectivity.rows();
        elementGLL = eigen::IMatXN_RM::Constant(nelem, spectral::nPEM, -1);
        
        // GLL tag start from 0
        int nGLL = 0;
        
        // loop over all elements
        for (int ielem = 0; ielem < nelem; ielem++) {
            // True for points not shared with any other
            // neighbouring element with a lower tag
            std::array<bool, spectral::nPEM> unshared;
            unshared.fill(true);
            
            // loop over all neighbouring elements with a smaller tag
            for (int inb = 0; inb < neighbours[ielem].rows(); inb++) {
                int inghb = neighbours[ielem](inb);
                if (inghb >= ielem) {
                    continue;
                }
                
                // shared node or edge
                int nodeOrEdgeMe = -1, nodeOrEdgeNb = -1;
                bool sharingEdge = false;
                sharedDOF(connectivity, ielem, inghb,
                          nodeOrEdgeMe, nodeOrEdgeNb, sharingEdge);
                
                // assign shared and mask unshared
                if (sharingEdge) {
                    for (int ip = 0; ip < spectral::nPED; ip++) {
                        int iq = spectral::nPol - ip; // reversed on other
                        elementGLL(ielem, gEdgeIPnt[nodeOrEdgeMe][ip]) =
                        elementGLL(inghb, gEdgeIPnt[nodeOrEdgeNb][iq]);
                        unshared[gEdgeIPnt[nodeOrEdgeMe][ip]] = false;
                    }
                } else {
                    elementGLL(ielem, gNodeIPnt[nodeOrEdgeMe]) =
                    elementGLL(inghb, gNodeIPnt[nodeOrEdgeNb]);
                    unshared[gNodeIPnt[nodeOrEdgeMe]] = false;
                }
            }
            
            // fill unshared
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                if (unshared[ipnt]) {
                    elementGLL(ielem, ipnt) = nGLL++;
                }
            }
        }
        return nGLL;
    }
    
    // build local to global GLL
    eigen::IColX buildL2G_GLL(int nLocalGLL,
                              const eigen::IColX &myL2G_Element,
                              const eigen::IMatXN_RM &myElementGLL,
                              const eigen::IMatXN_RM &globalElementGLL) {
        // allocate and fill with -1
        eigen::IColX myL2G_GLL = eigen::IColX::Constant(nLocalGLL, -1);
        
        // match local and global GLL tags by (element, ipnt)
        for (int letag = 0; letag < myL2G_Element.size(); letag++) {
            int getag = myL2G_Element[letag];
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                int lgtag = myElementGLL(letag, ipnt);
                myL2G_GLL[lgtag] = globalElementGLL(getag, ipnt);
            }
        }
        
        // return
        return myL2G_GLL;
    }
    
    // build mpi communication
    // std::vector as the output to use mpi::scatter
    // rank0 gll0_0 gll0_1 ... -1 rank1 gll1_0 gll1_1 ... -1 rank2 ...
    std::vector<int>
    buildCommMPI(const eigen::IColX &myL2G_Element,
                 const std::vector<eigen::IColX> &globalNeighbours,
                 const eigen::IColX &elemRank, int myRank,
                 const eigen::IMatX4_RM &globalConnectivity,
                 const eigen::IMatXN_RM &myElementGLL,
                 const eigen::IColX &myL2G_GLL) {
        using namespace constants;
        // structured output: std::map<otherRank, std::map<globalGLL, localGLL>>
        // A: Why each element is a map with key = globalGLL,
        //    as we only need localGLL for the final output?
        // Q: To have the elements sorted by globalGLL so that localGLL
        //    on the two communicating ranks have the same order.
        std::map<int, std::map<int, int>> mapGLLComm;
        for (int letag = 0; letag < myL2G_Element.size(); letag++) {
            // loop in global neighbours
            int getag = myL2G_Element[letag];
            for (int netag: globalNeighbours[getag]) {
                if (elemRank[netag] == myRank) {
                    // neighbour on the same rank
                    continue;
                }
                
                // insert an empty item with key = the other rank
                // nothing happens if exists
                mapGLLComm.insert({elemRank[netag], std::map<int, int>()});
                std::map<int, int> &g2lMap = mapGLLComm.at(elemRank[netag]);
                
                // find shared node or edge
                int nodeOrEdgeMe = -1, nodeOrEdgeNb = -1;
                bool sharingEdge = false;
                sharedDOF(globalConnectivity, getag, netag,
                          nodeOrEdgeMe, nodeOrEdgeNb, sharingEdge);
                
                // find shared GLL
                if (sharingEdge) {
                    for (int ip = 0; ip < spectral::nPED; ip++) {
                        int ipnt = gEdgeIPnt[nodeOrEdgeMe][ip];
                        int lGLL = myElementGLL(letag, ipnt);
                        g2lMap.insert({myL2G_GLL[lGLL], lGLL});
                    }
                } else {
                    int lGLL = myElementGLL(letag, gNodeIPnt[nodeOrEdgeMe]);
                    g2lMap.insert({myL2G_GLL[lGLL], lGLL});
                }
            }
        }
        
        // cast to a flattened vector for mpi::scatter
        std::vector<int> vecCommProcGLL;
        for (auto itr = mapGLLComm.begin(); itr != mapGLLComm.end(); itr++) {
            // itr->first is the other rank
            vecCommProcGLL.push_back(itr->first);
            const auto &g2lMap = itr->second;
            for (auto itg = g2lMap.begin(); itg != g2lMap.end(); itg++) {
                // itg->first is the global GLL tag only for sorting
                vecCommProcGLL.push_back(itg->second);
            }
            // separating ranks by -1
            vecCommProcGLL.push_back(-1);
        }
        return vecCommProcGLL;
    }
}
