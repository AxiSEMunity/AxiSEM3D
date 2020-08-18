//
//  mpi.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/6/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  MPI interfaces

#include "mpi.hpp"
// verbose
#include "bstring.hpp"
// io::cout
#include "io.hpp"
// grouping
#include "vector_tools.hpp"

namespace mpi {
    ////////////////////////////// internal //////////////////////////////
    namespace internal {
#ifndef _SERIAL_BUILD
        // group superior
        MPI_Comm iCommSuper = MPI_COMM_NULL;
        // group inferior
        MPI_Comm iCommInfer = MPI_COMM_NULL;
        // current MPI_Comm
        MPI_Comm iCommCurrent = MPI_COMM_NULL;
#endif
        // # proc per group
        int iNumProcPerGroup = 1;
        
        // min/max # proc on node
        int iMinNumProcOnNode = -1;
        int iMaxNumProcOnNode = -1;
        
        // number of groups
        int iNumGroups = -1;
        
        // am I a super proc in my group
        bool iSuper = true;
    }
    
    
    ////////////////////////////// basics //////////////////////////////
    // initialize mpi
    void initialize(int *argc, char ***argv) {
#ifndef _SERIAL_BUILD
        MPI_Init(argc, argv);
#endif
        // init current comm
        enterWorld();
        // by default, one proc per group (no grouping)
        setupGroup(1);
        // init io::cout
        io::cout.setMyWorldRank(rank());
    }
    
    // finalize mpi
    void finalize() {
        freeGroupComm();
#ifndef _SERIAL_BUILD
        MPI_Finalize();
#endif
    }
    
    // initialize group
    void setupGroup(int nprocPerGroup) {
        // free if exists
        freeGroupComm();
        
        // limit # proc per group
        nprocPerGroup = std::min(nprocPerGroup, nproc());
        nprocPerGroup = std::max(nprocPerGroup, 1);
        
        // get node id
        char charNID[MPI_MAX_PROCESSOR_NAME] = "Node";
        int lenNID = 4;
#ifndef _SERIAL_BUILD
        MPI_Get_processor_name(charNID, &lenNID);
#endif
        std::string nodeID(charNID);
        
        // assemble node id
        std::vector<std::string> nodeIDs;
        gather(nodeID, nodeIDs, MPI_ALL);
        
        // find nproc in each id and my index in my id
        int myIndexInID = 0;
        std::map<std::string, int> nodeID_NP;
        for (int irank = 0; irank < nproc(); irank++) {
            vector_tools::aggregate(nodeID_NP, nodeIDs[irank], 1);
            if (irank < rank() && nodeIDs[irank] == nodeID) {
                myIndexInID += 1;
            }
        }
        
        // find starting group index for each id
        int curGroupIndex = 0;
        internal::iMinNumProcOnNode = nproc() + 1;
        internal::iMaxNumProcOnNode = 0;
        std::map<std::string, int> nodeID_StartG;
        for (auto it = nodeID_NP.begin(); it != nodeID_NP.end(); ++it) {
            nodeID_StartG.insert({it->first, curGroupIndex});
            curGroupIndex += it->second / nprocPerGroup;
            if (it->second % nprocPerGroup > 0) {
                curGroupIndex += 1;
            }
            // find min/max number of procs on node
            internal::iMinNumProcOnNode =
            std::min(internal::iMinNumProcOnNode, it->second);
            internal::iMaxNumProcOnNode =
            std::max(internal::iMaxNumProcOnNode, it->second);
        }
        internal::iNumProcPerGroup = std::min(internal::iMaxNumProcOnNode,
                                              nprocPerGroup);
        
        // total number of groups
        internal::iNumGroups = curGroupIndex;
        
        // my group index
        int myGroupIndex = (myIndexInID / internal::iNumProcPerGroup +
                            nodeID_StartG.at(nodeID));
        int myIndexInGroup = myIndexInID % internal::iNumProcPerGroup;
        
        // super
        internal::iSuper = (myIndexInGroup == 0);
        
        // split MPI comm
#ifndef _SERIAL_BUILD
        // super: discontinuous in rank
        MPI_Comm_split(MPI_COMM_WORLD, myIndexInGroup, rank(),
                       &internal::iCommSuper);
        // infer: continuous in rank
        MPI_Comm_split(MPI_COMM_WORLD, myGroupIndex, rank(),
                       &internal::iCommInfer);
#endif
    }
    
    // finalize group
    void freeGroupComm() {
        // first go back to world
        enterWorld();
#ifndef _SERIAL_BUILD
        if (internal::iCommSuper != MPI_COMM_NULL) {
            MPI_Comm_free(&internal::iCommSuper);
        }
        if (internal::iCommInfer != MPI_COMM_NULL) {
            MPI_Comm_free(&internal::iCommInfer);
        }
#endif
    }
    
    // barrier
    void barrier() {
#ifndef _SERIAL_BUILD
        MPI_Barrier(internal::iCommCurrent);
#endif
    }
    
    // abort
    void abort(int err) {
#ifndef _SERIAL_BUILD
        MPI_Abort(MPI_COMM_WORLD, err);
#else
        exit(err);
#endif
    }
    
    
    ////////////////////////////// broadcast //////////////////////////////
    // specialization for string
    void bcast(std::string &str, int src) {
        // string to std::vector<char>
        std::vector<char> vchar;
        if (rank() == src) {
            vchar = std::vector<char>(str.begin(), str.end());
        }
        
        // bcast std::vector<char>
        bcast(vchar, src);
        
        // std::vector<char> to string
        if (rank() != src) {
            str = std::string(vchar.begin(), vchar.end());
        }
    }
    
    // specialization for vector of string
    void bcast(std::vector<std::string> &vecStr, int src) {
        // string to std::vector<char>
        std::vector<std::vector<char>> vvchar;
        vvchar.reserve(vecStr.size());
        if (rank() == src) {
            for (int istr = 0; istr < vecStr.size(); istr++) {
                vvchar.push_back(std::vector<char>(vecStr[istr].begin(),
                                                   vecStr[istr].end()));
            }
        }
        
        // bcast std::vector<char>
        bcast(vvchar, src);
        
        // std::vector<char> to string
        if (rank() != src) {
            vecStr.assign(vvchar.size(), "");
            for (int istr = 0; istr < vvchar.size(); istr++) {
                vecStr[istr] = std::string(vvchar[istr].begin(),
                                           vvchar[istr].end());
            }
        }
    }
    
    
    ////////////////////////////// isend / irecv //////////////////////////////
    // wait_all: must be implemented in .cpp
    void waitAll(std::vector<MPI_Request> &requests) {
#ifndef _SERIAL_BUILD
        MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
#endif
    }
    
    
    ////////////////////////////// gather //////////////////////////////
    // specialization for string
    void gather(const std::string &str,
                std::vector<std::string> &vecStr, int dest) {
        // string to std::vector<char>
        std::vector<char> vchar(str.begin(), str.end());
        
        // gather cstr
        std::vector<std::vector<char>> vvchar;
        gather(vchar, vvchar, dest);
        
        // std::vector<char> to string
        if (dest < 0 || dest == rank()) {
            vecStr.assign(nproc(), "");
            for (int iproc = 0; iproc < nproc(); iproc++) {
                vecStr[iproc] = std::string(vvchar[iproc].begin(),
                                            vvchar[iproc].end());
            }
        }
    }
    
    // specialization for vector of string
    void gather(const std::vector<std::string> &vecStr,
                std::vector<std::vector<std::string>> &vecVecStr, int dest) {
        // flattened
        std::string localFlattened = "";
        std::vector<int> localFlattenedLength(vecStr.size(), 0);
        for (int istr = 0; istr < vecStr.size(); istr++) {
            localFlattened += vecStr[istr];
            localFlattenedLength[istr] = (int)vecStr[istr].length();
        }
        
        // gather flattened
        std::vector<std::string> globalFlattened;
        std::vector<std::vector<int>> globalFlattenedLength;
        gather(localFlattened, globalFlattened, dest);
        gather(localFlattenedLength, globalFlattenedLength, dest);
        
        // cast back to vector
        if (dest < 0 || dest == rank()) {
            vecVecStr.clear();
            vecVecStr.reserve(nproc());
            for (int iproc = 0; iproc < nproc(); iproc++) {
                // each proc
                int nstr = (int)globalFlattenedLength[iproc].size();
                std::vector<std::string> vstr(nstr, "");
                int pos = 0;
                for (int jstr = 0; jstr < nstr; jstr++) {
                    // each string
                    int length = globalFlattenedLength[iproc][jstr];
                    vstr[jstr] = globalFlattened[iproc].substr(pos, length);
                    pos += length;
                }
                vecVecStr.push_back(vstr);
            }
        }
    }
    
    ////////////////////////////// external //////////////////////////////
    // verbose
    std::string verbose() {
        std::stringstream ss;
        ss << bstring::boxTitle("MPI");
#ifndef _SERIAL_BUILD
        ss << bstring::boxEquals(0, 24, "# processors", nproc());
        ss << bstring::boxEquals(0, 24, "# MPI groups", internal::iNumGroups);
        ss << bstring::boxEquals(0, 24, "# processors per group",
                                 internal::iNumProcPerGroup);
        ss << bstring::boxEquals(0, 24, "min # processors on node",
                                 internal::iMinNumProcOnNode);
        ss << bstring::boxEquals(0, 24, "max # processors on node",
                                 internal::iMaxNumProcOnNode);
#else
        ss << "  This is a serial build without MPI.\n";
#endif
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
}
