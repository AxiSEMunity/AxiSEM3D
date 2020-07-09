//
//  metis.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/17/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  dual graph partioning by Metis
//  compatible with both 32-bit and 64-bit builds of Metis

#include "metis.hpp"
#include "vector_tools.hpp"
#include <metis.h>

namespace eigen {
    // idx_t
    typedef Eigen::Matrix<idx_t, Eigen::Dynamic, 4, Eigen::RowMajor> IIMatX4_RM;
    typedef Eigen::Matrix<idx_t, Eigen::Dynamic, 1> IIColX;
}

namespace metis {
    /////////////////////// internal functions ///////////////////////
    // error handler
    void error(int retval, const std::string &func_name) {
        if (retval != METIS_OK) {
            throw std::runtime_error
            ("metis::error || Error in Metis function: " + func_name);
        }
    }
    
    // form adjacency
    void formAdjacency(const eigen::IIMatX4_RM &connectivity,
                       idx_t ncommon, idx_t *&xadj, idx_t *&adjncy) {
        // get unique node list
        idx_t nelem = (idx_t)connectivity.rows();
        std::vector<idx_t> nodes(connectivity.data(),
                                 connectivity.data() + nelem * 4);
        vector_tools::sortUnique(nodes);
        idx_t nnode = (idx_t)nodes.size();
        // form an inverse map from global to local for fast search
        std::map<idx_t, idx_t> nG2L;
        for (idx_t inl = 0; inl < nnode; inl++) {
            nG2L.insert({nodes[inl], inl});
        }
        
        // convert connectivity to CSR format
        eigen::IIColX eptr =
        eigen::IIColX::LinSpaced(nelem + 1, 0, nelem * 4);
        eigen::IIColX eind = eigen::IIColX(nelem * 4);
        for (idx_t ie = 0; ie < nelem; ie++) {
            for (idx_t ivt = 0; ivt < 4; ivt++) {
                eind(ie * 4 + ivt) = nG2L.at(connectivity(ie, ivt));
            }
        }
        
        // convert to dual graph
        idx_t numflag = 0;
        error(METIS_MeshToDual(&nelem, &nnode, eptr.data(), eind.data(),
                               &ncommon, &numflag, &xadj, &adjncy),
              "METIS_MeshToDual");
    }
    
    // free adjacency
    void freeAdjacency(idx_t *&xadj, idx_t *&adjncy) {
        // free CRS created by metis
        error(METIS_Free(xadj), "METIS_Free");
        error(METIS_Free(adjncy), "METIS_Free");
    }
    
    
    /////////////////////// external functions ///////////////////////
    // form neighbourhood of connectivity
    void formNeighbourhood(const eigen::IMatX4_RM &connectivity, int ncommon,
                           std::vector<eigen::IColX> &neighbours) {
        // form graph
        idx_t *xadj, *adjncy;
        formAdjacency(connectivity.cast<idx_t>(), (idx_t)ncommon, xadj, adjncy);
        
        // convert CRS format to neighbours
        neighbours.clear();
        neighbours.reserve(connectivity.rows());
        for (int ie = 0; ie < connectivity.rows(); ie++) {
            int nNeighb = (int)(xadj[ie + 1] - xadj[ie]);
            eigen::IColX neighb(nNeighb);
            for (int in = 0; in < nNeighb; in++) {
                neighb(in) = (int)adjncy[xadj[ie] + in];
            }
            neighbours.push_back(neighb);
        }
        
        // delete graph
        freeAdjacency(xadj, adjncy);
    }
    
    // domain decomposition
    double decompose(const eigen::IMatX4_RM &connectivity,
                     const eigen::DColX &weights, int npart, int rseed,
                     eigen::IColX &elemRank) {
        // this cause error in Metis
        if (npart == 1) {
            elemRank = eigen::IColX::Zero(connectivity.rows());
            return 0.;
        }
        
        //////////// input ////////////
        // trival
        idx_t nelem = (idx_t)connectivity.rows();
        idx_t ncon = 1;
        
        // form graph
        idx_t *xadj, *adjncy;
        formAdjacency(connectivity.cast<idx_t>(), 2, xadj, adjncy);
        
        // weights
        idx_t *vwgt = NULL;
        if (weights.rows() == nelem) {
            // first normalize then magnify to avoid over/underflow
            double imax = std::numeric_limits<idx_t>::max() * .9;
            eigen::IIColX iweights = (weights / weights.sum() *
                                      imax).array().round().cast<idx_t>();
            vwgt = iweights.data();
        }
        
        // trival
        idx_t npart_idx_t = (idx_t)npart;
        real_t ubvec = (real_t)1.001;
        
        // metis options
        idx_t metisOps[METIS_NOPTIONS];
        error(METIS_SetDefaultOptions(metisOps), "METIS_SetDefaultOptions");
        metisOps[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
        metisOps[METIS_OPTION_CONTIG] = 1;
        metisOps[METIS_OPTION_NCUTS] = 1;
        metisOps[METIS_OPTION_SEED] = rseed;
        
        //////////// output ////////////
        idx_t objval = 0;
        eigen::IIColX erank(nelem);
        
        //////////// run ////////////
        error(METIS_PartGraphKway(&nelem, &ncon, xadj, adjncy, vwgt, NULL,
                                  NULL, &npart_idx_t, NULL, &ubvec, metisOps,
                                  &objval, erank.data()),
              "METIS_PartGraphKway");
        
        // free graph
        freeAdjacency(xadj, adjncy);
        
        // cast result
        elemRank = erank.cast<int>();
        return objval * 1.;
    }
}
