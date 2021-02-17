//
//  mpi.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/6/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  MPI interfaces

#ifndef mpi_hpp
#define mpi_hpp

#include <vector>
#include <map>
#include "eigen.hpp"

#ifndef _SERIAL_BUILD
// mpi.h has extern "C" {} internally
#include <mpi.h>
#else
#define MPI_Request int
#define MPI_REQUEST_NULL 0
#define MPI_Op int
#define MPI_MIN 0
#define MPI_MAX 1
#define MPI_SUM 2
#define MPI_MAX_PROCESSOR_NAME 32
#endif

// flag for gathering operations
#define MPI_ALL -1


namespace mpi {
    ////////////////////////////// internal //////////////////////////////
    namespace internal {
#ifndef _SERIAL_BUILD
        // group superior
        extern MPI_Comm iCommSuper;
        // group inferior
        extern MPI_Comm iCommInfer;
        // current MPI_Comm
        extern MPI_Comm iCommCurrent;
#endif
        // # proc per group
        extern int iNumProcPerGroup;
        
        // number of groups
        extern int iNumGroups;
        
        // am I a super proc in my group
        extern bool iSuper;
        
        //////// data type ////////
#ifndef _SERIAL_BUILD
        // datatype covertion, e.g., int -> MPI_INT
        template <typename T>
        inline MPI_Datatype typeMPI() = delete;
        
        template <>
        inline MPI_Datatype typeMPI<int>() {
            return MPI_INT;
        }
        
        template <>
        inline MPI_Datatype typeMPI<long>() {
            return MPI_LONG;
        }
        
        template <>
        inline MPI_Datatype typeMPI<float>() {
            return MPI_FLOAT;
        }
        
        template <>
        inline MPI_Datatype typeMPI<double>() {
            return MPI_DOUBLE;
        }
        
        template <>
        inline MPI_Datatype typeMPI<std::complex<float>>() {
            return MPI_C_FLOAT_COMPLEX;
        }
        
        template <>
        inline MPI_Datatype typeMPI<std::complex<double>>() {
            return MPI_C_DOUBLE_COMPLEX;
        }
        
        template <>
        inline MPI_Datatype typeMPI<char>() {
            return MPI_CHAR;
        }
#endif
    }
    
    
    ////////////////////////////// basics //////////////////////////////
    // initialize MPI
    void initialize(int *argc, char ***argv);
    
    // finalize MPI
    void finalize();
    
    // setup group
    void setupGroup(int nprocPerGroup);
    
    // free group MPI_Comm
    void freeGroupComm();
    
    // number of processors
    inline int nproc() {
#ifndef _SERIAL_BUILD
        static int nproc;
        MPI_Comm_size(internal::iCommCurrent, &nproc);
        return nproc;
#else
        return 1;
#endif
    }
    
    // MPI rank
    inline int rank() {
#ifndef _SERIAL_BUILD
        static int rank;
        MPI_Comm_rank(internal::iCommCurrent, &rank);
        return rank;
#else
        return 0;
#endif
    }
    
    // number of processors in World
    inline int nprocWorld() {
#ifndef _SERIAL_BUILD
        static int nproc;
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        return nproc;
#else
        return 1;
#endif
    }
    
    // MPI rank in World
    inline int rankWorld() {
#ifndef _SERIAL_BUILD
        static int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
#else
        return 0;
#endif
    }
    
    // MPI rank string
    inline std::string strRank() {
        std::stringstream ss;
        ss << rank();
        return ss.str();
    }
    
    // enter world group level
    inline void enterWorld() {
#ifndef _SERIAL_BUILD
        internal::iCommCurrent = MPI_COMM_WORLD;
#endif
    }
    
    // enter super group level
    inline void enterSuper() {
#ifndef _SERIAL_BUILD
        internal::iCommCurrent = internal::iCommSuper;
#endif
    }
    
    // enter infer group level
    inline void enterInfer() {
#ifndef _SERIAL_BUILD
        internal::iCommCurrent = internal::iCommInfer;
#endif
    }
    
    // root
    inline bool root() {
        return rank() == 0;
    }
    
    // am I a super proc in my group
    inline bool super() {
        return internal::iSuper;
    }
    
    // barrier
    void barrier();
    
    // abort
    void abort(int err = 1);
    
    
    ////////////////////////////// broadcast //////////////////////////////
    // raw array, ALLOCATED
    template <typename T>
    void bcast(T *valptr, int size, int src = 0) {
#ifndef _SERIAL_BUILD
        MPI_Bcast(valptr, size, internal::typeMPI<T>(), src,
                  internal::iCommCurrent);
#endif
    }
    
    // single
    template <typename T>
    void bcast(T &value, int src = 0) {
        bcast(&value, 1, src);
    }
    
    // std::vector
    template <typename T>
    void bcast(std::vector<T> &vec, int src = 0) {
        // size
        int size = 0;
        if (rank() == src) {
            size = (int)vec.size();
        }
        bcast(size, src);
        // allocate
        if (rank() != src) {
            vec.resize(size);
        }
        // data
        bcast(vec.data(), size, src);
    }
    
    // std::vector<std::vector>
    template <typename T>
    void bcast(std::vector<std::vector<T>> &vecVec, int src = 0) {
        // concatenate before bcast for minimum communication
        std::vector<T> flattened;
        std::vector<int> lengths;
        if (rank() == src) {
            lengths.assign(vecVec.size(), 0);
            for (int ivec = 0; ivec < vecVec.size(); ivec++) {
                flattened.insert(flattened.end(), vecVec[ivec].begin(),
                                 vecVec[ivec].end());
                lengths[ivec] = (int)vecVec[ivec].size();
            }
        }
        
        // bcast
        bcast(flattened, src);
        bcast(lengths, src);
        
        // cast back
        if (rank() != src) {
            vecVec.assign(lengths.size(), std::vector<T>());
            int pos = 0;
            for (int ivec = 0; ivec < lengths.size(); ivec++) {
                vecVec[ivec].insert(vecVec[ivec].begin(),
                                    flattened.begin() + pos,
                                    flattened.begin() + pos + lengths[ivec]);
                pos += lengths[ivec];
            }
        }
    }
    
    // Eigen::Matrix
    template <typename EigenMat>
    void bcastEigen(EigenMat &mat, int src = 0) {
        // size
        int dim[2];
        if (rank() == src) {
            dim[0] = (int)mat.rows();
            dim[1] = (int)mat.cols();
        }
        bcast(dim, 2, src);
        // allocate
        if (rank() != src) {
            mat.resize(dim[0], dim[1]);
        }
        // data
        bcast(mat.data(), (int)mat.size(), src);
    }
    
    // specialization for string
    void bcast(std::string &str, int src = 0);
    
    // specialization for vector of string
    void bcast(std::vector<std::string> &vecStr, int src = 0);
    
    
    ////////////////////////////// isend / irecv //////////////////////////////
    // isend, for ALLOCATED Eigen::Matrix
    template <typename EigenMat>
    void isend(int dest, const EigenMat &mat, MPI_Request &request) {
#ifndef _SERIAL_BUILD
        MPI_Isend(mat.data(), (int)mat.size(),
                  internal::typeMPI<typename EigenMat::Scalar>(),
                  dest, dest, internal::iCommCurrent, &request);
#endif
    }
    
    // irecv, for ALLOCATED Eigen::Matrix
    template <typename EigenMat>
    void irecv(int source, EigenMat &mat, MPI_Request &request) {
#ifndef _SERIAL_BUILD
        MPI_Irecv(mat.data(), (int)mat.size(),
                  internal::typeMPI<typename EigenMat::Scalar>(),
                  source, rank(), internal::iCommCurrent, &request);
#endif
    }
    
    // wait_all: must be implemented in .cpp
    void waitAll(std::vector<MPI_Request> &requests);
    
    
    ////////////////////////////// send / recv //////////////////////////////
    template <typename EigenMat>
    void sendVecEigen(int dest, const std::vector<EigenMat> &vecMat, int tag) {
        if (dest == mpi::rank()) {
            return;
        }
        
        // structured to flattened
        int nmat = (int)vecMat.size();
        std::vector<int> rowCol(nmat * 2, 0);
        std::vector<typename EigenMat::Scalar> flattened;
        for (int imat = 0; imat < nmat; imat++) {
            rowCol[imat * 2 + 0] = (int)vecMat[imat].rows();
            rowCol[imat * 2 + 1] = (int)vecMat[imat].cols();
            flattened.insert(flattened.end(),
                             vecMat[imat].data(),
                             vecMat[imat].data() + vecMat[imat].size());
        }
        
        // check: size < max int
        if (flattened.size() > std::numeric_limits<int>::max()) {
            throw std::runtime_error
            ("mpi::sendVecEigen || "
             "Sending data with size larger than MAX_INT. || "
             "Consider using more processors.");
        }
        
        // size info
        std::array<int, 2> sizeInfo;
        sizeInfo[0] = nmat;
        sizeInfo[1] = (int)flattened.size();
        
        // send
#ifndef _SERIAL_BUILD
        // size info
        MPI_Send(sizeInfo.data(), 2, MPI_INT, dest,
                 tag * 3 + 0, internal::iCommCurrent);
        
        // row, col
        MPI_Send(rowCol.data(), nmat * 2, MPI_INT, dest,
                 tag * 3 + 1, internal::iCommCurrent);
        
        // data
        MPI_Send(flattened.data(), sizeInfo[1],
                 internal::typeMPI<typename EigenMat::Scalar>(), dest,
                 tag * 3 + 2, internal::iCommCurrent);
#endif
    }
    
    template <typename EigenMat>
    void recvVecEigen(int source, std::vector<EigenMat> &vecMat, int tag) {
        if (source == mpi::rank()) {
            return;
        }
        
        // size info
        std::array<int, 2> sizeInfo;
#ifndef _SERIAL_BUILD
        MPI_Recv(sizeInfo.data(), 2, MPI_INT, source,
                 tag * 3 + 0, internal::iCommCurrent, MPI_STATUS_IGNORE);
#endif
        int nmat = sizeInfo[0];
        int nflattened = sizeInfo[1];
        
        // row, col
        std::vector<int> rowCol(nmat * 2, 0);
#ifndef _SERIAL_BUILD
        MPI_Recv(rowCol.data(), nmat * 2, MPI_INT, source,
                 tag * 3 + 1, internal::iCommCurrent, MPI_STATUS_IGNORE);
#endif
        
        // data
        std::vector<typename EigenMat::Scalar> flattened(nflattened, 0);
#ifndef _SERIAL_BUILD
        MPI_Recv(flattened.data(), nflattened,
                 internal::typeMPI<typename EigenMat::Scalar>(), source,
                 tag * 3 + 2, internal::iCommCurrent, MPI_STATUS_IGNORE);
#endif
        
        // flattened to structured
        vecMat.clear();
        int pos = 0;
        for (int imat = 0; imat < nmat; imat++) {
            int row = rowCol[imat * 2 + 0];
            int col = rowCol[imat * 2 + 1];
            const EigenMat &mat =
            Eigen::Map<const EigenMat>(flattened.data() + pos, row, col);
            vecMat.push_back(mat);
            pos += mat.size();
        }
    }
    
    
    ////////////////////////////// reduce //////////////////////////////
    // min
    template <typename T>
    T min(const T &value) {
        T result = value;
#ifndef _SERIAL_BUILD
        MPI_Allreduce(&value, &result, 1, internal::typeMPI<T>(), MPI_MIN,
                      internal::iCommCurrent);
#endif
        return result;
    }
    
    // max
    template <typename T>
    T max(const T &value) {
        T result = value;
#ifndef _SERIAL_BUILD
        MPI_Allreduce(&value, &result, 1, internal::typeMPI<T>(), MPI_MAX,
                      internal::iCommCurrent);
#endif
        return result;
    }
    
    // sum
    template <typename T>
    T sum(const T &value) {
        T result = value;
#ifndef _SERIAL_BUILD
        MPI_Allreduce(&value, &result, 1, internal::typeMPI<T>(), MPI_SUM,
                      internal::iCommCurrent);
#endif
        return result;
    }
    
    // sum eigen in-place
    template <typename EigenMat>
    void sumEigen(EigenMat &value) {
        EigenMat total(value);
#ifndef _SERIAL_BUILD
        MPI_Allreduce(value.data(), total.data(), (int)value.size(),
                      internal::typeMPI<typename EigenMat::Scalar>(), MPI_SUM,
                      internal::iCommCurrent);
#endif
        value = total;
    }
    
    
    ////////////////////////////// gather //////////////////////////////
    // gather single
    template <typename T>
    void gather(T val, std::vector<T> &vecVal, int dest) {
        if (dest == MPI_ALL) {
            vecVal.resize(nproc());
#ifndef _SERIAL_BUILD
            MPI_Allgather(&val, 1, internal::typeMPI<T>(), vecVal.data(), 1,
                          internal::typeMPI<T>(), internal::iCommCurrent);
#else
            vecVal[0] = val;
#endif
        } else {
            if (rank() == dest) {
                vecVal.resize(nproc());
            }
#ifndef _SERIAL_BUILD
            MPI_Gather(&val, 1, internal::typeMPI<T>(), vecVal.data(), 1,
                       internal::typeMPI<T>(), dest, internal::iCommCurrent);
#else
            vecVal[0] = val;
#endif
        }
    }
    
    // vector to vector<vector>
    template <typename T>
    void gather(const std::vector<T> &vec,
                std::vector<std::vector<T>> &vecVec, int dest) {
        // size
        int size = (int)vec.size();
        std::vector<int> vSize;
        gather(size, vSize, dest);
        
        // allocate flattened
        std::vector<int> vDisp;
        std::vector<T> vecVecFlatten;
        if (dest == MPI_ALL || dest == rank()) {
            int totalSize = 0;
            vDisp.assign(nproc(), 0);
            for (int iproc = 0; iproc < nproc(); iproc++) {
                vDisp[iproc] = totalSize;
                totalSize += vSize[iproc];
            }
            vecVecFlatten.resize(totalSize);
        }
        
        // bcast flattened
        if (dest == MPI_ALL) {
#ifndef _SERIAL_BUILD
            MPI_Allgatherv(vec.data(), size, internal::typeMPI<T>(),
                           vecVecFlatten.data(), vSize.data(), vDisp.data(),
                           internal::typeMPI<T>(), internal::iCommCurrent);
#else
            vecVecFlatten = vec;
#endif
        } else {
#ifndef _SERIAL_BUILD
            MPI_Gatherv(vec.data(), size, internal::typeMPI<T>(),
                        vecVecFlatten.data(), vSize.data(), vDisp.data(),
                        internal::typeMPI<T>(), dest, internal::iCommCurrent);
#else
            vecVecFlatten = vec;
#endif
        }
        
        // cast back to nested vector
        if (dest == MPI_ALL || dest == rank()) {
            vecVec.clear();
            vecVec.reserve(nproc());
            for (int iproc = 0; iproc < nproc(); iproc++) {
                std::vector<T> vecRank(vecVecFlatten.begin() + vDisp[iproc],
                                       vecVecFlatten.begin() + vDisp[iproc] +
                                       vSize[iproc]);
                vecVec.push_back(vecRank);
            }
        }
    }
    
    // Eigen::Matrix to vector<Eigen::Matrix>
    template <typename EigenMat>
    void gatherEigen(const EigenMat &mat,
                     std::vector<EigenMat> &vecMat, int dest) {
        // eigen to vector
        std::vector<typename EigenMat::Scalar> vec(mat.data(),
                                                   mat.data() + mat.size());
        // need both row and col to handle a zero-sized matrix with
        // one of the dimensions known at compile time
        int row = (int)mat.rows();
        int col = (int)mat.cols();
        
        // gahter vector
        std::vector<std::vector<typename EigenMat::Scalar>> vecVec;
        std::vector<int> vecRow, vecCol;
        gather(vec, vecVec, dest);
        gather(row, vecRow, dest);
        gather(col, vecCol, dest);
        
        // vector to eigen
        if (dest == MPI_ALL || dest == rank()) {
            vecMat.reserve(nproc());
            for (int iproc = 0; iproc < nproc(); iproc++) {
                vecMat.push_back(Eigen::Map<EigenMat>(vecVec[iproc].data(),
                                                      vecRow[iproc],
                                                      vecCol[iproc]));
            }
        }
    }
    
    // specialization for string
    void gather(const std::string &str,
                std::vector<std::string> &vecStr, int dest);
    
    // specialization for vector of string
    void gather(const std::vector<std::string> &vecStr,
                std::vector<std::vector<std::string>> &vecVecStr, int dest);
    
    
    ////////////////////////////// scatter //////////////////////////////
    // single
    template <typename T>
    void scatter(const std::vector<T> &vecVal, T &val, int src) {
#ifndef _SERIAL_BUILD
        MPI_Scatter(vecVal.data(), 1, internal::typeMPI<T>(), &val, 1,
                    internal::typeMPI<T>(), src, internal::iCommCurrent);
#else
        val = vecVal[0];
#endif
    }
    
    // vector<vector> to vector
    template <typename T>
    void scatter(const std::vector<std::vector<T>> &vecVec,
                 std::vector<T> &vec, int src) {
        // to be formed on src
        std::vector<int> vSize;
        std::vector<int> vDisp;
        std::vector<T> vecVecFlatten;
        if (rank() == src) {
            int totalSize = 0;
            vSize.assign(nproc(), 0);
            vDisp.assign(nproc(), 0);
            for (int iproc = 0; iproc < nproc();iproc++) {
                vSize[iproc] = (int)vecVec[iproc].size();
                vDisp[iproc] = totalSize;
                totalSize += vSize[iproc];
                vecVecFlatten.insert(vecVecFlatten.end(),
                                     vecVec[iproc].begin(),
                                     vecVec[iproc].end());
            }
        }
        
        // scatter size
        int size = 0;
        scatter(vSize, size, src);
        vec.resize(size);
        
        // scatter flattened
#ifndef _SERIAL_BUILD
        MPI_Scatterv(vecVecFlatten.data(), vSize.data(), vDisp.data(),
                     internal::typeMPI<T>(), vec.data(), size,
                     internal::typeMPI<T>(), src, internal::iCommCurrent);
#else
        vec = vecVec[0];
#endif
    }
    
    // vector<Eigen::Matrix> to Eigen::Matrix
    template <typename EigenMat>
    void scatterEigen(const std::vector<EigenMat> &vecMat,
                      EigenMat &mat, int src) {
        // eigen to vector
        std::vector<std::vector<typename EigenMat::Scalar>> vecVec;
        std::vector<int> vecRow, vecCol;
        if (rank() == src) {
            vecVec.reserve(nproc());
            vecRow.resize(nproc());
            vecCol.resize(nproc());
            for (int iproc = 0; iproc < nproc(); iproc++) {
                vecVec.push_back(std::vector<typename EigenMat::Scalar>
                                 (vecMat[iproc].data(),
                                  vecMat[iproc].data() + vecMat[iproc].size()));
                vecRow[iproc] = (int)vecMat[iproc].rows();
                vecCol[iproc] = (int)vecMat[iproc].cols();
            }
        }
        
        // scatter vector
        std::vector<typename EigenMat::Scalar> vec;
        int row = 0, col = 0;
        scatter(vecVec, vec, src);
        scatter(vecRow, row, src);
        scatter(vecCol, col, src);
        
        // vector to eigen
        mat = Eigen::Map<EigenMat>(vec.data(), row, col);
    }
    
    
    ////////////////////////////// map //////////////////////////////
    // gather map
    template <typename KEY, typename VAL>
    void gather(const std::map<KEY, VAL> &map,
                std::vector<std::map<KEY, VAL>> &vecMap, int dest) {
        // seperate keys and vals
        std::vector<KEY> keys;
        std::vector<VAL> vals;
        for (auto it = map.begin(); it != map.end(); ++it) {
            keys.push_back(it->first);
            vals.push_back(it->second);
        }
        
        // gather keys and vals
        std::vector<std::vector<KEY>> allKeys;
        std::vector<std::vector<VAL>> allVals;
        gather(keys, allKeys, dest);
        gather(vals, allVals, dest);
        
        // back to map
        if (dest == MPI_ALL || dest == rank()) {
            vecMap.clear();
            for (int iproc = 0; iproc < nproc(); iproc++) {
                std::map<KEY, VAL> mapIP;
                for (int ikey = 0; ikey < allKeys[iproc].size(); ikey++) {
                    mapIP.insert({allKeys[iproc][ikey], allVals[iproc][ikey]});
                }
                vecMap.push_back(mapIP);
            }
        }
    }
    
    // aggregate map by sum
    template <typename KEY, typename VAL>
    void aggregate(std::map<KEY, VAL> &map, int dest, MPI_Op op) {
        // gather
        std::vector<std::map<KEY, VAL>> allMaps;
        gather(map, allMaps, dest);
        // aggregate all into one map
        if (dest == MPI_ALL || dest == rank()) {
            map.clear();
            for (int iproc = 0; iproc < nproc(); iproc++) {
                for (auto it = allMaps[iproc].begin();
                     it != allMaps[iproc].end(); ++it) {
                    const KEY &key = it->first;
                    if (op == MPI_SUM) {
                        map.insert({key, (VAL)0});
                        map.at(key) += it->second;
                    } else if (op == MPI_MAX) {
                        map.insert({key, std::numeric_limits<VAL>::lowest()});
                        map.at(key) = std::max(map.at(key), it->second);
                    } else if (op == MPI_MIN) {
                        map.insert({key, std::numeric_limits<VAL>::max()});
                        map.at(key) = std::min(map.at(key), it->second);
                    }
                }
            }
        }
    }
    
    
    ////////////////////////////// mesh //////////////////////////////
    // convert nodal field to elemental
    template <typename IndexMat, typename DataMat>
    DataMat nodalToElemental(const IndexMat &connectivity,
                             const DataMat &nodal, bool collective) {
        int nelem = (int)connectivity.rows();
        int nvetex = (int)connectivity.cols();
        DataMat elemental = DataMat::Zero(nelem, nodal.cols());
        for (int ieg = 0; ieg < nelem; ieg++) {
            if (ieg % nproc() == rank() || !collective) {
                for (int ivt = 0; ivt < nvetex; ivt++) {
                    elemental.row(ieg) += nodal.row(connectivity(ieg, ivt));
                }
            }
        }
        if (collective) {
            sumEigen(elemental);
        }
        elemental /= nvetex;
        return elemental;
    }
    
    
    ////////////////////////////// external //////////////////////////////
    // verbose
    std::string verbose();
}

#endif /* mpi_hpp */
