//
//  spectral.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  polynomial order of spectral elements

#ifndef spectral_hpp
#define spectral_hpp

#ifndef _NPOL
#define _NPOL 4
#endif

namespace spectral {
    // npol
    const int nPol = _NPOL;
    
    // number of points on an edge
    const int nPntEdge = nPol + 1;
    
    // number of points in an element
    const int nPntElem = nPntEdge * nPntEdge;
    
    // shortened alias
    const int nPED = nPntEdge;
    const int nPEM = nPntElem;
}

#endif /* spectral_hpp */
