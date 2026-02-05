//
//  fenv.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 6/7/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

// fenv.cpp

#ifndef _SKIP_DISABLE_SSE_DENORMS
  #if (defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64)) && defined(__SSE__)
    #include <xmmintrin.h>
    #define A3D_HAS_FTZ 1
  #else
    #define A3D_HAS_FTZ 0
  #endif
#else
  #define A3D_HAS_FTZ 0
#endif

void fenv_setup() {
#if A3D_HAS_FTZ
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
}
