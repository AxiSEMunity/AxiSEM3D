//
//  fenv.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 6/7/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

void fenv_setup() {
#ifndef _SKIP_DISABLE_SSE_DENORMS

  // Enable only if the compiler provides SSE intrinsics
  #if (defined(__i386__) || defined(__x86_64__)) && defined(__SSE__)
    #include <xmmintrin.h>
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  #endif

#endif
}
