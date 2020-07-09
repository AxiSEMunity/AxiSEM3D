//
//  fenv.c
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 6/7/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  setup floating-point environment

#ifndef _SKIP_DISABLE_SSE_DENORMS
#include <xmmintrin.h>
#endif

void fenv_setup() {
#ifndef _SKIP_DISABLE_SSE_DENORMS
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
}
