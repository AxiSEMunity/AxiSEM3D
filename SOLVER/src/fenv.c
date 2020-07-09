//
//  fenv.c
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 6/7/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  setup floating-point environment

#include <fenv.h>

void fenv_setup() {
#ifndef _SKIP_DISABLE_SSE_DENORMS
    // To activate flush to zero for denormal float handling
    // google: "Why does changing 0.1f to 0 slow down performance by 10x?"
    fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
#endif
}
