#ifndef DSFMT_PARAMS_H
#define DSFMT_PARAMS_H

#include "dSFMT.h"

/*----------------------
  the parameters of DSFMT
  following definitions are in dSFMT-paramsXXXX.h file.
  ----------------------*/
/** the pick up position of the array.
#define DSFMT_POS1 122 
*/

/** the parameter of shift left as four 32-bit registers.
#define DSFMT_SL1 18
 */

/** the parameter of shift right as four 32-bit registers.
#define DSFMT_SR1 12
*/

/** A bitmask, used in the recursion.  These parameters are introduced
 * to break symmetry of SIMD.
#define DSFMT_MSK1 (uint64_t)0xdfffffefULL
#define DSFMT_MSK2 (uint64_t)0xddfecb7fULL
*/

/** These definitions are part of a 128-bit period certification vector.
#define DSFMT_PCV1	UINT64_C(0x00000001)
#define DSFMT_PCV2	UINT64_C(0x00000000)
*/

#define DSFMT_LOW_MASK  UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_SR    12

/* for sse2 */
#if defined(HAVE_SSE2)
#define SSE2_SHUFF 0x1b
#elif defined(HAVE_ALTIVEC)
#if defined(__APPLE__)  /* For OSX */
#define ALTI_SR (vector unsigned char)(4)
#define ALTI_SR_PERM \
        (vector unsigned char)(15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14)
#define ALTI_SR_MSK \
        (vector unsigned int)(0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU)
#define ALTI_PERM \
        (vector unsigned char)(12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3)
#else
#define ALTI_SR      {4}
#define ALTI_SR_PERM {15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14}
#define ALTI_SR_MSK  {0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU}
#define ALTI_PERM    {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3}
#endif
#endif

#ifndef DSFMT_PARAMS19937_H
#define DSFMT_PARAMS19937_H

/* #define DSFMT_N	191 */
/* #define DSFMT_MAXDEGREE	19992 */
#define DSFMT_POS1    117
#define DSFMT_SL1    19
#define DSFMT_MSK1    UINT64_C(0x000ffafffffffb3f)
#define DSFMT_MSK2    UINT64_C(0x000ffdfffc90fffd)
#define DSFMT_MSK32_1    0x000ffaffU
#define DSFMT_MSK32_2    0xfffffb3fU
#define DSFMT_MSK32_3    0x000ffdffU
#define DSFMT_MSK32_4    0xfc90fffdU
#define DSFMT_FIX1    UINT64_C(0x90014964b32f4329)
#define DSFMT_FIX2    UINT64_C(0x3b8d12ac548a7c7a)
#define DSFMT_PCV1    UINT64_C(0x3d84e1ac0dc82880)
#define DSFMT_PCV2    UINT64_C(0x0000000000000001)
#define DSFMT_IDSTR    "dSFMT2-19937:117-19:ffafffffffb3f-ffdfffc90fffd"


/* PARAMETERS FOR ALTIVEC */
#if defined(__APPLE__)    /* For OSX */
#define ALTI_SL1 	(vector unsigned int)(3, 3, 3, 3)
#define ALTI_SL1_PERM \
    (vector unsigned char)(2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1)
#define ALTI_SL1_MSK \
    (vector unsigned int)(0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U)
#define ALTI_MSK	(vector unsigned int)(DSFMT_MSK32_1, \
            DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4)
#else	/* For OTHER OSs(Linux?) */
#define ALTI_SL1    {3, 3, 3, 3}
#define ALTI_SL1_PERM \
    {2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1}
#define ALTI_SL1_MSK \
    {0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U}
#define ALTI_MSK \
    {DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4}
#endif

#endif /* DSFMT_PARAMS19937_H */
#endif /* DSFMT_PARAMS_H */
