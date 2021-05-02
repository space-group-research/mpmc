/************************************************************************/
/* Cellular automata rule 30-derived pseudo random number generator	*/
/*									*/
/* This implementation is for both 32 and 64 bit architectures, but is	*/
/* ideally implemented for 64.  The code has been generalized so that	*/
/* class III systems other than rule 30 can be utilized for research	*/
/* purposes.  The scheme makes use a register array upon which the rule	*/
/* is imposed - the next step in the automata is generated in the set	*/
/* of output registers.  Each end of the register array is 'wrapped'	*/
/* into a circular register; while this reduces the periodic interval	*/
/* of rule 30, the literature notes that this probablistically occurs	*/
/* on the order of modern cryptographic systems.  To simplify the logic	*/
/* the rule is implemented directly into the output register array, and	*/
/* then shifted one bit to the right at the end in order for the output	*/
/* to align.  Notice that rule 30 could be more cheaply implemented	*/
/* (at the cost of generality) by (left XOR (middle OR right)) - that	*/
/* is not done here.  The designated center bit of each iteration is	*/
/* used to generate the mantissa of the double float returned, as per	*/
/* Wolfram, "A New Kind of Science".  It is known that Mathematica uses	*/
/* this exact method in it's implementation of Random[].		*/
/*									*/
/* Benchmark results:							*/
/*									*/
/* compile with:							*/
/*	gcc -o rule30rng -fomit-frame -funroll-loops -O3 rule30.rng.c	*/
/*									*/
/* @2005 Jonathan Belof							*/
/* Space Research Group							*/
/* Department of Chemistry						*/
/* University of South Florida						*/
/************************************************************************/

#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>

#include <time.h>

/* XXX JB */
#define OSX

#ifdef OSX
#include <limits.h>
#include <stdint.h>
#else
#include <bits/wordsize.h>
#include <values.h>
#endif /* OSX */

#define WORDSIZE __WORDSIZE

#if WORDSIZE == 64
/* 64-bit masks */
#define RULE30 0x000000000000001E  /* 0000000000000000000000000000000000000000000000000000000000011110 */
#define RULE110 0x000000000000006E /* 0000000000000000000000000000000000000000000000000000000001101110 */
#define RULE10 0x000000000000000A  /* 0000000000000000000000000000000000000000000000000000000000001010 */
#define RULE90 0x000000000000005A  /* 0000000000000000000000000000000000000000000000000000000001011010 */

#define CELL_MASK 0x0000000000000007      /* 0000000000000000000000000000000000000000000000000000000000000111 */
#define CENTER_MASK 0x0000000100000000    /* 0000000000000000000000000000000100000000000000000000000000000000 */
#define DELTA_CENTER 0x0000000000000020   /* 0000000000000000000000000000000000000000000000000000000000100000 */
#define RHS_ONE 0x0000000000000001        /* 0000000000000000000000000000000000000000000000000000000000000001 */
#define LHS_ONE 0x8000000000000000        /* 1000000000000000000000000000000000000000000000000000000000000000 */
#define LHS_ZERO 0x7FFFFFFFFFFFFFFF       /* 0111111111111111111111111111111111111111111111111111111111111111 */
#define INNER_COUNT 0x00000000000000FF    /* 0000000000000000000000000000000000000000000000000000000011111111 */
#define INNER_ONE 0x0000000000000001      /* 0000000000000000000000000000000000000000000000000000000000000001 */
#define INNER_ZERO 0xFFFFFFFFFFFFFF00     /* 1111111111111111111111111111111111111111111111111111111100000000 */
#define OUTER_COUNT 0x0000000000FFFF00    /* 0000000000000000000000000000000000000000111111111111111100000000 */
#define OUTER_ONE 0x0000000000000100      /* 0000000000000000000000000000000000000001000000000000000000000000 */
#define OUTER_ZERO 0xFFFFFFFFFF0000FF     /* 1111111111111111111111111111111111111111000000000000000011111111 */
#define DELTA_COUNT 0x0000000000000008    /* 0000000000000000000000000000000000000000000000000000000000001000 */
#define DELTA_MANTISSA 0x0000000000000034 /* 0000000000000000000000000000000000000000000000000000000000110100 */
#define MAX_MANTISSA 0x000FFFFFFFFFFFFF   /* 0000000000001111111111111111111111111111111111111111111111111111 */
#define SEED 0x38B1D098F2C40E5D
#else
/* 32-bit masks */
#define RULE30 0x0000001E  /* 00000000000000000000000000011110 */
#define RULE110 0x0000006E /* 00000000000000000000000001101110 */
#define RULE10 0x0000000A  /* 00000000000000000000000000001010 */
#define RULE90 0x0000005A  /* 00000000000000000000000001011010 */

#define CELL_MASK 0x00000007              /* 00000000000000000000000000000111 */
#define CENTER_MASK 0x00010000            /* 00000000000000010000000000000000 */
#define DELTA_CENTER 0x00000010           /* 00000000000000000000000000010000 */
#define RHS_ONE 0x00000001                /* 00000000000000000000000000000001 */
#define LHS_ONE 0x80000000                /* 10000000000000000000000000000000 */
#define LHS_ZERO 0x7FFFFFFF               /* 01111111111111111111111111111111 */
#define INNER_COUNT 0x000000FF            /* 00000000000000000000000011111111 */
#define INNER_ONE 0x00000001              /* 00000000000000000000000000000001 */
#define INNER_ZERO 0xFFFFFF00             /* 11111111111111111111111100000000 */
#define OUTER_COUNT 0x00FFFF00            /* 00000000111111111111111100000000 */
#define OUTER_ONE 0x00000100              /* 00000001000000000000000000000000 */
#define OUTER_ZERO 0xFF0000FF             /* 11111111000000000000000011111111 */
#define DELTA_COUNT 0x00000008            /* 00000000000000000000000000001000 */
#define DELTA_MANTISSA 0x0000000000000034 /* 0000000000000000000000000000000000000000000000000000000000110100 */
#define MAX_MANTISSA 0x00000000FFFFFFFF   /* 0000000000000000000000000000000011111111111111111111111111111111 */
#define SEED 0xF2C40E5D
#endif /* WORDSIZE == 64 */

#if WORDSIZE == 64 /* 64-bit */

double rule30_rng(unsigned long int seed) {
    register unsigned long int rule = RULE30; /* the rule to enforce */
    register unsigned long int in_reg1 = 0,   /* input registers */
        in_reg2 = 0,
                               in_reg3 = 0,
                               in_reg4 = 0,
                               in_reg5 = 0,
                               in_reg6 = 0,
                               in_reg7 = 0;
    register unsigned long int out_reg1 = 0, /* output registers */
        out_reg2 = 0,
                               out_reg3 = 0,
                               out_reg4 = 0,
                               out_reg5 = 0,
                               out_reg6 = 0,
                               out_reg7 = 0;
    register unsigned long int mp = 0;  /* multi-purpose register:					*/
                                        /* 	- the right-most 8 bits are for the inner loop counter	*/
                                        /* 	- the next 16 bits are for the outer loop counter	*/
                                        /* 	- the left-most bit is for carries			*/
    static unsigned long int last_reg1, /* static memory addrs to store results from the current run */
        last_reg2,
        last_reg3,
        last_reg4,
        last_reg5,
        last_reg6,
        last_reg7;
    double random_result = 0;                /* return a double from 0.0 to 1.0 */
    unsigned long int random_result_int = 0; /* integer version of the above for boolean ops */

    /* start with initial config */
    if (seed) {
        in_reg1 = in_reg2 = in_reg3 = in_reg4 = in_reg5 = in_reg6 = in_reg7 = seed;
    } else {
        /* already seeded - restore state from memory */
        in_reg1 = last_reg1;
        in_reg2 = last_reg2;
        in_reg3 = last_reg3;
        in_reg4 = last_reg4;
        in_reg5 = last_reg5;
        in_reg6 = last_reg6;
        in_reg7 = last_reg7;
    }

    for ((mp &= OUTER_ZERO); ((mp & OUTER_COUNT) >> DELTA_COUNT) < DELTA_MANTISSA; mp += OUTER_ONE) { /* <-- notice the fact that the increment here */
        for ((mp &= INNER_ZERO); (mp & INNER_COUNT) < WORDSIZE; mp += INNER_ONE) {                    /* will blow away the low-order bits doesn't matter */

            /* mask off first three bits and compare with rule */
            /* set the output register bit appropriately */
            out_reg1 |= ((rule >> (in_reg1 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg2 |= ((rule >> (in_reg2 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg3 |= ((rule >> (in_reg3 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg4 |= ((rule >> (in_reg4 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg5 |= ((rule >> (in_reg5 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg6 |= ((rule >> (in_reg6 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg7 |= ((rule >> (in_reg7 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);

            /* rotate all input registers one bit to the right, preserve carry */
            mp &= LHS_ZERO;                                /* clear the carry bit */
            mp |= ((in_reg7 & RHS_ONE) << (WORDSIZE - 1)); /* set carry bit if needed */
            in_reg7 >>= RHS_ONE;
            in_reg7 |= ((in_reg6 & RHS_ONE) << (WORDSIZE - 1));
            in_reg6 >>= RHS_ONE;
            in_reg6 |= ((in_reg5 & RHS_ONE) << (WORDSIZE - 1));
            in_reg5 >>= RHS_ONE;
            in_reg5 |= ((in_reg4 & RHS_ONE) << (WORDSIZE - 1));
            in_reg4 >>= RHS_ONE;
            in_reg4 |= ((in_reg3 & RHS_ONE) << (WORDSIZE - 1));
            in_reg3 >>= RHS_ONE;
            in_reg3 |= ((in_reg2 & RHS_ONE) << (WORDSIZE - 1));
            in_reg2 >>= RHS_ONE;
            in_reg2 |= ((in_reg1 & RHS_ONE) << (WORDSIZE - 1));
            in_reg1 >>= RHS_ONE;
            in_reg1 |= mp & LHS_ONE;
        }

        /* now must rotate output registers one bit to the left */
        mp &= LHS_ZERO;           /* clear the carry bit */
        mp |= out_reg1 & LHS_ONE; /* set the carry bit if needed */
        out_reg1 <<= RHS_ONE;
        out_reg1 |= ((out_reg2 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg2 <<= RHS_ONE;
        out_reg2 |= ((out_reg3 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg3 <<= RHS_ONE;
        out_reg3 |= ((out_reg4 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg4 <<= RHS_ONE;
        out_reg4 |= ((out_reg5 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg5 <<= RHS_ONE;
        out_reg5 |= ((out_reg6 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg6 <<= RHS_ONE;
        out_reg6 |= ((out_reg7 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg7 <<= RHS_ONE;
        out_reg7 |= ((mp & LHS_ONE) >> (WORDSIZE - 1));

        /* set output bits of random number */
        random_result_int |= ((out_reg4 & CENTER_MASK) >> DELTA_CENTER) << ((DELTA_MANTISSA - 1) - ((mp & OUTER_COUNT) >> DELTA_COUNT));

        /* swap the input and output registers */
        in_reg1 = out_reg1;
        in_reg2 = out_reg2;
        in_reg3 = out_reg3;
        in_reg4 = out_reg4;
        in_reg5 = out_reg5;
        in_reg6 = out_reg6;
        in_reg7 = out_reg7;

        /* clear output registers */
        out_reg1 = out_reg2 = out_reg3 = out_reg4 = out_reg5 = out_reg6 = out_reg7 = 0;
    }
    /* save last state point to static memory */
    last_reg1 = in_reg1;
    last_reg2 = in_reg2;
    last_reg3 = in_reg3;
    last_reg4 = in_reg4;
    last_reg5 = in_reg5;
    last_reg6 = in_reg6;
    last_reg7 = in_reg7;

    random_result = (double)random_result_int;
    random_result /= (double)MAX_MANTISSA; /* ensure that result is normalized from 0 to 1 */
    return (random_result);
}

#else  /* 32-bit */

double rule30_rng(unsigned long int seed) {
    register unsigned long int rule = RULE30; /* the rule to enforce */
    register unsigned long int in_reg1 = 0,   /* input registers */
        in_reg2 = 0,
                               in_reg3 = 0;
    register unsigned long int out_reg1 = 0, /* output registers */
        out_reg2 = 0,
                               out_reg3 = 0;
    register unsigned long int mp = 0;  /* multi-purpose register:					*/
                                        /* 	- the right-most 8 bits are for the inner loop counter	*/
                                        /* 	- the next 16 bits are for the outer loop counter	*/
                                        /* 	- the left-most bit is for carries			*/
    static unsigned long int last_reg1, /* static memory addrs to store results from the current run */
        last_reg2,
        last_reg3;

    double random_result = 0;                     /* return a double from 0.0 to 1.0 */
    unsigned long long int random_result_int = 0; /* integer version of the above for boolean ops */

    /* start with initial config */
    if (seed) {
        in_reg1 = in_reg2 = in_reg3 = seed;
    } else {
        /* already seeded - restore state from memory */
        in_reg1 = last_reg1;
        in_reg2 = last_reg2;
        in_reg3 = last_reg3;
    }

    for ((mp &= OUTER_ZERO); ((mp & OUTER_COUNT) >> DELTA_COUNT) < DELTA_MANTISSA; mp += OUTER_ONE) { /* <-- notice the fact that the increment here */
        for ((mp &= INNER_ZERO); (mp & INNER_COUNT) < WORDSIZE; mp += INNER_ONE) {                    /* will blow away the low-order bits doesn't matter */

            /* mask off first three bits and compare with rule */
            /* set the output register bit appropriately */
            out_reg1 |= ((rule >> (in_reg1 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg2 |= ((rule >> (in_reg2 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);
            out_reg3 |= ((rule >> (in_reg3 & CELL_MASK)) & RHS_ONE) << (mp & INNER_COUNT);

            /* rotate all input registers one bit to the right, preserve carry */
            mp &= LHS_ZERO;                                /* clear the carry bit */
            mp |= ((in_reg3 & RHS_ONE) << (WORDSIZE - 1)); /* set carry bit if needed */
            in_reg3 >>= RHS_ONE;
            in_reg3 |= ((in_reg2 & RHS_ONE) << (WORDSIZE - 1));
            in_reg2 >>= RHS_ONE;
            in_reg2 |= ((in_reg1 & RHS_ONE) << (WORDSIZE - 1));
            in_reg1 >>= RHS_ONE;
            in_reg1 |= mp & LHS_ONE;
        }

        /* now must rotate output registers one bit to the left */
        mp &= LHS_ZERO;           /* clear the carry bit */
        mp |= out_reg1 & LHS_ONE; /* set the carry bit if needed */
        out_reg1 <<= RHS_ONE;
        out_reg1 |= ((out_reg2 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg2 <<= RHS_ONE;
        out_reg2 |= ((out_reg3 & LHS_ONE) >> (WORDSIZE - 1));
        out_reg3 <<= RHS_ONE;
        out_reg3 |= ((mp & LHS_ONE) >> (WORDSIZE - 1));

        /* set output bits of random sequence */
        random_result_int |= ((out_reg2 & CENTER_MASK) >> DELTA_CENTER) << ((DELTA_MANTISSA - 1) - ((mp & OUTER_COUNT) >> DELTA_COUNT));

        /* swap the input and output registers */
        in_reg1 = out_reg1;
        in_reg2 = out_reg2;
        in_reg3 = out_reg3;

        /* clear the output registers */
        out_reg1 = out_reg2 = out_reg3 = 0;
    }
    /* save the last state point in static memory */
    last_reg1 = in_reg1;
    last_reg2 = in_reg2;
    last_reg3 = in_reg3;

    random_result = (double)random_result_int;
    random_result /= (double)MAX_MANTISSA; /* ensure that result is normalized from 0 to 1 */
    return (random_result);
}
#endif /* WORDSIZE == 64 */
