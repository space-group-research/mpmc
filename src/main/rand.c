#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mc.h>

/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

    To the extent possible under law, the author has dedicated all copyright
    and related and neighboring rights to this software to the public domain
    worldwide. This software is distributed without any warranty.

   see http://prng.di.unimi.it/ and
   http://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf

   This is xoroshiro128++ 1.0, one of our all-purpose, rock-solid,
   small-state generators. It is extremely (sub-ns) fast and it passes all
   tests we are aware of, but its state space is large enough only for
   mild parallelism.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */
static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

struct xoroshiro128pp_state {
    uint64_t s[2];
};

uint64_t xoroshiro128pp(struct xoroshiro128pp_state *state) {
    uint64_t *s = state->s;
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = rotl(s0 + s1, 17) + s0;

    s1 ^= s0;
    s[0] = rotl(s0, 49) ^ s1 ^ (s1 << 21);
    s[1] = rotl(s1, 28);

    return result;
}

/*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

    To the extent possible under law, the author has dedicated all copyright
    and related and neighboring rights to this software to the public domain
    worldwide. This software is distributed without any warranty.

   This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush; we're using it to seed
   our xoroshiro128++ state to avoid zeroland, etc */
struct splitmix64_state {
    uint64_t s;
};

uint64_t splitmix64(struct splitmix64_state *state) {
    uint64_t result = state->s;

    state->s = result + 0x9E3779B97f4A7C15;
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return result ^ (result >> 31);
}

/* Seed initial state with splitmix64 */
void xoroshiro128pp_init(struct xoroshiro128pp_state *state, uint64_t seed) {
    struct splitmix64_state smstate;
    smstate.s = seed;

    state->s[0] = splitmix64(&smstate);
    state->s[1] = splitmix64(&smstate);
}

/* Skips 2^96 ints, generates up to 2^32 non-overlapping sequences, we'll never need more than this */
void xoroshiro128pp_long_jump(struct xoroshiro128pp_state *state) {
    static const uint64_t LONG_JUMP[] = {0x360fd5f2cf8d5d99, 0x9c6e6877736c46e3};
    uint64_t *s = state->s;
    uint64_t s0 = 0;
    uint64_t s1 = 0;
    int i, b;
    for (i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
        for (b = 0; b < 64; b++) {
            if (LONG_JUMP[i] & UINT64_C(1) << b) {
                s0 ^= s[0];
                s1 ^= s[1];
            }
            xoroshiro128pp(state);
        }

    s[0] = s0;
    s[1] = s1;
}

double get_rand(system_t *system) {
    static struct xoroshiro128pp_state state;

    if (!system->rng_initialized) {
        system->rng_initialized = 1;
        if (system->preset_seeds_on) {
            xoroshiro128pp_init(&state, system->preset_seeds);
        } else {
            uint64_t seed;
            struct timespec timer;
            /* Try to seed using the nanosecond timer, use time(NULL) if that fails */
            if (clock_gettime(CLOCK_REALTIME, &timer) == 0)
                seed = timer.tv_nsec;
            else
                seed = (uint64_t)(time(NULL));
            xoroshiro128pp_init(&state, seed);
        }

        /* Avoid generating the same random numbers if we're using MPI */
        int i;
        for (i = 0; i < rank; i++)
            xoroshiro128pp_long_jump(&state);
    }

    return (xoroshiro128pp(&state) >> 11) * 0x1.0p-53; /* Get a double between 0 and 1 */
}
