#include <mc.h>

// the point of this is to store all the polar_wolf_alpha calculations in a big array, and then just look shit up
// that way we don't need to calculate erfc's and exp's over and over and over

#define OneOverSqrtPi 0.56418958354

double* polar_wolf_alpha_lookup_init(system_t* system) {
    double* rval;
    int i;
    double a = system->polar_wolf_alpha;
    double r, rr;

    system->polar_wolf_alpha_table_max = ceil(system->polar_wolf_alpha_lookup_cutoff) * 1000;
    rval = malloc(sizeof(double) * system->polar_wolf_alpha_table_max);

    for (i = 1; i < system->polar_wolf_alpha_table_max; i++) {
        r = (double)i / 1000.;
        rr = 1.0 / r;
        rval[i] = erfc(a * r) * rr * rr + 2.0 * a * OneOverSqrtPi * exp(-a * a * r * r) * rr;
    }
    //store the zeroth component without blowing up
    rval[0] = rval[1];

    return rval;
}

double polar_wolf_alpha_getval(system_t* system, double r) {
    int i = (int)(r * 1000);
    if (i >= system->polar_wolf_alpha_table_max) return 0.0;  //answer will be zero if cutoff is large enough
    return system->polar_wolf_alpha_table[i];
}
