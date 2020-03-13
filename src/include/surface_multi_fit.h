#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mc.h>

typedef struct {
    double *abInitioEnergy;
    double *fitEnergy, *lastFitEnergy;
    molecule_t **molecules;
    int nConfigs;
} multiConfigData_t;

typedef struct {
    char **atomtype;
    double *c6, *last_c6;
    double *c8, *last_c8;
    double *c10, *last_c10;
    double *epsilon, *last_epsilon;
    double *sigma, *last_sigma;
    double *omega, *last_omega;
    double *polarizability, *last_polarizability;
    int nParams;
} multiParamData_t;

void load_initial_multi_params(system_t *, multiParamData_t *);
void apply_multi_params(multiParamData_t *, multiConfigData_t *);
void accept_multi_params(multiParamData_t *, multiConfigData_t *);
void undo_multi_params(multiParamData_t *, multiConfigData_t *);

void read_multi_configs(system_t *, multiConfigData_t *, multiParamData_t *);
void perturb_multi_params(system_t *, multiParamData_t *);
double energy_multi_fit(system_t *);
double calc_multi_configurational_energy(system_t *);
double calc_multi_error(system_t *, multiConfigData_t *);
void output_multi_params(double, double, multiParamData_t *, system_t *, molecule_t *, multiConfigData_t *);
void output_best_config_energies(multiConfigData_t *);
int surface_multi_fit(system_t *);
