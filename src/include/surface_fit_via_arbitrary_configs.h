#include <stdio.h>
#include <stdlib.h>
#include <mc.h>

// atom (and/or site) type codes
static const unsigned int MOVEABLE = 1;
static const unsigned int FIXED = 0;

// A molecule-by-molecule, site-by-site Cartesian coord description of a system
// configuration, along with the ab initio Energy for that configuration, the model's
// current energy estimation for that configuration and the best fit (globally, so far)
// energy estimation for the configuration using the current set or any previous set of
// model-parameters.
//
// WHEN CREATING A moleculeCoords RECORD, ALWAYS set ALL pointers to zero, unless you are immediately
// populating them with data. Some functions (like storeSystemConfig()) overwrite existing data and,
// as such, will free any pointer found here in order to create a new array. If this is some number
// resulting from an uninitialized variable, a seg-fault will result.

typedef struct {
    int nSites;          // Number of atoms or sites in the molecule
    char **id;           // Array for the identifiers for each site
    double iCOM[3];      // initial Center of Mass (COM) for this molecule, this is the position about which perturbations will occur for moveable sites.
                         // Although the center of mass may change over time, this initial COM will be used for this entire run, as a fixed point, relative
                         // to which site positional perturbations can be directed along.
    unsigned int *type;  // Type code for the parameters that will be adjusted for each site
    double *x;           // an array of the x-coordinates for each site
    double *y;           // "             " y-coordinates "           "
    double *z;           // "             " z-coordinates "           "

} moleculeCoords;

typedef struct {
    char *id;
    double abInitioEnergy;     // Ab initio energy for this configuration
    double currentFitEnergy;   // Energy for this config using the current parameters
    double bestFitEnergy;      // Energy for this config using the best fit (so far) parameters
    int nMolecules;            // Number of molecules in this configuration
    moleculeCoords *molecule;  // An array of molecule descriptors for the molecules in this config
} configData_t;                // WHEN CREATING A configData_t RECORD, ALWAYS set this pointer to zero, unless you are immediately
                               // populating it with data. Some functions (like storeSystemConfig()) overwrite existing data and,
                               // as such, will free any pointer found here and create a new array. If this is some number
                               // resulting from an uninitialized variable, a seg-fault will result.

// calculates the current error for a surface energy versus the ab-initio curve
double error_calc_fr_arbitrary_configurations(system_t *system, int nConfigs, configData_t *configuration, double max_energy);

// apply_config() takes the information in a configData_t record and applies
// it to the system molecules. The configData_t records have an array of molecules,
// each of which has an array of x, y and z coordinates for every site in the molecule.
// This function loops through said molecules, and for each, loops through each site/atom
// transferring the Cartesian coords in the config array to the system atoms, which are
// located in a linked list, one list for each molecule, in a linked list of molecules.
void apply_configuration(system_t *system, configData_t *config);

// CaLculates the model-calculated energies associated with each an array of system configuration records
// The energy fields in said records are updated with the calculated values.
void get_configuration_energies(system_t *system, configData_t *configurations, int nConfigs, param_g *params, configData_t *orig);

// Reads a file of system configurations, and creates an array of configData_t records and stores the
// configuration information from the input file therein.
void read_config_fit_input_file(system_t *system, configData_t **c, int *nConfigs);

// fit potential energy parameters via simulated annealing
int surface_fit_via_arbitrary_configurations(system_t *system);

// store the current system configuration in a configData_t record.
void store_system_config(system_t *system, configData_t *originalConfig);

// Determine if an atom type is indicative of a FIXED type or a MOVEABLE type
int determine_mobility_code_from_atom_type(char *atomType);

void new_global_min_arbitrary_configs(system_t *, int, configData_t *);
