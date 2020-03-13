// .-._                                                   _,-,
//  `._`-._                                           _,-'_,'
//     `._ `-._                                   _,-' _,'
//        `._  `-._        __.-----.__        _,-'  _,'
//           `._   `#==="""           """===#'   _,'
//              `._/)  ._               _.  (\_,'
//               )*'     **.__     __.**     '*(
//               #  .==..__  ""   ""  __..==,  #
//               #   `"._(_).       .(_)_."'   #

/*

@2007, Jonathan Belof
Space Research Group
Department of Chemistry
University of South Florida

*/

/* herein lie sleeping demons */
/* i truly hate this fucking file, but dont see a better way JB */

#include <mc.h>
#include <quaternion.h>
#ifdef MPI
#include <mpi.h>
#endif

//re-enumerate atoms and molecules -> set atom and molecule id's which get messed up in UVT runs
void enumerate_particles(system_t *system) {
    molecule_t *mptr;
    atom_t *aptr;
    int aid, mid;
    aid = mid = 1;

    for (mptr = system->molecules; mptr; mptr = mptr->next) {
        mptr->id = mid++;
        for (aptr = mptr->atoms; aptr; aptr = aptr->next)
            aptr->id = aid++;
    }

    return;
}

// parallel tempering stuff
/* parallel tempering is treated differently than the other types of MC moves, because it is done
 * IN ADDITION to any other moves: i.e. some random move is performed and the energy is calculated
 * before and after. once that is all figured out, we then perform a temper step. The boltzmann factor
 * that is calculated will NOT be averaged into the BF quantity. It can be, but care must be taken so
 * there is no double-counting, and it really isn't so important to do so. */
void temper_system(system_t *system, double current_energy) {
#ifdef MPI
    // system->ptemp->index[j] is a mapping from core -> bath_index (bath_index 0 is lowest temperature, etc.)
    // bath2core maps from bath_index -> core
    // only half of the cores will be responsible for carrying out the calculations
    // partner_list maps from core -> swap partner, if core is a master. core -> -1, if core is a slave.
    int *bath2core, *partner_list, *is_master, master, slave, *update_index, new_index_val;
    double slave_energy, boltzmann_factor, *update_templist, new_templist_val;
    int i, j, lucky, accept_move;
    MPI_Status status;

    slave_energy = 0;

    //make the bath index
    bath2core = calloc(size, sizeof(int));
    memnullcheck(bath2core, size * sizeof(int), __LINE__ - 1, __FILE__);
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            if (system->ptemp->index[j] == i) bath2core[i] = j;

    //choose the lucky bath. it's not really lucky.. this is just the first bath that we consider for swapping
    if (!rank) lucky = floor(get_rand(system) * size);
    MPI_Bcast(&lucky, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //we will use this array to designate whether a particular core is a master or slave
    is_master = calloc(size, sizeof(int));
    memnullcheck(is_master, size * sizeof(int), __LINE__ - 1, __FILE__);

    //build the partner list
    partner_list = calloc(size, sizeof(int));
    memnullcheck(partner_list, size * sizeof(int), __LINE__ - 1, __FILE__);
    master = lucky;
    slave = (lucky + 1) % size;  //this assigns the slave bath index to the variable
    do {
        partner_list[bath2core[slave]] = bath2core[master];  //partner_list maps slave to master core
        partner_list[bath2core[master]] = bath2core[slave];  //partner_list maps master to slave
        is_master[bath2core[master]] = 1;                    //master flag is on
        is_master[bath2core[slave]] = 0;                     //master flag is off
        //now generate the next master
        master = (master + 2) % size;
        slave = (slave + 2) % size;
    } while ((master != lucky) && (slave != lucky));
    if (slave == lucky) {
        partner_list[bath2core[master]] = -1;  //if odd # of cores, we have one member unassigned
        is_master[bath2core[master]] = -1;     //neither master nor slave
    }
    //all cores (except 1 if size is odd) are mapped

    //communicate the energy to the master nodes
    if (!is_master[rank]) {
        MPI_Send(&current_energy, 1, MPI_DOUBLE, partner_list[rank], 0, MPI_COMM_WORLD);     //slave sends it's energy
        MPI_Recv(&accept_move, 1, MPI_INT, partner_list[rank], 1, MPI_COMM_WORLD, &status);  //receive result (accept/reject)
    } else if (is_master[rank] == 1) {
        //master receives energy from slave
        MPI_Recv(&slave_energy, 1, MPI_DOUBLE, partner_list[rank], 0, MPI_COMM_WORLD, &status);
        //calculate boltzmann factor exp(dE*dB)
        boltzmann_factor = exp((current_energy - slave_energy) *
                               (1.0 / system->ptemp->templist[rank] - 1.0 / system->ptemp->templist[partner_list[rank]]));
        if (get_rand(system) < boltzmann_factor)
            accept_move = 1;
        else
            accept_move = 0;
        //communicate the move result to slave
        MPI_Send(&accept_move, 1, MPI_INT, partner_list[rank], 1, MPI_COMM_WORLD);
    } else
        accept_move = 0;  //no partner

    if (accept_move) {
        //reassign local temperature
        system->temperature = system->ptemp->templist[partner_list[rank]];
        //update our temperature and index values to prepare for transmission to root
        new_templist_val = system->ptemp->templist[partner_list[rank]];
        new_index_val = system->ptemp->index[partner_list[rank]];
        //reassign fugacities
        if (!system->user_fugacities && system->fugacities) {
            MPI_Send(system->fugacities, 1, MPI_DOUBLE, partner_list[rank], 2, MPI_COMM_WORLD);           //send fugacity
            MPI_Recv(system->fugacities, 1, MPI_DOUBLE, partner_list[rank], 2, MPI_COMM_WORLD, &status);  //receive fugacity
        }
    } else {  //reject
        new_templist_val = system->ptemp->templist[rank];
        new_index_val = system->ptemp->index[rank];
    }

    //now we need to update the templist and index on each core
    update_templist = calloc(size, sizeof(double));
    update_index = calloc(size, sizeof(int));
    memnullcheck(update_templist, size * sizeof(double), __LINE__ - 1, __FILE__);
    memnullcheck(update_index, size * sizeof(int), __LINE__ - 1, __FILE__);
    // create updated arrays on head node
    MPI_Gather(&new_templist_val, 1, MPI_DOUBLE, update_templist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&new_index_val, 1, MPI_INT, update_index, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // build the updated array
    if (!rank) {
        for (i = 0; i < size; i++) {
            system->ptemp->templist[i] = update_templist[i];
            system->ptemp->index[i] = update_index[i];
        }
    }
    // transmit to each core
    MPI_Bcast(system->ptemp->templist, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(system->ptemp->index, size, MPI_INT, 0, MPI_COMM_WORLD);

    if (accept_move)
        system->nodestats->accept_ptemp++;
    else if (is_master[rank] != -1)
        system->nodestats->reject_ptemp++;

    free(update_templist);
    free(update_index);
    free(bath2core);
    free(is_master);
    free(partner_list);
#endif  //MPI
    return;
}

/* scale the volume : used in NPT ensemble */
void volume_change(system_t *system) {
    molecule_t *m;
    atom_t *a;
    double new_volume, log_new_volume, basis_scale_factor;
    double new_com[3], old_com[3], delta_pos[3];
    int i, j;

    // figure out what the new volume will be
    if (system->ensemble == ENSEMBLE_REPLAY)
        //if ensemble replay, then we're just trying to calculate the pressure via dV change
        new_volume = system->pbc->volume + system->calc_pressure_dv;
    else {
        log_new_volume = log(system->pbc->volume) + (get_rand(system) - 0.5) * system->volume_change_factor;
        new_volume = exp(log_new_volume);
    }

    //scale basis
    basis_scale_factor = pow(new_volume / system->pbc->volume, 1.0 / 3.0);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            system->pbc->basis[i][j] *= basis_scale_factor;

    //recalculate PBC stuff (volume/reciprocal basis/cutoff)
    pbc(system);
    system->observables->volume = system->pbc->volume;

    //scale molecule positions
    for (m = system->molecules; m; m = m->next) {
        for (i = 0; i < 3; i++) {
            old_com[i] = m->com[i];  //molecule ptr's com will be udated during pairs routine
            new_com[i] = m->com[i] * basis_scale_factor;
            delta_pos[i] = new_com[i] - old_com[i];
        }
        for (a = m->atoms; a; a = a->next) {  //calculate new atomic positions based on new com
            for (i = 0; i < 3; i++) {
                a->pos[i] += delta_pos[i];
                a->wrapped_pos[i] += delta_pos[i];
            }
        }
    }

    return;
}

/* revert changes to volume */
void revert_volume_change(system_t *system) {
    double basis_scale_factor, new_volume;
    double old_com[3], new_com[3], delta_pos[3];
    molecule_t *m;
    atom_t *a;
    int i, j;
    new_volume = system->checkpoint->observables->volume;

    //un-scale basis
    basis_scale_factor = pow(new_volume / system->pbc->volume, 1.0 / 3.0);  //system->pbc->volume still has rejected value until pbc() is called
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            system->pbc->basis[i][j] *= basis_scale_factor;

    //recalc PBC stuff
    pbc(system);
    system->observables->volume = system->pbc->volume;

    //scale molecule positions
    for (m = system->molecules; m; m = m->next) {
        for (i = 0; i < 3; i++) {
            old_com[i] = m->com[i];  //molecule ptr's com will ned to be updated since this is a revert (no energy call following)
            new_com[i] = m->com[i] * basis_scale_factor;
            m->com[i] = new_com[i];
            delta_pos[i] = new_com[i] - old_com[i];
        }
        for (a = m->atoms; a; a = a->next) {  //calculate new atomic positions based on new com
            for (i = 0; i < 3; i++) {
                a->pos[i] += delta_pos[i];
                a->wrapped_pos[i] += delta_pos[i];
            }
        }
    }
    wrapall(system->molecules, system->pbc);
    return;
}

/* make an exact copy of src */
molecule_t *copy_molecule(system_t *system, molecule_t *src) {
    molecule_t *dst;
    atom_t *atom_dst_ptr, *prev_atom_dst_ptr, *atom_src_ptr;
    pair_t *pair_dst_ptr, *prev_pair_dst_ptr, *pair_src_ptr;

    /* allocate the start of the new lists */
    dst = calloc(1, sizeof(molecule_t));
    memnullcheck(dst, sizeof(molecule_t), __LINE__ - 1, __FILE__);
    /* copy molecule attributes */
    dst->id = src->id;
    strcpy(dst->moleculetype, src->moleculetype);
    dst->mass = src->mass;
    dst->frozen = src->frozen;
    dst->adiabatic = src->adiabatic;
    dst->spectre = src->spectre;
    dst->target = src->target;
    dst->nuclear_spin = src->nuclear_spin;
    dst->rot_partfunc_g = src->rot_partfunc_g;
    dst->rot_partfunc_u = src->rot_partfunc_u;
    dst->rot_partfunc = src->rot_partfunc;
    memcpy(dst->com, src->com, 3 * sizeof(double));
    memcpy(dst->wrapped_com, src->wrapped_com, 3 * sizeof(double));

#ifdef QM_ROTATION
    int i, j;
    if (system->quantum_rotation) {
        dst->quantum_rotational_energies = calloc(system->quantum_rotation_level_max, sizeof(double));
        memnullcheck(dst->quantum_rotational_energies, system->quantum_rotation_level_max * sizeof(double), __LINE__ - 1, __FILE__);
        dst->quantum_rotational_eigenvectors = calloc(system->quantum_rotation_level_max, sizeof(complex_t *));
        memnullcheck(dst->quantum_rotational_eigenvectors, system->quantum_rotation_level_max * sizeof(complex_t *), __LINE__ - 1, __FILE__);
        for (i = 0; i < system->quantum_rotation_level_max; i++) {
            dst->quantum_rotational_eigenvectors[i] = calloc((system->quantum_rotation_l_max + 1) * (system->quantum_rotation_l_max + 1), sizeof(complex_t));
            memnullcheck(dst->quantum_rotational_eigenvectors[i], (system->quantum_rotation_l_max + 1) * (system->quantum_rotation_l_max + 1) * sizeof(complex_t), __LINE__ - 1, __FILE__);
        }
        dst->quantum_rotational_eigensymmetry = calloc(system->quantum_rotation_level_max, sizeof(int));
        memnullcheck(dst->quantum_rotational_eigensymmetry, system->quantum_rotation_level_max * sizeof(int), __LINE__ - 1, __FILE__);

        memcpy(dst->quantum_rotational_energies, src->quantum_rotational_energies, system->quantum_rotation_level_max * sizeof(double));
        for (i = 0; i < system->quantum_rotation_level_max; i++) {
            for (j = 0; j < (system->quantum_rotation_l_max + 1) * (system->quantum_rotation_l_max + 1); j++) {
                dst->quantum_rotational_eigenvectors[i][j].real = src->quantum_rotational_eigenvectors[i][j].real;
                dst->quantum_rotational_eigenvectors[i][j].imaginary = src->quantum_rotational_eigenvectors[i][j].imaginary;
            }
        }
        memcpy(dst->quantum_rotational_eigensymmetry, src->quantum_rotational_eigensymmetry, system->quantum_rotation_level_max * sizeof(int));
    }
#endif /* QM_ROTATION */

    dst->next = NULL;

    /* new atoms list */
    dst->atoms = calloc(1, sizeof(atom_t));
    memnullcheck(dst->atoms, sizeof(atom_t), __LINE__ - 1, __FILE__);
    prev_atom_dst_ptr = dst->atoms;

    for (atom_dst_ptr = dst->atoms, atom_src_ptr = src->atoms; atom_src_ptr; atom_dst_ptr = atom_dst_ptr->next, atom_src_ptr = atom_src_ptr->next) {
        atom_dst_ptr->id = atom_src_ptr->id;
        strcpy(atom_dst_ptr->atomtype, atom_src_ptr->atomtype);
        atom_dst_ptr->frozen = atom_src_ptr->frozen;
        atom_dst_ptr->adiabatic = atom_src_ptr->adiabatic;
        atom_dst_ptr->spectre = atom_src_ptr->spectre;
        atom_dst_ptr->target = atom_src_ptr->target;
        atom_dst_ptr->mass = atom_src_ptr->mass;
        atom_dst_ptr->charge = atom_src_ptr->charge;
        atom_dst_ptr->gwp_alpha = atom_src_ptr->gwp_alpha;
        atom_dst_ptr->gwp_spin = atom_src_ptr->gwp_spin;
        atom_dst_ptr->polarizability = atom_src_ptr->polarizability;
        atom_dst_ptr->omega = atom_src_ptr->omega;
        atom_dst_ptr->epsilon = atom_src_ptr->epsilon;
        atom_dst_ptr->sigma = atom_src_ptr->sigma;
        atom_dst_ptr->gwp_spin = atom_src_ptr->gwp_spin;
        atom_dst_ptr->c6 = atom_src_ptr->c6;
        atom_dst_ptr->c8 = atom_src_ptr->c8;
        atom_dst_ptr->c10 = atom_src_ptr->c10;
        atom_dst_ptr->c9 = atom_src_ptr->c9;

        memcpy(atom_dst_ptr->pos, atom_src_ptr->pos, 3 * sizeof(double));
        memcpy(atom_dst_ptr->wrapped_pos, atom_src_ptr->wrapped_pos, 3 * sizeof(double));
        memcpy(atom_dst_ptr->ef_static, atom_src_ptr->ef_static, 3 * sizeof(double));
        memcpy(atom_dst_ptr->ef_induced, atom_src_ptr->ef_induced, 3 * sizeof(double));
        memcpy(atom_dst_ptr->mu, atom_src_ptr->mu, 3 * sizeof(double));
        memcpy(atom_dst_ptr->old_mu, atom_src_ptr->old_mu, 3 * sizeof(double));
        memcpy(atom_dst_ptr->new_mu, atom_src_ptr->new_mu, 3 * sizeof(double));

        atom_dst_ptr->pairs = calloc(1, sizeof(pair_t));
        memnullcheck(atom_dst_ptr->pairs, sizeof(pair_t), __LINE__ - 1, __FILE__);
        pair_dst_ptr = atom_dst_ptr->pairs;
        prev_pair_dst_ptr = pair_dst_ptr;
        for (pair_src_ptr = atom_src_ptr->pairs; pair_src_ptr; pair_src_ptr = pair_src_ptr->next) {
            pair_dst_ptr->rd_energy = pair_src_ptr->rd_energy;
            pair_dst_ptr->es_real_energy = pair_src_ptr->es_real_energy;
            pair_dst_ptr->es_self_intra_energy = pair_src_ptr->es_self_intra_energy;

            pair_dst_ptr->frozen = pair_src_ptr->frozen;
            pair_dst_ptr->rd_excluded = pair_src_ptr->rd_excluded;
            pair_dst_ptr->es_excluded = pair_src_ptr->es_excluded;
            //			pair_dst_ptr->charge = pair_src_ptr->charge;
            //			pair_dst_ptr->gwp_alpha = pair_src_ptr->gwp_alpha;
            //			pair_dst_ptr->gwp_spin = pair_src_ptr->gwp_spin;
            pair_dst_ptr->epsilon = pair_src_ptr->epsilon;
            pair_dst_ptr->lrc = pair_src_ptr->lrc;
            pair_dst_ptr->sigma = pair_src_ptr->sigma;
            pair_dst_ptr->r = pair_src_ptr->r;
            pair_dst_ptr->rimg = pair_src_ptr->rimg;

            pair_dst_ptr->next = calloc(1, sizeof(pair_t));
            memnullcheck(pair_dst_ptr->next, sizeof(pair_t), __LINE__ - 1, __FILE__);
            prev_pair_dst_ptr = pair_dst_ptr;
            pair_dst_ptr = pair_dst_ptr->next;
        }
        prev_pair_dst_ptr->next = NULL;
        free(pair_dst_ptr);
        /* handle an empty list */
        if (!atom_src_ptr->pairs) atom_dst_ptr->pairs = NULL;

        prev_atom_dst_ptr = atom_dst_ptr;
        atom_dst_ptr->next = calloc(1, sizeof(atom_t));
        memnullcheck(atom_dst_ptr->next, sizeof(atom_t), __LINE__ - 1, __FILE__);
    }

    prev_atom_dst_ptr->next = NULL;
    free(atom_dst_ptr);

    return (dst);
}

/* perform a general random translation */
void translate(system_t *system, molecule_t *molecule, pbc_t *pbc, double scale) {
    atom_t *atom_ptr;
    double trans_x, trans_y, trans_z;

    trans_x = scale * get_rand(system) * pbc->cutoff;
    trans_y = scale * get_rand(system) * pbc->cutoff;
    trans_z = scale * get_rand(system) * pbc->cutoff;
    if (get_rand(system) < 0.5) trans_x *= -1.0;
    if (get_rand(system) < 0.5) trans_y *= -1.0;
    if (get_rand(system) < 0.5) trans_z *= -1.0;

    molecule->com[0] += trans_x;
    molecule->com[1] += trans_y;
    molecule->com[2] += trans_z;

    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        atom_ptr->pos[0] += trans_x;
        atom_ptr->pos[1] += trans_y;
        atom_ptr->pos[2] += trans_z;
    }
}

/* perform a general random rotation */
/* now with quaternions AH */
void rotate(system_t *system, molecule_t *molecule, pbc_t *pbc, double scale) {
    atom_t *atom_ptr;
    double u1, u2, u3;
    double com[3];
    int i, ii, n;
    double *new_coord_array;

    struct quaternion position_vector, rnd_rotation, rnd_rotation_conjugate, answer;

    // create a random rotation
    // http://planning.cs.uiuc.edu/node198.html
    u1 = get_rand(system);
    u2 = get_rand(system);
    u3 = get_rand(system);

    // simple linear interpolation for adjusting the scale
    quaternion_construct_xyzw(&rnd_rotation, scale * sqrt(1 - u1) * sin(2 * M_PI * u2), scale * sqrt(1 - u1) * cos(2 * M_PI * u2), scale * sqrt(u1) * sin(2 * M_PI * u3), 1 - scale + scale * sqrt(u1) * cos(2 * M_PI * u3)); /* make a random quaternion */
    quaternion_normalize(&rnd_rotation);
    quaternion_conjugate(&rnd_rotation, &rnd_rotation_conjugate); /* needed to transform coordinates */

    /* count the number of atoms in a molecule, and allocate new coords array */
    for (atom_ptr = molecule->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next)
        ++n;
    new_coord_array = calloc(n * 3, sizeof(double));
    memnullcheck(new_coord_array, n * 3 * sizeof(double), __LINE__ - 1, __FILE__);

    /* save the com coordinate */
    com[0] = molecule->com[0];
    com[1] = molecule->com[1];
    com[2] = molecule->com[2];

    /* translate the molecule to the origin */
    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        atom_ptr->pos[0] -= com[0];
        atom_ptr->pos[1] -= com[1];
        atom_ptr->pos[2] -= com[2];
    }

    /* quaternion multiply */
    for (atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {
        ii = i * 3;

        //position_vector = position
        quaternion_construct_xyzw(&position_vector, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2], 0.);

        //answer = rnd_rotation*(position*rnd_rotation_conjugate)
        quaternion_multiplication(&position_vector, &rnd_rotation_conjugate, &answer);
        quaternion_multiplication(&rnd_rotation, &answer, &answer);

        //set the new coords
        new_coord_array[ii + 0] = answer.x;
        new_coord_array[ii + 1] = answer.y;
        new_coord_array[ii + 2] = answer.z;
    }

    /* set the new coordinates and then translate back from the origin */
    for (atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {
        ii = i * 3;
        atom_ptr->pos[0] = new_coord_array[ii + 0];
        atom_ptr->pos[1] = new_coord_array[ii + 1];
        atom_ptr->pos[2] = new_coord_array[ii + 2];

        atom_ptr->pos[0] += com[0];
        atom_ptr->pos[1] += com[1];
        atom_ptr->pos[2] += com[2];
    }

    /* free our temporary array */
    free(new_coord_array);
}

/* perform a 1D translation without periodic boundaries */
void displace_1D(system_t *system, molecule_t *molecule, double scale) {
    atom_t *atom_ptr;
    double trans;

    trans = scale * get_rand(system);
    if (get_rand(system) < 0.5) trans *= -1.0;
    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next)
        atom_ptr->pos[0] += trans;

    molecule->com[0] += trans;
}

/* perform a random translation/rotation of a molecule */
void displace(system_t *system, molecule_t *molecule, pbc_t *pbc, double trans_scale, double rot_scale) {
    translate(system, molecule, pbc, trans_scale);
    rotate(system, molecule, pbc, rot_scale);
}

/* perform a perturbation to a gaussian width */
void displace_gwp(system_t *system, molecule_t *molecule, double scale) {
    atom_t *atom_ptr;
    double perturb;

    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        if (atom_ptr->gwp_spin) {
            perturb = scale * (get_rand(system) - 0.5);
            atom_ptr->gwp_alpha += perturb;
            /* make sure the width remains positive */
            atom_ptr->gwp_alpha = fabs(atom_ptr->gwp_alpha);
        }
    }
}

/* perform a random translation/rotation of a molecule for SPECTRE algorithm */
void spectre_displace(system_t *system, molecule_t *molecule, double trans_scale, double max_charge, double max_target) {
    atom_t *atom_ptr;
    int p;
    double trans[3];
    double delta_charge;

    /* randomly translate */
    for (p = 0; p < 3; p++) {
        trans[p] = trans_scale * get_rand(system) * max_target;
        if (get_rand(system) < 0.5)
            trans[p] *= -1.0;
    }

    for (atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
        /* position reassignment */
        for (p = 0; p < 3; p++)
            atom_ptr->pos[p] += trans[p];

        /* charge reassignment */
        do {
            delta_charge = get_rand(system);
            if (get_rand(system) < 0.5) delta_charge *= -1.0;
        } while (fabs(atom_ptr->charge + delta_charge) > max_charge);
        atom_ptr->charge += delta_charge;
    }

    /* restrict the SPECTRE charges to the domain */
    spectre_wrapall(system);

    /* renormalize all SPECTRE charges */
    spectre_charge_renormalize(system);
}

/* ensure neutrality for the SPECTRE system */
void spectre_charge_renormalize(system_t *system) {
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    int num_spectre;
    double residual_charge, frac_charge;

    /* calculate the net charge */
    num_spectre = 0;
    residual_charge = 0;
    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
            if (atom_ptr->spectre) {
                ++num_spectre;
                residual_charge += atom_ptr->charge;
            }
        }
    }

    /* spread that charge out among the other SPECTRE sites*/
    frac_charge = -1.0 * residual_charge / ((double)num_spectre);

    for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
        for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
            if (atom_ptr->spectre) atom_ptr->charge += frac_charge;
}

/* apply what was already determined in checkpointing */
void make_move(system_t *system) {
    int i, j, k, p, q;
    cavity_t *cavities_array;
    int cavities_array_counter, random_index;
    double com[3], rand[3];
    molecule_t *molecule_ptr;
    atom_t *atom_ptr;
    pair_t *pair_ptr;

    /* update the cavity grid prior to making a move */
    if (system->cavity_bias) {
        cavity_update_grid(system);
        system->checkpoint->biased_move = 0;
    }

    switch (system->checkpoint->movetype) {
        case MOVETYPE_INSERT: /* insert a molecule at a random pos and orientation */
            /* umbrella sampling */
            if (system->cavity_bias && system->cavities_open) {
                /* doing a biased move - this flag lets mc.c know about it */
                system->checkpoint->biased_move = 1;
                /* make an array of possible insertion points */
                cavities_array = calloc(system->cavities_open, sizeof(cavity_t));
                memnullcheck(cavities_array, system->cavities_open * sizeof(cavity_t), __LINE__ - 1, __FILE__);
                for (i = 0, cavities_array_counter = 0; i < system->cavity_grid_size; i++) {
                    for (j = 0; j < system->cavity_grid_size; j++) {
                        for (k = 0; k < system->cavity_grid_size; k++) {
                            if (!system->cavity_grid[i][j][k].occupancy) {
                                for (p = 0; p < 3; p++)
                                    cavities_array[cavities_array_counter].pos[p] = system->cavity_grid[i][j][k].pos[p];
                                ++cavities_array_counter;
                            }
                        } /* end k */
                    }     /* end j */
                }         /* end i */
                /* insert randomly at one of the free cavity points */
                random_index = (system->cavities_open - 1) - (int)rint(((double)(system->cavities_open - 1)) * get_rand(system));
                for (p = 0; p < 3; p++)
                    com[p] = cavities_array[random_index].pos[p];
                /* free the insertion array */
                free(cavities_array);
            }  // end umbrella

            else {
                /* insert the molecule to a random location within the unit cell */
                for (p = 0; p < 3; p++)
                    rand[p] = 0.5 - get_rand(system);
                for (p = 0; p < 3; p++)
                    for (q = 0, com[p] = 0; q < 3; q++)
                        com[p] += system->pbc->basis[q][p] * rand[q];
            }

            /* process the inserted molecule */
            for (atom_ptr = system->checkpoint->molecule_backup->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                /* move the molecule back to the origin and then assign it to com */
                for (p = 0; p < 3; p++)
                    atom_ptr->pos[p] += com[p] - system->checkpoint->molecule_backup->com[p];
            }

            /* update the molecular com */
            for (p = 0; p < 3; p++)
                system->checkpoint->molecule_backup->com[p] = com[p];
            /* give it a random orientation */
            rotate(system, system->checkpoint->molecule_backup, system->pbc, 1.0);

            // insert into the list
            if (system->num_insertion_molecules) {
                // If inserting a molecule from an insertion list, we will always insert at the end
                system->checkpoint->head->next = system->checkpoint->molecule_backup;
                system->checkpoint->molecule_backup->next = NULL;
            } else {
                if (!system->checkpoint->head) {  // if we're at the start of the list:
                    system->molecules = system->checkpoint->molecule_backup;
                } else {
                    system->checkpoint->head->next = system->checkpoint->molecule_backup;
                }
                system->checkpoint->molecule_backup->next = system->checkpoint->molecule_altered;
            }

            /* set new altered and tail to reflect the insertion */
            system->checkpoint->molecule_altered = system->checkpoint->molecule_backup;
            system->checkpoint->tail = system->checkpoint->molecule_altered->next;
            system->checkpoint->molecule_backup = NULL;

            if (system->num_insertion_molecules) {  //multi sorbate
                // Free all pair memory in the list
                for (molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
                    for (atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
                        pair_ptr = atom_ptr->pairs;
                        while (pair_ptr) {
                            pair_t *temp = pair_ptr;
                            pair_ptr = pair_ptr->next;
                            free(temp);
                        }
                    }
                }
                // Generate new pairs lists for all atoms in system
                setup_pairs(system);
            }  // only one sorbate
            else
                update_pairs_insert(system);

            //reset atom and molecule id's
            enumerate_particles(system);

            break;
        case MOVETYPE_REMOVE: /* remove a randomly chosen molecule */

            if (system->cavity_bias) {
                if (get_rand(system) < pow((1.0 - system->avg_observables->cavity_bias_probability),
                                           ((double)system->cavity_grid_size * system->cavity_grid_size * system->cavity_grid_size)))
                    system->checkpoint->biased_move = 0;
                else
                    system->checkpoint->biased_move = 1;
            }

            /* remove 'altered' from the list */
            if (!system->checkpoint->head) { /* handle the case where we're removing from the start of the list */
                system->checkpoint->molecule_altered = system->molecules;
                system->molecules = system->molecules->next;
            } else {
                system->checkpoint->head->next = system->checkpoint->tail;
            }
            free_molecule(system, system->checkpoint->molecule_altered);
            system->checkpoint->molecule_altered = NULL; /* Insurance against memory errors */
            update_pairs_remove(system);

            //reset atom and molecule id's
            enumerate_particles(system);

            break;
        case MOVETYPE_DISPLACE:

            /* change coords of 'altered' */
            if (system->rd_anharmonic)
                displace_1D(system, system->checkpoint->molecule_altered, system->move_factor);
            else if (system->spectre)
                spectre_displace(system, system->checkpoint->molecule_altered, system->move_factor,
                                 system->spectre_max_charge, system->spectre_max_target);
            else if (system->gwp) {
                if (system->checkpoint->molecule_altered->atoms->gwp_spin) {
                    displace(system, system->checkpoint->molecule_altered, system->pbc, system->gwp_probability, system->rot_factor);
                    displace_gwp(system, system->checkpoint->molecule_altered, system->gwp_probability);
                } else
                    displace(system, system->checkpoint->molecule_altered, system->pbc, system->move_factor, system->rot_factor);
            } else
                displace(system, system->checkpoint->molecule_altered, system->pbc, system->move_factor, system->rot_factor);

            break;
        case MOVETYPE_ADIABATIC:
            /* change coords of 'altered' */
            displace(system, system->checkpoint->molecule_altered, system->pbc, system->adiabatic_probability, 1.0);

            break;
        case MOVETYPE_SPINFLIP:

            if (get_rand(system) < 0.5)
                system->checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_PARA;
            else
                system->checkpoint->molecule_altered->nuclear_spin = NUCLEAR_SPIN_ORTHO;

            break;
        case MOVETYPE_VOLUME:

            volume_change(system);  // I don't want to contribute to the god damned mess -- kmclaugh

            break;
        default:
            error(
                "MC_MOVES: invalid mc move\n");
            die(-1);
    }

    return;
}

/* this function (a) undoes what make_move() did and (b) determines the next move sequence by calling checkpoint() */
void restore(system_t *system) {
    // restore the remaining observables
    memcpy(system->observables, system->checkpoint->observables, sizeof(observables_t));

    /* restore state by undoing the steps of make_move() */
    switch (system->checkpoint->movetype) {
        case MOVETYPE_INSERT:

            /* take altered out of the list */
            if (!system->checkpoint->head) { /* handle the case where we inserted at the head of the list */
                system->molecules = system->molecules->next;
            } else {
                system->checkpoint->head->next = system->checkpoint->tail;
            }
            unupdate_pairs_insert(system);
            free_molecule(system, system->checkpoint->molecule_altered);
            system->checkpoint->molecule_altered = NULL; /* Insurance against memory errors */

            //reset atom and molecule id's
            enumerate_particles(system);

            break;
        case MOVETYPE_REMOVE:

            /* put backup back into the list */
            if (!system->checkpoint->head) {
                system->molecules = system->checkpoint->molecule_backup;
            } else {
                system->checkpoint->head->next = system->checkpoint->molecule_backup;
            }
            system->checkpoint->molecule_backup->next = system->checkpoint->tail;
            unupdate_pairs_remove(system);
            system->checkpoint->molecule_backup = NULL;

            //reset atom and molecule id's
            enumerate_particles(system);

            break;
        case MOVETYPE_VOLUME:

            revert_volume_change(system);

            break;
        default:

            /* link the backup into the working list again */
            if (!system->checkpoint->head)
                system->molecules = system->checkpoint->molecule_backup;
            else
                system->checkpoint->head->next = system->checkpoint->molecule_backup;
            system->checkpoint->molecule_backup->next = system->checkpoint->tail;
            free_molecule(system, system->checkpoint->molecule_altered);
            system->checkpoint->molecule_altered = NULL; /* Insurance against memory errors */
            system->checkpoint->molecule_backup = NULL;
    }

    /* renormalize charges */
    if (system->spectre) spectre_charge_renormalize(system);

    /* establish the previous checkpoint again */
    checkpoint(system);

    return;
}
