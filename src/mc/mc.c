/* 

Space Research Group
Department of Chemistry
University of South Florida

*/

/* note on multiple sorbate systems:

N1 -> N1 + 1
prob(N1,N2) = exp(Beta*(mu1*N1+mu2*N2))V^N f(E) / (LAMBDA^3N * N1! * N2!)
bias = 1/2 (when inserting, 50% it's an N1)

N1' + 1 -> N1'
prob(N1'+1,N2) = exp(Beta*(mu1*(N1'+1)+mu2*N2))V^(N'+1) f(E) / (LAMBDA^3(N'+1) * (N1'+1)! * N2!)
bias = (N1' + 1)/(N'+1) (we pick a molecule randomly from the system, then delete)

following receipe in Frenkel p 129.
insertion: 
substitute N = N'-1
acc_rate = 2.0 * exp(beta*mu1)V/(LAMBDA^3 * N)
deletion:
substitute N' = N+1
acc_rate = 0.5 * exp(-beta*mu1) * (LAMBDA^3 * (N+1))/V

This result was verified via ideal gas (non-interacting system). See tools/ideal_GCMC.c.
--kmclaugh 2012-10-17
*/

#include <mc.h>
#ifdef MPI
#include <mpi.h>
#endif

/* the prime quantity of interest */
void boltzmann_factor(system_t *system, double initial_energy, double final_energy, double rot_partfunc) {
    double delta_energy;
    double v_new, v_old;
    double fugacity;

    delta_energy = final_energy - initial_energy;

    switch (system->ensemble) {
        case ENSEMBLE_UVT:
            //obtain the correct fugacity value
            if (system->h2_fugacity || system->co2_fugacity || system->ch4_fugacity || system->n2_fugacity)
                fugacity = system->fugacities[0];
            else if (system->user_fugacities)
                fugacity = system->fugacities[system->sorbateInsert];
            else
                fugacity = system->pressure;
            /* if biased_move not set, no cavity available so do normal evaluation below */
            if (system->cavity_bias && system->checkpoint->biased_move) {
                /* modified metropolis function */
                switch (system->checkpoint->movetype) {
                    case MOVETYPE_INSERT:
                        system->nodestats->boltzmann_factor =
                            (system->cavity_volume * system->avg_nodestats->cavity_bias_probability * fugacity * ATM2REDUCED /
                             (system->temperature * system->observables->N)) *
                            exp(-delta_energy / system->temperature) * system->sorbateCount;
                        break;
                    case MOVETYPE_REMOVE:
                        system->nodestats->boltzmann_factor =
                            (system->temperature * (system->observables->N + 1.0)) /
                            (system->cavity_volume * system->avg_nodestats->cavity_bias_probability * fugacity * ATM2REDUCED) *
                            exp(-delta_energy / system->temperature) / system->sorbateCount;
                        break;
                    case MOVETYPE_DISPLACE:
                        system->nodestats->boltzmann_factor = exp(-delta_energy / system->temperature);
                        break;
                    default:
                        error(
                            "MC: invalid mc move (not implemented for biased moves?)\n");
                        die(-1);
                }  //mc move switch
            }      //end biased move
            else {
                switch (system->checkpoint->movetype) {
                    case MOVETYPE_INSERT:
                        system->nodestats->boltzmann_factor =
                            system->pbc->volume * fugacity * ATM2REDUCED / (system->temperature * (double)(system->observables->N)) *
                            exp(-delta_energy / system->temperature) *
                            (double)(system->sorbateCount);  //for add/delete bias
                        break;
                    case MOVETYPE_REMOVE:
                        system->nodestats->boltzmann_factor =
                            system->temperature * ((double)(system->observables->N) + 1.0) / (system->pbc->volume * fugacity * ATM2REDUCED) *
                            exp(-delta_energy / system->temperature) /
                            (double)(system->sorbateCount);  //for add/delete bias
                        break;
                    case MOVETYPE_DISPLACE:
                        system->nodestats->boltzmann_factor = exp(-delta_energy / system->temperature);
                        break;
                    case MOVETYPE_SPINFLIP:
                        system->nodestats->boltzmann_factor = rot_partfunc * exp(-delta_energy / system->temperature);
                        break;
                    default:
                        error(
                            "MC: invalid mc move (not implemented for binary mixtures?)\n");
                        die(-1);
                }   //MC move switch
            }       //end biased or not?
            break;  //end UVT

        case ENSEMBLE_NVT:
            switch (system->checkpoint->movetype) {
                case MOVETYPE_SPINFLIP:
                    system->nodestats->boltzmann_factor = rot_partfunc * exp(-delta_energy / system->temperature);
                    break;
                default: /*DISPLACE*/
                    system->nodestats->boltzmann_factor = exp(-delta_energy / system->temperature);
            }
            break;  //end NVT

        case ENSEMBLE_NPT:
            //either displace or change volume
            switch (system->checkpoint->movetype) {
                case MOVETYPE_VOLUME:

                    v_old = system->checkpoint->observables->volume;
                    v_new = system->observables->volume;
                    system->nodestats->boltzmann_factor =
                        exp(-(delta_energy + system->pressure * ATM2REDUCED * (v_new - v_old) - (system->observables->N + 1) * system->temperature * log(v_new / v_old)) / system->temperature);
                    break;
                default: /*displace*/
                    system->nodestats->boltzmann_factor = exp(-delta_energy / system->temperature);
            }
            break;

        case ENSEMBLE_NVE:
            system->nodestats->boltzmann_factor = pow((system->total_energy - final_energy), 3.0 * system->N / 2.0);
            system->nodestats->boltzmann_factor /= pow((system->total_energy - initial_energy), 3.0 * system->N / 2.0);
            break;

        default:
            error(
                "MC: invalid ensemble. aborting.\n");
            die(-1);
    }

    return;
}

/* keep track of which specific moves were accepted */
void register_accept(system_t *system) {
    ++system->nodestats->accept;
    switch (system->checkpoint->movetype) {
        case MOVETYPE_INSERT:
            ++system->nodestats->accept_insert;
            break;
        case MOVETYPE_REMOVE:
            ++system->nodestats->accept_remove;
            break;
        case MOVETYPE_DISPLACE:
            ++system->nodestats->accept_displace;
            break;
        case MOVETYPE_ADIABATIC:
            ++system->nodestats->accept_adiabatic;
            break;
        case MOVETYPE_SPINFLIP:
            ++system->nodestats->accept_spinflip;
            break;
        case MOVETYPE_VOLUME:
            ++system->nodestats->accept_volume;
            break;
    }
}

/* keep track of which specific moves were rejected */
void register_reject(system_t *system) {
    ++system->nodestats->reject;
    switch (system->checkpoint->movetype) {
        case MOVETYPE_INSERT:
            ++system->nodestats->reject_insert;
            break;
        case MOVETYPE_REMOVE:
            ++system->nodestats->reject_remove;
            break;
        case MOVETYPE_DISPLACE:
            ++system->nodestats->reject_displace;
            break;
        case MOVETYPE_ADIABATIC:
            ++system->nodestats->reject_adiabatic;
            break;
        case MOVETYPE_SPINFLIP:
            ++system->nodestats->reject_spinflip;
            break;
        case MOVETYPE_VOLUME:
            ++system->nodestats->reject_volume;
            break;
    }
}

/* implements the Markov chain */
int mc(system_t *system) {
    int j, msgsize;
    double initial_energy, final_energy, current_energy;
    double rot_partfunc;
    observables_t *observables_mpi;
    avg_nodestats_t *avg_nodestats_mpi;
    sorbateInfo_t *sinfo_mpi = 0;
    double *temperature_mpi = 0;
    char *snd_strct = 0, *rcv_strct = 0;
    system->count_autorejects = 0;  // initialize counter for skipped close (unphysical) contacts
                                    // char linebuf[MAXLINE];  (unused variable)
#ifdef MPI
    MPI_Datatype msgtype;
#endif /* MPI */

    /* allocate the statistics structures */
    observables_mpi = calloc(1, sizeof(observables_t));
    memnullcheck(observables_mpi, sizeof(observables_t), __LINE__ - 1, __FILE__);
    avg_nodestats_mpi = calloc(1, sizeof(avg_nodestats_t));
    memnullcheck(avg_nodestats_mpi, sizeof(avg_nodestats_t), __LINE__ - 1, __FILE__);
    // if multiple-sorbates, allocate sorb statistics struct
    if (system->sorbateCount > 1) {
        sinfo_mpi = calloc(system->sorbateCount, sizeof(sorbateInfo_t));
        memnullcheck(sinfo_mpi, sizeof(sorbateInfo_t), __LINE__ - 1, __FILE__);
        system->sorbateGlobal = calloc(system->sorbateCount, sizeof(sorbateAverages_t));
        memnullcheck(system->sorbateGlobal, sizeof(sorbateAverages_t), __LINE__ - 1, __FILE__);
    }

    // compute message size
    msgsize = sizeof(observables_t) + sizeof(avg_nodestats_t);
    if (system->calc_hist) msgsize += system->n_histogram_bins * sizeof(int);
    if (system->sorbateCount > 1) msgsize += system->sorbateCount * sizeof(sorbateInfo_t);

#ifdef MPI
    MPI_Type_contiguous(msgsize, MPI_BYTE, &msgtype);
    MPI_Type_commit(&msgtype);
#endif /* MPI */

    /* allocate MPI structures */
    snd_strct = calloc(msgsize, 1);
    memnullcheck(snd_strct, sizeof(msgsize), __LINE__ - 1, __FILE__);
    if (!rank) {
        rcv_strct = calloc(size, msgsize);
        memnullcheck(rcv_strct, size * sizeof(msgsize), __LINE__ - 1, __FILE__);
        temperature_mpi = calloc(size, sizeof(double));  //temperature list for parallel tempering
        memnullcheck(temperature_mpi, size * sizeof(double), __LINE__ - 1, __FILE__);
    }

    /* update the grid for the first time */
    if (system->cavity_bias) cavity_update_grid(system);

    /* set volume observable */
    system->observables->volume = system->pbc->volume;

    /* get the initial energy of the system */
    initial_energy = energy(system);

#ifdef QM_ROTATION
    /* solve for the rotational energy levels */
    if (system->quantum_rotation) quantum_system_rotational_energies(system);
#endif /* QM_ROTATION */

    /* be a bit forgiving of the initial state */
    if (!isfinite(initial_energy))
        initial_energy = system->observables->energy = MAXVALUE;

    /* if root, open necessary output files */
    if (!rank)
        if (open_files(system) < 0) {
            error(
                "MC: could not open files\n");
            return (-1);
        }

    // write initial observables to stdout and logs
    if (!rank) {
        calc_system_mass(system);
        // average in the initial values once  (we don't want to double-count the initial state when using MPI)
        update_root_averages(system, system->observables, system->avg_observables);
        // average in the initial sorbate values
        if (system->sorbateCount > 1) {
            update_sorbate_info(system);                             //local update
            update_root_sorb_averages(system, system->sorbateInfo);  //global update
        }
        // write initial observables exactly once
        if (system->file_pointers.fp_energy)
            write_observables(system->file_pointers.fp_energy, system, system->observables, system->temperature);
        if (system->file_pointers.fp_energy_csv)
            write_observables_csv(system->file_pointers.fp_energy_csv, system, system->observables, system->temperature);
        output(
            "MC: initial values:\n");
        write_averages(system);
    }

    /* save the initial state */
    checkpoint(system);

    /* main MC loop */
    for (system->step = 1; system->step <= system->numsteps; (system->step)++) {
        /* restore the last accepted energy */
        initial_energy = system->observables->energy;

        /* perturb the system */
        make_move(system);

        /* calculate the energy change */
        final_energy = energy(system);

#ifdef QM_ROTATION
        /* solve for the rotational energy levels */
        if (system->quantum_rotation && (system->checkpoint->movetype == MOVETYPE_SPINFLIP))
            quantum_system_rotational_energies(system);
#endif /* QM_ROTATION */
        if (system->checkpoint->movetype != MOVETYPE_REMOVE)
            rot_partfunc = system->checkpoint->molecule_altered->rot_partfunc;
        else
            rot_partfunc = system->checkpoint->molecule_backup->rot_partfunc;

        /* treat a bad contact as a reject */
        if (!isfinite(final_energy)) {
            system->observables->energy = MAXVALUE;
            system->nodestats->boltzmann_factor = 0;
        } else
            boltzmann_factor(system, initial_energy, final_energy, rot_partfunc);

        /* Metropolis function */
        if ((get_rand(system) < system->nodestats->boltzmann_factor) && (system->iter_success == 0)) {
            /////////// ACCEPT

            current_energy = final_energy;

            /* checkpoint */
            checkpoint(system);
            register_accept(system);

            /* SA */
            if (system->simulated_annealing) {
                if (system->simulated_annealing_linear == 1) {
                    system->temperature = system->temperature + (system->simulated_annealing_target - system->temperature) / (system->numsteps - system->step);
                    if (system->numsteps - system->step == 0)
                        system->temperature = system->simulated_annealing_target;
                } else
                    system->temperature = system->simulated_annealing_target + (system->temperature - system->simulated_annealing_target) * system->simulated_annealing_schedule;
            }

        } else {
            /////////////// REJECT

            current_energy = initial_energy;  //used in parallel tempering

            //reset the polar iterative failure flag
            system->iter_success = 0;

            /* restore from last checkpoint */
            restore(system);
            register_reject(system);

        }  // END REJECT

        // perform parallel_tempering
        if ((system->parallel_tempering) && (system->step % system->ptemp_freq == 0))
            temper_system(system, current_energy);

        /* track the acceptance_rate */
        track_ar(system->nodestats);

        /* each node calculates its stats */
        update_nodestats(system->nodestats, system->avg_nodestats);

        /* do this every correlation time, and at the very end */
        if (!(system->step % system->corrtime) || (system->step == system->numsteps)) {
            /* copy observables and avgs to the mpi send buffer */
            /* histogram array is at the end of the message */
            if (system->calc_hist) {
                zero_grid(system->grids->histogram->grid, system);
                population_histogram(system);
            }

            // update frozen and total system mass
            calc_system_mass(system);

            // update sorbate info on each node
            if (system->sorbateCount > 1)
                update_sorbate_info(system);

/*write trajectory files for each node -> one at a time to avoid disk congestion*/
#ifdef MPI
            for (j = 0; j < size; j++) {
                MPI_Barrier(MPI_COMM_WORLD);
                if (j == rank) write_states(system);
            }
#else
            write_states(system);
#endif

            /*restart files for each node -> one at a time to avoid disk congestion*/
            if (write_molecules_wrapper(system, system->pqr_restart) < 0) {
                error(
                    "MC: could not write restart state to disk\n");
                return (-1);
            }

/*dipole/field data for each node -> one at a time to avoid disk congestion*/
#ifdef MPI
            if (system->polarization) {
                for (j = 0; j < size; j++) {
                    MPI_Barrier(MPI_COMM_WORLD);
                    if (j == rank) {
                        write_dipole(system);
                        write_field(system);
                    }
                }
            }
#else
            if (system->polarization) {
                write_dipole(system);
                write_field(system);
            }
#endif

            /* zero the send buffer */
            memset(snd_strct, 0, msgsize);
            memcpy(snd_strct, system->observables, sizeof(observables_t));
            memcpy(snd_strct + sizeof(observables_t), system->avg_nodestats, sizeof(avg_nodestats_t));
            if (system->calc_hist)
                mpi_copy_histogram_to_sendbuffer(snd_strct + sizeof(observables_t) + sizeof(avg_nodestats_t),
                                                 system->grids->histogram->grid, system);
            if (system->sorbateCount > 1)
                memcpy(snd_strct + sizeof(observables_t) + sizeof(avg_nodestats_t) + (system->calc_hist) * system->n_histogram_bins * sizeof(int),  //compensate for the size of hist data, if neccessary
                       system->sorbateInfo,
                       system->sorbateCount * sizeof(sorbateInfo_t));

            if (!rank) memset(rcv_strct, 0, size * msgsize);

#ifdef MPI
            MPI_Gather(snd_strct, 1, msgtype, rcv_strct, 1, msgtype, 0, MPI_COMM_WORLD);
            MPI_Gather(&(system->temperature), 1, MPI_DOUBLE, temperature_mpi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//need to gather shit for sorbate stats also
#else
            memcpy(rcv_strct, snd_strct, msgsize);
            temperature_mpi[0] = system->temperature;
#endif /* MPI */
            /* head node collects all observables and averages */
            if (!rank) {
                /* clear avg_nodestats to avoid double-counting */
                clear_avg_nodestats(system);
                //loop for each core -> shift data into variable_mpi, then average into avg_observables
                for (j = 0; j < size; j++) {
                    /* copy from the mpi buffer */
                    memcpy(observables_mpi, rcv_strct + j * msgsize, sizeof(observables_t));
                    memcpy(avg_nodestats_mpi, rcv_strct + j * msgsize + sizeof(observables_t), sizeof(avg_nodestats_t));
                    if (system->calc_hist)
                        mpi_copy_rcv_histogram_to_data(rcv_strct + j * msgsize + sizeof(observables_t) + sizeof(avg_nodestats_t), system->grids->histogram->grid, system);
                    if (system->sorbateCount > 1)
                        memcpy(sinfo_mpi,
                               rcv_strct + j * msgsize + sizeof(observables_t) + sizeof(avg_nodestats_t) + (system->calc_hist) * system->n_histogram_bins * sizeof(int),  //compensate for the size of hist data, if neccessary
                               system->sorbateCount * sizeof(sorbateInfo_t));

                    /* write observables */
                    if (system->file_pointers.fp_energy)
                        write_observables(system->file_pointers.fp_energy, system, observables_mpi, temperature_mpi[j]);
                    if (system->file_pointers.fp_energy_csv)
                        write_observables_csv(system->file_pointers.fp_energy_csv, system, observables_mpi, temperature_mpi[j]);
                    if (system->file_pointers.fp_xyz) {
                        write_molecules_xyz(system, system->file_pointers.fp_xyz);  //L
                    }
                    /* collect the averages */
                    /* if parallel tempering, we will collect obserables from the coldest bath. this can't be done for
					 * nodestats though, since nodestats are averaged over each corrtime, rather than based on a single 
					 * taken at the corrtime */
                    update_root_nodestats(system, avg_nodestats_mpi, system->avg_observables);
                    if (!system->parallel_tempering) {
                        update_root_averages(system, observables_mpi, system->avg_observables);
                        if (system->calc_hist) update_root_histogram(system);
                        if (system->sorbateCount > 1) update_root_sorb_averages(system, sinfo_mpi);
                    } else if (system->ptemp->index[j] == 0) {
                        update_root_averages(system, observables_mpi, system->avg_observables);
                        if (system->calc_hist) update_root_histogram(system);
                        if (system->sorbateCount > 1) update_root_sorb_averages(system, sinfo_mpi);
                    }
                }

                /* write the averages to stdout */
                if (system->file_pointers.fp_histogram)
                    write_histogram(system->file_pointers.fp_histogram, system->grids->avg_histogram->grid, system);

                if (write_performance(system->step, system) < 0) {
                    error(
                        "MC: could not write performance data to stdout\n");
                    return (-1);
                }
                if (write_averages(system) < 0) {
                    error(
                        "MC: could not write statistics to stdout\n");
                    return (-1);
                }

            } /* !rank */
        }     /* corrtime */
    }         /* main loop */

    /* write output, close any open files */
    free(snd_strct);

    // restart files for each node
    if (write_molecules_wrapper(system, system->pqr_output) < 0) {
        error(
            "MC: could not write final state to disk\n");
        return (-1);
    }

    if (!rank) {
        close_files(system);
        free(rcv_strct);
        free(temperature_mpi);
    }

    if (system->sorbateCount > 1) {
        free(system->sorbateGlobal);
        free(sinfo_mpi);
    }

    free(observables_mpi);
    free(avg_nodestats_mpi);

    printf(
        "MC: Total auto-rejected moves: %i\n", system->count_autorejects);

    return (0);
}
