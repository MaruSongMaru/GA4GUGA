import numpy as np
import sys
import os
from GA_mod import population as pop
from GA_mod import process_df
from GA_mod import sampling
from GA_mod import measure_fitness
from GA_mod import crossover as co
from FCIDUMP_tools import IntegralClass
import subprocess

def perform_GA(fitness_function, num_chroms, elite_size, mutation_rates, 
               restricted_ordering_len, generations, co_function, fcidump, norb,
               on_site_permutation=(1,), num_prefix=0, num_suffix=0, 
               restart_filename=None, **kwargs):
    """
    restricted_ordering_len is the length of a "restricted" ordering and the GA
    is performed on chromosomes with this length. But fcidump, and norb
    are for the full ordering.
    
    Output is written to stdout and can be redirected: python script.py > output.log
    
    Checkpoint feature: Create a file named by 'checkpoint_trigger' (default: 'WRITE_CHECKPOINT')
    during the GA run to trigger writing an FCIDUMP with the current best ordering.
    The trigger file will be deleted after the checkpoint is written.
    """

    expected_norb = num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation)
    if expected_norb != int(norb):
        raise ValueError(
            "Inconsistent norb: num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation) "
            f"= {expected_norb}, but norb = {norb}"
        )

    pop_filename = kwargs.get('pop_file_name', 'current_pop.log')
    checkpoint_trigger = kwargs.get('checkpoint_trigger', 'WRITE_CHECKPOINT')
    checkpoint_prefix = kwargs.get('checkpoint_prefix', 'FCIDUMP_checkpoint')
    git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    sms_ref_csf = kwargs.get('sms_ref_csf', None)
    sms_ref_ordering = kwargs.get('sms_ref_ordering', None)
    cluster_period = kwargs.get('cluster_period', 5)
    stagnation_limit = kwargs.get('stagnation_limit', max(100, generations // 100))

    print("Genetic Algoritm simulation started", file=sys.stdout)
    print(f"GIT hash: {git_hash}", file=sys.stdout)
    print("", file=sys.stdout)
    if restart_filename is not None:
        print(f"Restarting from population file: {restart_filename}\n", file=sys.stdout)
    print(f"- Number of chromosomes: {num_chroms}", file=sys.stdout)
    print(f"- Elite size: {elite_size}", file=sys.stdout)
    print(f"- Mutation rates: {mutation_rates}", file=sys.stdout)
    print(f"- Generations: {generations}", file=sys.stdout)
    print(f"- Crossover function: {co_function.__name__}", file=sys.stdout)
    print(f"- Fitness function: {fitness_function}", file=sys.stdout)
    print(f"- Clustering period: {cluster_period}", file=sys.stdout)
    print(f"- Stagnation limit: {stagnation_limit}", file=sys.stdout)
    if sms_ref_csf is not None and sms_ref_ordering is not None:
        print(f"- S-Ms mapping reference CSF: {sms_ref_csf}", file=sys.stdout)
        print(f"- S-Ms mapping reference ordering: {sms_ref_ordering}", file=sys.stdout)
    print("", file=sys.stdout)
    print(f"Checkpoint trigger: Create file '{checkpoint_trigger}' to write current best ordering", file=sys.stdout)
    print("", file=sys.stdout)

    POPClass = pop.Population(num_chroms, restricted_ordering_len, elite_size,
                               sms_ref_csf, sms_ref_ordering,
                               restart_filename=restart_filename)

    FCIDUMPClass = IntegralClass.FCIDUMPReader(fcidump)

    # Track best fitness and stagnation
    stagnation_counter = 0
    mutation_rate_index = 0
    current_mutation_rate = mutation_rates[mutation_rate_index]

    # 0th generation
    reduced_fitness_dict = \
        measure_fitness.calculate_fitness(fitness_function, POPClass, FCIDUMPClass,
                                          norb, on_site_permutation, num_prefix,
                                          num_suffix, **kwargs)
    bestchrom = max(reduced_fitness_dict, key=reduced_fitness_dict.get)

    best_fitness = reduced_fitness_dict[bestchrom]
    print("# Generation  Ordering  Fitness", file=sys.stdout)
    print(f"0  {bestchrom}  {best_fitness}", file=sys.stdout)

    # Subsequent generations
    for i in range(1, generations + 1):
        best_fitness_prev = best_fitness
        if i % cluster_period == 0:
            _co_function = co.shuffle_cluster
        else:
            _co_function = co_function
        POPClass.next_generation(reduced_fitness_dict, _co_function,
                                  sampling.roullette_wheel_sampling,
                                  current_mutation_rate)

        reduced_fitness_dict = \
            measure_fitness.calculate_fitness(fitness_function, POPClass,
                                              FCIDUMPClass, norb, 
                                              on_site_permutation, num_prefix,
                                              num_suffix, **kwargs)
        bestchrom = max(reduced_fitness_dict, key=reduced_fitness_dict.get)
        best_fitness = reduced_fitness_dict[bestchrom]

        if abs(best_fitness - best_fitness_prev) < 1e-6:
            stagnation_counter += 1
            if stagnation_counter == stagnation_limit:
                mutation_rate_index = (mutation_rate_index + 1) % len(mutation_rates)
                current_mutation_rate = mutation_rates[mutation_rate_index]
                print(f"# Stagnation detected. Change mutation rate to {current_mutation_rate}", file=sys.stdout)
                stagnation_counter = 0
        else:
            stagnation_counter = 0

        print(f"{i}  {bestchrom}  {best_fitness}", file=sys.stdout)

        # Check for checkpoint trigger file
        if os.path.exists(checkpoint_trigger):
            checkpoint_filename = f"{checkpoint_prefix}_gen{i}"
            print(f"\n# Checkpoint trigger detected at generation {i}", file=sys.stdout)
            print(f"# Writing FCIDUMP with current best ordering to '{checkpoint_filename}'", file=sys.stdout)
            FCIDUMPClass.dump_integrals(checkpoint_filename, bestchrom)
            print(f"# Checkpoint written. Removing trigger file and continuing...\n", file=sys.stdout)
            try:
                os.remove(checkpoint_trigger)
            except OSError as e:
                print(f"# Warning: Could not remove trigger file: {e}", file=sys.stdout)

        with open(pop_filename, 'w') as log_file:
            log_file.write(f"# Chromosomes in the {i}th generation and their fitnesses\n")
            for chrom in POPClass.current_pop:
                log_file.write(f"{chrom} {reduced_fitness_dict[chrom]}\n")

    # Generate an FCIDUMP file with the best ordering.
    print("\n\nGenetic Algorithm simulation completed.", file=sys.stdout)
    print(f"Best ordering found: {bestchrom}", file=sys.stdout)
    print("Writing reordered FCIDUMP to 'FCIDUMP_bestordering'", file=sys.stdout)
    FCIDUMPClass.dump_integrals('FCIDUMP_bestordering', bestchrom)


#------------------------------------------------------------------------------#
def perform_GA_test_use_df(num_chroms, ordering_len, elite_size, mutation_rates,
                          generations, co_function, fitness_dict,
                          reference_dict, **kwargs):
    """
    Function to perform the Genetic Algorithm test where the full data is provided.
 
    - reference_dict: Dictionary that contains the reference values to be used
                      check the target fitness (fitness_dict) is a good metric.
                      Mostly L4Norm is used as reference.
    """
    POPClass = pop.Population(num_chroms, ordering_len, elite_size)
    fitness = np.zeros(generations + 1)
    reference_fitness = np.zeros(generations + 1)
    cluster_period = kwargs.get('cluster_period', 5)
 
    # 0th generation
    reduced_fitness_dict = \
        process_df.gen_pop_fitness_ht(POPClass.current_pop, fitness_dict)
    fitness[0] = max(reduced_fitness_dict.values())
 
    reduced_reference_dict = \
        process_df.gen_pop_fitness_ht(POPClass.current_pop, reference_dict)
    reference_fitness[0] = max(reduced_reference_dict.values())
 
    # Subsequent generations
    for i in range(1, generations + 1):
        if i % cluster_period == 0:
            _co_function = co.shuffle_cluster
        else:
            _co_function = co_function
        POPClass.next_generation(reduced_fitness_dict, _co_function,
                                  sampling.roullette_wheel_sampling,
                                  mutation_rates)
 
        reduced_fitness_dict = \
            process_df.gen_pop_fitness_ht(POPClass.current_pop, fitness_dict)
        fitness[i] = max(reduced_fitness_dict.values())
 
        reduced_reference_dict = \
            process_df.gen_pop_fitness_ht(POPClass.current_pop, reference_dict)
        reference_fitness[i] = max(reduced_reference_dict.values())
 
    return fitness, reference_fitness

def GA_ensembles_read_df(num_ensembles, pop_size, restricted_ordering_len,
                         elite_size, mutation_rates, generations, co_function,
                         fitness_dict, reference_dict, **kwargs):

    fitness_ensemble_array = np.zeros((num_ensembles, generations + 1))
    reference_ensemble_array = np.zeros((num_ensembles, generations + 1))

    for ensemble in range(num_ensembles):
        fitness_ensemble_array[ensemble], reference_ensemble_array[ensemble] = \
            perform_GA_test_use_df(pop_size, restricted_ordering_len, elite_size,
                                   mutation_rates, generations, co_function,
                                   fitness_dict, reference_dict, **kwargs)

    return fitness_ensemble_array, reference_ensemble_array
