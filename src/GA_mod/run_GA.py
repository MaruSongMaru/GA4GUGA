import numpy as np
from GA_mod import population as pop
from GA_mod import process_df
from GA_mod import sampling
from GA_mod import measure_fitness
from GA_mod import crossover as co
from FCIDUMP_tools import IntegralClass
import subprocess

def perform_GA(fitness_function, num_chroms, restricted_ordering_len, elite_size,
               mutation_rates, generations, co_function, fcidump, s, nel, norb,
               restart_filename=None, **kwargs):
    """
    restricted_ordering_len is the length of a "restricted" ordering and the GA
    is performed on chromosomes with this length. But fcidump, s, nel, and norb
    are for the full ordering.
    """

    log_file_name = kwargs.get('log_file_name', 'progress.log')
    pop_filename = kwargs.get('pop_file_name', 'current_pop.log')
    git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    sms_ref_csf = kwargs.get('sms_ref_csf', None)
    sms_ref_ordering = kwargs.get('sms_ref_ordering', None)
    cluster_period = kwargs.get('cluster_period', 5)
    stagnation_limit = kwargs.get('stagnation_limit', max(100, generations // 100))

    with open(log_file_name, 'w') as log_file:
        log_file.write("Genetic Algoritm simulation started\n")
        log_file.write("GIT hash: {}\n".format(git_hash))
        log_file.write("\n")
        if restart_filename is not None:
            log_file.write("Restarting from population file: {}\n\n"
                           .format(restart_filename))
        log_file.write("- Number of chromosomes: {}\n".format(num_chroms))
        log_file.write("- Elite size: {}\n".format(elite_size))
        log_file.write("- Mutation rates: {}\n".format(mutation_rates))
        log_file.write("- Generations: {}\n".format(generations))
        log_file.write("- Crossover function: {}\n".format(co_function.__name__))
        log_file.write("- Fitness function: {}\n".format(fitness_function))
        log_file.write("- Clustering period: {}\n".format(cluster_period))
        log_file.write("- Stagnation limit: {}\n".format(stagnation_limit))
        log_file.write("\n")

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
                                          s, nel, norb, **kwargs)
    bestchrom = max(reduced_fitness_dict, key=reduced_fitness_dict.get)

    best_fitness = reduced_fitness_dict[bestchrom]
    with open(log_file_name, 'a') as log_file:
        log_file.write("# generation  chromosome  fitness\n")
        log_file.write("0  {}  {}\n".format(bestchrom, best_fitness))

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
                                              FCIDUMPClass, s, nel, norb,
                                              **kwargs)
        bestchrom = max(reduced_fitness_dict, key=reduced_fitness_dict.get)
        best_fitness = reduced_fitness_dict[bestchrom]

        if abs(best_fitness - best_fitness_prev) < 1e-6:
            stagnation_counter += 1
            if stagnation_counter == stagnation_limit:
                mutation_rate_index = (mutation_rate_index + 1) % len(mutation_rates)
                current_mutation_rate = mutation_rates[mutation_rate_index]
                with open(log_file_name, 'a') as log_file:
                    log_file.write("# Stagnation detected. Change mutation rate to {}\n".format(current_mutation_rate))
                stagnation_counter = 0
        else:
            stagnation_counter = 0

        with open(log_file_name, 'a') as log_file:
            log_file.write("{}  {}  {}\n".format(i, bestchrom, best_fitness))

        with open(pop_filename, 'w') as log_file:
            log_file.write(f"# Chromosomes in the {i}th generation and their fitnesses\n")
            for chrom in POPClass.current_pop:
                log_file.write(f"{chrom} {reduced_fitness_dict[chrom]}\n")


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
