import random

class Population:
    def __init__(self, num_chroms, ordering_len, elite_size, sms_ref_csf=None,
                 sms_ref_ordering=None, restart_filename=None):
        self.num_chroms = num_chroms
        self.ordering_len = ordering_len
        self.elite_size = elite_size
        # sms_ref_csf assumes the stepvector notation
        self.sms_ref_provided = sms_ref_csf is not None and sms_ref_ordering is not None
        if self.sms_ref_provided:
            self.sms_mapping_dict = {k:v for k,v in zip(sms_ref_ordering, sms_ref_csf)}
        else:
            self.sms_mapping_dict = None
        if restart_filename is not None:
            self.current_pop = self.read_population(restart_filename)
        else:
            self.current_pop = self.gen_pop_random()

    def gen_pop_random(self):
        """
        Return a list of num_chroms random orderings.
        """
        population = []
        while len(population) < self.num_chroms:
            chrom = tuple(random.sample(range(1, self.ordering_len + 1),
                                        self.ordering_len))
            if ((self.sms_ref_provided and is_csf_valid(chrom, self.sms_mapping_dict))
                or not self.sms_ref_provided):
                population.append(chrom)

        return population

    def mutate(self, mutation_rate):
# -------------------- Original mutation --------------------------------------#
        # mutate only non-elite chromosomes
        for i in range(self.elite_size, self.num_chroms):
            for gene1 in range(0, self.ordering_len):
                if random.random() < mutation_rate:
                    mutated_chrom = list(self.current_pop[i])
                    gene2 = random.sample(range(0, self.ordering_len), 1)[0]
                    mutated_chrom[gene1], mutated_chrom[gene2]\
                  = mutated_chrom[gene2], mutated_chrom[gene1]

                    # prevent invalid CSF
                    if self.sms_ref_provided:
                        while not is_csf_valid(mutated_chrom, self.sms_mapping_dict):
                            mutated_chrom = list(self.current_pop[i])
                            gene2 = random.sample(range(0, self.ordering_len), 1)[0]
                            mutated_chrom[gene1], mutated_chrom[gene2]\
                          = mutated_chrom[gene2], mutated_chrom[gene1]
                    self.current_pop[i] = tuple(mutated_chrom)

# -------------------- Cluster mutation --------------------------------------#
#        for i in range(self.elite_size, self.num_chroms):
#            chrom = list(self.current_pop[i])
#
#            num_clusters = random.randint(2,self.ordering_len)
#
#            breakpoints = [0]
#            breakpoints.extend(sorted(random.sample(range(1, self.ordering_len), num_clusters - 1)))
#            breakpoints.append(self.ordering_len)
#
#            result = []
#            clusters = []
#
#            for j in range(len(breakpoints) - 1):
#                start, end = breakpoints[j], breakpoints[j + 1]
#
#                cluster = chrom[start:end]
#                if random.random() < 0.10:
#                    cluster = cluster[::-1]
#                clusters.append(cluster)
#
#            for cluster_idx1 in range(len(clusters)):
#                if random.random() < mutation_rate:
#                    cluster_idx2 = random.sample(range(len(clusters)), 1)[0]
#                    clusters[cluster_idx1], clusters[cluster_idx2]\
#                  = clusters[cluster_idx2], clusters[cluster_idx1]
#                    
#            for cluster in clusters:
#                result.extend(cluster)
#
#            if self.sms_ref_provided:
#                while not is_csf_valid(tuple(result), self.sms_mapping_dict):
#                    result = []
#                    clusters = []
#
#                    for j in range(len(breakpoints) - 1):
#                        start, end = breakpoints[j], breakpoints[j + 1]
#
#                        cluster = chrom[start:end]
#                        if random.random() < 0.10:
#                            cluster = cluster[::-1]
#                        clusters.append(cluster)
#
#                    for cluster_idx1 in range(len(clusters)):
#                        if random.random() < mutation_rate:
#                            cluster_idx2 = random.sample(range(len(clusters)), 1)[0]
#                            clusters[cluster_idx1], clusters[cluster_idx2]\
#                          = clusters[cluster_idx2], clusters[cluster_idx1]
#                            
#                    for cluster in clusters:
#                        result.extend(cluster)
#
#            self.current_pop[i] = tuple(result)
# ----------------------------------------------------------------------------#

    def get_elite_chromosomes(self, reduced_fitness):
        """
        From the current population, get a list of chromosomes with the highest
        fitness.
        """
        elites = []
        sorted_chromosomes = sorted(self.current_pop,
                key=lambda x: reduced_fitness[x], reverse=True)

        self.current_pop = sorted_chromosomes

        for i in range(len(self.current_pop)):
            if self.current_pop[i] not in elites:
                elites.append(self.current_pop[i])
                if len(elites) == self.elite_size:
                    break
        return elites
    
    def gen_mating_pool(self, reduced_fitness, sampling_function):
        """
        Generate a pool of chromosomes for crossover.
        """
        pool = []
        for _ in range(0, self.num_chroms - self.elite_size):
            pool.append(sampling_function(self.current_pop, reduced_fitness))
        return pool
    
    def next_generation(self, reduced_fitness, crossover_function,
                        sampling_function, mutation_rate):
        elites = self.get_elite_chromosomes(reduced_fitness)
        pool = self.gen_mating_pool(reduced_fitness, sampling_function)
        offspring = self.perform_crossover(pool, crossover_function)
        self.current_pop = elites + offspring
        self.mutate(mutation_rate)

    def perform_crossover(self, pool, crossover_function):
        """
        Perform crossover on the given population (list of tuples).
        """
        random.shuffle(pool) # Is it redundant?

        pool = [list(individual) for individual in pool]

        offsprings = []

        for i in range(0, len(pool)):
            crossover_counter = 0
            max_attempts = 10
            while crossover_counter < max_attempts:
                offspring = crossover_function(pool[i], pool[len(pool)-i-1])
                if ((self.sms_ref_provided
                     and is_csf_valid(offspring, self.sms_mapping_dict))
                    or not self.sms_ref_provided):
                    offsprings.append(tuple(offspring))
                    break
                crossover_counter += 1
            if crossover_counter == max_attempts:
#                print('Crossover failed after {} attempts'.format(max_attempts))
#                print('Create a new offspring randomly')
                offspring = random.sample(range(1, self.ordering_len + 1),
                                          self.ordering_len)
                while not is_csf_valid(offspring, self.sms_mapping_dict):
                    offspring = random.sample(range(1, self.ordering_len + 1),
                                              self.ordering_len)
                offsprings.append(tuple(offspring))

        return offsprings

    def read_population(self, filename):
        """
        Initialize population from a file.
        
        File format:
        - Lines beginning with '#' are ignored as comments
        - Each line contains a chromosome in tuple format e.g. (1,3,2,5,...)
        - The last column contain the fitness score (which is ignored)
        """
        chromosomes = []
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                    
                # Extract tuple part from the line
                tuple_part = line.split('(')[1].split(')')[0] if '(' in line else line
                
                # Parse chromosome values, handling various formats
                if ',' in tuple_part:
                    # Handle tuple format (1,2,3,...)
                    values = [int(x.strip()) for x in tuple_part.split(',')]
                else:
                    # Handle space-separated format
                    values = [int(x) for x in tuple_part.split()]

                chromosome = tuple(values)
                    
                # If the last value is significantly different, it might be the fitness score
                if len(chromosome) != self.ordering_len:
                    raise ValueError(f"Chromosome length {len(chromosome)} "
                                     f"doesn't match expected length {self.ordering_len}")
                if self.sms_ref_provided and not is_csf_valid(chromosome, self.sms_mapping_dict):
                    raise ValueError(f"Chromosome {chromosome} is not CSF compatible")
                    
                chromosomes.append(chromosome)
        
        # Sanity checks
        if len(chromosomes) != self.num_chroms:
            raise ValueError(f"Number of chromosomes in file ({len(chromosomes)}) "
                             f"doesn't match expected number ({self.num_chroms})")
        
        return chromosomes

#-------------------------------------------------------------------------------

def is_csf_valid(ordering, sms_mapping_dict):
    cumul_spin = 0
    for i in ordering:
        if sms_mapping_dict[i] == 1:
            cumul_spin += 1
        elif sms_mapping_dict[i] == 2:
            cumul_spin -= 1

        if cumul_spin < 0:
            return False
    return True
