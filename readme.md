# GA4GUGA

A genetic algorithm (GA) designed to find an optimal ordering for a compact wave function representation in GUGA bases.
The code reads in an FCIDUMP file, perform a GA simulation using diagonal element information (see Section 3 of Reference [1] for details of available fitness functions), and creates a reordered FCIDUMP file using the best ordering found from the simulation.

## Installation

```bash
cd GA4GUGA
pip install -e .
```

This installs the package in editable mode.

### Verify installation

```bash
python test_installation.py
```

## Quick Start

```python
from GA_mod import run_GA
from GA_mod import crossover as co
from GA_mod.measure_fitness import FitnessFunction

# Run genetic algorithm
fitness = run_GA.perform_GA(
    fitness_function=FitnessFunction.NEEL_FAST_DIAG_MIN_OSONLY,
    num_chroms=200,
    restricted_ordering_len=20,
    elite_size=20,
    mutation_rates=[0.05],
    generations=1000,
    co_function=co.order_co,
    fcidump="path/to/FCIDUMP",
    norb=20,
)
```

### Checkpoint feature

```bash
touch WRITE_CHECKPOINT
```

creates a trigger file during a run to create an FCIDUMP file with the current best ordering.


See `examples/` directory for complete working examples.

### Input options

- `fitness_function`: Choose the fitness function for the GA simulation
    - `...`: Reference [1]
- `num_chroms`: The number of chromosomes in the population
- `elite_size`: The number of elite chromosoms

#### Restricting orderings
GA simulations can be performed on restricted orderings. 
- `restricted_ordering_len`: The effective ordering length for GA
- `on_site_permutation`: Each label in the effective ordering is extended with this information.
- `num_prefix`: The number of first ordering sequence that are not changed during the GA simulation.
- `num_suffix`: The number of last ordering sequence that are not changed during the GA simulation.
- `norb`: The number of orbitals of the system. This is used for a sanity check: norb = num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation)
For example, if `norb = 10`, `num_prefix = 3`, `num_suffix = 2`, `restricted_ordering_len = 5`, and `on_site_permutation=(1,)` (default), GA is only perforemd for {4,5,6,7,8}. Thus the orderings keep this form  (1,2,3,{4,5,6,7,8},9,10).
E.g., (1,2,3,4,6,7,5,8,9,10), (1,2,3,6,8,7,4,5,9,10).
If `norb = 10`, `num_prefix = 2`, `num_suffix = 2`, `restricted_ordering_len = 3`, and `on_site_permutation=(2,1)`, then the orderings keep (1,2,{(4,3),(6,5),(8,7)},9,10) form.
E.g., (1,2,5,6,4,3,8,7,9,10).
work in progress..


## Citation

If you cite this package, please cite:

[1] M. Song and G. Li Manni, A Genetic Algorithm Approach for Compact Wave Function Representations in Spin-Adapted Bases, J. Chem. Theory Comput. **XXX**, XXXX (2025), DOI: [10.1021/acs.jctc.5c01264](https://doi.org/10.1021/acs.jctc.5c01264)
