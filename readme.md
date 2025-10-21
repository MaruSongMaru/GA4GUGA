# GA4GUGA

A genetic algorithm designed to find an optimal ordering for a compact wave function representation in GUGA bases.

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

See `examples/` directory for complete working examples.

### Input options

- `fitness_function`: Choose the fitness function for the GA simulation
    - `...`: Reference [1]
- `num_chroms`: Number of chromosomes in the population
- `elite_size`: Number of elite chromosoms
work in progress..

### Checkpoint feature

```bash
touch WRITE_CHECKPOINT
```

creates a trigger file during a run to create an FCIDUMP file with the current best ordering.


## Citation

If you cite this package, please cite:

[1] M. Song and G. Li Manni, A Genetic Algorithm Approach for Compact Wave Function Representations in Spin-Adapted Bases, J. Chem. Theory Comput. **XXX**, XXXX (2025), DOI: [10.1021/acs.jctc.5c01264](https://doi.org/10.1021/acs.jctc.5c01264)

## License

MIT License
