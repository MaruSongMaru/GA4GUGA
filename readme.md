# GA4GUGA

A genetic algorithm (GA) for finding optimal orbital orderings that yield compact wave function representations in GUGA bases.
The code reads an FCIDUMP file, performs a GA simulation using diagonal-element information (see Section 3 of Reference [1] for details of the available fitness functions), and writes a reordered FCIDUMP file using the best ordering found.

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

Once installed, you can directly run the example scripts located in the
`examples/` directory. Below is a brief description of each:

- `fast_20site_Heisenberg_chain.py`
    This example runs a GA simulation for a 20-site NN Heisenberg chain
    (Subsection 4.1 of Reference [1]) using the simplified fitness function
    with the S-Ms mapping. Only ordering-dependent terms of the CSF energy
    are evaluated for fitness measure.

- `20site_Heisenberg_chain.py`
    This example is the same as `fast_20site_Heisenberg_chain`, but evaluates
    the full CSF energies, thus slower than `fast_20site_Heisenberg_chain`.

- `P-cluster_48i40.py`
    This example runs a GA simulation for the CAS(48e,40o) P-cluster consisting
    of 8 Fe sites (Subsection 5.2.1 of Reference [1]). The fitness function 
    evaluates the CSF energives over the 14 collinear VVS CSFs and maximizes the
    highest and the lowest CSF energies.

### Checkpoint feature

```bash
touch WRITE_CHECKPOINT
```

During a GA run, creating this file triggers writing an FCIDUMP file with the
current best ordering.


## Input arguments

Required arguments (black) and optional arguments (blue) are listed below. CSFs
are assumed to be in the step-vector notation where $0$, $1$, $2$, and $3$
represent empty, increasing intermediate spin by 1/2, decreasing by 1/2,
and doubly occupied orbital couplings, respectively.

- fitness_function `measure_fitness.FitnessFunction.FITNESS_FUNCTION_NAME`
    Choose the fitness function for the GA simulation. See Section 3 of Reference [1]. Options include: DIAG_ELEM_SMS_MAPPING, MIN_MAX_DIFF, MAX_DIAG_ELEM, MIN_DIAG_ELEM, NEEL_FAST_DIAG, NEEL_FAST_DIAG_MIN, NEEL_FAST_DIAG_MAX, MAX_FAST_DIAG, MIN_FAST_DIAG, NEEL_FAST_DIAG_MIN_OSONLY.

- co_function `crossover.CO_FUNCTION_NAME`
    Crossover operator to use, e.g., `crossover.CO_FUNCION_NAME`.

- num_chroms $n$(int)
    Population size (number of chromosomes).

- elite_size $n$(int)
    Number of elite chromosomes retained each generation.

- mutation_rates $mlist$(list[float])
    One or more mutation rates. If multiple rates are provided, the algorithm can cycle through them when progress stalls.

- restricted_ordering_len $n$(int)
    Length of the effective (restricted) ordering that the GA optimizes before expansion.

- generations $n$(int)
    Number of generations to run.

- fcidump `PATH_TO_FCIDUMP`
    Path to the source FCIDUMP file.

- norb $n$(int)
    Total number of orbitals. Must satisfy the consistency check:
    `norb = num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation)`

#### Fitness-function specific arguments

- <span style="color:blue">csf_list</span> $csf_list$(list[list[int]])
    List of CSF step-vectors to evaluate. Required for certain fitness functions (e.g., MIN_MAX_DIFF, MAX_DIAG_ELEM, MIN_DIAG_ELEM, MAX_FAST_DIAG, MIN_FAST_DIAG).

- <span style="color:blue">tMinimize</span> $l$(bool)
    For DIAG_ELEM_SMS_MAPPING only: if True, uses the absolute value to target minima.

- <span style="color:blue">SMS</span> $l$(bool)

#### Restricting orderings

If `norb = 10`, `num_prefix = 3`, `num_suffix = 2`, `restricted_ordering_len = 5`, and `on_site_permutation = (1,)` (default), the GA is performed only for the middle set `{4,5,6,7,8}`. Orderings maintain the structure `(1,2,3, {4,5,6,7,8}, 9,10)`, e.g., `(1,2,3,4,6,7,5,8,9,10)` or `(1,2,3,6,8,7,4,5,9,10)`.

If `norb = 10`, `num_prefix = 2`, `num_suffix = 2`, `restricted_ordering_len = 3`, and `on_site_permutation = (2,1)`, then orderings maintain the structure `(1,2, {(4,3), (6,5), (8,7)}, 9,10)`, e.g., `(1,2,5,6,4,3,8,7,9,10)`.

- <span style="color:blue">on_site_permutation</span> $perm$(tuple[int])
    Pattern applied to each gene when expanding the restricted ordering to the full ordering. Default: `(1,)`.

- <span style="color:blue">num_prefix</span> $n$(int)
    Number of leading orbitals that remain fixed (not optimized). Default: `0`.

- <span style="color:blue">num_suffix</span> $n$(int)
    Number of trailing orbitals that remain fixed (not optimized). Default: `0`.

#### Restart

- <span style="color:blue">restart_filename</span> `pop_filename`
    If provided, restart the GA from a saved state.

#### Logging arguments

- <span style="color:blue">checkpoint_trigger</span> n(str)
    File name whose presence triggers writing a checkpoint FCIDUMP of the current best ordering during the run. Default: `"WRITE_CHECKPOINT"`.

- <span style="color:blue">checkpoint_prefix</span> n(str)
    Prefix for the checkpoint FCIDUMP file name. Default: `"FCIDUMP_checkpoint"`.

- <span style="color:blue">pop_file_name</span> n(str)
    Where to log the current population each generation. Default: `"current_pop.log"`.


## Citation

If you cite this package, please cite:

[1] M. Song and G. Li Manni, A Genetic Algorithm Approach for Compact Wave Function Representations in Spin-Adapted Bases, J. Chem. Theory Comput. **XXX**, XXXX (2025), DOI: [10.1021/acs.jctc.5c01264](https://doi.org/10.1021/acs.jctc.5c01264)
