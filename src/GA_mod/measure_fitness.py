# Wrapper to measure fitness scores.
# The larger the score, the better the fitness.
# The sign of the fitness doesn't matter, as it is rescaled in sampling.py
# But when targetting a minimum, -1 has to be multipled at the end to give the
# minimum a larger score.
from GA_mod import read_fcidump
from GA_mod import population as pop
from GA_mod import gen_ref_dicts
from mypytoolbox import convert_reps
from enum import Enum, auto
from GA_mod import extend_ordering
import numpy as np
from GA_mod import GUGA_diag

class FitnessFunction(Enum):
    REF_DIAGELEM = auto()
    MIN_MAX_DIFF = auto()
    MAX_DIAG_ELEM = auto()
    MIN_DIAG_ELEM = auto()
    NEEL_FAST_DIAG = auto()
    NEEL_FAST_DIAG_MIN = auto()
    NEEL_FAST_DIAG_MAX = auto()
    MIN_FAST_DIAG = auto()
    MAX_FAST_DIAG = auto()
    # Only ud CSFs
    NEEL_FAST_DIAG_MIN_OSONLY = auto()

def _gen_ref_diagelem_fitness(population, extended_pop, ref_dict, FCIDUMPClass, s, nel, norb, tHeisenberg):
    """Internal function that calculates fitness based on reference diagonal elements."""
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        # permute CSF
        CSF_ref = [ref_dict[orb] for orb in extended_chrom]

        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)

        diagelem = GUGAClass.calc_diag_elem(CSF_ref, add_core=True)

        if tHeisenberg:
            reduced_fitness[chrom] = abs(diagelem)
        else:
            reduced_fitness[chrom] = diagelem
    return reduced_fitness

def _min_max_diff_fitness(population, extended_pop, FCIDUMPClass, s, nel, norb, csf_list):
    """Internal function that calculates fitness based on min-max difference."""
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)
        min_val = float('inf')
        max_val = float('-inf')
        for csf_stepvec in csf_list:
            diagelem = GUGAClass.calc_diag_elem(csf_stepvec, add_core=True)
            min_val = min(min_val, diagelem)
            max_val = max(max_val, diagelem)
        reduced_fitness[chrom] = max_val - min_val
    return reduced_fitness

def _max_diagelem(population, extended_pop, FCIDUMPClass, s, nel, norb, csf_list):
    """
    Internal function that calculates fitness based the max diagonal element
    using csf_list.
    """
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)
        max_val = float('-inf')
        for csf_stepvec in csf_list:
            diagelem = GUGAClass.calc_diag_elem(csf_stepvec, add_core=True)
            max_val = max(max_val, diagelem)
        reduced_fitness[chrom] = max_val
    return reduced_fitness

def _min_diagelem(population, extended_pop, FCIDUMPClass, s, nel, norb, csf_list):
    """
    Internal function that calculates fitness based the min diagonal element
    using csf_list.
    """
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)
        min_val = float('inf')
        for csf_stepvec in csf_list:
            diagelem = GUGAClass.calc_diag_elem(csf_stepvec, add_core=True)
            min_val = min(min_val, diagelem)
        reduced_fitness[chrom] = min_val * -1
    return reduced_fitness

def _neel_fast_diag(population, extended_pop, J, ref_dict, tMin):
    reduced_fitness = {}
    J = np.array(J)
    for chrom, extended_chrom  in zip(population, extended_pop):
        sorting_arr = np.array(extended_chrom) - 1
        J_reordered = J[np.ix_(sorting_arr, sorting_arr)]

        csf_stepvec = [ref_dict[i] for i in extended_chrom]
        X = X_matrix(csf_stepvec)
        N = N_matrix(csf_stepvec)

        if tMin:
            reduced_fitness[chrom] = -np.sum(N + J_reordered * X)
        else:
            reduced_fitness[chrom] = np.sum(N + J_reordered * X)
    return reduced_fitness

def _max_fast_diag(population, extended_pop, J, csf_list):
    reduced_fitness = {}
    J = np.array(J)
    for chrom, extended_chrom  in zip(population, extended_pop):
        sorting_arr = np.array(extended_chrom) - 1
        J_reordered = J[np.ix_(sorting_arr, sorting_arr)]
        max_val = float('-inf')
        for csf_stepvec in csf_list:

            X = X_matrix(csf_stepvec)
            N = N_matrix(csf_stepvec)
            diagelem = -np.sum(N + J_reordered * X)
            print(chrom, diagelem, csf_stepvec)
            max_val = max(max_val, diagelem)

        reduced_fitness[chrom] = max_val

    return reduced_fitness

def _min_fast_diag(population, extended_pop, J, csf_list):
    reduced_fitness = {}
    J = np.array(J)
    for chrom, extended_chrom  in zip(population, extended_pop):
        sorting_arr = np.array(extended_chrom) - 1
        J_reordered = J[np.ix_(sorting_arr, sorting_arr)]
        min_val = float('inf')
        for csf_stepvec in csf_list:

            X = X_matrix(csf_stepvec)
            N = N_matrix(csf_stepvec)
            diagelem = -np.sum(N + J_reordered * X)
            min_val = min(min_val, diagelem)

        if min_val > 0.0:
            # Since we **maximize** the fitness, we change the sign of min_val
            # assuming the original fitness score is always negative.
            raise ValueError(f'The method works for negative diagE. diagE = {min_val}')

        reduced_fitness[chrom] = -min_val

    return reduced_fitness

def _neel_fast_diag_min_openshell(population, extended_pop, J, ref_dict):
    reduced_fitness = {}
    J = np.array(J)
    for chrom, extended_chrom  in zip(population, extended_pop):
        sorting_arr = np.array(extended_chrom) - 1
        J_reordered = J[np.ix_(sorting_arr, sorting_arr)]

        csf_stepvec = [ref_dict[i] for i in extended_chrom]
        X = X_matrix_openshell_only(csf_stepvec)

        reduced_fitness[chrom] = np.sum(J_reordered * X) * -1
    return reduced_fitness

#------------------------------------------------------------------------------#
def J_mat_from_fcidump(fcidump_file, norb):
    """
    Generate the J matrix from the FCIDUMP file.
    We onlny use exchange integrals.

    Args:
        fcidump_file (str): Path to the FCIDUMP file.

    Returns:
        np.ndarray: J matrix.
    """
    reader = read_fcidump.FCIDUMPReader(fcidump_file)
    J = np.zeros((norb, norb))
    for i in range(norb):
        for j in range(i + 1, norb):
            J[i, j] = reader.get_integral(j + 1, i + 1, j + 1, i + 1)
            J[j, i] = reader.get_integral(j + 1, i + 1, j + 1, i + 1)

    return J

def N_matrix(d_vec):
    n = len(d_vec)
    N = np.zeros((n, n))
    f = lambda idx: d_vec[idx] - d_vec[idx] // 2
    for i, j in zip(*np.triu_indices(n, k=1)):
        N[i, j] = f(i) * f(j)
        N[j, i] = N[i, j]

    return N

def X_matrix(d_vec):
    N = len(d_vec)

    db = lambda dv: dv - 3 * (dv // 2)
    A2 = lambda b, x, y: (b + x) / (b + y)

    d_vec = np.array(d_vec)
    b_vec = np.cumsum(db(d_vec))
    X = np.zeros((N, N))

    def f(b, d):
        if d == 1:
            return A2(b, 2, 0) * A2(b, -1, 1)
        elif d == 2:
            return A2(b, 0, 2) * A2(b, 3, 1)
        else:
            return 1

    for i, j in zip(*np.triu_indices(N, k=1)):
        di, dj = d_vec[i], d_vec[j]
        bi, bj = b_vec[i], b_vec[j]
        if di == 1 and dj == 1:
            Xij = A2(bi, 2, 0) * A2(bj, -1, 1)
        elif di == 1 and dj == 2:
            Xij = A2(bi, 2, 0) * A2(bj, 3, 1)
        elif di == 2 and dj == 1:
            Xij = A2(bi, 0, 2) * A2(bj, -1, 1)
        elif di == 2 and dj == 2:
            Xij = A2(bi, 0, 2) * A2(bj, 3, 1)
        else:
            Xij = 0

        if Xij != 0:
            for k in range(i + 1, j):
                Xij *= f(b_vec[k], d_vec[k])
            Xij = np.sqrt(Xij)

        if di * dj == 2:
            Xij *= -1.0

        X[i, j] = Xij
        X[j, i] = Xij

    np.fill_diagonal(X, 0)

    return X

def X_matrix_openshell_only(d_vec):
    """
    Compute the matrix X_ij for a given step-vector d.
    See equation (A.10) in the appendix (A.2) of Werner Dobrautz's phd thesis.
    The exact form used here will be found in our paper (update this once the
    paper is published).
    **Note that this assumes the step-vector only contains 1's and 2's.**

    Args:
        d_vec (list): step-vector not including 0 and 3

    Returns:
        N x N array where X[i][j] corresponds to the computed value for i < j.
    """

    N = len(d_vec)
    d_vec = np.array(d_vec)
    s_vec = 2 * (d_vec == 1) - 1
    b_vec = np.cumsum(s_vec)

    A = (b_vec - 2 * d_vec + 4) / (b_vec + 2 * d_vec - 2)
    B = (b_vec + 4 * d_vec - 5) / (b_vec + 1)
    X = np.ones((N, N))

    for i, j in zip(*np.triu_indices(N, k=1)):
        X[i, j] = X[i, j - 1] * A[j - 1] * B[j] 
    X = np.sqrt(X)
    np.fill_diagonal(X, 0)
    for i, j in zip(*np.triu_indices(N, k=1)):
        X[i, j] *= s_vec[i] * s_vec[j]
        X[j, i] = X[i, j]

    return X
#------------------------------------------------------------------------------#

def calculate_fitness(method: FitnessFunction, POPClass, FCIDUMPClass, s, nel, norb, **kwargs):
    """
    Wrapper function to calculate fitness using specified method.

    Args:
        method (FitnessFunction): The method to use for fitness score evaluation
        POPClass: Population class instance
        FCIDUMPClass: FCIDUMP class instance
        s: Spin
        nel: Number of electrons
        norb: Number of orbitals
        **kwargs: Optional arguments depending on method:
            - For MIN_MAX_DIFF: 
                csf_list (list): List of CSFs to measure
            - For future methods:
                Add parameters as needed

    Returns:
        dict: Dictionary mapping chromosomes to their fitness values

    Raises:
        ValueError: If required parameters are missing for chosen method
    """
    on_site_permutation = kwargs.get('on_site_permutation', None)
    num_prefix = kwargs.get('num_prefix', 0)
    num_suffix = kwargs.get('num_suffix', 0)
    ref_dict = POPClass.ref_dict
    csf_list = kwargs.get('csf_list', None)
    tHeisenberg = kwargs.get('tHeisenberg', False)
    J = kwargs.get('J', None)
    if on_site_permutation is None and (num_prefix != 0 or num_suffix != 0):
        raise ValueError("on_site_permutation is required for num_prefix and"
        " num_suffix. Use (1,) for on_site_permutation if no permutation is "
        "needed.")

    if on_site_permutation is not None:
        extended_pop = [extend_ordering.extend_ordering(ordering, 
            on_site_permutation, num_prefix, num_suffix) 
            for ordering in POPClass.current_pop]
    else:
        extended_pop = POPClass.current_pop

    if method == FitnessFunction.MIN_MAX_DIFF:
        if csf_list is None:
            raise ValueError("csf_list is required for MIN_MAX_DIFF method")
        fitness_ht = _min_max_diff_fitness(POPClass.current_pop, extended_pop, FCIDUMPClass, s, nel, norb, csf_list)
    elif method == FitnessFunction.REF_DIAGELEM:
        fitness_ht =  _gen_ref_diagelem_fitness(POPClass.current_pop, extended_pop, ref_dict, FCIDUMPClass, s, nel, norb, tHeisenberg)
    elif method == FitnessFunction.MAX_DIAG_ELEM:
        fitness_ht = _max_diagelem(POPClass.current_pop, extended_pop, FCIDUMPClass, s, nel, norb, csf_list)
    elif method == FitnessFunction.MIN_DIAG_ELEM:
        fitness_ht = _min_diagelem(POPClass.current_pop, extended_pop, FCIDUMPClass, s, nel, norb, csf_list)
    elif method == FitnessFunction.NEEL_FAST_DIAG:
        if J is None:
            raise ValueError("J matrix is required for NEEL_FAST_DIAG method")
        if ref_dict is None:
            raise ValueError("ref_dict is required for NEEL_FAST_DIAG method")
        fitness_ht = _neel_fast_diag(POPClass.current_pop, extended_pop, J, ref_dict, True)
    elif method == FitnessFunction.NEEL_FAST_DIAG_MIN:
        if J is None:
            raise ValueError("J matrix is required for NEEL_FAST_DIAG_MIN method")
        if ref_dict is None:
            raise ValueError("ref_dict is required for NEEL_FAST_DIAG_MIN method")
        fitness_ht = _neel_fast_diag(POPClass.current_pop, extended_pop, J, ref_dict, True)
    elif method == FitnessFunction.NEEL_FAST_DIAG_MAX:
        if J is None:
            raise ValueError("J matrix is required for NEEL_FAST_DIAG_MAX method")
        if ref_dict is None:
            raise ValueError("ref_dict is required for NEEL_FAST_DIAG_MAX method")
        fitness_ht = _neel_fast_diag(POPClass.current_pop, extended_pop, J, ref_dict, False)
    elif method == FitnessFunction.MAX_FAST_DIAG:
        if J is None:
            raise ValueError("J matrix is required for MAX_FAST_DIAG method")
        if csf_list is None:
            raise ValueError("csf_list is required for MAX_FAST_DIAG method")
        fitness_ht = _max_fast_diag(POPClass.current_pop, extended_pop, J, csf_list)
    elif method == FitnessFunction.MIN_FAST_DIAG:
        if J is None:
            raise ValueError("J matrix is required for MIN_FAST_DIAG method")
        if csf_list is None:
            raise ValueError("csf_list is required for MIN_FAST_DIAG method")
        fitness_ht = _min_fast_diag(POPClass.current_pop, extended_pop, J, csf_list)
    elif method == FitnessFunction.NEEL_FAST_DIAG_MIN_OSONLY:
        if J is None:
            raise ValueError("J matrix is required for NEEL_FAST_DIAG_MIN_OSONLY method")
        if ref_dict is None:
            raise ValueError("ref_dict is required for NEEL_FAST_DIAG_MIN_OSONLY method")
        fitness_ht = _neel_fast_diag_min_openshell(POPClass.current_pop, extended_pop, J, ref_dict)
    else:
        raise ValueError(f"Unknown fitness method: {method}")

    return fitness_ht
