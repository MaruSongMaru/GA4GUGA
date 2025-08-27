"""
FCIDUMP permutation utilities for orbital reordering.

Copyright (c) 2025 Maru Song
"""

import itertools

def permgen(norb: int, iterator: bool = False):
    """
    Generate all permutations of orbitals.

    Args:
        norb (int): Number of orbitals.
        iterator (bool): If True, return an iterator; if False, return a list.

    Returns:
        iterator or list: All permutations of range(1, norb + 1).
    """
    if iterator:
        return itertools.permutations(range(1, norb + 1))
    else:
        return [p for p in itertools.permutations(range(1, norb + 1))]

def expand_perm_multielec(perm: tuple[int], multiplier: int):
    """
    Expand permutation to account for multiple electrons per site.

    Args:
        perm (tuple[int]): Original permutation.
        multiplier (int): Number of electrons per site.

    Returns:
        tuple[int]: Expanded permutation with multiplier electrons on each site.

    Example:
        >>> expand_perm_multielec((1, 4, 3, 2), 2)
        (1, 2, 7, 8, 5, 6, 3, 4)
    """
    return tuple(x * multiplier - i \
                 for x in perm for i in range(multiplier - 1, -1, -1))

def convert_perm_rep(perm: tuple[int]) -> tuple[int]:
    """
    Convert a permutation between "orbital representation" and "index representation".

    Args:
        perm (tuple[int]): Permutation to convert.

    Returns:
        tuple[int]: Converted permutation.

    Note:
        See `convert_fcidump_idx` for details on representations.
    """
    return(tuple([perm.index(i + 1) + 1 for i in range(len(perm))]))

def convert_fcidump_idx(fcidump_in: str, fcidump_out: str, permutation: tuple[int], orb_rep: bool = True):
    """
    Convert orbital indices in FCIDUMP file according to given permutation.

    Args:
        fcidump_in (str): Input FCIDUMP filename.
        fcidump_out (str): Output FCIDUMP filename.
        permutation (tuple[int]): Permutation to apply to orbital indices.
        orb_rep (bool): If True, interpret permutation as "orbital representation";
                       if False, interpret as "index representation". Defaults to True.

    Note:
        At the moment, the header is not converted and only the integral part is changed.

        Orbital representation: Numbers in permutation represent orbital labels in the given FCIDUMP,
        positions of the numbers represent the indices in FCIDUMP. For example, (4, 1, 2, 3) converts original indices
        {1, 2, 3, 4} to {2, 3, 4, 1}.

        Index representation: Numbers in permutation represent indices in the given FCIDUMP,
        positions represent orbitals. The same permutation (4, 1, 2, 3) converts original 
        indices {1, 2, 3, 4} to {4, 1, 2, 3}.

    Raises:
        Exception: If permutation length doesn't match number of orbitals in FCIDUMP.
        Exception: If there's an error in processing integral lines.
    """
    if not orb_rep:
        permutation = (0, ) + convert_perm_rep(permutation)
    else:
        permutation = (0, ) + permutation

    with open(fcidump_out, "w") as f:
        source = iter(open(fcidump_in, "r"))

        # Copy the header and take the number of orbitals
        for line in source:
            f.write(line)
            if 'NORB' in line:
                l = line.replace(',', ' ').replace('=', ' ').split()
                norb = int(l[l.index('NORB') + 1])
                if norb != len(permutation) - 1:
                    raise Exception("The given permutation length ({})is not the same as \
                                     the number of orbitals of the FCIDUMP.".format(norb))
            if ('/'  in line) or ('END' in line):
                break

        # Convert integral indices
        for line in source:
            x = line.split()[0]
            try:
                i = str(permutation.index(int(line.split()[1])))
                j = str(permutation.index(int(line.split()[2])))
                k = str(permutation.index(int(line.split()[3])))
                l = str(permutation.index(int(line.split()[4])))
            except ValueError:
                raise Exception("Something's wrong in", line)

            # Making the output molprolike
            if int(i) < int(j):
                i, j = j, i
            if int(k) < int(l):
                k, l = l, k
            if int(i) < int(k) or (int(i) == int(k) and int(j) < int(l)):
                i, j, k, l = k, l, i, j

            new_line = (x, i, j, k, l)

            # Just to pritify the output file
            if x[0] == '-':
                f.write(' ')
            else:
                f.write('  ')
            f.write('   '.join(new_line))
            f.write('\n')

        source.close()
