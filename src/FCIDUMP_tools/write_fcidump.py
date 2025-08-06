"""
collection of functions for writing an FCIDUMP file
"""

def main(fcidumpname, nel, norb, ms, integral_list):
    "integral list has to be [[x1 i1 j1 k1 l1], [x2 i2 j2 k2 l2], ...]"

    header = gen_header(nel, norb, ms)

    with open(fcidumpname, 'w') as f:
        f.write(header)
        for ints in integral_list:
            f.write(gen_int_manual(*ints))

def gen_header(nel, norb, ms):
    """
    Generates the header for a manual fcidump file
    Currently, it is assumed that there is no point group symmetry
    """

    header = f""" &FCI NORB=  {norb}, NELEC=  {nel}, MS2=  {ms},
  ORBSYM= {', '.join(['1'] * norb)}
  ISYM=0
 &END\n"""

    return header
    
def gen_int_manual(x, i, j, k, l, digits=12):
    integral = f"""    {x:.{digits}f}    {i}    {j}    {k}    {l}\n"""
    return integral


if __name__ == '__main__':
    nel = 6
    norb = 6
    ms = 0
    fcidumpname = 'FCIDUMP'
    integral_list = [[1.0000, 1, 2, 2, 1], [1.0000, 1, 3, 3, 1],
                     [2.0000 , 1, 4, 4, 1], [2.0000, 2, 3, 3, 2]]

    main(fcidumpname, nel, norb, ms, integral_list)