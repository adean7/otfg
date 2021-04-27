import numba

import input_output as io


def number_of_nodes(func, nmin, nmax):
    if nmin   < 1        : io.abort('Error: number_of_nodes: nmin < 1')
    if nmin+1 > len(func): io.abort('Error: number_of_nodes: nmin+1>len(func)')
    if nmax-1 < 1        : io.abort('Error: number_of_nodes: nmax-1 < 1')
    if nmax   > len(func): io.abort('Error: number_of_nodes: nmax>size(func)')

    return number_of_nodes_calc(func, nmin, nmax)

@numba.njit()
def number_of_nodes_calc(func, nmin, nmax):
    num_nodes = 0

    for i in range(nmin+1, nmax):
        if func[i-1] * func[i] < 0.0:
            num_nodes += 1

    return num_nodes