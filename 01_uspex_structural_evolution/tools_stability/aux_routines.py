import numpy as np
import sys

sys.path.append("..")
from genutils.format_tools import stoich2formula


def bin_paths(stoich):
    """
    :param stoich: stoichiometry list e.g. [3, 6] for C3H6 etc.
    :return: (left_fragments_ndarray, right_fragments_ndarray)
    """
    if np.sum(stoich) == 0:
        return [list(stoich)], [list(stoich)]
    elif np.sum(stoich) == 1:
        return [list(stoich)], [list(np.zeros(stoich.shape))]
    variable_bases = [x + 1 for x in stoich]
    moduli = [np.prod(variable_bases[k:]) for k in range(len(stoich))]
    moduli = moduli[1:] + [1]

    total = int(np.ceil(np.prod(variable_bases) / 2))
    odometer = []

    for k in range(1, total):
        q, p = divmod(k, moduli[0])
        r = [q]
        for i in range(1, len(stoich)):
            q, p = divmod(p, moduli[i])
            r.append(q)
        odometer.append(r)
    lhs_stoich = np.array(odometer)
    rhs_stoich = np.array(stoich) * np.ones(lhs_stoich.shape) - lhs_stoich

    return [list(x) for x in lhs_stoich], [list(x) for x in rhs_stoich]


def table2list_fmt(v_array, h_array, en_table):
    '''
    Turns 2D table into the list-formatted data. List-formatted data here is the data, represented as
    [[v0, h0, val(i0, j0)],
     [v1, h1, val(i1, j1)],
     ...]
    The first two indices will serve as indices for 2D matrix of values. The missing values
    will be replaced by placeholder.
    @param v_array: array of values corresponding to rows of matrix (e.g. x-coordinate)
    @param h_array: array of values corresponding to columns of matrix (e.g. y-coordinate)
    @param en_table: 2D array of data
    @return: list of lists as described above
    '''
    en_table = np.array(en_table)
    num_total = len(v_array) * len(h_array)
    assert num_total == len(en_table.flatten())
    list_fmt_out = []
    # output_dummy = np.zeros((num_total, 3))
    for v_idx in v_array:
        for h_idx in h_array:
            if (v_idx > 0) or (h_idx > 0):
                list_fmt_out.append([v_idx, h_idx, en_table[v_array == v_idx, h_array == h_idx][0]])
    return np.array(list_fmt_out)


def list_fmt2table(list_fmt_data, outfile='', placeholder=0.0):
    '''
    Turns list-formatted data to table. List-formatted data here is the data, represented as
    [[i0, j0, val(i0, j0)],
     [i1, j1, val(i1, j1)],
     ...]
    The first two indices will serve as indices for 2D matrix of values. The missing values
    will be replaced by placeholder.
    @param list_fmt_data:
    @param outfile:
    @param placeholder:
    @return:
    '''
    npdata = np.array(list_fmt_data)
    v_array = [int(x) for x in npdata[:, 0]]
    h_array = [int(x) for x in npdata[:, 1]]
    max_v = np.max(v_array[:])
    min_v = np.min(v_array[:])
    max_h = np.max(h_array[:])
    min_h = np.min(h_array[:])
    out_table = np.full((max_v - min_v + 1, max_h - min_h + 1), placeholder)
    for entry in npdata:
        v_idx = int(entry[0] - min_v)
        h_idx = int(entry[1] - min_h)
        out_table[v_idx, h_idx] = float(entry[2])
    if len(outfile) > 0:
        np.savetxt(outfile, out_table, fmt='%.6e')
    outv = np.array([int(x) for x in np.unique(v_array)])
    outh = np.array([int(x) for x in np.unique(h_array)])
    return out_table, outv, outh
