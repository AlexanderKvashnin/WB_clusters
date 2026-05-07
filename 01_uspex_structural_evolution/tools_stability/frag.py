import numpy as np
import sys
sys.path.append("..")
from genutils.format_tools import stoich2formula
from .aux_routines import bin_paths


def fragmentation(fmt_data, atomlabels, outfile):
    avail_stoich = fmt_data[:,:-1].tolist()
    all_en = fmt_data[:,-1]
    fmt_frag = []
    with open(outfile, 'w') as out_fid:
        for entry in fmt_data:
            current_stoich = entry[:-1]
            current_formula = stoich2formula(current_stoich, atomlabels)
            cl_en = entry[-1]
            l_frags, r_frags = bin_paths(current_stoich)
            reactives = []
            react_ens = []
            for l_frag, r_frag in zip(l_frags, r_frags):
                if (l_frag in avail_stoich) and (r_frag in avail_stoich):
                    reactives.append({'l':l_frag, 'r':r_frag})
                    l_idx = avail_stoich.index(l_frag)
                    r_idx = avail_stoich.index(r_frag)
                    e_react = all_en[l_idx] + all_en[r_idx] - cl_en # E_parts - E_whole
                    react_ens.append([e_react])
            if len(react_ens) > 0:
                e_frag_min = np.min(react_ens)
                min_idx = react_ens.index(e_frag_min)
                l_formula = stoich2formula(reactives[min_idx]['l'], atomlabels)
                r_formula = stoich2formula(reactives[min_idx]['r'], atomlabels)
                out_fid.write('%s -> %s + %s, E_frag = %6.3f\n'%(current_formula, l_formula, r_formula, e_frag_min) )
                fmt_frag.append(list(current_stoich) + [e_frag_min])
    return fmt_frag





