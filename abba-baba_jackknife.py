#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy as np

try:
    range = xrange
except NameError:
    pass

def parse_tsv(fn):
    counts = {}
    quadlist = []
    with open(fn) as f:
        next(f) # skip header
        for line in f:
            fields = line.split()
            pops = tuple(fields[2:6])
            if pops not in counts:
                counts[pops] = []
                quadlist.append(pops)
            counts[pops].append(map(int, fields[6:]))
    return quadlist, counts

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} abba-baba.txt".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    txt_fn = sys.argv[1]
    quadlist, counts = parse_tsv(txt_fn)

    AAAA=0
    AAAB=1
    AABA=2
    ABAA=3
    BAAA=4
    BBAA=5
    ABBA=6
    BABA=7
    nsites=8

    print("P1\tP2\tP3\tP4\tD\tSTD\tZ\tDbc\tSTDbc\tZbc\tABBA\tBABA\tNSITES\tBC")

    for quad in quadlist:
        m = counts[quad]
        sums = np.sum(m, axis=0)
        D_jack = []
        D_jack_bc = [] # branch corrected
        bc_jack = []
        weights = []
        for row in m:
            jrow = [(sums[i]-row[i]) for i in range(9)]
            numer = jrow[ABBA] - jrow[BABA]
            denom = jrow[ABBA] + jrow[BABA]
            if denom == 0:
                continue

            if jrow[AAAA]:
                # branch length correction, see Green et al. SOM pg138-139.
                bc = (jrow[AABA]-jrow[AAAB])*(jrow[ABAA]-jrow[BAAA]) / float(jrow[AAAA])
            else:
                bc = 0

            weight = float(jrow[nsites]) / sums[nsites]
            weights.append(weight)
            D_jack.append(float(numer)/denom)
            D_jack_bc.append(float(numer-bc)/denom)
            bc_jack.append(bc)

        n = len(D_jack)
        D_mean = np.sum(w*d for d,w in zip(D_jack,weights)) / float(n-1)
        D_var = np.sum(w*(d-D_mean)**2 for d,w in zip(D_jack,weights)) * float(n-1) / n
        D_std = np.sqrt(D_var)
        if D_std == 0:
            Z = float('nan')
        else:
            Z = D_mean/D_std

        D_mean_bc = np.sum(w*d for d,w in zip(D_jack_bc,weights)) / float(n-1)
        D_var_bc = np.sum(w*(d-D_mean_bc)**2 for d,w in zip(D_jack_bc,weights)) * float(n-1) / n
        D_std_bc = np.sqrt(D_var_bc)
        if D_std_bc == 0:
            Z_bc = float('nan')
        else:
            Z_bc = D_mean_bc/D_std_bc

        bc = np.sum(w*e for e,w in zip(bc_jack,weights)) / float(n-1)

        print(quad[0], quad[1], quad[2], quad[3], D_mean, D_std, Z, D_mean_bc, D_std_bc, Z_bc, sums[ABBA], sums[BABA], sums[nsites], bc, sep="\t")
