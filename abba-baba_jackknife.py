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
            counts[pops].append(map(float, fields[6:]))
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
    F4sum=9
    Ddensum=10
    F4bc=11

    print("P1\tP2\tP3\tP4\tD\tD_STD\tD_Z\tDbc\tDbc_STD\tDbc_Z\tD2\tD2_STD\tD2_Z\tD2bc\tD2bc_STD\tD2bc_Z\tF4\tF4_STD\tF4_Z\tF4bc\tF4bc_STD\tF4bc_Z\tABBA\tBABA\tNSITES")

    def do_jack(xx, weights):
        n = len(xx)
        mean = np.sum(w*x for x,w in zip(xx,weights)) / float(n-1)
        var = np.sum(w*(x-mean)**2 for x,w in zip(xx,weights)) * float(n-1) / n
        std = np.sqrt(var)
        if std == 0:
            Z = float('nan')
        else:
            Z = mean/std
        return mean,std,Z

    for quad in quadlist:
        m = counts[quad]
        sums = np.sum(m, axis=0)
        D_jack = []
        Dbc_jack = [] # branch corrected
        bc_jack = []
        f4_jack = []
        f4bc_jack = []
        D2_jack = []
        D2bc_jack = []
        weights = []
        for row in m:
            jrow = [(sums[i]-row[i]) for i in range(F4bc+1)]
            numer = jrow[ABBA] - jrow[BABA]
            denom = jrow[ABBA] + jrow[BABA]
            if denom == 0 or jrow[Ddensum] == 0:
                continue

            if jrow[AAAA]:
                # ABBA-BABA branch length correction, see Green et al. SOM pg138-139.
                bc = (jrow[AABA]-jrow[AAAB])*(jrow[ABAA]-jrow[BAAA]) / float(jrow[AAAA])
                bc_abba = (jrow[AABA]*jrow[ABAA] + jrow[AAAB]*jrow[BAAA]) / float(jrow[AAAA])
                bc_baba = (jrow[AABA]*jrow[BAAA] + jrow[AAAB]*jrow[ABAA]) / float(jrow[AAAA])
            else:
                bc = 0

            weight = float(jrow[nsites]) / sums[nsites]
            weights.append(weight)

            # ABBA-BABA (randomly sample allele from population)
            D_jack.append(float(numer)/denom)
            Dbc_jack.append(float(numer-bc)/(denom-bc_abba-bc_baba))

            # ABBA-BABA branch correction
            bc_jack.append(bc)

            # D statistic (using allele frequencies)
            D2_jack.append(jrow[F4sum] / jrow[Ddensum])
            D2bc_jack.append((jrow[F4sum] -jrow[F4bc]) / jrow[Ddensum])

            # F4
            f4_jack.append(jrow[F4sum] / jrow[nsites])
            f4bc_jack.append((jrow[F4sum] -jrow[F4bc]) / jrow[nsites])

        D_mean, D_std, D_Z = do_jack(D_jack, weights)
        Dbc_mean, Dbc_std, Dbc_Z = do_jack(Dbc_jack, weights)
        D2_mean, D2_std, D2_Z = do_jack(D2_jack, weights)
        D2bc_mean, D2bc_std, D2bc_Z = do_jack(D2bc_jack, weights)
        f4_mean, f4_std, f4_Z = do_jack(f4_jack, weights)
        f4bc_mean, f4bc_std, f4bc_Z = do_jack(f4bc_jack, weights)

        #bc,_,_,_ = do_jack(bc_jack, weights)

        print(quad[0], quad[1], quad[2], quad[3],
                D_mean, D_std, D_Z,
                Dbc_mean, Dbc_std, Dbc_Z,
                D2_mean, D2_std, D2_Z,
                D2bc_mean, D2bc_std, D2bc_Z,
                f4_mean, f4_std, f4_Z,
                f4bc_mean, f4bc_std, f4bc_Z,
                int(sums[ABBA]), int(sums[BABA]), int(sums[nsites]),
                sep="\t")
