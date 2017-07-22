#!/usr/bin/env python
import os
import numpy as np


def main():
    """
    Run all inla tmb comaprisons
    """
    rhos = np.arange(.01, 1, .1)
    sigmas = np.arange(.01, 2, .2)
    ranges = np.arange(.01, 1, .1)
    Ns = np.array([50, 100, 500, 1000, 10000])
    file_dir = os.path.dirname(os.path.realpath(__file__))
    execf_ = os.path.join(file_dir, "compare_2D.R")
    qsub_template = (
            "qsub -b y -pe multi_slot 10 -now no -P proj_forecasting "
            "-N {r}{s}{p}{n} {exec_file} --rho {p} -N {n} --sigma {s} --range {r}"
            )
    for p in rhos:
        for s in sigmas:
            for r in ranges:
                for n in Ns:
                    qsub = qsub_template.format(n=n, r=r, s=s, p=p)
                    print qsub
                    os.popen(qsub)

if __name__ == "__main__":
    main()
