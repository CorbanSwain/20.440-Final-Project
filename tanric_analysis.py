#!/bin/python3
# utils.py
# Corban Swain, 2018

import os
import numpy as np
import matplotlib.pyplot as plt
from utils import *
from scipy.stats import ttest_ind

MIN_NORMAL_SAMPLES = 5


def import_all_data():
    tanric_dir = os.path.join('data', 'tanric_data')
    names = []
    with open(os.path.join(tanric_dir, 'names.txt')) as names_file:
        for line in names_file.readlines():
            if not line.find('TCGA') is -1:
                names.append(line.strip(' \n'))

    datasets = []
    for name in names:  # TODO - this can be parallelized
        metafile = os.path.join(tanric_dir, name + '-rnaexpr-META.tsv')
        tds = TanricDataset(file2dict(metafile))

        # criteria for including datasets
        if tds.n_normal_samples < MIN_NORMAL_SAMPLES:
            continue

        cachefile = os.path.join(tanric_dir, 'np_cache', name + '.npy')
        try:
            data = np.load(cachefile)
        except FileNotFoundError:
            datafile = os.path.join(tanric_dir, name + '-rnaexpr.tsv')
            # this takes a while
            data = np.genfromtxt(datafile, delimiter='\t', dtype=None,
                                 encoding=None, names=True, deletechars='')
            np.save(cachefile, data)
        tds.parse_exprdata(data)
        print('\tLoaded in dataset for %s (%3d, %3d)'
              % (tds.cancer_type, tds.n_normal_samples, tds.n_tumor_samples))
        datasets.append(tds)
    return datasets


def perform_t_test(datasets, expr_cutoff=0.3, procedure='crit', **kwargs):
    # FIXME - May need to do some memory management here
    find_signif = mh_tests[procedure]
    for ds in datasets:
        print('\t' + ds.cancer_type)
        valid = np.any(ds.exprdata > expr_cutoff, 1)
        norm_valid = ds.normal_samples[valid]
        tumor_valid = ds.tumor_samples[valid]
        t, p = ttest_ind(norm_valid, tumor_valid, axis=1)
        signif = find_signif(p, **kwargs)
        print('\t\tnum signif = %d' % len(signif))


if __name__ == "__main__":
    print('1-Beginning Data Import')
    datasets = import_all_data()

    print('2-Performing t-tests')
    perform_t_test(datasets, procedure='bonferoni')

